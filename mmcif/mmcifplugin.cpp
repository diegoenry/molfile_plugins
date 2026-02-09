/*
 * mmCIF Reader Plugin for VMD
 *
 * Modern read-only plugin for PDBx/mmCIF files using a built-in CIF parser.
 * Parses _atom_site, _cell, and _struct_conn categories to provide atom
 * structure, coordinates, unit cell parameters, and explicit bonds (disulfides,
 * covalent links, metal coordination).
 *
 * Multi-model support: NMR ensembles loaded as trajectory frames.
 * Chain filtering: via filename query param (?chain=A,B) or MMCIF_CHAIN env var.
 *
 * Author: Diego Enry Barreto Gomes
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <tuple>

#include "cifparse.h"
#include "periodic_table.h"
#include "molfile_plugin.h"

/* ========================================================================
 * Plugin handle - stores all parsed data
 * ======================================================================== */

struct mmcif_data {
    char *filepath;            /* actual file path (query params stripped) */
    int natoms;
    int nmodels;
    int current_model;         /* timestep counter for read_next_timestep */

    /* Per-atom structure data (from model 1) */
    std::vector<std::string> atom_names;     /* _atom_site.label_atom_id    */
    std::vector<std::string> type_symbols;   /* _atom_site.type_symbol      */
    std::vector<std::string> resnames;       /* _atom_site.label_comp_id    */
    std::vector<int>         resids;         /* _atom_site.auth_seq_id      */
    std::vector<std::string> chains;         /* _atom_site.auth_asym_id     */
    std::vector<std::string> segids;         /* _atom_site.label_asym_id    */
    std::vector<std::string> altlocs;        /* _atom_site.label_alt_id     */
    std::vector<std::string> insertions;     /* _atom_site.pdbx_PDB_ins_code */
    std::vector<float>       occupancies;
    std::vector<float>       bfactors;
    std::vector<float>       charges;
    std::vector<int>         atomicnumbers;

    /* Coordinates: model 1 stored separately, models 2..N in vector */
    std::vector<float>                first_coords;   /* natoms*3 */
    std::vector<std::vector<float>>   model_coords;   /* models 2..N */

    /* Bonds from _struct_conn */
    int nbonds;
    int *bond_from;            /* 1-based */
    int *bond_to;              /* 1-based */
    float *bond_order;

    /* Unit cell from _cell */
    float cell[6];             /* A, B, C, alpha, beta, gamma */
    bool has_cell;

    /* Chain filter for subset loading */
    std::set<std::string> filter_chains;

    /* Atom lookup for bond resolution: (chain, resid, atom_name) -> 0-based */
    std::map<std::tuple<std::string,int,std::string>, int> atom_lookup;

    mmcif_data() : filepath(nullptr), natoms(0), nmodels(0), current_model(0),
                   nbonds(0), bond_from(nullptr), bond_to(nullptr),
                   bond_order(nullptr), has_cell(false) {
        cell[0] = cell[1] = cell[2] = 0.0f;
        cell[3] = cell[4] = cell[5] = 90.0f;
    }

    ~mmcif_data() {
        if (filepath) free(filepath);
        if (bond_from) free(bond_from);
        if (bond_to) free(bond_to);
        if (bond_order) free(bond_order);
    }
};


/* ========================================================================
 * Column access helpers
 * ======================================================================== */

/* Return column value as string, treating "." and "?" as empty */
static std::string col_str(const std::vector<std::string> &row, int col) {
    if (col < 0 || col >= (int)row.size()) return std::string();
    const std::string &v = row[col];
    if (v == "." || v == "?") return std::string();
    return v;
}

static float col_float(const std::vector<std::string> &row, int col) {
    if (col < 0 || col >= (int)row.size()) return 0.0f;
    const std::string &v = row[col];
    if (v == "." || v == "?") return 0.0f;
    return strtof(v.c_str(), nullptr);
}

static int col_int(const std::vector<std::string> &row, int col) {
    if (col < 0 || col >= (int)row.size()) return 0;
    const std::string &v = row[col];
    if (v == "." || v == "?") return 0;
    return atoi(v.c_str());
}


/* ========================================================================
 * CIF table processing functions
 * ======================================================================== */

static void process_cell(const cif::Table *table, mmcif_data *data) {
    if (!table || table->rows.empty()) return;

    int c_a     = table->find_column("length_a");
    int c_b     = table->find_column("length_b");
    int c_c     = table->find_column("length_c");
    int c_alpha = table->find_column("angle_alpha");
    int c_beta  = table->find_column("angle_beta");
    int c_gamma = table->find_column("angle_gamma");

    const auto &row = table->rows[0];
    float a = col_float(row, c_a);
    float b = col_float(row, c_b);
    float c = col_float(row, c_c);
    float alpha = (c_alpha >= 0) ? col_float(row, c_alpha) : 90.0f;
    float beta  = (c_beta >= 0)  ? col_float(row, c_beta)  : 90.0f;
    float gamma = (c_gamma >= 0) ? col_float(row, c_gamma) : 90.0f;

    if (a > 0 && b > 0 && c > 0) {
        data->cell[0] = a;
        data->cell[1] = b;
        data->cell[2] = c;
        data->cell[3] = alpha;
        data->cell[4] = beta;
        data->cell[5] = gamma;
        data->has_cell = true;
    }
}


static void process_atom_site(const cif::Table *table, mmcif_data *data) {
    if (!table || table->rows.empty()) return;

    /* Required columns */
    int c_x = table->find_column("Cartn_x");
    int c_y = table->find_column("Cartn_y");
    int c_z = table->find_column("Cartn_z");
    if (c_x < 0 || c_y < 0 || c_z < 0) return;

    /* Optional columns */
    int c_symbol    = table->find_column("type_symbol");
    int c_atom_name = table->find_column("label_atom_id");
    int c_comp_id   = table->find_column("label_comp_id");
    int c_asym_id   = table->find_column("label_asym_id");
    int c_auth_asym = table->find_column("auth_asym_id");
    int c_auth_seq  = table->find_column("auth_seq_id");
    int c_label_seq = table->find_column("label_seq_id");
    int c_altloc    = table->find_column("label_alt_id");
    int c_ins_code  = table->find_column("pdbx_PDB_ins_code");
    int c_occ       = table->find_column("occupancy");
    int c_bfactor   = table->find_column("B_iso_or_equiv");
    int c_charge    = table->find_column("pdbx_formal_charge");
    int c_model     = table->find_column("pdbx_PDB_model_num");

    /* auth_seq_id -> atom->resid (author-assigned residue number)
       also used for bond lookup (struct_conn uses auth identifiers) */
    int c_seq = (c_auth_seq >= 0) ? c_auth_seq : c_label_seq;
    bool filtering = !data->filter_chains.empty();

    int current_model_num = -1;
    std::vector<float> extra_coords;

    for (const auto &row : table->rows) {
        float x = col_float(row, c_x);
        float y = col_float(row, c_y);
        float z = col_float(row, c_z);

        std::string symbol    = col_str(row, c_symbol);
        std::string atom_name = col_str(row, c_atom_name);
        std::string comp_id   = col_str(row, c_comp_id);
        std::string asym_id   = col_str(row, c_asym_id);
        std::string auth_asym = col_str(row, c_auth_asym);
        std::string altloc    = col_str(row, c_altloc);
        std::string ins_code  = col_str(row, c_ins_code);
        int seq_id    = col_int(row, c_seq);
        int model_num = (c_model >= 0) ? col_int(row, c_model) : 1;
        float occ = (c_occ >= 0)     ? col_float(row, c_occ)     : 1.0f;
        float bf  = (c_bfactor >= 0) ? col_float(row, c_bfactor) : 0.0f;
        float chg = (c_charge >= 0)  ? col_float(row, c_charge)  : 0.0f;

        /* Detect model transitions */
        if (current_model_num < 0)
            current_model_num = model_num;
        if (model_num != current_model_num) {
            if (!extra_coords.empty()) {
                data->model_coords.push_back(std::move(extra_coords));
                extra_coords.clear();
            }
            current_model_num = model_num;
        }

        /* Chain filtering */
        const std::string &chain_val = (c_auth_asym >= 0) ? auth_asym :
                                       (c_asym_id >= 0)  ? asym_id :
                                       auth_asym;
        if (filtering && data->filter_chains.find(chain_val) == data->filter_chains.end())
            continue;

        /* Skip alternate locations beyond the first */
        if (!altloc.empty() && altloc != "A" && altloc != ".")
            continue;

        /* Model 1: record structure + coords */
        if (data->model_coords.empty() && extra_coords.empty()) {
            data->first_coords.push_back(x);
            data->first_coords.push_back(y);
            data->first_coords.push_back(z);

            data->atom_names.push_back(atom_name);
            data->type_symbols.push_back(symbol);
            data->resnames.push_back(comp_id);
            data->resids.push_back(seq_id);
            data->chains.push_back((c_auth_asym >= 0) ? auth_asym : asym_id);
            data->segids.push_back((c_asym_id >= 0) ? asym_id : std::string());
            data->altlocs.push_back(altloc);
            data->insertions.push_back(ins_code);
            data->occupancies.push_back(occ);
            data->bfactors.push_back(bf);
            data->charges.push_back(chg);

            int anum = get_pte_idx_from_string(symbol.c_str());
            data->atomicnumbers.push_back(anum);

            auto key = std::make_tuple(chain_val, seq_id, atom_name);
            int idx = (int)data->atom_names.size() - 1;
            data->atom_lookup[key] = idx;
        } else {
            extra_coords.push_back(x);
            extra_coords.push_back(y);
            extra_coords.push_back(z);
        }
    }

    /* Flush last extra model if any */
    if (!extra_coords.empty())
        data->model_coords.push_back(std::move(extra_coords));

    data->natoms = (int)data->atom_names.size();
    data->nmodels = 1 + (int)data->model_coords.size();
}


static void process_struct_conn(const cif::Table *table, mmcif_data *data) {
    if (!table || table->rows.empty()) return;

    int c_type   = table->find_column("conn_type_id");
    int c_chain1 = table->find_column("ptnr1_auth_asym_id");
    int c_seq1   = table->find_column("ptnr1_auth_seq_id");
    int c_atom1  = table->find_column("ptnr1_label_atom_id");
    int c_chain2 = table->find_column("ptnr2_auth_asym_id");
    int c_seq2   = table->find_column("ptnr2_auth_seq_id");
    int c_atom2  = table->find_column("ptnr2_label_atom_id");
    int c_order  = table->find_column("pdbx_value_order");

    /* Fall back to label_ columns if auth_ not available */
    if (c_chain1 < 0) c_chain1 = table->find_column("ptnr1_label_asym_id");
    if (c_seq1 < 0)   c_seq1   = table->find_column("ptnr1_label_seq_id");
    if (c_chain2 < 0) c_chain2 = table->find_column("ptnr2_label_asym_id");
    if (c_seq2 < 0)   c_seq2   = table->find_column("ptnr2_label_seq_id");

    if (c_atom1 < 0 || c_atom2 < 0) return;

    std::vector<int> from_list, to_list;
    std::vector<float> order_list;

    for (const auto &row : table->rows) {
        /* Only include relevant connection types */
        if (c_type >= 0) {
            std::string ct = col_str(row, c_type);
            for (auto &c : ct) c = tolower((unsigned char)c);
            if (!ct.empty() &&
                ct != "disulf" && ct != "covale" && ct != "covale_base" &&
                ct != "covale_phosph" && ct != "covale_sugar" && ct != "metalc")
                continue;
        }

        std::string ch1 = col_str(row, c_chain1);
        int s1          = col_int(row, c_seq1);
        std::string at1 = col_str(row, c_atom1);
        std::string ch2 = col_str(row, c_chain2);
        int s2          = col_int(row, c_seq2);
        std::string at2 = col_str(row, c_atom2);

        auto key1 = std::make_tuple(ch1, s1, at1);
        auto key2 = std::make_tuple(ch2, s2, at2);
        auto it1 = data->atom_lookup.find(key1);
        auto it2 = data->atom_lookup.find(key2);

        if (it1 != data->atom_lookup.end() && it2 != data->atom_lookup.end()) {
            from_list.push_back(it1->second + 1);  /* 1-based */
            to_list.push_back(it2->second + 1);

            float bo = 1.0f;
            std::string os = col_str(row, c_order);
            if (!os.empty()) {
                for (auto &c : os) c = toupper((unsigned char)c);
                if (os == "DOUB") bo = 2.0f;
                else if (os == "TRIP") bo = 3.0f;
            }
            order_list.push_back(bo);
        }
    }

    data->nbonds = (int)from_list.size();
    if (data->nbonds > 0) {
        data->bond_from  = (int *)malloc(data->nbonds * sizeof(int));
        data->bond_to    = (int *)malloc(data->nbonds * sizeof(int));
        data->bond_order = (float *)malloc(data->nbonds * sizeof(float));
        for (int i = 0; i < data->nbonds; i++) {
            data->bond_from[i]  = from_list[i];
            data->bond_to[i]    = to_list[i];
            data->bond_order[i] = order_list[i];
        }
    }
}


/* ========================================================================
 * Helper: parse chain filter from filepath query params or env var
 * ======================================================================== */

static void parse_chain_filter(mmcif_data *data, const char *filepath) {
    /* Check for ?chain=A,B in filepath */
    const char *q = strchr(filepath, '?');
    if (q) {
        const char *cp = strstr(q, "chain=");
        if (cp) {
            cp += 6;
            std::string chains_str;
            while (*cp && *cp != '&') {
                chains_str += *cp++;
            }
            /* Split by comma */
            size_t start = 0;
            for (size_t i = 0; i <= chains_str.size(); i++) {
                if (i == chains_str.size() || chains_str[i] == ',') {
                    if (i > start) {
                        data->filter_chains.insert(chains_str.substr(start, i - start));
                    }
                    start = i + 1;
                }
            }
        }
        /* Store filepath without query params */
        data->filepath = (char *)malloc(q - filepath + 1);
        memcpy(data->filepath, filepath, q - filepath);
        data->filepath[q - filepath] = '\0';
    } else {
        data->filepath = strdup(filepath);
    }

    /* Check MMCIF_CHAIN environment variable (only if no query param filter) */
    if (data->filter_chains.empty()) {
        const char *env = getenv("MMCIF_CHAIN");
        if (env && *env) {
            std::string chains_str(env);
            size_t start = 0;
            for (size_t i = 0; i <= chains_str.size(); i++) {
                if (i == chains_str.size() || chains_str[i] == ',') {
                    if (i > start) {
                        std::string c = chains_str.substr(start, i - start);
                        /* Trim whitespace */
                        while (!c.empty() && isspace(c.front())) c.erase(0, 1);
                        while (!c.empty() && isspace(c.back())) c.pop_back();
                        if (!c.empty())
                            data->filter_chains.insert(c);
                    }
                    start = i + 1;
                }
            }
        }
    }
}


/* ========================================================================
 * VMD Plugin Callbacks
 * ======================================================================== */

static void *open_mmcif_read(const char *filepath, const char *filetype,
                             int *natoms) {
    mmcif_data *data = new mmcif_data();

    /* Parse chain filter and set filepath */
    parse_chain_filter(data, filepath);

    /* Parse the CIF file */
    try {
        cif::DataBlock block = cif::parse_file(data->filepath);
        process_cell(block.find_table("cell"), data);
        process_atom_site(block.find_table("atom_site"), data);
        process_struct_conn(block.find_table("struct_conn"), data);
    } catch (std::exception &e) {
        fprintf(stderr, "mmcifplugin) Error parsing %s: %s\n",
                data->filepath, e.what());
        delete data;
        return NULL;
    }

    if (data->natoms <= 0) {
        fprintf(stderr, "mmcifplugin) No atoms found in %s\n", data->filepath);
        delete data;
        return NULL;
    }

    *natoms = data->natoms;

    printf("mmcifplugin) Opened %s: %d atoms, %d model(s)",
           data->filepath, data->natoms, data->nmodels);
    if (data->nbonds > 0)
        printf(", %d explicit bonds", data->nbonds);
    if (data->has_cell)
        printf(", unit cell: %.1f %.1f %.1f",
               data->cell[0], data->cell[1], data->cell[2]);
    if (!data->filter_chains.empty()) {
        printf(", chains:");
        for (auto &c : data->filter_chains)
            printf(" %s", c.c_str());
    }
    printf("\n");

    return data;
}


static int read_mmcif_structure(void *v, int *optflags, molfile_atom_t *atoms) {
    mmcif_data *data = (mmcif_data *)v;

    *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS |
                MOLFILE_OCCUPANCY | MOLFILE_BFACTOR | MOLFILE_CHARGE |
                MOLFILE_ALTLOC | MOLFILE_INSERTION | MOLFILE_BONDSSPECIAL;

    for (int i = 0; i < data->natoms; i++) {
        molfile_atom_t *atom = &atoms[i];

        /* Name */
        strncpy(atom->name, data->atom_names[i].c_str(), 15);
        atom->name[15] = '\0';

        /* Type = atom name (matches pdbx behaviour; element is in atomicnumber) */
        strncpy(atom->type, data->atom_names[i].c_str(), 15);
        atom->type[15] = '\0';

        /* Residue name */
        strncpy(atom->resname, data->resnames[i].c_str(), 7);
        atom->resname[7] = '\0';

        /* Residue ID */
        atom->resid = data->resids[i];

        /* Chain (max 1 char in current ABI) */
        if (!data->chains[i].empty()) {
            atom->chain[0] = data->chains[i][0];
            atom->chain[1] = '\0';
        } else {
            atom->chain[0] = '\0';
        }

        /* Segment ID */
        strncpy(atom->segid, data->segids[i].c_str(), 7);
        atom->segid[7] = '\0';

        /* Alt location */
        if (!data->altlocs[i].empty()) {
            atom->altloc[0] = data->altlocs[i][0];
            atom->altloc[1] = '\0';
        } else {
            atom->altloc[0] = '\0';
        }

        /* Insertion code */
        if (!data->insertions[i].empty()) {
            atom->insertion[0] = data->insertions[i][0];
            atom->insertion[1] = '\0';
        } else {
            atom->insertion[0] = '\0';
        }

        /* Numeric fields */
        atom->occupancy = data->occupancies[i];
        atom->bfactor = data->bfactors[i];
        atom->charge = data->charges[i];
        atom->atomicnumber = data->atomicnumbers[i];
        atom->mass = get_pte_mass(data->atomicnumbers[i]);
        atom->radius = get_pte_vdw_radius(data->atomicnumbers[i]);
    }

    return MOLFILE_SUCCESS;
}


static int read_mmcif_bonds(void *v, int *nbonds, int **fromptr, int **toptr,
                            float **bondorderptr, int **bondtype,
                            int *nbondtypes, char ***bondtypename) {
    mmcif_data *data = (mmcif_data *)v;

    if (data->nbonds > 0) {
        *nbonds = data->nbonds;
        *fromptr = data->bond_from;
        *toptr = data->bond_to;
        *bondorderptr = data->bond_order;
    } else {
        *nbonds = 0;
        *fromptr = NULL;
        *toptr = NULL;
        *bondorderptr = NULL;
    }

    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;

    return MOLFILE_SUCCESS;
}


static int read_mmcif_timestep(void *v, int natoms, molfile_timestep_t *ts) {
    mmcif_data *data = (mmcif_data *)v;

    if (data->current_model >= data->nmodels)
        return MOLFILE_EOF;

    if (ts) {
        const float *src;

        if (data->current_model == 0) {
            src = data->first_coords.data();
        } else {
            int mi = data->current_model - 1;
            if (mi >= (int)data->model_coords.size())
                return MOLFILE_EOF;
            /* Verify coord count matches */
            if ((int)data->model_coords[mi].size() != natoms * 3) {
                fprintf(stderr, "mmcifplugin) Model %d has %d coords, expected %d\n",
                        data->current_model + 1,
                        (int)data->model_coords[mi].size(), natoms * 3);
                return MOLFILE_EOF;
            }
            src = data->model_coords[mi].data();
        }

        memcpy(ts->coords, src, natoms * 3 * sizeof(float));

        /* Unit cell */
        if (data->has_cell) {
            ts->A     = data->cell[0];
            ts->B     = data->cell[1];
            ts->C     = data->cell[2];
            ts->alpha = data->cell[3];
            ts->beta  = data->cell[4];
            ts->gamma = data->cell[5];
        }
    }

    data->current_model++;
    return MOLFILE_SUCCESS;
}


static void close_mmcif_read(void *v) {
    mmcif_data *data = (mmcif_data *)v;
    delete data;
}


/* ========================================================================
 * Plugin registration
 * ======================================================================== */

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
    memset(&plugin, 0, sizeof(molfile_plugin_t));
    plugin.abiversion = vmdplugin_ABIVERSION;
    plugin.type = MOLFILE_PLUGIN_TYPE;
    plugin.name = "mmcif";
    plugin.prettyname = "mmCIF Reader (Modern)";
    plugin.author = "Diego Enry Barreto Gomes";
    plugin.majorv = 1;
    plugin.minorv = 0;
    plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
    plugin.filename_extension = "mmcif,cif";

    plugin.open_file_read = open_mmcif_read;
    plugin.read_structure = read_mmcif_structure;
    plugin.read_bonds = read_mmcif_bonds;
    plugin.read_next_timestep = read_mmcif_timestep;
    plugin.close_file_read = close_mmcif_read;

    return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
    (*cb)(v, (vmdplugin_t *)&plugin);
    return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
    return VMDPLUGIN_SUCCESS;
}
