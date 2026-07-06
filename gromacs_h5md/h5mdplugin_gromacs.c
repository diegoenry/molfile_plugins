/*
 * GROMACS H5MD Plugin for VMD
 * 
 * University of Illinois Open Source License
 *
 * Reads H5MD trajectory files produced by GROMACS 2026.0+ with the
 * gromacs_topology module extension.
 *
 * Author: Diego E. B. Gomes - Auburn University, Alabama, AL - USA 
 * E-mail: dgomes@auburn.edu 
 * Version: 0.1
 *
 * File Format:
 *   H5MD (HDF5-based Molecular Data) with GROMACS topology extension
 *   - Positions in /particles/system/position/value[frames][atoms][3]
 *   - Topology in /h5md/modules/gromacs_topology/
 *   - Units: nm (converted to Angstroms)
 *
 * Features:
 *   - Reads atomic structure from GROMACS topology module
 *   - Loads trajectory frames with positions
 *   - Optional velocities
 *   - Bond connectivity from /connectivity/bonds
 *   - Periodic boundary conditions (orthogonal and triclinic boxes)
 *
 * Limitations:
 *   - Thread-unsafe (HDF5 library limitation)
 *   - Forces are not read
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <hdf5.h>

#include "molfile_plugin.h"

/* Plugin version (experimental; matches the 0.1 module spec) */
#define H5MD_GROMACS_VERSION_MAJOR 0
#define H5MD_GROMACS_VERSION_MINOR 1

/* Unit conversion: nm to Angstroms */
#define NM_TO_ANGSTROM 10.0f

/* Error codes */
#define H5MD_SUCCESS 0
#define H5MD_ERROR -1
#define H5MD_BADFORMAT -2

/* Maximum string lengths */
#define MAX_NAME_LEN 256
#define MAX_MOLTYPE_NAME 64
#define MAX_ATOM_NAME 16
#define MAX_RESIDUE_NAME 8
#define MAX_SEGMENT_NAME 8

/*
 * Molecule Type Structure
 *
 * Represents one type of molecule in the GROMACS topology.
 * The full system is built by replicating molecule types
 * according to molecule_counts.
 */
typedef struct {
    char name[MAX_MOLTYPE_NAME];    /* Molecule type name */
    int particle_count;              /* Atoms per molecule */
    int residue_count;               /* Residues per molecule */

    /* Atomic properties (per-molecule-type) */
    int64_t *ids;                    /* Atom IDs */
    float *masses;                   /* Atomic masses */
    float *charges;                  /* Partial charges */
    int32_t *species;                /* Atomic numbers (Z) */
    int32_t *residue_ids;            /* Residue IDs within molecule */

    /* Decoded names (from lookup tables), one entry per particle */
    char **particle_names;           /* Atom names (per particle) */
    char **atom_residue_names;       /* Residue name resolved per particle */
} molecule_type_t;

/*
 * Main Plugin Data Structure
 */
typedef struct {
    /* HDF5 file handles */
    hid_t file_id;
    hid_t position_dataset;
    hid_t velocity_dataset;
    hid_t box_dataset;
    hid_t time_dataset;

    /* Atom and frame counts */
    int natoms;                      /* After selection */
    int natoms_file;                 /* Before selection */
    int nframes;
    int current_frame;

    /* Topology data */
    molecule_type_t *molecule_types; /* Array of molecule types */
    int *molecule_counts;            /* Copies of each type */
    int num_molecule_types;

    /* Optional data flags */
    int has_velocities;

    /* Bond connectivity */
    int nbonds;
    int *bond_from;                  /* Bond start atom indices */
    int *bond_to;                    /* Bond end atom indices */

    /* Cached file path for error messages */
    char filename[MAX_NAME_LEN];
} h5md_gromacs_t;

/*
 * Helper Functions
 */

/* Safe string copy with null termination */
static void safe_strcpy(char *dest, const char *src, size_t dest_size) {
    if (dest_size > 0) {
        strncpy(dest, src, dest_size - 1);
        dest[dest_size - 1] = '\0';
    }
}

/* Check if HDF5 dataset exists */
static int h5_dataset_exists(hid_t loc_id, const char *name) {
    htri_t exists = H5Lexists(loc_id, name, H5P_DEFAULT);
    return (exists > 0) ? 1 : 0;
}

/* Read HDF5 integer attribute */
static int h5_read_int_attribute(hid_t obj_id, const char *attr_name, int *value) {
    hid_t attr_id;
    herr_t status;

    if (!H5Aexists(obj_id, attr_name)) {
        return H5MD_ERROR;
    }

    attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT);
    if (attr_id < 0) return H5MD_ERROR;

    status = H5Aread(attr_id, H5T_NATIVE_INT, value);
    H5Aclose(attr_id);

    return (status >= 0) ? H5MD_SUCCESS : H5MD_ERROR;
}

/* Read string dataset (handles both variable-length and fixed-length strings) */
static char **h5_read_vlen_string_dataset(hid_t dataset_id, int *num_strings) {
    hid_t space_id, file_dtype_id, mem_dtype_id;
    hsize_t dims[1];
    char **strings = NULL;
    char *buffer = NULL;
    size_t str_size;
    int i;
    htri_t is_vlen;

    /* Get dataset dimensions */
    space_id = H5Dget_space(dataset_id);
    if (H5Sget_simple_extent_ndims(space_id) != 1) {
        H5Sclose(space_id);
        return NULL;
    }

    H5Sget_simple_extent_dims(space_id, dims, NULL);
    *num_strings = (int)dims[0];

    /* Get file datatype */
    file_dtype_id = H5Dget_type(dataset_id);
    is_vlen = H5Tis_variable_str(file_dtype_id);

    /* Allocate string array */
    strings = (char **)malloc(*num_strings * sizeof(char *));
    for (i = 0; i < *num_strings; i++) {
        strings[i] = (char *)malloc(MAX_ATOM_NAME);
        strings[i][0] = '\0';
    }

    if (is_vlen > 0) {
        /* Variable-length strings */
        hvl_t *rdata = (hvl_t *)malloc(*num_strings * sizeof(hvl_t));
        mem_dtype_id = H5Tcopy(H5T_C_S1);
        H5Tset_size(mem_dtype_id, H5T_VARIABLE);

        if (H5Dread(dataset_id, mem_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata) >= 0) {
            for (i = 0; i < *num_strings; i++) {
                if (rdata[i].len > 0) {
                    safe_strcpy(strings[i], (char *)rdata[i].p, MAX_ATOM_NAME);
                }
            }
            H5Dvlen_reclaim(mem_dtype_id, space_id, H5P_DEFAULT, rdata);
        }
        free(rdata);
        H5Tclose(mem_dtype_id);
    } else {
        /* Fixed-length strings. The read buffer MUST be sized by the file's real
         * element width (str_size); only the per-element copy into strings[] is
         * capped at MAX_ATOM_NAME. Clamping the buffer stride would overflow it. */
        str_size = H5Tget_size(file_dtype_id);
        if (str_size == 0) {
            str_size = MAX_ATOM_NAME;
        }

        buffer = (char *)malloc((size_t)(*num_strings) * str_size);
        if (buffer && H5Dread(dataset_id, file_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) >= 0) {
            for (i = 0; i < *num_strings; i++) {
                /* Fixed-length HDF5 strings are not guaranteed NUL-terminated when
                 * they fill the field, so bound the copy by both str_size and the
                 * destination capacity. */
                size_t n = str_size < MAX_ATOM_NAME ? str_size : (size_t)(MAX_ATOM_NAME - 1);
                memcpy(strings[i], buffer + (size_t)i * str_size, n);
                strings[i][n] = '\0';
            }
        }
        free(buffer);
    }

    H5Tclose(file_dtype_id);
    H5Sclose(space_id);

    return strings;
}

/* Free molecule type structure */
static void free_molecule_type(molecule_type_t *moltype) {
    int i;

    if (moltype->ids) free(moltype->ids);
    if (moltype->masses) free(moltype->masses);
    if (moltype->charges) free(moltype->charges);
    if (moltype->species) free(moltype->species);
    if (moltype->residue_ids) free(moltype->residue_ids);

    if (moltype->particle_names) {
        for (i = 0; i < moltype->particle_count; i++) {
            if (moltype->particle_names[i]) free(moltype->particle_names[i]);
        }
        free(moltype->particle_names);
    }

    if (moltype->atom_residue_names) {
        for (i = 0; i < moltype->particle_count; i++) {
            if (moltype->atom_residue_names[i]) free(moltype->atom_residue_names[i]);
        }
        free(moltype->atom_residue_names);
    }
}

/* Angle (degrees) between two 3-vectors; returns 90 for degenerate input. */
static double vec_angle_deg(const double u[3], const double v[3]) {
    double lu = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    double lv = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (lu <= 0.0 || lv <= 0.0) return 90.0;
    double c = (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]) / (lu * lv);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return acos(c) * (180.0 / M_PI);
}

/*
 * Decode GROMACS Molecule Type
 *
 * Reads a molecule type from the GROMACS topology module, including:
 * - Direct arrays: id, mass, charge, species, residue_id
 * - Lookup-table-encoded: particle_name, residue_name
 */
static int decode_molecule_type(hid_t group_id, molecule_type_t *moltype) {
    hid_t dataset_id;
    int32_t *particle_name_indices = NULL;
    int32_t *residue_name_indices = NULL;
    char **particle_name_table = NULL;
    char **residue_name_table = NULL;
    int num_particle_names = 0, num_residue_names = 0;
    int i;

    /* Read attributes */
    if (h5_read_int_attribute(group_id, "particle_count", &moltype->particle_count) != H5MD_SUCCESS) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading particle_count attribute\n");
        return H5MD_ERROR;
    }

    if (h5_read_int_attribute(group_id, "residue_count", &moltype->residue_count) != H5MD_SUCCESS) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading residue_count attribute\n");
        return H5MD_ERROR;
    }

    if (moltype->particle_count <= 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Invalid particle_count (%d)\n", moltype->particle_count);
        return H5MD_ERROR;
    }

    /* Allocate arrays. The per-particle pointer arrays are calloc'd so that a
     * partially-populated molecule type can be freed safely on error. */
    moltype->ids = (int64_t *)malloc(moltype->particle_count * sizeof(int64_t));
    moltype->masses = (float *)malloc(moltype->particle_count * sizeof(float));
    moltype->charges = (float *)malloc(moltype->particle_count * sizeof(float));
    moltype->species = (int32_t *)malloc(moltype->particle_count * sizeof(int32_t));
    moltype->residue_ids = (int32_t *)malloc(moltype->particle_count * sizeof(int32_t));
    moltype->particle_names = (char **)calloc(moltype->particle_count, sizeof(char *));
    moltype->atom_residue_names = (char **)calloc(moltype->particle_count, sizeof(char *));

    if (!moltype->ids || !moltype->masses || !moltype->charges ||
        !moltype->species || !moltype->residue_ids ||
        !moltype->particle_names || !moltype->atom_residue_names) {
        fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
        return H5MD_ERROR;
    }

    /* Read direct arrays */
    dataset_id = H5Dopen(group_id, "id", H5P_DEFAULT);
    if (dataset_id < 0 || H5Dread(dataset_id, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, moltype->ids) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading id dataset\n");
        if (dataset_id >= 0) H5Dclose(dataset_id);
        return H5MD_ERROR;
    }
    H5Dclose(dataset_id);

    dataset_id = H5Dopen(group_id, "mass", H5P_DEFAULT);
    if (dataset_id < 0 || H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, moltype->masses) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading mass dataset\n");
        if (dataset_id >= 0) H5Dclose(dataset_id);
        return H5MD_ERROR;
    }
    H5Dclose(dataset_id);

    dataset_id = H5Dopen(group_id, "charge", H5P_DEFAULT);
    if (dataset_id < 0 || H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, moltype->charges) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading charge dataset\n");
        if (dataset_id >= 0) H5Dclose(dataset_id);
        return H5MD_ERROR;
    }
    H5Dclose(dataset_id);

    dataset_id = H5Dopen(group_id, "species", H5P_DEFAULT);
    if (dataset_id < 0 || H5Dread(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, moltype->species) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading species dataset\n");
        if (dataset_id >= 0) H5Dclose(dataset_id);
        return H5MD_ERROR;
    }
    H5Dclose(dataset_id);

    dataset_id = H5Dopen(group_id, "residue_id", H5P_DEFAULT);
    if (dataset_id < 0 || H5Dread(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, moltype->residue_ids) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading residue_id dataset\n");
        if (dataset_id >= 0) H5Dclose(dataset_id);
        return H5MD_ERROR;
    }
    H5Dclose(dataset_id);

    /* Read particle_name lookup table */
    dataset_id = H5Dopen(group_id, "particle_name_table", H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error opening particle_name_table\n");
        return H5MD_ERROR;
    }
    particle_name_table = h5_read_vlen_string_dataset(dataset_id, &num_particle_names);
    H5Dclose(dataset_id);

    if (!particle_name_table) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading particle_name_table\n");
        return H5MD_ERROR;
    }

    /* Read particle_name indices */
    particle_name_indices = (int32_t *)malloc(moltype->particle_count * sizeof(int32_t));
    if (!particle_name_indices) {
        fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
        goto decode_error;
    }
    dataset_id = H5Dopen(group_id, "particle_name", H5P_DEFAULT);
    if (dataset_id < 0 || H5Dread(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, particle_name_indices) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading particle_name indices\n");
        if (dataset_id >= 0) H5Dclose(dataset_id);
        goto decode_error;
    }
    H5Dclose(dataset_id);

    /* Decode particle names (one per particle) */
    for (i = 0; i < moltype->particle_count; i++) {
        int idx = particle_name_indices[i];
        moltype->particle_names[i] = (char *)malloc(MAX_ATOM_NAME);
        if (!moltype->particle_names[i]) {
            fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
            goto decode_error;
        }
        if (idx >= 0 && idx < num_particle_names) {
            safe_strcpy(moltype->particle_names[i], particle_name_table[idx], MAX_ATOM_NAME);
        } else {
            safe_strcpy(moltype->particle_names[i], "X", MAX_ATOM_NAME);
        }
    }

    /* Read residue_name lookup table */
    dataset_id = H5Dopen(group_id, "residue_name_table", H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error opening residue_name_table\n");
        goto decode_error;
    }
    residue_name_table = h5_read_vlen_string_dataset(dataset_id, &num_residue_names);
    H5Dclose(dataset_id);

    if (!residue_name_table) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading residue_name_table\n");
        goto decode_error;
    }

    /* Read residue_name indices (one per particle) */
    residue_name_indices = (int32_t *)malloc(moltype->particle_count * sizeof(int32_t));
    if (!residue_name_indices) {
        fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
        goto decode_error;
    }
    dataset_id = H5Dopen(group_id, "residue_name", H5P_DEFAULT);
    if (dataset_id < 0 || H5Dread(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, residue_name_indices) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading residue_name indices\n");
        if (dataset_id >= 0) H5Dclose(dataset_id);
        goto decode_error;
    }
    H5Dclose(dataset_id);

    /* Resolve the residue name for each particle directly from the per-particle
     * residue_name index. NOTE: residue_id is intentionally NOT used as a table
     * index. GROMACS writes residue_id = resind + 1 (values run 2..residue_count+1,
     * not 1..residue_count), so indexing a per-type residue table by residue_id
     * mislabels every single-residue molecule (water, ions) and the last residue
     * of each chain as "UNK". The per-particle residue_name index is authoritative. */
    for (i = 0; i < moltype->particle_count; i++) {
        int idx = residue_name_indices[i];
        moltype->atom_residue_names[i] = (char *)malloc(MAX_RESIDUE_NAME);
        if (!moltype->atom_residue_names[i]) {
            fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
            goto decode_error;
        }
        if (idx >= 0 && idx < num_residue_names) {
            safe_strcpy(moltype->atom_residue_names[i], residue_name_table[idx], MAX_RESIDUE_NAME);
        } else {
            safe_strcpy(moltype->atom_residue_names[i], "UNK", MAX_RESIDUE_NAME);
        }
    }

    /* Cleanup (success). moltype's own arrays are owned by the caller. */
    free(particle_name_indices);
    free(residue_name_indices);
    for (i = 0; i < num_particle_names; i++) free(particle_name_table[i]);
    free(particle_name_table);
    for (i = 0; i < num_residue_names; i++) free(residue_name_table[i]);
    free(residue_name_table);

    return H5MD_SUCCESS;

decode_error:
    /* Free only the local scratch buffers; moltype's arrays are released by the
     * caller via free_molecule_type(). */
    if (particle_name_indices) free(particle_name_indices);
    if (residue_name_indices) free(residue_name_indices);
    if (particle_name_table) {
        for (i = 0; i < num_particle_names; i++) free(particle_name_table[i]);
        free(particle_name_table);
    }
    if (residue_name_table) {
        for (i = 0; i < num_residue_names; i++) free(residue_name_table[i]);
        free(residue_name_table);
    }
    return H5MD_ERROR;
}

/*
 * VMD Plugin Interface Functions
 */

static void *open_h5md_gromacs_read(const char *filepath, const char *filetype, int *natoms) {
    h5md_gromacs_t *h5md = NULL;
    hid_t topology_group = -1;
    hid_t moltype_group = -1;
    hid_t position_group = -1;
    hid_t attr_id = -1;
    hid_t aspace_id = -1;
    hid_t dtype_id = -1;
    hsize_t attr_dims[1];
    char **moltype_names = NULL;
    hvl_t *rdata = NULL;
    int i, total_atoms;
    (void)filetype;  /* Plugin API requires this parameter */

    /* Allocate plugin handle */
    h5md = (h5md_gromacs_t *)calloc(1, sizeof(h5md_gromacs_t));
    if (!h5md) {
        fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
        return NULL;
    }

    safe_strcpy(h5md->filename, filepath, MAX_NAME_LEN);

    /* calloc leaves hid_t fields at 0, which is not a valid handle. Initialize
     * the optional trajectory handles to -1 so "< 0" / ">= 0" guards behave. */
    h5md->position_dataset = -1;
    h5md->velocity_dataset = -1;
    h5md->box_dataset = -1;
    h5md->time_dataset = -1;

    /* Disable HDF5 error messages to stderr */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    /* Open HDF5 file */
    h5md->file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h5md->file_id < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Unable to open file: %s\n", filepath);
        free(h5md);
        return NULL;
    }

    /* Verify H5MD structure - check for /h5md/modules/gromacs_topology */
    if (!h5_dataset_exists(h5md->file_id, "/h5md/modules/gromacs_topology")) {
        fprintf(stderr, "h5mdplugin_gromacs) File does not contain GROMACS topology module\n");
        H5Fclose(h5md->file_id);
        free(h5md);
        return NULL;
    }

    /* Open topology group */
    topology_group = H5Gopen(h5md->file_id, "/h5md/modules/gromacs_topology", H5P_DEFAULT);
    if (topology_group < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error opening topology group\n");
        H5Fclose(h5md->file_id);
        free(h5md);
        return NULL;
    }

    /* Read molecule_block_names attribute */
    if (!H5Aexists(topology_group, "molecule_block_names")) {
        fprintf(stderr, "h5mdplugin_gromacs) Missing molecule_block_names attribute\n");
        H5Gclose(topology_group);
        H5Fclose(h5md->file_id);
        free(h5md);
        return NULL;
    }

    attr_id = H5Aopen(topology_group, "molecule_block_names", H5P_DEFAULT);
    aspace_id = H5Aget_space(attr_id);
    attr_dims[0] = 0;
    if (attr_id < 0 || aspace_id < 0 ||
        H5Sget_simple_extent_dims(aspace_id, attr_dims, NULL) < 0 ||
        (int)attr_dims[0] <= 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Invalid molecule_block_names attribute\n");
        if (aspace_id >= 0) H5Sclose(aspace_id);
        if (attr_id >= 0) H5Aclose(attr_id);
        H5Gclose(topology_group);
        H5Fclose(h5md->file_id);
        free(h5md);
        return NULL;
    }
    h5md->num_molecule_types = (int)attr_dims[0];

    /* Get attribute datatype to check if variable-length or fixed-length */
    hid_t file_dtype_id = H5Aget_type(attr_id);
    htri_t is_vlen = H5Tis_variable_str(file_dtype_id);

    /* Allocate molecule type names array (calloc so the error path can free
     * partially-populated slots safely). */
    moltype_names = (char **)calloc(h5md->num_molecule_types, sizeof(char *));
    if (!moltype_names) {
        fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
        H5Tclose(file_dtype_id);
        H5Sclose(aspace_id);
        H5Aclose(attr_id);
        H5Gclose(topology_group);
        H5Fclose(h5md->file_id);
        free(h5md);
        return NULL;
    }
    for (i = 0; i < h5md->num_molecule_types; i++) {
        moltype_names[i] = (char *)malloc(MAX_MOLTYPE_NAME);
        if (!moltype_names[i]) {
            fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
            H5Tclose(file_dtype_id);
            H5Sclose(aspace_id);
            H5Aclose(attr_id);
            goto error;
        }
        moltype_names[i][0] = '\0';
    }

    if (is_vlen > 0) {
        /* Variable-length strings */
        dtype_id = H5Tcopy(H5T_C_S1);
        H5Tset_size(dtype_id, H5T_VARIABLE);

        rdata = (hvl_t *)malloc(h5md->num_molecule_types * sizeof(hvl_t));
        if (H5Aread(attr_id, dtype_id, rdata) >= 0) {
            for (i = 0; i < h5md->num_molecule_types; i++) {
                if (rdata[i].len > 0) {
                    safe_strcpy(moltype_names[i], (char *)rdata[i].p, MAX_MOLTYPE_NAME);
                }
            }
            H5Dvlen_reclaim(dtype_id, aspace_id, H5P_DEFAULT, rdata);
        }
        free(rdata);
        H5Tclose(dtype_id);
    } else {
        /* Fixed-length strings: size the buffer by the real element width, cap
         * only the per-element copy (see h5_read_vlen_string_dataset). */
        size_t str_size = H5Tget_size(file_dtype_id);
        if (str_size == 0) {
            str_size = MAX_MOLTYPE_NAME;
        }

        char *buffer = (char *)malloc((size_t)h5md->num_molecule_types * str_size);
        if (buffer && H5Aread(attr_id, file_dtype_id, buffer) >= 0) {
            for (i = 0; i < h5md->num_molecule_types; i++) {
                size_t n = str_size < MAX_MOLTYPE_NAME ? str_size : (size_t)(MAX_MOLTYPE_NAME - 1);
                memcpy(moltype_names[i], buffer + (size_t)i * str_size, n);
                moltype_names[i][n] = '\0';
            }
        }
        free(buffer);
    }

    H5Tclose(file_dtype_id);
    H5Sclose(aspace_id);
    H5Aclose(attr_id);

    /* Read molecule_block_counts attribute */
    h5md->molecule_counts = (int *)malloc(h5md->num_molecule_types * sizeof(int));
    if (!h5md->molecule_counts) {
        fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
        goto error;
    }
    attr_id = H5Aopen(topology_group, "molecule_block_counts", H5P_DEFAULT);
    if (attr_id < 0 || H5Aread(attr_id, H5T_NATIVE_INT, h5md->molecule_counts) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading molecule_block_counts attribute\n");
        if (attr_id >= 0) H5Aclose(attr_id);
        goto error;
    }
    H5Aclose(attr_id);

    /* Allocate molecule types array */
    h5md->molecule_types = (molecule_type_t *)calloc(h5md->num_molecule_types, sizeof(molecule_type_t));
    if (!h5md->molecule_types) {
        fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
        goto error;
    }

    /* Decode each molecule type */
    total_atoms = 0;
    for (i = 0; i < h5md->num_molecule_types; i++) {
        safe_strcpy(h5md->molecule_types[i].name, moltype_names[i], MAX_MOLTYPE_NAME);

        moltype_group = H5Gopen(topology_group, moltype_names[i], H5P_DEFAULT);
        if (moltype_group < 0) {
            fprintf(stderr, "h5mdplugin_gromacs) Error opening molecule type: %s\n", moltype_names[i]);
            goto error;
        }

        if (decode_molecule_type(moltype_group, &h5md->molecule_types[i]) != H5MD_SUCCESS) {
            fprintf(stderr, "h5mdplugin_gromacs) Error decoding molecule type: %s\n", moltype_names[i]);
            H5Gclose(moltype_group);
            goto error;
        }

        H5Gclose(moltype_group);

        if (h5md->molecule_counts[i] < 0) {
            fprintf(stderr, "h5mdplugin_gromacs) Invalid molecule count %d for %s\n",
                    h5md->molecule_counts[i], moltype_names[i]);
            goto error;
        }

        total_atoms += h5md->molecule_types[i].particle_count * h5md->molecule_counts[i];
        free(moltype_names[i]);
        moltype_names[i] = NULL;  /* prevent double-free on the error path */
    }
    free(moltype_names);
    moltype_names = NULL;

    H5Gclose(topology_group);
    topology_group = -1;

    h5md->natoms_file = total_atoms;
    h5md->natoms = total_atoms;

    /* Open position dataset to get frame count */
    position_group = H5Gopen(h5md->file_id, "/particles/system/position", H5P_DEFAULT);
    if (position_group < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error opening position group\n");
        goto error;
    }

    h5md->position_dataset = H5Dopen(position_group, "value", H5P_DEFAULT);
    if (h5md->position_dataset < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error opening position/value dataset\n");
        H5Gclose(position_group);
        goto error;
    }

    /* Get dataset dimensions [nframes, natoms, 3] */
    hid_t space_id = H5Dget_space(h5md->position_dataset);
    hsize_t dims[3] = {0, 0, 0};
    if (H5Sget_simple_extent_ndims(space_id) != 3 ||
        H5Sget_simple_extent_dims(space_id, dims, NULL) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error: position/value must be a 3D [frames,atoms,3] dataset\n");
        H5Sclose(space_id);
        H5Gclose(position_group);
        goto error;
    }
    H5Sclose(space_id);

    h5md->nframes = (int)dims[0];

    /* The plugin reconstructs the full system from the topology and reads each
     * frame with a hyperslab sized natoms_file, so the counts must agree. */
    if ((int)dims[1] != h5md->natoms_file) {
        fprintf(stderr, "h5mdplugin_gromacs) Error: position atom count (%d) != topology atom count (%d)\n",
                (int)dims[1], h5md->natoms_file);
        H5Gclose(position_group);
        goto error;
    }

    /* Open time dataset */
    h5md->time_dataset = H5Dopen(position_group, "time", H5P_DEFAULT);
    H5Gclose(position_group);

    /* Check for velocities */
    h5md->has_velocities = h5_dataset_exists(h5md->file_id, "/particles/system/velocity/value");
    if (h5md->has_velocities) {
        hid_t vel_group = H5Gopen(h5md->file_id, "/particles/system/velocity", H5P_DEFAULT);
        if (vel_group >= 0) {
            h5md->velocity_dataset = H5Dopen(vel_group, "value", H5P_DEFAULT);
            H5Gclose(vel_group);
        }
    }

    /* Check for box dimensions */
    if (h5_dataset_exists(h5md->file_id, "/particles/system/box/edges/value")) {
        hid_t box_group = H5Gopen(h5md->file_id, "/particles/system/box/edges", H5P_DEFAULT);
        if (box_group >= 0) {
            h5md->box_dataset = H5Dopen(box_group, "value", H5P_DEFAULT);
            H5Gclose(box_group);
        }
    }

    *natoms = h5md->natoms;
    h5md->current_frame = 0;

    return h5md;

error:
    if (moltype_names) {
        for (i = 0; i < h5md->num_molecule_types; i++) {
            if (moltype_names[i]) free(moltype_names[i]);
        }
        free(moltype_names);
    }

    if (h5md->molecule_types) {
        for (i = 0; i < h5md->num_molecule_types; i++) {
            free_molecule_type(&h5md->molecule_types[i]);
        }
        free(h5md->molecule_types);
    }

    if (h5md->molecule_counts) free(h5md->molecule_counts);

    if (topology_group >= 0) H5Gclose(topology_group);
    if (h5md->position_dataset >= 0) H5Dclose(h5md->position_dataset);
    if (h5md->velocity_dataset >= 0) H5Dclose(h5md->velocity_dataset);
    if (h5md->box_dataset >= 0) H5Dclose(h5md->box_dataset);
    if (h5md->time_dataset >= 0) H5Dclose(h5md->time_dataset);
    if (h5md->file_id >= 0) H5Fclose(h5md->file_id);

    free(h5md);
    return NULL;
}

static int read_h5md_gromacs_structure(void *mydata, int *optflags, molfile_atom_t *atoms) {
    h5md_gromacs_t *h5md = (h5md_gromacs_t *)mydata;
    int vmd_atom_idx = 0;
    int moltype_idx, mol_copy, atom_idx;

    /* Set optional flags - we provide mass, charge, and atomic number */
    *optflags = MOLFILE_MASS | MOLFILE_CHARGE | MOLFILE_ATOMICNUMBER;

    /* Reconstruct full system from molecule types */
    for (moltype_idx = 0; moltype_idx < h5md->num_molecule_types; moltype_idx++) {
        molecule_type_t *moltype = &h5md->molecule_types[moltype_idx];
        int num_copies = h5md->molecule_counts[moltype_idx];

        for (mol_copy = 0; mol_copy < num_copies; mol_copy++) {
            for (atom_idx = 0; atom_idx < moltype->particle_count; atom_idx++) {
                molfile_atom_t *atom = &atoms[vmd_atom_idx];

                /* REQUIRED FIELDS */
                safe_strcpy(atom->name, moltype->particle_names[atom_idx], sizeof(atom->name));
                safe_strcpy(atom->type, moltype->particle_names[atom_idx], sizeof(atom->type));

                /* Residue name is resolved per particle from the lookup table
                 * (see decode_molecule_type); residue_id must not be used as a
                 * table index because GROMACS offsets it by one. */
                safe_strcpy(atom->resname, moltype->atom_residue_names[atom_idx], sizeof(atom->resname));

                /* residue_id is already a unique per-residue key; offset by the
                 * molecule copy so distinct copies stay distinct in VMD. */
                atom->resid = moltype->residue_ids[atom_idx] + (mol_copy * 10000);

                safe_strcpy(atom->segid, moltype->name, sizeof(atom->segid));
                atom->chain[0] = '\0';
                atom->chain[1] = '\0';

                /* OPTIONAL FIELDS */
                atom->mass = moltype->masses[atom_idx];
                atom->charge = moltype->charges[atom_idx];
                atom->atomicnumber = moltype->species[atom_idx];

                vmd_atom_idx++;
            }
        }
    }

    h5md->natoms = vmd_atom_idx;

    return MOLFILE_SUCCESS;
}

static int read_h5md_gromacs_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
    h5md_gromacs_t *h5md = (h5md_gromacs_t *)mydata;
    hid_t mem_space = -1, file_space = -1;
    hsize_t mem_dims[2], file_start[3], file_count[3];
    float *temp_coords = NULL;
    int i;
    (void)natoms;  /* We use h5md->natoms_file, sized from the topology */

    /* Check EOF */
    if (h5md->current_frame >= h5md->nframes) {
        return MOLFILE_EOF;
    }

    /* Skip if ts is NULL (for seeking) */
    if (ts == NULL) {
        h5md->current_frame++;
        return MOLFILE_SUCCESS;
    }

    /* Define hyperslab for this frame */
    file_start[0] = h5md->current_frame;
    file_start[1] = 0;
    file_start[2] = 0;

    file_count[0] = 1;
    file_count[1] = h5md->natoms_file;
    file_count[2] = 3;

    file_space = H5Dget_space(h5md->position_dataset);
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

    /* Memory space for coordinate data */
    mem_dims[0] = h5md->natoms_file;
    mem_dims[1] = 3;
    mem_space = H5Screate_simple(2, mem_dims, NULL);

    /* Allocate temporary buffer */
    temp_coords = (float *)malloc(h5md->natoms_file * 3 * sizeof(float));
    if (!temp_coords) {
        fprintf(stderr, "h5mdplugin_gromacs) Memory allocation failed\n");
        H5Sclose(mem_space);
        H5Sclose(file_space);
        return MOLFILE_ERROR;
    }

    /* Read coordinates */
    if (H5Dread(h5md->position_dataset, H5T_NATIVE_FLOAT, mem_space, file_space, H5P_DEFAULT, temp_coords) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error reading coordinates for frame %d\n", h5md->current_frame);
        free(temp_coords);
        H5Sclose(mem_space);
        H5Sclose(file_space);
        return MOLFILE_ERROR;
    }

    /* Convert units: nm to Angstroms and copy to VMD format */
    for (i = 0; i < h5md->natoms_file; i++) {
        ts->coords[3 * i + 0] = temp_coords[3 * i + 0] * NM_TO_ANGSTROM;
        ts->coords[3 * i + 1] = temp_coords[3 * i + 1] * NM_TO_ANGSTROM;
        ts->coords[3 * i + 2] = temp_coords[3 * i + 2] * NM_TO_ANGSTROM;
    }

    free(temp_coords);
    H5Sclose(mem_space);
    H5Sclose(file_space);

    /* Read velocities if available and requested */
    if (h5md->has_velocities && h5md->velocity_dataset >= 0 && ts->velocities) {
        float *temp_vels = (float *)malloc(h5md->natoms_file * 3 * sizeof(float));
        if (temp_vels) {
            file_space = H5Dget_space(h5md->velocity_dataset);
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

            mem_space = H5Screate_simple(2, mem_dims, NULL);

            if (H5Dread(h5md->velocity_dataset, H5T_NATIVE_FLOAT, mem_space, file_space, H5P_DEFAULT, temp_vels) >= 0) {
                /* Convert (nm/ps -> Angstrom/ps) and copy velocities */
                for (i = 0; i < h5md->natoms_file; i++) {
                    ts->velocities[3 * i + 0] = temp_vels[3 * i + 0] * NM_TO_ANGSTROM;
                    ts->velocities[3 * i + 1] = temp_vels[3 * i + 1] * NM_TO_ANGSTROM;
                    ts->velocities[3 * i + 2] = temp_vels[3 * i + 2] * NM_TO_ANGSTROM;
                }
            }

            free(temp_vels);
            H5Sclose(mem_space);
            H5Sclose(file_space);
        }
    }

    /* Read box dimensions if available */
    if (h5md->box_dataset >= 0) {
        float box[9];  /* 3x3 matrix */
        hsize_t box_start[3] = {h5md->current_frame, 0, 0};
        hsize_t box_count[3] = {1, 3, 3};

        file_space = H5Dget_space(h5md->box_dataset);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, box_start, NULL, box_count, NULL);

        hsize_t box_dims[2] = {3, 3};
        mem_space = H5Screate_simple(2, box_dims, NULL);

        if (H5Dread(h5md->box_dataset, H5T_NATIVE_FLOAT, mem_space, file_space, H5P_DEFAULT, box) >= 0) {
            /* box is a row-major 3x3 matrix whose rows are the cell vectors
             * a, b, c (GROMACS triclinic convention). Derive VMD's A/B/C lengths
             * and alpha/beta/gamma angles so non-orthogonal boxes are correct. */
            double a[3] = { box[0], box[1], box[2] };
            double b[3] = { box[3], box[4], box[5] };
            double c[3] = { box[6], box[7], box[8] };

            ts->A = (float)(sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) * NM_TO_ANGSTROM);
            ts->B = (float)(sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]) * NM_TO_ANGSTROM);
            ts->C = (float)(sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) * NM_TO_ANGSTROM);
            ts->alpha = (float)vec_angle_deg(b, c);  /* angle between b and c */
            ts->beta  = (float)vec_angle_deg(a, c);  /* angle between a and c */
            ts->gamma = (float)vec_angle_deg(a, b);  /* angle between a and b */
        }

        H5Sclose(mem_space);
        H5Sclose(file_space);
    }

    /* Read physical time */
    if (h5md->time_dataset >= 0) {
        double time_val;
        hsize_t time_start[1] = {h5md->current_frame};
        hsize_t time_count[1] = {1};

        file_space = H5Dget_space(h5md->time_dataset);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, time_start, NULL, time_count, NULL);

        hsize_t time_dim[1] = {1};
        mem_space = H5Screate_simple(1, time_dim, NULL);

        if (H5Dread(h5md->time_dataset, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &time_val) >= 0) {
            ts->physical_time = time_val;
        }

        H5Sclose(mem_space);
        H5Sclose(file_space);
    }

    h5md->current_frame++;
    return MOLFILE_SUCCESS;
}

/*
 * Read bond connectivity from /connectivity/bonds dataset
 *
 * GROMACS H5MD files store bonds as [nbonds, 2] array of atom indices.
 * This function reads and translates bonds accounting for atom selection.
 */
static int read_h5md_gromacs_bonds(void *mydata, int *nbonds, int **fromptr, int **toptr,
                                    float **bondorder, int **bondtype,
                                    int *nbondtypes, char ***bondtypename) {
    h5md_gromacs_t *h5md = (h5md_gromacs_t *)mydata;
    hid_t bonds_dataset = -1;
    hid_t space_id = -1;
    hsize_t dims[2];
    int64_t *bonds_buffer = NULL;
    int i, valid_bonds;

    /* Initialize outputs */
    *nbonds = 0;
    *fromptr = NULL;
    *toptr = NULL;
    *bondorder = NULL;
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;

    /* Check if /connectivity/bonds dataset exists */
    if (!h5_dataset_exists(h5md->file_id, "/connectivity/bonds")) {
        /* No bonds in file - this is not an error */
        return MOLFILE_SUCCESS;
    }

    /* Open bonds dataset */
    bonds_dataset = H5Dopen(h5md->file_id, "/connectivity/bonds", H5P_DEFAULT);
    if (bonds_dataset < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Warning: Cannot open /connectivity/bonds\n");
        return MOLFILE_SUCCESS;  /* Not a fatal error */
    }

    /* Get dataset dimensions */
    space_id = H5Dget_space(bonds_dataset);
    if (H5Sget_simple_extent_ndims(space_id) != 2) {
        fprintf(stderr, "h5mdplugin_gromacs) Error: /connectivity/bonds must be 2D\n");
        H5Sclose(space_id);
        H5Dclose(bonds_dataset);
        return MOLFILE_ERROR;
    }

    H5Sget_simple_extent_dims(space_id, dims, NULL);
    if (dims[1] != 2) {
        fprintf(stderr, "h5mdplugin_gromacs) Error: /connectivity/bonds must have shape [N, 2]\n");
        H5Sclose(space_id);
        H5Dclose(bonds_dataset);
        return MOLFILE_ERROR;
    }

    h5md->nbonds = (int)dims[0];

    /* Allocate temporary buffer for reading file bonds */
    bonds_buffer = (int64_t *)malloc(h5md->nbonds * 2 * sizeof(int64_t));
    if (!bonds_buffer) {
        fprintf(stderr, "h5mdplugin_gromacs) Error: Cannot allocate bond buffer\n");
        H5Sclose(space_id);
        H5Dclose(bonds_dataset);
        return MOLFILE_ERROR;
    }

    /* Read bonds from file */
    if (H5Dread(bonds_dataset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, bonds_buffer) < 0) {
        fprintf(stderr, "h5mdplugin_gromacs) Error: Cannot read bond data\n");
        free(bonds_buffer);
        H5Sclose(space_id);
        H5Dclose(bonds_dataset);
        return MOLFILE_ERROR;
    }

    H5Sclose(space_id);
    H5Dclose(bonds_dataset);

    /* Allocate bond arrays (initially for all bonds) */
    h5md->bond_from = (int *)malloc(h5md->nbonds * sizeof(int));
    h5md->bond_to = (int *)malloc(h5md->nbonds * sizeof(int));

    if (!h5md->bond_from || !h5md->bond_to) {
        fprintf(stderr, "h5mdplugin_gromacs) Error: Cannot allocate bond arrays\n");
        free(bonds_buffer);
        return MOLFILE_ERROR;
    }

    /*
     * GROMACS stores bonds as 0-based global atom indices; VMD expects 1-based,
     * so add 1. Skip out-of-range and self-bonds defensively.
     */
    valid_bonds = 0;
    for (i = 0; i < h5md->nbonds; i++) {
        int file_from = (int)bonds_buffer[2*i + 0];
        int file_to = (int)bonds_buffer[2*i + 1];
        int vmd_from, vmd_to;

        /* Check bounds */
        if (file_from < 0 || file_from >= h5md->natoms_file ||
            file_to < 0 || file_to >= h5md->natoms_file) {
            fprintf(stderr, "h5mdplugin_gromacs) Warning: Bond %d has out-of-bounds indices (%d, %d)\n",
                    i, file_from, file_to);
            continue;
        }

        /* Skip self-bonds */
        if (file_from == file_to) {
            continue;
        }

        vmd_from = file_from;
        vmd_to = file_to;

        /* Store valid bond with 1-based indices (VMD requirement) */
        h5md->bond_from[valid_bonds] = vmd_from + 1;
        h5md->bond_to[valid_bonds] = vmd_to + 1;
        valid_bonds++;
    }

    free(bonds_buffer);

    /* Update bond count to valid bonds only */
    h5md->nbonds = valid_bonds;

    /* Set output pointers */
    *nbonds = h5md->nbonds;
    *fromptr = h5md->bond_from;
    *toptr = h5md->bond_to;

    /* GROMACS H5MD doesn't provide bond order or type information */
    *bondorder = NULL;
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;

    return MOLFILE_SUCCESS;
}

static void close_h5md_gromacs_read(void *mydata) {
    h5md_gromacs_t *h5md = (h5md_gromacs_t *)mydata;
    int i;

    if (!h5md) return;

    /* Close HDF5 datasets */
    if (h5md->position_dataset >= 0) H5Dclose(h5md->position_dataset);
    if (h5md->velocity_dataset >= 0) H5Dclose(h5md->velocity_dataset);
    if (h5md->box_dataset >= 0) H5Dclose(h5md->box_dataset);
    if (h5md->time_dataset >= 0) H5Dclose(h5md->time_dataset);

    /* Close HDF5 file */
    if (h5md->file_id >= 0) H5Fclose(h5md->file_id);

    /* Free topology data */
    if (h5md->molecule_types) {
        for (i = 0; i < h5md->num_molecule_types; i++) {
            free_molecule_type(&h5md->molecule_types[i]);
        }
        free(h5md->molecule_types);
    }

    if (h5md->molecule_counts) free(h5md->molecule_counts);

    /* Free bond data */
    if (h5md->bond_from) free(h5md->bond_from);
    if (h5md->bond_to) free(h5md->bond_to);

    free(h5md);
}

static int read_h5md_gromacs_metadata(void *mydata, molfile_timestep_metadata_t *meta) {
    h5md_gromacs_t *h5md = (h5md_gromacs_t *)mydata;

    meta->count = h5md->nframes;
    meta->has_velocities = h5md->has_velocities;

    return MOLFILE_SUCCESS;
}

/*
 * Plugin Registration
 */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
    memset(&plugin, 0, sizeof(molfile_plugin_t));

    plugin.abiversion = vmdplugin_ABIVERSION;
    plugin.type = MOLFILE_PLUGIN_TYPE;
    plugin.name = "h5md_gromacs";
    plugin.prettyname = "GROMACS H5MD";
    plugin.author = "Diego E. B. Gomes";
    plugin.majorv = H5MD_GROMACS_VERSION_MAJOR;
    plugin.minorv = H5MD_GROMACS_VERSION_MINOR;
    plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;  /* HDF5 limitation */
    plugin.filename_extension = "h5md";

    /* Set function pointers */
    plugin.open_file_read = open_h5md_gromacs_read;
    plugin.read_structure = read_h5md_gromacs_structure;
    plugin.read_bonds = read_h5md_gromacs_bonds;
    plugin.read_next_timestep = read_h5md_gromacs_timestep;
    plugin.close_file_read = close_h5md_gromacs_read;
    plugin.read_timestep_metadata = read_h5md_gromacs_metadata;

    return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
    (*cb)(v, (vmdplugin_t *)&plugin);
    return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
    return VMDPLUGIN_SUCCESS;
}
