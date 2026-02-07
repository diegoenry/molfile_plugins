/*
 * SDF/MOL V2000 plugin for VMD
 *
 * Native reader/writer for MDL SDF (Structure Data File) and MOL formats.
 * Supports V2000 connection tables with 3D coordinates, bonds with bond
 * orders, formal charges (both atom block charge codes and M CHG lines),
 * and multi-molecule SDF files as trajectory frames.
 *
 * Author: Diego Enry Barreto Gomes
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "molfile_plugin.h"

#define SDF_LINESIZE 256

/* ========================================================================
 * Element lookup table
 * ======================================================================== */

typedef struct {
  const char *symbol;
  int atomicnumber;
  float mass;
} sdf_element_t;

static const sdf_element_t element_table[] = {
  {"H",   1,   1.008f},
  {"He",  2,   4.003f},
  {"Li",  3,   6.941f},
  {"Be",  4,   9.012f},
  {"B",   5,  10.811f},
  {"C",   6,  12.011f},
  {"N",   7,  14.007f},
  {"O",   8,  15.999f},
  {"F",   9,  18.998f},
  {"Ne", 10,  20.180f},
  {"Na", 11,  22.990f},
  {"Mg", 12,  24.305f},
  {"Al", 13,  26.982f},
  {"Si", 14,  28.086f},
  {"P",  15,  30.974f},
  {"S",  16,  32.065f},
  {"Cl", 17,  35.453f},
  {"Ar", 18,  39.948f},
  {"K",  19,  39.098f},
  {"Ca", 20,  40.078f},
  {"Ti", 22,  47.867f},
  {"V",  23,  50.942f},
  {"Cr", 24,  51.996f},
  {"Mn", 25,  54.938f},
  {"Fe", 26,  55.845f},
  {"Co", 27,  58.933f},
  {"Ni", 28,  58.693f},
  {"Cu", 29,  63.546f},
  {"Zn", 30,  65.380f},
  {"Ga", 31,  69.723f},
  {"Ge", 32,  72.640f},
  {"As", 33,  74.922f},
  {"Se", 34,  78.960f},
  {"Br", 35,  79.904f},
  {"Kr", 36,  83.798f},
  {"Rb", 37,  85.468f},
  {"Sr", 38,  87.620f},
  {"Y",  39,  88.906f},
  {"Zr", 40,  91.224f},
  {"Mo", 42,  95.960f},
  {"Ru", 44, 101.070f},
  {"Rh", 45, 102.906f},
  {"Pd", 46, 106.420f},
  {"Ag", 47, 107.868f},
  {"Cd", 48, 112.411f},
  {"In", 49, 114.818f},
  {"Sn", 50, 118.710f},
  {"Sb", 51, 121.760f},
  {"Te", 52, 127.600f},
  {"I",  53, 126.904f},
  {"Xe", 54, 131.293f},
  {"Cs", 55, 132.905f},
  {"Ba", 56, 137.327f},
  {"Pt", 78, 195.078f},
  {"Au", 79, 196.967f},
  {"Hg", 80, 200.590f},
  {"Tl", 81, 204.383f},
  {"Pb", 82, 207.200f},
  {"Bi", 83, 208.980f},
  {"U",  92, 238.029f},
  {NULL,  0,   0.000f}
};

#define NUM_ELEMENTS (sizeof(element_table) / sizeof(element_table[0]) - 1)

/* Look up element by symbol; returns pointer or NULL */
static const sdf_element_t *lookup_element_by_symbol(const char *sym) {
  int i;
  for (i = 0; element_table[i].symbol != NULL; i++) {
    if (strcasecmp(sym, element_table[i].symbol) == 0)
      return &element_table[i];
  }
  return NULL;
}

/* Look up element by atomic number; returns pointer or NULL */
static const sdf_element_t *lookup_element_by_number(int z) {
  int i;
  for (i = 0; element_table[i].symbol != NULL; i++) {
    if (element_table[i].atomicnumber == z)
      return &element_table[i];
  }
  return NULL;
}


/* ========================================================================
 * V2000 charge code conversion
 *   Code 0 = 0 (uncharged)
 *   Code 1 = +3, Code 2 = +2, Code 3 = +1
 *   Code 4 = doublet radical (treat as 0)
 *   Code 5 = -1, Code 6 = -2, Code 7 = -3
 * ======================================================================== */

static float charge_code_to_formal(int code) {
  switch (code) {
    case 1: return  3.0f;
    case 2: return  2.0f;
    case 3: return  1.0f;
    case 4: return  0.0f;  /* doublet radical */
    case 5: return -1.0f;
    case 6: return -2.0f;
    case 7: return -3.0f;
    default: return 0.0f;
  }
}

static int formal_to_charge_code(float charge) {
  int q = (int)charge;
  switch (q) {
    case  3: return 1;
    case  2: return 2;
    case  1: return 3;
    case -1: return 5;
    case -2: return 6;
    case -3: return 7;
    default: return 0;
  }
}


/* ========================================================================
 * Plugin handle
 * ======================================================================== */

typedef struct {
  FILE *file;
  char *filepath;
  int natoms;
  int nbonds;
  int coords_read;              /* have we read the first frame yet? */
  float *first_coords;          /* cached coords from read_structure parse */

  /* Bond data (reader) */
  int *from, *to;               /* 1-based bond indices */
  float *bondorder;

  /* Atom data for element lookup */
  int *atomicnumbers;
  float *masses;
  float *charges;

  /* Writer cached data */
  molfile_atom_t *atomlist;     /* cached from write_structure */
  int optflags;
  int *write_from, *write_to;
  float *write_bondorder;
  int write_nbonds;
} sdf_handle;


/* ========================================================================
 * Helper: read a complete MOL block (header through M  END) and parse it.
 * Stores atom names, bonds, charges, and coords into the handle.
 * Returns 0 on success, -1 on error/EOF.
 *
 * If 'atoms' is non-NULL, fills in the molfile_atom_t array (structure read).
 * If 'ts' is non-NULL, fills in coordinates.
 * If 'atoms' is NULL and 'ts' is NULL, we just parse to count atoms/bonds
 *   (used by open_sdf_read).
 * ======================================================================== */

static int parse_mol_block(sdf_handle *sdf, molfile_atom_t *atoms,
                           molfile_timestep_t *ts, int *out_natoms,
                           int *out_nbonds, int first_pass) {
  char line[SDF_LINESIZE];
  int natoms, nbonds, i;
  char version[8];

  /* Line 1: molecule name */
  if (fgets(line, SDF_LINESIZE, sdf->file) == NULL)
    return -1;

  /* Line 2: program/timestamp */
  if (fgets(line, SDF_LINESIZE, sdf->file) == NULL)
    return -1;

  /* Line 3: comment */
  if (fgets(line, SDF_LINESIZE, sdf->file) == NULL)
    return -1;

  /* Line 4: counts line
   * Format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
   * aaa = natoms (3 chars), bbb = nbonds (3 chars), rest follows */
  if (fgets(line, SDF_LINESIZE, sdf->file) == NULL)
    return -1;

  /* Ensure line is long enough */
  if (strlen(line) < 6) {
    fprintf(stderr, "sdfplugin) Counts line too short\n");
    return -1;
  }

  version[0] = '\0';
  /* Parse fixed-width fields */
  {
    char buf_a[4], buf_b[4];
    strncpy(buf_a, line, 3); buf_a[3] = '\0';
    strncpy(buf_b, line + 3, 3); buf_b[3] = '\0';
    natoms = atoi(buf_a);
    nbonds = atoi(buf_b);

    /* Check for V2000/V3000 tag */
    if (strlen(line) >= 39) {
      strncpy(version, line + 33, 6);
      version[6] = '\0';
    }
  }

  if (natoms <= 0) {
    fprintf(stderr, "sdfplugin) Invalid atom count: %d\n", natoms);
    return -1;
  }

  /* Check for V3000 - not supported */
  if (strstr(version, "V3000") != NULL) {
    fprintf(stderr, "sdfplugin) V3000 format not supported\n");
    return -1;
  }

  if (out_natoms) *out_natoms = natoms;
  if (out_nbonds) *out_nbonds = nbonds;

  /* ---- Atom block ---- */
  for (i = 0; i < natoms; i++) {
    float x, y, z;
    char symbol[4];
    int chg_code = 0;

    if (fgets(line, SDF_LINESIZE, sdf->file) == NULL) {
      fprintf(stderr, "sdfplugin) Unexpected EOF in atom block at atom %d\n", i+1);
      return -1;
    }

    /* V2000 atom line format (fixed-width):
     * xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
     * Positions: x(0-9), y(10-19), z(20-29), space(30), symbol(31-33),
     *            mass_diff(34-35), charge(36-38) */
    if (strlen(line) < 34) {
      fprintf(stderr, "sdfplugin) Atom line too short at atom %d\n", i+1);
      return -1;
    }

    /* Parse coordinates */
    {
      char buf[11];
      strncpy(buf, line, 10); buf[10] = '\0';
      x = (float)atof(buf);
      strncpy(buf, line + 10, 10); buf[10] = '\0';
      y = (float)atof(buf);
      strncpy(buf, line + 20, 10); buf[10] = '\0';
      z = (float)atof(buf);
    }

    /* Parse element symbol (columns 31-33, 3 chars, left-justified with space padding) */
    {
      char raw[4];
      int k, j;
      strncpy(raw, line + 31, 3);
      raw[3] = '\0';
      /* Strip leading/trailing whitespace */
      j = 0;
      for (k = 0; k < 3 && raw[k] != '\0'; k++) {
        if (!isspace((unsigned char)raw[k]))
          symbol[j++] = raw[k];
      }
      symbol[j] = '\0';
    }

    /* Parse charge code (columns 36-38) */
    if (strlen(line) >= 39) {
      char buf[4];
      strncpy(buf, line + 36, 3); buf[3] = '\0';
      chg_code = atoi(buf);
    }

    /* Store coordinates */
    if (first_pass && sdf->first_coords) {
      sdf->first_coords[3*i    ] = x;
      sdf->first_coords[3*i + 1] = y;
      sdf->first_coords[3*i + 2] = z;
    }
    if (ts) {
      ts->coords[3*i    ] = x;
      ts->coords[3*i + 1] = y;
      ts->coords[3*i + 2] = z;
    }

    /* Store atom data if doing structure read */
    if (atoms) {
      molfile_atom_t *atom = &atoms[i];
      const sdf_element_t *elem;

      strncpy(atom->name, symbol, 15);
      atom->name[15] = '\0';
      strncpy(atom->type, symbol, 15);
      atom->type[15] = '\0';

      strcpy(atom->resname, "MOL");
      atom->resid = 1;
      atom->chain[0] = '\0';
      atom->segid[0] = '\0';

      elem = lookup_element_by_symbol(symbol);
      if (elem) {
        atom->atomicnumber = elem->atomicnumber;
        atom->mass = elem->mass;
      } else {
        atom->atomicnumber = 0;
        atom->mass = 0.0f;
      }

      /* Initial charge from atom block charge code */
      atom->charge = charge_code_to_formal(chg_code);
    }

    /* Store per-atom data in handle for later use */
    if (first_pass) {
      const sdf_element_t *elem = lookup_element_by_symbol(symbol);
      if (sdf->atomicnumbers) {
        sdf->atomicnumbers[i] = elem ? elem->atomicnumber : 0;
      }
      if (sdf->masses) {
        sdf->masses[i] = elem ? elem->mass : 0.0f;
      }
      if (sdf->charges) {
        sdf->charges[i] = charge_code_to_formal(chg_code);
      }
    }
  }

  /* ---- Bond block ---- */
  for (i = 0; i < nbonds; i++) {
    int a1, a2, btype;

    if (fgets(line, SDF_LINESIZE, sdf->file) == NULL) {
      fprintf(stderr, "sdfplugin) Unexpected EOF in bond block at bond %d\n", i+1);
      return -1;
    }

    /* V2000 bond line format (fixed-width):
     * 111222tttsssxxxrrrccc
     * atom1(0-2), atom2(3-5), type(6-8) */
    if (strlen(line) < 9) {
      fprintf(stderr, "sdfplugin) Bond line too short at bond %d\n", i+1);
      return -1;
    }

    {
      char buf[4];
      strncpy(buf, line, 3); buf[3] = '\0';
      a1 = atoi(buf);
      strncpy(buf, line + 3, 3); buf[3] = '\0';
      a2 = atoi(buf);
      strncpy(buf, line + 6, 3); buf[3] = '\0';
      btype = atoi(buf);
    }

    if (first_pass && sdf->from && sdf->to && sdf->bondorder) {
      sdf->from[i] = a1;       /* 1-based as required by VMD */
      sdf->to[i] = a2;
      switch (btype) {
        case 1: sdf->bondorder[i] = 1.0f; break;
        case 2: sdf->bondorder[i] = 2.0f; break;
        case 3: sdf->bondorder[i] = 3.0f; break;
        case 4: sdf->bondorder[i] = 1.5f; break;  /* aromatic */
        default: sdf->bondorder[i] = 1.0f; break;
      }
    }
  }

  /* ---- Properties block and M lines ---- */
  while (fgets(line, SDF_LINESIZE, sdf->file) != NULL) {
    /* M  END terminates the connection table */
    if (strncmp(line, "M  END", 6) == 0)
      break;

    /* M  CHG line: formal charges that override atom block charges
     * Format: M  CHG  n  aaa vvv  aaa vvv ...
     * n = number of entries, aaa = atom number (1-based), vvv = charge */
    if (strncmp(line, "M  CHG", 6) == 0) {
      int nentries, j;
      char *ptr = line + 6;
      nentries = atoi(ptr);
      ptr += 3;  /* skip count field (3 chars) */

      for (j = 0; j < nentries; j++) {
        int atom_idx;
        float chg_val;

        if (strlen(ptr) < 8) break;
        {
          char buf_idx[5], buf_val[5];
          /* Each entry is: space + aaa(4 chars) + space + vvv(4 chars) = 9 chars */
          strncpy(buf_idx, ptr, 4); buf_idx[4] = '\0';
          atom_idx = atoi(buf_idx);
          strncpy(buf_val, ptr + 4, 4); buf_val[4] = '\0';
          chg_val = (float)atoi(buf_val);
          ptr += 8;
        }

        if (atom_idx >= 1 && atom_idx <= natoms) {
          if (first_pass && sdf->charges) {
            sdf->charges[atom_idx - 1] = chg_val;
          }
          if (atoms) {
            atoms[atom_idx - 1].charge = chg_val;
          }
        }
      }
    }
    /* Ignore other M lines */
  }

  return 0;
}


/* Helper: skip to end of current molecule (past $$$$) */
static int skip_to_next_molecule(FILE *f) {
  char line[SDF_LINESIZE];
  while (fgets(line, SDF_LINESIZE, f) != NULL) {
    if (strncmp(line, "$$$$", 4) == 0)
      return 0;
  }
  return -1; /* EOF */
}


/* ========================================================================
 * Reader functions
 * ======================================================================== */

static void *open_sdf_read(const char *path, const char *filetype, int *natoms) {
  FILE *fd;
  sdf_handle *sdf;
  int na = 0, nb = 0;
  char line[SDF_LINESIZE];

  fd = fopen(path, "r");
  if (!fd) {
    fprintf(stderr, "sdfplugin) Cannot open file: %s\n", path);
    return NULL;
  }

  /* Quick parse: read header (3 lines) + counts line to get natoms/nbonds */
  /* Line 1: name */
  if (fgets(line, SDF_LINESIZE, fd) == NULL) {
    fprintf(stderr, "sdfplugin) File is empty\n");
    fclose(fd);
    return NULL;
  }
  /* Line 2: program/date */
  if (fgets(line, SDF_LINESIZE, fd) == NULL) {
    fprintf(stderr, "sdfplugin) Unexpected EOF at header line 2\n");
    fclose(fd);
    return NULL;
  }
  /* Line 3: comment */
  if (fgets(line, SDF_LINESIZE, fd) == NULL) {
    fprintf(stderr, "sdfplugin) Unexpected EOF at header line 3\n");
    fclose(fd);
    return NULL;
  }
  /* Line 4: counts */
  if (fgets(line, SDF_LINESIZE, fd) == NULL) {
    fprintf(stderr, "sdfplugin) Unexpected EOF at counts line\n");
    fclose(fd);
    return NULL;
  }

  if (strlen(line) < 6) {
    fprintf(stderr, "sdfplugin) Counts line too short\n");
    fclose(fd);
    return NULL;
  }

  {
    char buf_a[4], buf_b[4];
    strncpy(buf_a, line, 3); buf_a[3] = '\0';
    strncpy(buf_b, line + 3, 3); buf_b[3] = '\0';
    na = atoi(buf_a);
    nb = atoi(buf_b);
  }

  if (na <= 0) {
    fprintf(stderr, "sdfplugin) Invalid atom count in counts line: %d\n", na);
    fclose(fd);
    return NULL;
  }

  /* Check for V3000 */
  if (strlen(line) >= 39 && strstr(line + 33, "V3000") != NULL) {
    fprintf(stderr, "sdfplugin) V3000 format not supported\n");
    fclose(fd);
    return NULL;
  }

  /* Allocate handle */
  sdf = (sdf_handle *)calloc(1, sizeof(sdf_handle));
  if (!sdf) {
    fclose(fd);
    return NULL;
  }

  sdf->file = fd;
  sdf->filepath = strdup(path);
  sdf->natoms = na;
  sdf->nbonds = nb;
  sdf->coords_read = 0;

  /* Allocate arrays for first-pass data */
  sdf->first_coords = (float *)calloc(na * 3, sizeof(float));
  sdf->from = (int *)calloc(nb > 0 ? nb : 1, sizeof(int));
  sdf->to = (int *)calloc(nb > 0 ? nb : 1, sizeof(int));
  sdf->bondorder = (float *)calloc(nb > 0 ? nb : 1, sizeof(float));
  sdf->atomicnumbers = (int *)calloc(na, sizeof(int));
  sdf->masses = (float *)calloc(na, sizeof(float));
  sdf->charges = (float *)calloc(na, sizeof(float));

  /* Now rewind and do full parse of first molecule */
  rewind(sdf->file);
  if (parse_mol_block(sdf, NULL, NULL, NULL, NULL, 1) != 0) {
    fprintf(stderr, "sdfplugin) Failed to parse first molecule\n");
    free(sdf->first_coords);
    free(sdf->from);
    free(sdf->to);
    free(sdf->bondorder);
    free(sdf->atomicnumbers);
    free(sdf->masses);
    free(sdf->charges);
    free(sdf->filepath);
    fclose(fd);
    free(sdf);
    return NULL;
  }

  /* Skip past data fields and $$$$ to position for trajectory reading */
  skip_to_next_molecule(sdf->file);

  *natoms = na;
  printf("sdfplugin) Opened %s: %d atoms, %d bonds\n", path, na, nb);

  return sdf;
}


static int read_sdf_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  sdf_handle *sdf = (sdf_handle *)v;
  int i;

  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_CHARGE;

  for (i = 0; i < sdf->natoms; i++) {
    molfile_atom_t *atom = &atoms[i];
    const sdf_element_t *elem;

    /* Use cached data from first parse */
    elem = lookup_element_by_number(sdf->atomicnumbers[i]);

    if (elem) {
      strncpy(atom->name, elem->symbol, 15);
      atom->name[15] = '\0';
      strncpy(atom->type, elem->symbol, 15);
      atom->type[15] = '\0';
    } else {
      strcpy(atom->name, "X");
      strcpy(atom->type, "X");
    }

    strcpy(atom->resname, "MOL");
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';

    atom->atomicnumber = sdf->atomicnumbers[i];
    atom->mass = sdf->masses[i];
    atom->charge = sdf->charges[i];
  }

  return MOLFILE_SUCCESS;
}


static int read_sdf_bonds(void *v, int *nbonds, int **fromptr, int **toptr,
                           float **bondorderptr, int **bondtype,
                           int *nbondtypes, char ***bondtypename) {
  sdf_handle *sdf = (sdf_handle *)v;

  if (sdf->nbonds > 0) {
    *nbonds = sdf->nbonds;
    *fromptr = sdf->from;
    *toptr = sdf->to;
    *bondorderptr = sdf->bondorder;
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


static int read_sdf_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  sdf_handle *sdf = (sdf_handle *)v;

  /* Frame 0: return cached coords from initial parse */
  if (sdf->coords_read == 0) {
    sdf->coords_read = 1;
    if (ts && sdf->first_coords) {
      memcpy(ts->coords, sdf->first_coords, natoms * 3 * sizeof(float));
    }
    return MOLFILE_SUCCESS;
  }

  /* Frame 1+: try to read next molecule from multi-SDF file */
  {
    int next_natoms = 0, next_nbonds = 0;
    int rc;

    /* Check for EOF before attempting to parse */
    if (feof(sdf->file))
      return MOLFILE_EOF;

    /* Peek at next character to check for EOF */
    {
      int c = fgetc(sdf->file);
      if (c == EOF) return MOLFILE_EOF;
      ungetc(c, sdf->file);
    }

    /* Parse next MOL block */
    rc = parse_mol_block(sdf, NULL, ts, &next_natoms, &next_nbonds, 0);
    if (rc != 0)
      return MOLFILE_EOF;

    /* Skip data fields and $$$$ separator */
    skip_to_next_molecule(sdf->file);

    /* Check if atom count matches - if so, this is a valid trajectory frame */
    if (next_natoms == sdf->natoms) {
      return MOLFILE_SUCCESS;
    }

    /* Atom count mismatch - this is a multi-molecule SDF, not a trajectory */
    fprintf(stderr, "sdfplugin) Multi-molecule SDF detected: molecule with %d atoms "
            "differs from first molecule (%d atoms).\n", next_natoms, sdf->natoms);
    fprintf(stderr, "sdfplugin) VMD loads each molecule as a separate entry. "
            "Use 'mol new file.sdf' for each molecule, or split the SDF file.\n");
    return MOLFILE_EOF;
  }
}


static void close_sdf_read(void *v) {
  sdf_handle *sdf = (sdf_handle *)v;
  if (sdf) {
    if (sdf->file) fclose(sdf->file);
    if (sdf->filepath) free(sdf->filepath);
    if (sdf->first_coords) free(sdf->first_coords);
    if (sdf->from) free(sdf->from);
    if (sdf->to) free(sdf->to);
    if (sdf->bondorder) free(sdf->bondorder);
    if (sdf->atomicnumbers) free(sdf->atomicnumbers);
    if (sdf->masses) free(sdf->masses);
    if (sdf->charges) free(sdf->charges);
    free(sdf);
  }
}


/* ========================================================================
 * Writer functions
 * ======================================================================== */

static void *open_sdf_write(const char *filename, const char *filetype,
                             int natoms) {
  FILE *fd;
  sdf_handle *sdf;

  fd = fopen(filename, "w");
  if (!fd) {
    fprintf(stderr, "sdfplugin) Error: unable to open file %s for writing\n",
            filename);
    return NULL;
  }

  sdf = (sdf_handle *)calloc(1, sizeof(sdf_handle));
  if (!sdf) {
    fclose(fd);
    return NULL;
  }

  sdf->file = fd;
  sdf->natoms = natoms;
  sdf->write_nbonds = 0;

  return sdf;
}


static int write_sdf_bonds(void *v, int nbonds, int *fromptr, int *toptr,
                            float *bondorderptr, int *bondtype,
                            int nbondtypes, char **bondtypename) {
  sdf_handle *sdf = (sdf_handle *)v;
  int i;

  sdf->write_nbonds = nbonds;

  if (nbonds > 0) {
    sdf->write_from = (int *)malloc(nbonds * sizeof(int));
    sdf->write_to = (int *)malloc(nbonds * sizeof(int));
    for (i = 0; i < nbonds; i++) {
      sdf->write_from[i] = fromptr[i];
      sdf->write_to[i] = toptr[i];
    }

    if (bondorderptr != NULL) {
      sdf->write_bondorder = (float *)malloc(nbonds * sizeof(float));
      for (i = 0; i < nbonds; i++) {
        sdf->write_bondorder[i] = bondorderptr[i];
      }
    }
  }

  return MOLFILE_SUCCESS;
}


static int write_sdf_structure(void *v, int optflags,
                                const molfile_atom_t *atoms) {
  sdf_handle *sdf = (sdf_handle *)v;
  sdf->optflags = optflags;
  sdf->atomlist = (molfile_atom_t *)malloc(sdf->natoms * sizeof(molfile_atom_t));
  memcpy(sdf->atomlist, atoms, sdf->natoms * sizeof(molfile_atom_t));
  return MOLFILE_SUCCESS;
}


static int write_sdf_timestep(void *v, const molfile_timestep_t *ts) {
  sdf_handle *sdf = (sdf_handle *)v;
  const molfile_atom_t *atom;
  const float *pos;
  int i;
  int has_mchg;
  int mchg_count;

  /* Line 1: molecule name */
  fprintf(sdf->file, "VMD generated molecule\n");

  /* Line 2: program/timestamp (IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR) */
  fprintf(sdf->file, "  VMD           3D\n");

  /* Line 3: comment */
  fprintf(sdf->file, "\n");

  /* Line 4: counts line  aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv */
  fprintf(sdf->file, "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n",
          sdf->natoms, sdf->write_nbonds);

  /* Atom block */
  atom = sdf->atomlist;
  pos = ts->coords;
  for (i = 0; i < sdf->natoms; i++) {
    const char *sym;
    int chg_code;

    /* Determine element symbol from atom name/type */
    sym = atom->name;
    if (sym[0] == '\0') sym = atom->type;
    if (sym[0] == '\0') sym = "X";

    /* Convert formal charge to V2000 charge code */
    chg_code = 0;
    if (sdf->optflags & MOLFILE_CHARGE) {
      chg_code = formal_to_charge_code(atom->charge);
    }

    /* V2000 atom line: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee */
    fprintf(sdf->file, "%10.4f%10.4f%10.4f %-3s 0%3d  0  0  0  0  0  0  0  0  0  0\n",
            pos[0], pos[1], pos[2], sym, chg_code);

    ++atom;
    pos += 3;
  }

  /* Bond block */
  for (i = 0; i < sdf->write_nbonds; i++) {
    int btype = 1;
    if (sdf->write_bondorder != NULL) {
      float bo = sdf->write_bondorder[i];
      if (bo >= 2.75f) btype = 3;
      else if (bo >= 1.75f) btype = 2;
      else if (bo >= 1.25f) btype = 4;  /* aromatic */
      else btype = 1;
    }
    fprintf(sdf->file, "%3d%3d%3d  0  0  0  0\n",
            sdf->write_from[i], sdf->write_to[i], btype);
  }

  /* M  CHG lines for non-zero charges */
  has_mchg = 0;
  mchg_count = 0;
  if (sdf->optflags & MOLFILE_CHARGE) {
    for (i = 0; i < sdf->natoms; i++) {
      if (sdf->atomlist[i].charge != 0.0f) {
        has_mchg = 1;
        mchg_count++;
      }
    }
  }

  if (has_mchg) {
    /* Write M  CHG lines (max 8 entries per line per spec) */
    int atom_idx = 0;
    while (mchg_count > 0) {
      int line_count = mchg_count > 8 ? 8 : mchg_count;
      int entries[8][2];
      int j, n = 0;

      /* Collect next batch of charged atoms */
      while (n < line_count && atom_idx < sdf->natoms) {
        if (sdf->atomlist[atom_idx].charge != 0.0f) {
          entries[n][0] = atom_idx + 1;  /* 1-based */
          entries[n][1] = (int)sdf->atomlist[atom_idx].charge;
          n++;
        }
        atom_idx++;
      }

      fprintf(sdf->file, "M  CHG%3d", n);
      for (j = 0; j < n; j++) {
        fprintf(sdf->file, " %3d %3d", entries[j][0], entries[j][1]);
      }
      fprintf(sdf->file, "\n");
      mchg_count -= n;
    }
  }

  /* M  END */
  fprintf(sdf->file, "M  END\n");

  /* $$$$ separator */
  fprintf(sdf->file, "$$$$\n");

  return MOLFILE_SUCCESS;
}


static void close_sdf_write(void *v) {
  sdf_handle *sdf = (sdf_handle *)v;
  if (sdf) {
    if (sdf->file) fclose(sdf->file);
    if (sdf->atomlist) free(sdf->atomlist);
    if (sdf->write_from) free(sdf->write_from);
    if (sdf->write_to) free(sdf->write_to);
    if (sdf->write_bondorder) free(sdf->write_bondorder);
    free(sdf);
  }
}


/* ========================================================================
 * Plugin registration
 * ======================================================================== */

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "sdf";
  plugin.prettyname = "SDF/MOL (V2000)";
  plugin.author = "Diego Enry Barreto Gomes";
  plugin.majorv = 1;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "sdf,mol";

  plugin.open_file_read = open_sdf_read;
  plugin.read_structure = read_sdf_structure;
  plugin.read_bonds = read_sdf_bonds;
  plugin.read_next_timestep = read_sdf_timestep;
  plugin.close_file_read = close_sdf_read;

  plugin.open_file_write = open_sdf_write;
  plugin.write_structure = write_sdf_structure;
  plugin.write_timestep = write_sdf_timestep;
  plugin.close_file_write = close_sdf_write;
  plugin.write_bonds = write_sdf_bonds;

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}
