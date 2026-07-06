/*
 * SDF/MOL V2000 plugin for VMD
 *
 * Native reader/writer for MDL SDF (Structure Data File) and MOL formats.
 * Supports V2000 connection tables with 3D coordinates, bonds with bond
 * orders, formal charges (both atom block charge codes and M CHG lines),
 * and multi-molecule SDF files as trajectory frames.
 *
 * Author: Diego Enry Barreto Gomes
 *
 * Hardening notes (v1.1):
 *   The reader is written to be safe against malformed, truncated, adversarial,
 *   or merely unusual SDF/MOL input.  In particular it never writes past the
 *   caller-owned coordinate buffer, validates all bond indices before handing
 *   them to VMD, tolerates over-long lines without desynchronizing, and checks
 *   every allocation.  The writer refuses to emit structures that the fixed-width
 *   V2000 format cannot represent (>999 atoms/bonds) instead of producing a
 *   silently corrupt file, and formats charges per the CTfile specification.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <locale.h>

#include "molfile_plugin.h"

/* Locale-independent, ASCII-only case-insensitive compare.  Avoids strcasecmp
 * (POSIX <strings.h>, absent on MSVC) and the locale sensitivity of tolower();
 * element symbols are pure ASCII so this is always correct. */
static int sdf_ascii_casecmp(const char *a, const char *b) {
  for (;; a++, b++) {
    int ca = (unsigned char)*a, cb = (unsigned char)*b;
    if (ca >= 'A' && ca <= 'Z') ca += 32;
    if (cb >= 'A' && cb <= 'Z') cb += 32;
    if (ca != cb) return ca - cb;
    if (ca == 0) return 0;
  }
}

/* SDF/V2000 mandates '.' as the decimal separator regardless of host locale.
 * localeconv() reports the active decimal point so we can translate to/from it;
 * this keeps numeric I/O correct even when VMD runs under a comma-decimal
 * locale (e.g. de_DE, fr_FR). */
static char locale_decimal_point(void) {
  struct lconv *lc = localeconv();
  return (lc && lc->decimal_point && lc->decimal_point[0]) ? lc->decimal_point[0] : '.';
}

#define SDF_LINESIZE 256

/* V2000 fixed-width count/index fields are 3 columns => max representable value */
#define V2000_MAX_FIELD 999

/* parse_mol_block return codes */
#define PARSE_OK        0
#define PARSE_EOF      -1   /* EOF or hard error */
#define PARSE_MISMATCH -2   /* atom count differs from expected (multi-molecule) */

/* ========================================================================
 * Element lookup table (Z = 1..118, standard atomic weights)
 * ======================================================================== */

typedef struct {
  const char *symbol;
  int atomicnumber;
  float mass;
} sdf_element_t;

static const sdf_element_t element_table[] = {
  {"H",   1,   1.008f},   {"He",  2,   4.003f},   {"Li",  3,   6.941f},
  {"Be",  4,   9.012f},   {"B",   5,  10.811f},   {"C",   6,  12.011f},
  {"N",   7,  14.007f},   {"O",   8,  15.999f},   {"F",   9,  18.998f},
  {"Ne", 10,  20.180f},   {"Na", 11,  22.990f},   {"Mg", 12,  24.305f},
  {"Al", 13,  26.982f},   {"Si", 14,  28.086f},   {"P",  15,  30.974f},
  {"S",  16,  32.065f},   {"Cl", 17,  35.453f},   {"Ar", 18,  39.948f},
  {"K",  19,  39.098f},   {"Ca", 20,  40.078f},   {"Sc", 21,  44.956f},
  {"Ti", 22,  47.867f},   {"V",  23,  50.942f},   {"Cr", 24,  51.996f},
  {"Mn", 25,  54.938f},   {"Fe", 26,  55.845f},   {"Co", 27,  58.933f},
  {"Ni", 28,  58.693f},   {"Cu", 29,  63.546f},   {"Zn", 30,  65.380f},
  {"Ga", 31,  69.723f},   {"Ge", 32,  72.640f},   {"As", 33,  74.922f},
  {"Se", 34,  78.960f},   {"Br", 35,  79.904f},   {"Kr", 36,  83.798f},
  {"Rb", 37,  85.468f},   {"Sr", 38,  87.620f},   {"Y",  39,  88.906f},
  {"Zr", 40,  91.224f},   {"Nb", 41,  92.906f},   {"Mo", 42,  95.960f},
  {"Tc", 43,  98.000f},   {"Ru", 44, 101.070f},   {"Rh", 45, 102.906f},
  {"Pd", 46, 106.420f},   {"Ag", 47, 107.868f},   {"Cd", 48, 112.411f},
  {"In", 49, 114.818f},   {"Sn", 50, 118.710f},   {"Sb", 51, 121.760f},
  {"Te", 52, 127.600f},   {"I",  53, 126.904f},   {"Xe", 54, 131.293f},
  {"Cs", 55, 132.905f},   {"Ba", 56, 137.327f},   {"La", 57, 138.905f},
  {"Ce", 58, 140.116f},   {"Pr", 59, 140.908f},   {"Nd", 60, 144.242f},
  {"Pm", 61, 145.000f},   {"Sm", 62, 150.360f},   {"Eu", 63, 151.964f},
  {"Gd", 64, 157.250f},   {"Tb", 65, 158.925f},   {"Dy", 66, 162.500f},
  {"Ho", 67, 164.930f},   {"Er", 68, 167.259f},   {"Tm", 69, 168.934f},
  {"Yb", 70, 173.054f},   {"Lu", 71, 174.967f},   {"Hf", 72, 178.490f},
  {"Ta", 73, 180.948f},   {"W",  74, 183.840f},   {"Re", 75, 186.207f},
  {"Os", 76, 190.230f},   {"Ir", 77, 192.217f},   {"Pt", 78, 195.078f},
  {"Au", 79, 196.967f},   {"Hg", 80, 200.590f},   {"Tl", 81, 204.383f},
  {"Pb", 82, 207.200f},   {"Bi", 83, 208.980f},   {"Po", 84, 209.000f},
  {"At", 85, 210.000f},   {"Rn", 86, 222.000f},   {"Fr", 87, 223.000f},
  {"Ra", 88, 226.000f},   {"Ac", 89, 227.000f},   {"Th", 90, 232.038f},
  {"Pa", 91, 231.036f},   {"U",  92, 238.029f},   {"Np", 93, 237.000f},
  {"Pu", 94, 244.000f},   {"Am", 95, 243.000f},   {"Cm", 96, 247.000f},
  {"Bk", 97, 247.000f},   {"Cf", 98, 251.000f},   {"Es", 99, 252.000f},
  {"Fm", 100, 257.000f},  {"Md", 101, 258.000f},  {"No", 102, 259.000f},
  {"Lr", 103, 262.000f},  {"Rf", 104, 267.000f},  {"Db", 105, 268.000f},
  {"Sg", 106, 271.000f},  {"Bh", 107, 272.000f},  {"Hs", 108, 270.000f},
  {"Mt", 109, 276.000f},  {"Ds", 110, 281.000f},  {"Rg", 111, 280.000f},
  {"Cn", 112, 285.000f},  {"Nh", 113, 286.000f},  {"Fl", 114, 289.000f},
  {"Mc", 115, 290.000f},  {"Lv", 116, 293.000f},  {"Ts", 117, 294.000f},
  {"Og", 118, 294.000f},
  {NULL,  0,   0.000f}
};

#define NUM_ELEMENTS (sizeof(element_table) / sizeof(element_table[0]) - 1)

/* Look up element by symbol; returns pointer or NULL */
static const sdf_element_t *lookup_element_by_symbol(const char *sym) {
  int i;
  if (sym == NULL || sym[0] == '\0')
    return NULL;
  for (i = 0; element_table[i].symbol != NULL; i++) {
    if (sdf_ascii_casecmp(sym, element_table[i].symbol) == 0)
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

/* Return the integer formal charge to write for a value, or 0 if it is not a
 * representable non-zero formal charge.  SDF/V2000 only encodes INTEGER formal
 * charges; a fractional force-field partial charge (e.g. -0.834) is NOT a
 * formal charge and must not be coerced into one, so only values that lie very
 * close to a non-zero integer are emitted. */
static int formal_charge_of(float c) {
  long r = lroundf(c);
  if (r == 0)
    return 0;
  if (fabsf(c - (float)r) > 0.05f)
    return 0;               /* genuine partial charge -> not representable */
  if (r >  99) r =  99;     /* keep the 3-column M CHG value field well-formed */
  if (r < -99) r = -99;
  return (int)r;
}


/* ========================================================================
 * Line I/O helper
 *
 * read_line() behaves like fgets() but ALWAYS consumes the whole physical
 * line: if the line is longer than the buffer, the tail is read and discarded
 * so the next read starts on the following line.  Without this, an over-long
 * name/comment/atom line (legal in SDF) would leave a fragment in the stream
 * that the next fgets() would misinterpret, desynchronizing the parser.
 * Returns buf, or NULL on EOF with no data read.
 * ======================================================================== */

static char *read_line(char *buf, int size, FILE *f) {
  size_t len;
  if (fgets(buf, size, f) == NULL)
    return NULL;
  len = strlen(buf);
  /* The line was truncated only if fgets completely filled the buffer without
   * capturing a newline.  Testing "buffer full" (len == size-1) rather than
   * "last char isn't newline" is robust to embedded NUL bytes, which would
   * otherwise make strlen under-report and spuriously drain the next line. */
  if (len == (size_t)(size - 1) && buf[len - 1] != '\n') {
    int c;
    while ((c = fgetc(f)) != EOF && c != '\n')
      ; /* discard remainder of the over-long line */
  }
  return buf;
}


/* ========================================================================
 * Plugin handle
 * ======================================================================== */

typedef struct {
  FILE *file;
  char *filepath;
  int natoms;
  int nbonds;                   /* number of VALID bonds actually stored */
  int coords_read;              /* have we returned the first frame yet? */
  float *first_coords;          /* cached coords from the initial parse */

  /* Bond data (reader) */
  int *from, *to;               /* 1-based bond indices */
  float *bondorder;

  /* Per-atom data cached at open time */
  int *atomicnumbers;
  float *masses;
  float *charges;
  char (*symbols)[8];           /* original element symbols, for fidelity */

  /* Writer cached data */
  molfile_atom_t *atomlist;     /* cached from write_structure */
  int optflags;
  int *write_from, *write_to;
  float *write_bondorder;
  int write_nbonds;
  int write_is_mol;             /* target is a plain .mol (no $$$$ separator) */
} sdf_handle;


/* ========================================================================
 * Helper: parse a signed integer from a fixed-width field of the line.
 * Copies up to width chars starting at offset into a temp buffer, NUL
 * terminates, and atoi()s it.  Never reads past the string terminator.
 * ======================================================================== */

static int parse_int_field(const char *line, int offset, int width) {
  char buf[16];
  int i, n = 0;
  size_t len = strlen(line);
  if ((size_t)offset >= len)
    return 0;
  if (width > (int)sizeof(buf) - 1)
    width = (int)sizeof(buf) - 1;
  for (i = 0; i < width && line[offset + i] != '\0'; i++)
    buf[n++] = line[offset + i];
  buf[n] = '\0';
  return atoi(buf);
}

static float parse_float_field(const char *line, int offset, int width) {
  char buf[32];
  char dp;
  int i, n = 0;
  size_t len = strlen(line);
  if ((size_t)offset >= len)
    return 0.0f;
  if (width > (int)sizeof(buf) - 1)
    width = (int)sizeof(buf) - 1;
  for (i = 0; i < width && line[offset + i] != '\0'; i++)
    buf[n++] = line[offset + i];
  buf[n] = '\0';
  /* The file always uses '.'; if the locale expects a different decimal point,
   * translate so atof() parses the value correctly. */
  dp = locale_decimal_point();
  if (dp != '.') {
    char *p;
    for (p = buf; *p; p++)
      if (*p == '.') *p = dp;
  }
  return (float)atof(buf);
}


/* ========================================================================
 * Helper: read and parse a complete MOL block (header through M END).
 *
 * coord_capacity : max number of atoms that ts->coords / first_coords can hold.
 *                  Coordinate writes past this index are suppressed so the
 *                  caller-owned buffer can never be overflowed.
 * expected_natoms: if >= 0 and the parsed atom count differs, the function
 *                  returns PARSE_MISMATCH BEFORE reading the atom block, so no
 *                  coordinates are written (prevents heap overflow on a
 *                  multi-molecule SDF whose later record is larger).
 * first_pass     : if non-zero, cache per-atom data and bonds into the handle.
 * saw_separator  : set to 1 if a "$$$$" record separator was consumed while
 *                  scanning the properties block (i.e. no M END was present).
 *
 * Returns PARSE_OK, PARSE_EOF, or PARSE_MISMATCH.
 * ======================================================================== */

static int parse_mol_block(sdf_handle *sdf, molfile_atom_t *atoms,
                           molfile_timestep_t *ts, int coord_capacity,
                           int expected_natoms, int first_pass,
                           int *out_natoms, int *out_nbonds,
                           int *saw_separator) {
  char line[SDF_LINESIZE];
  char version[8];
  int natoms, nbonds, i;
  int nvalid_bonds = 0;
  int mchg_reset = 0;

  if (saw_separator) *saw_separator = 0;

  /* Lines 1-3: name, program/timestamp, comment */
  if (read_line(line, SDF_LINESIZE, sdf->file) == NULL) return PARSE_EOF;
  if (read_line(line, SDF_LINESIZE, sdf->file) == NULL) return PARSE_EOF;
  if (read_line(line, SDF_LINESIZE, sdf->file) == NULL) return PARSE_EOF;

  /* Line 4: counts line "aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv" */
  if (read_line(line, SDF_LINESIZE, sdf->file) == NULL) return PARSE_EOF;
  if (strlen(line) < 6) {
    fprintf(stderr, "sdfplugin) Counts line too short\n");
    return PARSE_EOF;
  }

  natoms = parse_int_field(line, 0, 3);
  nbonds = parse_int_field(line, 3, 3);

  version[0] = '\0';
  if (strlen(line) >= 39) {
    strncpy(version, line + 33, 6);
    version[6] = '\0';
  }

  if (natoms <= 0) {
    fprintf(stderr, "sdfplugin) Invalid atom count: %d\n", natoms);
    return PARSE_EOF;
  }
  if (nbonds < 0)
    nbonds = 0;

  if (strstr(version, "V3000") != NULL) {
    fprintf(stderr, "sdfplugin) V3000 format not supported\n");
    return PARSE_EOF;
  }

  /* Multi-molecule guard: if this record's atom count differs from what the
   * caller's coordinate buffer was sized for, stop BEFORE writing any coords. */
  if (expected_natoms >= 0 && natoms != expected_natoms) {
    if (out_natoms) *out_natoms = natoms;
    if (out_nbonds) *out_nbonds = nbonds;
    return PARSE_MISMATCH;
  }

  if (out_natoms) *out_natoms = natoms;
  if (out_nbonds) *out_nbonds = nbonds;

  /* ---- Atom block ---- */
  for (i = 0; i < natoms; i++) {
    float x, y, z;
    char symbol[8];
    int chg_code = 0;

    if (read_line(line, SDF_LINESIZE, sdf->file) == NULL) {
      fprintf(stderr, "sdfplugin) Unexpected EOF in atom block at atom %d\n", i + 1);
      return PARSE_EOF;
    }

    /* V2000 atom line: x(0-9) y(10-19) z(20-29) space symbol(31-33) ... charge(36-38) */
    if (strlen(line) < 34) {
      fprintf(stderr, "sdfplugin) Atom line too short at atom %d\n", i + 1);
      return PARSE_EOF;
    }

    x = parse_float_field(line, 0, 10);
    y = parse_float_field(line, 10, 10);
    z = parse_float_field(line, 20, 10);

    /* Element symbol: columns 31-33, whitespace-stripped */
    {
      int k, j = 0;
      for (k = 31; k < 34 && line[k] != '\0'; k++) {
        if (!isspace((unsigned char)line[k]))
          symbol[j++] = line[k];
      }
      symbol[j] = '\0';
    }

    /* Charge code in columns 36-38 (if present) */
    if (strlen(line) >= 37)
      chg_code = parse_int_field(line, 36, 3);

    /* Store coordinates, guarded against the caller's buffer size. */
    if (i < coord_capacity) {
      if (first_pass && sdf->first_coords) {
        sdf->first_coords[3 * i]     = x;
        sdf->first_coords[3 * i + 1] = y;
        sdf->first_coords[3 * i + 2] = z;
      }
      if (ts && ts->coords) {
        ts->coords[3 * i]     = x;
        ts->coords[3 * i + 1] = y;
        ts->coords[3 * i + 2] = z;
      }
    }

    if (atoms) {
      molfile_atom_t *atom = &atoms[i];
      const sdf_element_t *elem;

      if (symbol[0] != '\0') {
        strncpy(atom->name, symbol, sizeof(atom->name) - 1);
        atom->name[sizeof(atom->name) - 1] = '\0';
        strncpy(atom->type, symbol, sizeof(atom->type) - 1);
        atom->type[sizeof(atom->type) - 1] = '\0';
      } else {
        strcpy(atom->name, "X");
        strcpy(atom->type, "X");
      }

      strcpy(atom->resname, "MOL");
      atom->resid = 1;
      atom->chain[0] = '\0';
      atom->segid[0] = '\0';

      elem = lookup_element_by_symbol(symbol);
      atom->atomicnumber = elem ? elem->atomicnumber : 0;
      atom->mass         = elem ? elem->mass : 0.0f;
      atom->charge       = charge_code_to_formal(chg_code);
    }

    if (first_pass) {
      const sdf_element_t *elem = lookup_element_by_symbol(symbol);
      if (sdf->atomicnumbers) sdf->atomicnumbers[i] = elem ? elem->atomicnumber : 0;
      if (sdf->masses)        sdf->masses[i]        = elem ? elem->mass : 0.0f;
      if (sdf->charges)       sdf->charges[i]       = charge_code_to_formal(chg_code);
      if (sdf->symbols) {
        strncpy(sdf->symbols[i], symbol, sizeof(sdf->symbols[i]) - 1);
        sdf->symbols[i][sizeof(sdf->symbols[i]) - 1] = '\0';
      }
    }
  }

  /* ---- Bond block ---- */
  {
  int pending = 0;   /* set when the bond loop over-reads into a property line */
  for (i = 0; i < nbonds; i++) {
    int a1, a2, btype;
    const char *p;

    if (read_line(line, SDF_LINESIZE, sdf->file) == NULL) {
      fprintf(stderr, "sdfplugin) Unexpected EOF in bond block at bond %d\n", i + 1);
      /* Keep whatever valid bonds we already have. */
      break;
    }

    /* If the counts line over-declared nbonds, the bond block ends early and
     * this line is actually a property/terminator line (M END, $$$$, M CHG,
     * an atom-alias block, ...).  A real bond record starts, after optional
     * leading spaces, with a digit.  Anything else means the bond block is
     * over; hand this line to the properties phase instead of mis-reading
     * "M  END"/"$$$$" as skippable short bond lines (which would swallow the
     * record terminator and desynchronize multi-molecule parsing). */
    p = line;
    while (*p == ' ') p++;
    if (!isdigit((unsigned char)*p)) {
      pending = 1;
      break;
    }

    if (strlen(line) < 9) {
      fprintf(stderr, "sdfplugin) Bond line too short at bond %d (skipped)\n", i + 1);
      continue;
    }

    a1    = parse_int_field(line, 0, 3);
    a2    = parse_int_field(line, 3, 3);
    btype = parse_int_field(line, 6, 3);

    /* Validate: indices must reference real atoms and not be self-bonds.
     * VMD indexes its atom array with these directly, so out-of-range values
     * are a memory-safety hazard downstream. */
    if (a1 < 1 || a1 > natoms || a2 < 1 || a2 > natoms || a1 == a2) {
      fprintf(stderr, "sdfplugin) Ignoring invalid bond %d: %d-%d (natoms=%d)\n",
              i + 1, a1, a2, natoms);
      continue;
    }

    if (first_pass && sdf->from && sdf->to && sdf->bondorder) {
      sdf->from[nvalid_bonds] = a1;
      sdf->to[nvalid_bonds]   = a2;
      switch (btype) {
        case 1:  sdf->bondorder[nvalid_bonds] = 1.0f; break;
        case 2:  sdf->bondorder[nvalid_bonds] = 2.0f; break;
        case 3:  sdf->bondorder[nvalid_bonds] = 3.0f; break;
        case 4:  sdf->bondorder[nvalid_bonds] = 1.5f; break;  /* aromatic */
        default: sdf->bondorder[nvalid_bonds] = 1.0f; break;
      }
    }
    nvalid_bonds++;
  }

  if (first_pass)
    sdf->nbonds = nvalid_bonds;
  if (out_nbonds)
    *out_nbonds = nvalid_bonds;

  /* ---- Properties block ----
   * If the bond loop over-read a property line, process that pending line
   * first instead of reading a new one. */
  for (;;) {
    if (pending) {
      pending = 0;
    } else if (read_line(line, SDF_LINESIZE, sdf->file) == NULL) {
      break;
    }

    /* A record separator here means the record had no M END; stop and let the
     * caller know the separator was already consumed. */
    if (strncmp(line, "$$$$", 4) == 0) {
      if (saw_separator) *saw_separator = 1;
      break;
    }
    /* M END terminates the connection table. */
    if (strncmp(line, "M  END", 6) == 0)
      break;

    /* Per the CTfile spec, once any M CHG OR M RAD line is present, the
     * atom-block charge column is void for ALL atoms; on the first such line
     * we clear every previously-derived charge before applying overrides. */
    if (strncmp(line, "M  CHG", 6) == 0 || strncmp(line, "M  RAD", 6) == 0) {
      if (!mchg_reset) {
        int k;
        if (first_pass && sdf->charges)
          for (k = 0; k < natoms; k++) sdf->charges[k] = 0.0f;
        if (atoms)
          for (k = 0; k < natoms; k++) atoms[k].charge = 0.0f;
        mchg_reset = 1;
      }
    }

    /* M CHG entries: "M  CHGnnn aaa vvv aaa vvv ...", 8-char entries after the
     * count.  (M RAD sets radical state we do not model, but it still voids the
     * block charges above.) */
    if (strncmp(line, "M  CHG", 6) == 0) {
      int nentries = parse_int_field(line, 6, 3);
      int j;
      const char *ptr = line + 9;  /* skip "M  CHG" + 3-char count */

      for (j = 0; j < nentries; j++) {
        int atom_idx, chg_val;
        if (strlen(ptr) < 8) break;
        atom_idx = parse_int_field(ptr, 0, 4);
        chg_val  = parse_int_field(ptr, 4, 4);
        ptr += 8;

        if (atom_idx >= 1 && atom_idx <= natoms) {
          if (first_pass && sdf->charges) sdf->charges[atom_idx - 1] = (float)chg_val;
          if (atoms)                      atoms[atom_idx - 1].charge = (float)chg_val;
        }
      }
    }
    /* Ignore all other M / property lines. */
  }
  }  /* end bond+properties scope */

  return PARSE_OK;
}


/* Helper: advance past the current record's "$$$$" separator. */
static int skip_to_next_molecule(FILE *f) {
  char line[SDF_LINESIZE];
  while (read_line(line, SDF_LINESIZE, f) != NULL) {
    if (strncmp(line, "$$$$", 4) == 0)
      return 0;
  }
  return -1; /* EOF */
}


/* ========================================================================
 * Reader functions
 * ======================================================================== */

static void free_read_handle(sdf_handle *sdf) {
  if (!sdf) return;
  if (sdf->file)          fclose(sdf->file);
  if (sdf->filepath)      free(sdf->filepath);
  if (sdf->first_coords)  free(sdf->first_coords);
  if (sdf->from)          free(sdf->from);
  if (sdf->to)            free(sdf->to);
  if (sdf->bondorder)     free(sdf->bondorder);
  if (sdf->atomicnumbers) free(sdf->atomicnumbers);
  if (sdf->masses)        free(sdf->masses);
  if (sdf->charges)       free(sdf->charges);
  if (sdf->symbols)       free(sdf->symbols);
  free(sdf);
}

static void *open_sdf_read(const char *path, const char *filetype, int *natoms) {
  FILE *fd;
  sdf_handle *sdf;
  int na = 0, nb = 0;
  int saw_sep = 0;
  char line[SDF_LINESIZE];

  (void)filetype;

  fd = fopen(path, "r");
  if (!fd) {
    fprintf(stderr, "sdfplugin) Cannot open file: %s\n", path);
    return NULL;
  }

  /* Quick scan of the header to size allocations. */
  if (read_line(line, SDF_LINESIZE, fd) == NULL) {
    fprintf(stderr, "sdfplugin) File is empty\n");
    fclose(fd);
    return NULL;
  }
  if (read_line(line, SDF_LINESIZE, fd) == NULL ||
      read_line(line, SDF_LINESIZE, fd) == NULL ||
      read_line(line, SDF_LINESIZE, fd) == NULL) {
    fprintf(stderr, "sdfplugin) Unexpected EOF in header\n");
    fclose(fd);
    return NULL;
  }
  if (strlen(line) < 6) {
    fprintf(stderr, "sdfplugin) Counts line too short\n");
    fclose(fd);
    return NULL;
  }

  na = parse_int_field(line, 0, 3);
  nb = parse_int_field(line, 3, 3);

  if (na <= 0) {
    fprintf(stderr, "sdfplugin) Invalid atom count in counts line: %d\n", na);
    fclose(fd);
    return NULL;
  }
  if (nb < 0) nb = 0;

  if (strlen(line) >= 39 && strstr(line + 33, "V3000") != NULL) {
    fprintf(stderr, "sdfplugin) V3000 format not supported\n");
    fclose(fd);
    return NULL;
  }

  sdf = (sdf_handle *)calloc(1, sizeof(sdf_handle));
  if (!sdf) {
    fclose(fd);
    return NULL;
  }

  sdf->file        = fd;
  sdf->filepath    = strdup(path);
  sdf->natoms      = na;
  sdf->nbonds      = nb;
  sdf->coords_read = 0;

  sdf->first_coords  = (float *)calloc((size_t)na * 3, sizeof(float));
  sdf->from          = (int *)calloc(nb > 0 ? nb : 1, sizeof(int));
  sdf->to            = (int *)calloc(nb > 0 ? nb : 1, sizeof(int));
  sdf->bondorder     = (float *)calloc(nb > 0 ? nb : 1, sizeof(float));
  sdf->atomicnumbers = (int *)calloc(na, sizeof(int));
  sdf->masses        = (float *)calloc(na, sizeof(float));
  sdf->charges       = (float *)calloc(na, sizeof(float));
  sdf->symbols       = (char (*)[8])calloc(na, sizeof(*sdf->symbols));

  if (!sdf->filepath || !sdf->first_coords || !sdf->from || !sdf->to ||
      !sdf->bondorder || !sdf->atomicnumbers || !sdf->masses ||
      !sdf->charges || !sdf->symbols) {
    fprintf(stderr, "sdfplugin) Out of memory opening %s\n", path);
    free_read_handle(sdf);   /* also closes fd */
    return NULL;
  }

  /* Rewind and fully parse the first record. */
  rewind(sdf->file);
  if (parse_mol_block(sdf, NULL, NULL, na, -1, 1, NULL, NULL, &saw_sep) != PARSE_OK) {
    fprintf(stderr, "sdfplugin) Failed to parse first molecule\n");
    free_read_handle(sdf);
    return NULL;
  }

  /* Position at the start of the next record for trajectory reads. */
  if (!saw_sep)
    skip_to_next_molecule(sdf->file);

  *natoms = na;
  printf("sdfplugin) Opened %s: %d atoms, %d bonds\n", path, na, sdf->nbonds);

  return sdf;
}


static int read_sdf_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  sdf_handle *sdf = (sdf_handle *)v;
  int i;

  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_CHARGE;

  for (i = 0; i < sdf->natoms; i++) {
    molfile_atom_t *atom = &atoms[i];
    const char *sym = sdf->symbols[i];

    if (sym[0] != '\0') {
      strncpy(atom->name, sym, sizeof(atom->name) - 1);
      atom->name[sizeof(atom->name) - 1] = '\0';
      strncpy(atom->type, sym, sizeof(atom->type) - 1);
      atom->type[sizeof(atom->type) - 1] = '\0';
    } else {
      strcpy(atom->name, "X");
      strcpy(atom->type, "X");
    }

    strcpy(atom->resname, "MOL");
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';

    atom->atomicnumber = sdf->atomicnumbers[i];
    atom->mass         = sdf->masses[i];
    atom->charge       = sdf->charges[i];
  }

  return MOLFILE_SUCCESS;
}


static int read_sdf_bonds(void *v, int *nbonds, int **fromptr, int **toptr,
                           float **bondorderptr, int **bondtype,
                           int *nbondtypes, char ***bondtypename) {
  sdf_handle *sdf = (sdf_handle *)v;

  if (sdf->nbonds > 0) {
    *nbonds       = sdf->nbonds;
    *fromptr      = sdf->from;
    *toptr        = sdf->to;
    *bondorderptr = sdf->bondorder;
  } else {
    *nbonds       = 0;
    *fromptr      = NULL;
    *toptr        = NULL;
    *bondorderptr = NULL;
  }

  *bondtype     = NULL;
  *nbondtypes   = 0;
  *bondtypename = NULL;

  return MOLFILE_SUCCESS;
}


static int read_sdf_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  sdf_handle *sdf = (sdf_handle *)v;

  /* Frame 0: return the coordinates cached during open. */
  if (sdf->coords_read == 0) {
    sdf->coords_read = 1;
    if (ts && ts->coords && sdf->first_coords) {
      int n = natoms < sdf->natoms ? natoms : sdf->natoms;  /* defensive clamp */
      memcpy(ts->coords, sdf->first_coords, (size_t)n * 3 * sizeof(float));
    }
    return MOLFILE_SUCCESS;
  }

  /* Frame 1+: read the next record from a multi-frame/multi-molecule SDF. */
  {
    int next_natoms = 0, next_nbonds = 0, saw_sep = 0, rc;
    int c;

    if (feof(sdf->file))
      return MOLFILE_EOF;

    c = fgetc(sdf->file);
    if (c == EOF)
      return MOLFILE_EOF;
    ungetc(c, sdf->file);

    /* Only accept a record with EXACTLY the established atom count as a frame.
     * A differing count returns PARSE_MISMATCH before any coordinate is written,
     * so the caller's ts->coords buffer can never be overrun. */
    rc = parse_mol_block(sdf, NULL, ts, natoms, sdf->natoms, 0,
                         &next_natoms, &next_nbonds, &saw_sep);

    if (rc == PARSE_MISMATCH) {
      fprintf(stderr, "sdfplugin) Multi-molecule SDF: record with %d atoms differs "
              "from the first molecule (%d atoms).\n", next_natoms, sdf->natoms);
      fprintf(stderr, "sdfplugin) VMD loads only the first molecule's trajectory; "
              "load additional molecules separately or split the SDF file.\n");
      return MOLFILE_EOF;
    }
    if (rc != PARSE_OK)
      return MOLFILE_EOF;

    if (!saw_sep)
      skip_to_next_molecule(sdf->file);

    return MOLFILE_SUCCESS;
  }
}


static void close_sdf_read(void *v) {
  free_read_handle((sdf_handle *)v);
}


/* ========================================================================
 * Writer functions
 * ======================================================================== */

/* Format a coordinate into a fixed 10-column field.  Non-finite values are
 * zeroed and out-of-range magnitudes are clamped, then the decimal precision
 * is reduced as needed so the result always fits exactly 10 columns. */
static void format_coord(char *out, size_t outsize, float v) {
  int prec;
  char dp;
  if (!isfinite(v)) v = 0.0f;
  if (v >  1.0e8f) v =  1.0e8f;
  if (v < -1.0e8f) v = -1.0e8f;
  out[0] = '\0';
  for (prec = 4; prec >= 0; prec--) {
    char tmp[32];
    snprintf(tmp, sizeof(tmp), "%.*f", prec, (double)v);
    if (strlen(tmp) <= 10) {
      snprintf(out, outsize, "%10.*f", prec, (double)v);
      break;
    }
  }
  if (out[0] == '\0')
    snprintf(out, outsize, "%10.0f", (double)v);
  /* Force '.' as the decimal separator regardless of the host locale. */
  dp = locale_decimal_point();
  if (dp != '.') {
    char *p;
    for (p = out; *p; p++)
      if (*p == dp) *p = '.';
  }
}

/* Choose the 1-3 char element symbol to write for an atom, preferring the
 * atomic number (reliable) over the atom name (which may be a PDB-style label
 * like "CA" that is NOT the calcium symbol). */
static void writer_symbol(const molfile_atom_t *atom, int optflags,
                          char *out, size_t outsize) {
  const sdf_element_t *elem = NULL;

  if ((optflags & MOLFILE_ATOMICNUMBER) && atom->atomicnumber > 0)
    elem = lookup_element_by_number(atom->atomicnumber);

  if (elem) {
    strncpy(out, elem->symbol, outsize - 1);
    out[outsize - 1] = '\0';
    return;
  }

  /* Fall back to the atom name/type, cleaned and capped to 3 chars. */
  {
    const char *src = (atom->name[0] != '\0') ? atom->name :
                      (atom->type[0] != '\0') ? atom->type : "X";
    size_t j = 0, k;
    for (k = 0; src[k] != '\0' && j < 3 && j < outsize - 1; k++) {
      if (!isspace((unsigned char)src[k]))
        out[j++] = src[k];
    }
    out[j] = '\0';
    if (out[0] == '\0') { out[0] = 'X'; out[1] = '\0'; }
  }
}

/* Case-insensitive test for a trailing extension (e.g. ".mol"). */
static int has_extension(const char *path, const char *ext) {
  size_t lp, le;
  if (!path || !ext) return 0;
  lp = strlen(path);
  le = strlen(ext);
  if (lp < le) return 0;
  return sdf_ascii_casecmp(path + (lp - le), ext) == 0;
}

static void *open_sdf_write(const char *filename, const char *filetype,
                             int natoms) {
  FILE *fd;
  sdf_handle *sdf;

  (void)filetype;

  if (natoms < 1) {
    fprintf(stderr, "sdfplugin) Refusing to write %d atoms\n", natoms);
    return NULL;
  }
  if (natoms > V2000_MAX_FIELD) {
    fprintf(stderr, "sdfplugin) V2000 cannot represent %d atoms (max %d); "
            "use a V3000-capable format.\n", natoms, V2000_MAX_FIELD);
    return NULL;
  }

  fd = fopen(filename, "w");
  if (!fd) {
    fprintf(stderr, "sdfplugin) Error: unable to open file %s for writing\n", filename);
    return NULL;
  }

  sdf = (sdf_handle *)calloc(1, sizeof(sdf_handle));
  if (!sdf) {
    fclose(fd);
    return NULL;
  }

  sdf->file        = fd;
  sdf->natoms      = natoms;
  sdf->write_nbonds = 0;
  /* A plain .mol file ends at "M  END" with no "$$$$" record separator. */
  sdf->write_is_mol = has_extension(filename, ".mol");

  return sdf;
}


static int write_sdf_bonds(void *v, int nbonds, int *fromptr, int *toptr,
                            float *bondorderptr, int *bondtype,
                            int nbondtypes, char **bondtypename) {
  sdf_handle *sdf = (sdf_handle *)v;
  int i;

  (void)bondtype; (void)nbondtypes; (void)bondtypename;

  if (nbonds < 0) nbonds = 0;
  sdf->write_nbonds = nbonds;

  if (nbonds > 0) {
    sdf->write_from = (int *)malloc(nbonds * sizeof(int));
    sdf->write_to   = (int *)malloc(nbonds * sizeof(int));
    if (!sdf->write_from || !sdf->write_to) {
      fprintf(stderr, "sdfplugin) Out of memory caching bonds\n");
      return MOLFILE_ERROR;
    }
    for (i = 0; i < nbonds; i++) {
      sdf->write_from[i] = fromptr[i];
      sdf->write_to[i]   = toptr[i];
    }
    if (bondorderptr != NULL) {
      sdf->write_bondorder = (float *)malloc(nbonds * sizeof(float));
      if (!sdf->write_bondorder) {
        fprintf(stderr, "sdfplugin) Out of memory caching bond orders\n");
        return MOLFILE_ERROR;
      }
      for (i = 0; i < nbonds; i++)
        sdf->write_bondorder[i] = bondorderptr[i];
    }
  }

  return MOLFILE_SUCCESS;
}


static int write_sdf_structure(void *v, int optflags,
                                const molfile_atom_t *atoms) {
  sdf_handle *sdf = (sdf_handle *)v;
  sdf->optflags = optflags;
  sdf->atomlist = (molfile_atom_t *)malloc(sdf->natoms * sizeof(molfile_atom_t));
  if (!sdf->atomlist) {
    fprintf(stderr, "sdfplugin) Out of memory caching structure\n");
    return MOLFILE_ERROR;
  }
  memcpy(sdf->atomlist, atoms, sdf->natoms * sizeof(molfile_atom_t));
  return MOLFILE_SUCCESS;
}


static int write_sdf_timestep(void *v, const molfile_timestep_t *ts) {
  sdf_handle *sdf = (sdf_handle *)v;
  const molfile_atom_t *atom;
  const float *pos;
  int i, mchg_count = 0;

  if (!sdf->atomlist) {
    fprintf(stderr, "sdfplugin) write_timestep called before write_structure\n");
    return MOLFILE_ERROR;
  }
  if (sdf->write_nbonds > V2000_MAX_FIELD) {
    fprintf(stderr, "sdfplugin) V2000 cannot represent %d bonds (max %d)\n",
            sdf->write_nbonds, V2000_MAX_FIELD);
    return MOLFILE_ERROR;
  }

  /* Header (3 lines) + counts line. */
  fprintf(sdf->file, "VMD generated molecule\n");
  fprintf(sdf->file, "  VMD           3D\n");
  fprintf(sdf->file, "\n");
  fprintf(sdf->file, "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n",
          sdf->natoms, sdf->write_nbonds);

  /* Atom block.  The charge column is always 0: formal charges are emitted via
   * M CHG lines below (per the CTfile spec, that voids the atom-block column,
   * so writing them here as well would double-specify the charge). */
  atom = sdf->atomlist;
  pos  = ts->coords;
  for (i = 0; i < sdf->natoms; i++) {
    char sym[8];
    char sx[16], sy[16], sz[16];

    writer_symbol(atom, sdf->optflags, sym, sizeof(sym));
    format_coord(sx, sizeof(sx), pos[0]);
    format_coord(sy, sizeof(sy), pos[1]);
    format_coord(sz, sizeof(sz), pos[2]);

    fprintf(sdf->file, "%s%s%s %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n",
            sx, sy, sz, sym);
    ++atom;
    pos += 3;
  }

  /* Bond block. */
  for (i = 0; i < sdf->write_nbonds; i++) {
    int btype = 1;
    if (sdf->write_bondorder != NULL) {
      float bo = sdf->write_bondorder[i];
      if (bo >= 2.75f)      btype = 3;
      else if (bo >= 1.75f) btype = 2;
      else if (bo >= 1.25f) btype = 4;  /* aromatic */
      else                  btype = 1;
    }
    fprintf(sdf->file, "%3d%3d%3d  0  0  0  0\n",
            sdf->write_from[i], sdf->write_to[i], btype);
  }

  /* M CHG lines for non-zero integer formal charges (max 8 entries per line). */
  if (sdf->optflags & MOLFILE_CHARGE) {
    for (i = 0; i < sdf->natoms; i++) {
      if (formal_charge_of(sdf->atomlist[i].charge) != 0) mchg_count++;
    }
  }

  if (mchg_count > 0) {
    int atom_idx = 0;
    while (mchg_count > 0) {
      int line_count = mchg_count > 8 ? 8 : mchg_count;
      int entries[8][2];
      int j, n = 0;

      while (n < line_count && atom_idx < sdf->natoms) {
        int fq = formal_charge_of(sdf->atomlist[atom_idx].charge);
        if (fq != 0) {
          entries[n][0] = atom_idx + 1;  /* 1-based */
          entries[n][1] = fq;
          n++;
        }
        atom_idx++;
      }

      fprintf(sdf->file, "M  CHG%3d", n);
      for (j = 0; j < n; j++)
        fprintf(sdf->file, " %3d %3d", entries[j][0], entries[j][1]);
      fprintf(sdf->file, "\n");
      mchg_count -= n;
    }
  }

  fprintf(sdf->file, "M  END\n");
  /* The "$$$$" separator is SDF-specific; a plain MOL file ends at M END. */
  if (!sdf->write_is_mol)
    fprintf(sdf->file, "$$$$\n");

  return MOLFILE_SUCCESS;
}


static void close_sdf_write(void *v) {
  sdf_handle *sdf = (sdf_handle *)v;
  if (sdf) {
    if (sdf->file)           fclose(sdf->file);
    if (sdf->atomlist)       free(sdf->atomlist);
    if (sdf->write_from)     free(sdf->write_from);
    if (sdf->write_to)       free(sdf->write_to);
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
  plugin.minorv = 1;
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
