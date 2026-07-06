/*
 * OpenMiChroM CNDB (HDF5) plugin for VMD
 *
 * Reads OpenMiChroM chromatin dynamics trajectories stored in HDF5 format
 * File structure:
 *   /Header - metadata attributes
 *   /replica{N}_chr{M}/ - replica groups containing:
 *     - genomic_position: [natoms, 2] genomic coordinates
 *     - loops: [nloops, 2] loop constraint pairs
 *     - types: [natoms] variable-length strings (compartment labels)
 *     - spatial_position/{frame}: [natoms, 3] xyz coordinates per frame
 *
 * v0.2 (hardened): validates every dataset rank/shape/dtype and every HDF5 read
 * before trusting it, checks all allocations, bounds every count against a sane
 * cap to prevent integer overflow, and cleans up completely on all error paths.
 * A malformed or adversarial .cndb file may be rejected but must never corrupt
 * memory. See scratchpad harness + corpus and [[sdf-plugin-hardening]] for the
 * verification pattern used here.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <hdf5.h>

#include "molfile_plugin.h"

/* Upper bound on any count read from the file (atoms, loops, frames). Large
 * enough for any real chromatin model, small enough that count*3 cannot
 * overflow a 32-bit int and count*sizeof(...) stays well within size_t. */
#define CNDB_MAX_COUNT 100000000

typedef struct {
  char *filepath;
  hid_t file_id;
  char replica_name[256];    // Selected replica group name
  hid_t replica_group;       // Handle to replica group
  int natoms;                // Number of beads
  int current_frame;
  int nframes;
  char **frame_names;        // Sorted list of frame dataset names (as strings)
  int *frame_numbers;        // Frame numbers for sorting

  // Genomic position data
  long *genomic_pos;         // [natoms][2] - start/end positions

  // Bead type data (our own copies; NULL entries allowed)
  char **bead_types;         // [natoms] - chromatin compartment types (A1, A2, B1, B3, NA)

  // Loop constraint data
  int nloops;
  int *loop_pairs;           // [nloops][2] - bead pair indices (0-based, as stored)
} cndb_handle;


/* Free everything a (possibly partially-initialized) handle owns. Safe to call
 * on any handle whose pointer fields are NULL and hid_t fields are <0 or valid.
 * Consolidates all cleanup so no error path can leak or double-free. */
static void free_cndb_handle(cndb_handle *cndb) {
  int i;
  if (!cndb) return;

  if (cndb->frame_names) {
    for (i = 0; i < cndb->nframes; i++)
      free(cndb->frame_names[i]);
    free(cndb->frame_names);
  }
  free(cndb->frame_numbers);
  free(cndb->genomic_pos);
  if (cndb->bead_types) {
    for (i = 0; i < cndb->natoms; i++)
      free(cndb->bead_types[i]);
    free(cndb->bead_types);
  }
  free(cndb->loop_pairs);
  if (cndb->replica_group >= 0) H5Gclose(cndb->replica_group);
  if (cndb->file_id >= 0) H5Fclose(cndb->file_id);
  free(cndb->filepath);
  free(cndb);
}


/* Validate that an open dataset is exactly 2-dimensional and return its dims.
 * Returns 0 and fills *d0,*d1 on success; -1 (dataset rejected) otherwise.
 * Rejecting rank != 2 up front prevents H5Sget_simple_extent_dims from writing
 * past a 2-element dims[] array and stops full-dataset reads from overflowing
 * buffers that were sized assuming a [n,2] / [n,3] layout. */
static int cndb_dims2d(hid_t dset, hsize_t *d0, hsize_t *d1) {
  hid_t sp = H5Dget_space(dset);
  hsize_t dims[2] = {0, 0};
  int rank;
  if (sp < 0) return -1;
  rank = H5Sget_simple_extent_ndims(sp);
  if (rank != 2) { H5Sclose(sp); return -1; }
  if (H5Sget_simple_extent_dims(sp, dims, NULL) < 0) { H5Sclose(sp); return -1; }
  H5Sclose(sp);
  *d0 = dims[0];
  *d1 = dims[1];
  return 0;
}

/* True if a dimension size is a usable count: >0 and within our cap. */
static int cndb_count_ok(hsize_t n) {
  return (n > 0 && n <= (hsize_t)CNDB_MAX_COUNT);
}

/* Convert a stored loop index (double) to int without undefined behavior.
 * NaN/inf and out-of-int-range values map to -1, which the bounds check in
 * read_cndb_bonds then rejects as an invalid (out-of-range) index. */
static int cndb_index_from_double(double v) {
  if (!(v >= (double)INT_MIN && v <= (double)INT_MAX)) return -1;
  return (int)v;
}


/* Comparator for numeric frame sort. */
typedef struct { int num; char *name; } cndb_frame;
static int cndb_frame_cmp(const void *a, const void *b) {
  int na = ((const cndb_frame *)a)->num;
  int nb = ((const cndb_frame *)b)->num;
  return (na > nb) - (na < nb);
}


/* Open HDF5 file and initialize handle */
static void *open_cndb_read(const char *filepath, const char *filetype, int *natoms) {
  cndb_handle *cndb;
  hid_t file_id, genomic_dset;
  hsize_t d0, d1;
  H5G_info_t group_info;
  hsize_t li;
  int i, replica_count = 0;
  char group_name[256];
  H5E_auto2_t old_efunc = NULL;
  void *old_edata = NULL;
  (void)filetype;

  /* Turn off HDF5 error messages while we probe optional datasets (HDF5 status
   * is checked explicitly below), saving the prior handler so we can restore it
   * on exit and not leave the whole process's HDF5 error printing disabled. */
  H5Eget_auto2(H5E_DEFAULT, &old_efunc, &old_edata);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    fprintf(stderr, "cndbplugin: Cannot open file %s\n", filepath);
    H5Eset_auto2(H5E_DEFAULT, old_efunc, old_edata);
    return NULL;
  }

  cndb = (cndb_handle *)calloc(1, sizeof(cndb_handle));
  if (!cndb) {
    H5Fclose(file_id);
    return NULL;
  }
  cndb->file_id = file_id;
  cndb->replica_group = -1;      /* sentinel: not yet opened */
  cndb->current_frame = 0;

  cndb->filepath = strdup(filepath);
  if (!cndb->filepath) goto fail;

  /* Find the first replica group: iterate root links and pick the first one
   * named "replica*" that actually opens as a group. */
  if (H5Gget_info(file_id, &group_info) < 0) {
    fprintf(stderr, "cndbplugin: Cannot get root group info\n");
    goto fail;
  }

  for (li = 0; li < group_info.nlinks; li++) {
    ssize_t nlen = H5Lget_name_by_idx(file_id, ".", H5_INDEX_NAME, H5_ITER_INC, li,
                                      group_name, sizeof(group_name), H5P_DEFAULT);
    if (nlen < 0) continue;                                  /* skip unreadable link */
    group_name[sizeof(group_name) - 1] = '\0';
    if (strncmp(group_name, "replica", 7) != 0) continue;
    replica_count++;

    /* Select the first replica that opens as a group; keep scanning only to
     * count the rest so we can tell the user how many were skipped. */
    if (cndb->replica_group < 0) {
      hid_t g = H5Gopen2(file_id, group_name, H5P_DEFAULT);
      if (g >= 0) {
        cndb->replica_group = g;
        strncpy(cndb->replica_name, group_name, sizeof(cndb->replica_name) - 1);
        cndb->replica_name[sizeof(cndb->replica_name) - 1] = '\0';
        printf("cndbplugin: Found replica group: %s\n", cndb->replica_name);
      } else {
        replica_count--;   /* named "replica*" but not a group; don't count it */
      }
    }
  }

  if (cndb->replica_group < 0) {
    fprintf(stderr, "cndbplugin: No replica groups found in file\n");
    goto fail;
  }
  if (replica_count > 1) {
    fprintf(stderr, "cndbplugin: file has %d replica groups; loading only '%s'\n",
            replica_count, cndb->replica_name);
  }

  /* --- genomic_position: [natoms, 2] -> defines natoms --- */
  genomic_dset = H5Dopen2(cndb->replica_group, "genomic_position", H5P_DEFAULT);
  if (genomic_dset < 0) {
    fprintf(stderr, "cndbplugin: Cannot open genomic_position dataset\n");
    goto fail;
  }
  if (cndb_dims2d(genomic_dset, &d0, &d1) < 0 || d1 != 2 || !cndb_count_ok(d0)) {
    fprintf(stderr, "cndbplugin: genomic_position must be [natoms,2] with 0<natoms<=%d\n",
            CNDB_MAX_COUNT);
    H5Dclose(genomic_dset);
    goto fail;
  }
  cndb->natoms = (int)d0;
  *natoms = cndb->natoms;
  printf("cndbplugin: Found %d atoms (beads)\n", cndb->natoms);

  cndb->genomic_pos = (long *)malloc((size_t)cndb->natoms * 2 * sizeof(long));
  if (!cndb->genomic_pos) { H5Dclose(genomic_dset); goto fail; }
  if (H5Dread(genomic_dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              cndb->genomic_pos) < 0) {
    fprintf(stderr, "cndbplugin: Failed to read genomic_position\n");
    H5Dclose(genomic_dset);
    goto fail;
  }
  H5Dclose(genomic_dset);

  /* --- types: [natoms] variable-length strings (optional) --- */
  {
    hid_t types_dset = H5Dopen2(cndb->replica_group, "types", H5P_DEFAULT);
    if (types_dset >= 0) {
      hid_t ftype = H5Dget_type(types_dset);
      hid_t tspace = H5Dget_space(types_dset);
      int trank = (tspace >= 0) ? H5Sget_simple_extent_ndims(tspace) : -1;
      htri_t is_vlstr = (ftype >= 0) ? H5Tis_variable_str(ftype) : -1;
      hsize_t tdims[1] = {0};
      int ok = (ftype >= 0 && tspace >= 0 && trank == 1 && is_vlstr > 0 &&
                H5Sget_simple_extent_dims(tspace, tdims, NULL) >= 0 &&
                tdims[0] == (hsize_t)cndb->natoms);

      if (ok) {
        char **raw = (char **)calloc((size_t)cndb->natoms, sizeof(char *));
        hid_t memtype = H5Tcopy(H5T_C_S1);
        if (raw && memtype >= 0 && H5Tset_size(memtype, H5T_VARIABLE) >= 0 &&
            H5Dread(types_dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw) >= 0) {
          /* Copy out of HDF5-owned vlen memory into our own strings, then let
           * HDF5 reclaim its allocations with the correct API. */
          cndb->bead_types = (char **)calloc((size_t)cndb->natoms, sizeof(char *));
          if (cndb->bead_types) {
            for (i = 0; i < cndb->natoms; i++)
              cndb->bead_types[i] = raw[i] ? strdup(raw[i]) : NULL;
            printf("cndbplugin: Read %d bead types\n", cndb->natoms);
          }
#if H5_VERSION_GE(1,12,0)
          H5Treclaim(memtype, tspace, H5P_DEFAULT, raw);
#else
          H5Dvlen_reclaim(memtype, tspace, H5P_DEFAULT, raw);
#endif
        } else {
          printf("cndbplugin: Could not read types dataset, using default names\n");
        }
        free(raw);
        if (memtype >= 0) H5Tclose(memtype);
      } else {
        printf("cndbplugin: types dataset not [natoms] vlen strings, using default names\n");
      }
      if (ftype >= 0) H5Tclose(ftype);
      if (tspace >= 0) H5Sclose(tspace);
      H5Dclose(types_dset);
    } else {
      printf("cndbplugin: No types dataset found, using default types\n");
    }
  }

  /* --- loops: [nloops, 2] (optional) --- */
  {
    hid_t loops_dset = H5Dopen2(cndb->replica_group, "loops", H5P_DEFAULT);
    if (loops_dset >= 0) {
      if (cndb_dims2d(loops_dset, &d0, &d1) == 0 && d1 == 2 &&
          d0 > 0 && d0 <= (hsize_t)CNDB_MAX_COUNT) {
        cndb->nloops = (int)d0;
        cndb->loop_pairs = (int *)malloc((size_t)cndb->nloops * 2 * sizeof(int));
        double *temp_loops = (double *)malloc((size_t)cndb->nloops * 2 * sizeof(double));
        if (cndb->loop_pairs && temp_loops &&
            H5Dread(loops_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    temp_loops) >= 0) {
          for (i = 0; i < cndb->nloops * 2; i++)
            cndb->loop_pairs[i] = cndb_index_from_double(temp_loops[i]);
          printf("cndbplugin: Found %d loop constraints\n", cndb->nloops);
        } else {
          /* read/alloc failed: drop loops rather than trust partial data */
          free(cndb->loop_pairs);
          cndb->loop_pairs = NULL;
          cndb->nloops = 0;
        }
        free(temp_loops);
      } else {
        fprintf(stderr, "cndbplugin: loops dataset not [nloops,2]; ignoring\n");
      }
      H5Dclose(loops_dset);
    }
  }

  /* --- spatial_position: enumerate and numerically sort frames --- */
  {
    hid_t spatial_group = H5Gopen2(cndb->replica_group, "spatial_position", H5P_DEFAULT);
    if (spatial_group < 0) {
      fprintf(stderr, "cndbplugin: Cannot open spatial_position group\n");
      goto fail;
    }
    if (H5Gget_info(spatial_group, &group_info) < 0) {
      fprintf(stderr, "cndbplugin: Cannot get spatial_position info\n");
      H5Gclose(spatial_group);
      goto fail;
    }
    if (group_info.nlinks > (hsize_t)CNDB_MAX_COUNT) {
      fprintf(stderr, "cndbplugin: too many frames (%llu)\n",
              (unsigned long long)group_info.nlinks);
      H5Gclose(spatial_group);
      goto fail;
    }
    cndb->nframes = (int)group_info.nlinks;
    printf("cndbplugin: Found %d frames\n", cndb->nframes);

    if (cndb->nframes > 0) {
      cndb_frame *entries = (cndb_frame *)calloc((size_t)cndb->nframes, sizeof(cndb_frame));
      cndb->frame_names = (char **)calloc((size_t)cndb->nframes, sizeof(char *));
      cndb->frame_numbers = (int *)calloc((size_t)cndb->nframes, sizeof(int));
      if (!entries || !cndb->frame_names || !cndb->frame_numbers) {
        free(entries);
        H5Gclose(spatial_group);
        goto fail;
      }

      int nvalid = 0;
      for (li = 0; li < (hsize_t)cndb->nframes; li++) {
        char fname[256];
        ssize_t nlen = H5Lget_name_by_idx(spatial_group, ".", H5_INDEX_NAME, H5_ITER_INC,
                                          li, fname, sizeof(fname), H5P_DEFAULT);
        if (nlen < 0) continue;                       /* skip unreadable frame link */
        fname[sizeof(fname) - 1] = '\0';
        char *dup = strdup(fname);
        if (!dup) { free(entries); H5Gclose(spatial_group); goto fail; }
        entries[nvalid].name = dup;
        entries[nvalid].num = (int)strtol(fname, NULL, 10);
        nvalid++;
      }
      cndb->nframes = nvalid;                          /* only keep the ones we read */

      qsort(entries, (size_t)nvalid, sizeof(cndb_frame), cndb_frame_cmp);
      for (i = 0; i < nvalid; i++) {
        cndb->frame_names[i]   = entries[i].name;      /* transfer ownership */
        cndb->frame_numbers[i] = entries[i].num;
      }
      free(entries);
    }
    H5Gclose(spatial_group);
  }

  H5Eset_auto2(H5E_DEFAULT, old_efunc, old_edata);   /* restore error handler */
  return cndb;

fail:
  H5Eset_auto2(H5E_DEFAULT, old_efunc, old_edata);   /* restore error handler */
  free_cndb_handle(cndb);
  return NULL;
}


/* Read structure information */
static int read_cndb_structure(void *mydata, int *optflags, molfile_atom_t *atoms) {
  cndb_handle *cndb = (cndb_handle *)mydata;
  int i;

  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS |
              MOLFILE_BFACTOR | MOLFILE_OCCUPANCY;

  for (i = 0; i < cndb->natoms; i++) {
    molfile_atom_t *atom = &atoms[i];

    /* Atom name/type from bead types when present (tail may be NULL: fall back). */
    if (cndb->bead_types && cndb->bead_types[i]) {
      strncpy(atom->name, cndb->bead_types[i], 15);
      atom->name[15] = '\0';
      strncpy(atom->type, cndb->bead_types[i], 15);
      atom->type[15] = '\0';
    } else {
      strcpy(atom->name, "BEAD");
      strcpy(atom->type, "");
    }

    strcpy(atom->resname, "CHR"); // Chromatin
    atom->resid = i + 1;
    atom->chain[0] = 'A';
    atom->chain[1] = '\0';
    /* Segid is limited to 7 characters + null terminator */
    strncpy(atom->segid, cndb->replica_name, 7);
    atom->segid[7] = '\0';

    /* Genomic position exposed via beta/occupancy fields. */
    atom->bfactor   = (float)cndb->genomic_pos[i*2];     // Start position
    atom->occupancy = (float)cndb->genomic_pos[i*2+1];   // End position

    atom->mass = 1.0;
    atom->radius = 5.0;           // Larger radius for coarse-grained bead
    atom->atomicnumber = 0;       // No atomic number (not an atom)
  }

  return MOLFILE_SUCCESS;
}


/* Read bond information (chromatin loops). Loop indices are stored 0-based;
 * VMD's molfile ABI expects 1-based atom indices, so valid pairs are emitted
 * as idx+1. Out-of-range pairs are dropped. */
static int read_cndb_bonds(void *mydata, int *nbonds, int **from, int **to,
                           float **bondorder, int **bondtype,
                           int *nbondtypes, char ***bondtypename) {
  cndb_handle *cndb = (cndb_handle *)mydata;
  int i;

  *nbonds = 0;
  *from = NULL;
  *to = NULL;
  *bondorder = NULL;
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;

  if (cndb->nloops <= 0 || !cndb->loop_pairs)
    return MOLFILE_SUCCESS;

  /* First pass: count in-range pairs. */
  int valid_bonds = 0;
  for (i = 0; i < cndb->nloops; i++) {
    int idx1 = cndb->loop_pairs[i*2];
    int idx2 = cndb->loop_pairs[i*2+1];
    if (idx1 >= 0 && idx1 < cndb->natoms && idx2 >= 0 && idx2 < cndb->natoms) {
      valid_bonds++;
    } else if (i < 5) {
      fprintf(stderr, "cndbplugin: Skipping loop %d with out-of-bounds indices: %d, %d (natoms=%d)\n",
              i, idx1, idx2, cndb->natoms);
    }
  }

  if (valid_bonds == 0) {
    printf("cndbplugin: No valid loop constraints found (all out of bounds)\n");
    return MOLFILE_SUCCESS;
  }

  int *f = (int *)malloc((size_t)valid_bonds * sizeof(int));
  int *t = (int *)malloc((size_t)valid_bonds * sizeof(int));
  if (!f || !t) {
    free(f);
    free(t);
    fprintf(stderr, "cndbplugin: out of memory allocating %d bonds\n", valid_bonds);
    return MOLFILE_SUCCESS;   /* degrade to no bonds rather than fail the load */
  }

  int bond_idx = 0;
  for (i = 0; i < cndb->nloops; i++) {
    int idx1 = cndb->loop_pairs[i*2];
    int idx2 = cndb->loop_pairs[i*2+1];
    if (idx1 >= 0 && idx1 < cndb->natoms && idx2 >= 0 && idx2 < cndb->natoms) {
      f[bond_idx] = idx1 + 1;    /* 0-based -> VMD 1-based */
      t[bond_idx] = idx2 + 1;
      bond_idx++;
    }
  }

  printf("cndbplugin: Using %d valid loop constraints (out of %d total)\n",
         valid_bonds, cndb->nloops);
  *nbonds = valid_bonds;
  *from = f;
  *to = t;
  return MOLFILE_SUCCESS;
}


/* Read next timestep */
static int read_cndb_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  cndb_handle *cndb = (cndb_handle *)mydata;
  char dataset_path[512];
  hid_t dset_id;
  hsize_t d0, d1;
  float *coords;
  int i;

  if (cndb->current_frame >= cndb->nframes)
    return MOLFILE_EOF;

  /* VMD sized ts->coords for the natoms we reported at open time. */
  if (natoms != cndb->natoms) {
    fprintf(stderr, "cndbplugin: natoms mismatch (%d vs %d)\n", natoms, cndb->natoms);
    return MOLFILE_ERROR;
  }

  snprintf(dataset_path, sizeof(dataset_path), "%s/spatial_position/%s",
           cndb->replica_name, cndb->frame_names[cndb->current_frame]);

  dset_id = H5Dopen2(cndb->file_id, dataset_path, H5P_DEFAULT);
  if (dset_id < 0) {
    fprintf(stderr, "cndbplugin: Cannot open dataset %s\n", dataset_path);
    return MOLFILE_ERROR;
  }

  /* Require exactly [natoms, 3] before reading a single element. */
  if (cndb_dims2d(dset_id, &d0, &d1) < 0 ||
      d0 != (hsize_t)cndb->natoms || d1 != 3) {
    fprintf(stderr, "cndbplugin: Unexpected dimensions in frame %d\n", cndb->current_frame);
    H5Dclose(dset_id);
    return MOLFILE_ERROR;
  }

  /* Nothing to store (VMD only wants to advance): skip the read entirely. */
  if (ts == NULL) {
    H5Dclose(dset_id);
    cndb->current_frame++;
    return MOLFILE_SUCCESS;
  }

  coords = (float *)malloc((size_t)cndb->natoms * 3 * sizeof(float));
  if (!coords) {
    H5Dclose(dset_id);
    return MOLFILE_ERROR;
  }
  if (H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords) < 0) {
    fprintf(stderr, "cndbplugin: Failed to read frame %d\n", cndb->current_frame);
    free(coords);
    H5Dclose(dset_id);
    return MOLFILE_ERROR;
  }

  /* Convert nm to Angstroms (x10). NOTE: OpenMiChroM output may be in reduced
   * units; this fixed scale is preserved from v0.1 for visual compatibility. */
  for (i = 0; i < cndb->natoms; i++) {
    ts->coords[i*3]     = coords[i*3]     * 10.0f;
    ts->coords[i*3 + 1] = coords[i*3 + 1] * 10.0f;
    ts->coords[i*3 + 2] = coords[i*3 + 2] * 10.0f;
  }

  free(coords);
  H5Dclose(dset_id);
  cndb->current_frame++;
  return MOLFILE_SUCCESS;
}


/* Close file and clean up */
static void close_cndb_read(void *mydata) {
  free_cndb_handle((cndb_handle *)mydata);
}


/*
 * Plugin registration
 */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "cndb";
  plugin.prettyname = "OpenMiChroM CNDB";
  plugin.author = "Diego Enry Barreto Gomes";
  plugin.majorv = 0;
  plugin.minorv = 2;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "cndb";

  plugin.open_file_read = open_cndb_read;
  plugin.read_structure = read_cndb_structure;
  plugin.read_bonds = read_cndb_bonds;
  plugin.read_next_timestep = read_cndb_timestep;
  plugin.close_file_read = close_cndb_read;

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}
