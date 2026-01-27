/*
 * OpenMiChroM CNDB (HDF5) plugin for VMD
 *
 * Reads OpenMiChroM chromatin dynamics trajectories stored in HDF5 format
 * File structure:
 *   /Header - metadata attributes
 *   /replica{N}_chr{M}/ - replica groups containing:
 *     - genomic_position: [natoms, 2] genomic coordinates
 *     - loops: [nloops, 2] loop constraint pairs
 *     - spatial_position/{frame}: [natoms, 3] xyz coordinates per frame
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <hdf5.h>

#include "molfile_plugin.h"

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

  // Bead type data
  char **bead_types;         // [natoms] - chromatin compartment types (A1, A2, B1, B3, NA)

  // Loop constraint data
  int nloops;
  int *loop_pairs;           // [nloops][2] - bead pair indices
} cndb_handle;


/* Open HDF5 file and initialize handle */
static void *open_cndb_read(const char *filepath, const char *filetype, int *natoms) {
  cndb_handle *cndb;
  hid_t file_id, space_id;
  hsize_t dims[2];
  H5G_info_t group_info;
  herr_t status;
  int i;
  char group_name[256];

  /* Turn off HDF5 error messages */
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  /* Open HDF5 file */
  file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    fprintf(stderr, "cndbplugin: Cannot open file %s\n", filepath);
    return NULL;
  }

  /* Allocate handle */
  cndb = (cndb_handle *)calloc(1, sizeof(cndb_handle));
  if (!cndb) {
    H5Fclose(file_id);
    return NULL;
  }

  cndb->filepath = strdup(filepath);
  cndb->file_id = file_id;
  cndb->current_frame = 0;

  /* Find first replica group - look for groups matching replica*_chr* pattern */
  status = H5Gget_info(file_id, &group_info);
  if (status < 0) {
    fprintf(stderr, "cndbplugin: Cannot get root group info\n");
    free(cndb->filepath);
    free(cndb);
    H5Fclose(file_id);
    return NULL;
  }

  /* Iterate through root groups to find first replica */
  int found_replica = 0;
  for (i = 0; i < group_info.nlinks; i++) {
    H5Lget_name_by_idx(file_id, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                       group_name, sizeof(group_name), H5P_DEFAULT);

    /* Check if group name starts with "replica" */
    if (strncmp(group_name, "replica", 7) == 0) {
      strncpy(cndb->replica_name, group_name, sizeof(cndb->replica_name) - 1);
      found_replica = 1;
      printf("cndbplugin: Found replica group: %s\n", group_name);
      break;
    }
  }

  if (!found_replica) {
    fprintf(stderr, "cndbplugin: No replica groups found in file\n");
    free(cndb->filepath);
    free(cndb);
    H5Fclose(file_id);
    return NULL;
  }

  /* Open replica group */
  cndb->replica_group = H5Gopen2(file_id, cndb->replica_name, H5P_DEFAULT);
  if (cndb->replica_group < 0) {
    fprintf(stderr, "cndbplugin: Cannot open replica group %s\n", cndb->replica_name);
    free(cndb->filepath);
    free(cndb);
    H5Fclose(file_id);
    return NULL;
  }

  /* Read genomic_position dataset to get number of atoms */
  hid_t genomic_dset = H5Dopen2(cndb->replica_group, "genomic_position", H5P_DEFAULT);
  if (genomic_dset < 0) {
    fprintf(stderr, "cndbplugin: Cannot open genomic_position dataset\n");
    H5Gclose(cndb->replica_group);
    free(cndb->filepath);
    free(cndb);
    H5Fclose(file_id);
    return NULL;
  }

  space_id = H5Dget_space(genomic_dset);
  H5Sget_simple_extent_dims(space_id, dims, NULL);
  cndb->natoms = dims[0];
  *natoms = cndb->natoms;

  printf("cndbplugin: Found %d atoms (beads)\n", cndb->natoms);

  /* Allocate and read genomic positions */
  cndb->genomic_pos = (long *)malloc(cndb->natoms * 2 * sizeof(long));
  H5Dread(genomic_dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, cndb->genomic_pos);
  H5Sclose(space_id);
  H5Dclose(genomic_dset);

  /* Read types dataset for bead types */
  hid_t types_dset = H5Dopen2(cndb->replica_group, "types", H5P_DEFAULT);
  if (types_dset >= 0) {
    hid_t str_type = H5Dget_type(types_dset);

    /* Allocate array of string pointers */
    cndb->bead_types = (char **)malloc(cndb->natoms * sizeof(char *));

    /* Read variable-length strings */
    H5Dread(types_dset, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, cndb->bead_types);

    printf("cndbplugin: Read %d bead types\n", cndb->natoms);

    H5Tclose(str_type);
    H5Dclose(types_dset);
  } else {
    /* If no types dataset, set to NULL */
    cndb->bead_types = NULL;
    printf("cndbplugin: No types dataset found, using default types\n");
  }

  /* Read loops dataset */
  hid_t loops_dset = H5Dopen2(cndb->replica_group, "loops", H5P_DEFAULT);
  if (loops_dset >= 0) {
    space_id = H5Dget_space(loops_dset);
    H5Sget_simple_extent_dims(space_id, dims, NULL);
    cndb->nloops = dims[0];

    /* Allocate and read loop pairs */
    cndb->loop_pairs = (int *)malloc(cndb->nloops * 2 * sizeof(int));
    double *temp_loops = (double *)malloc(cndb->nloops * 2 * sizeof(double));
    H5Dread(loops_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_loops);

    /* Convert to integers */
    for (i = 0; i < cndb->nloops * 2; i++) {
      cndb->loop_pairs[i] = (int)temp_loops[i];
    }
    free(temp_loops);

    printf("cndbplugin: Found %d loop constraints\n", cndb->nloops);
    H5Sclose(space_id);
    H5Dclose(loops_dset);
  }

  /* Open spatial_position group and enumerate frames */
  hid_t spatial_group = H5Gopen2(cndb->replica_group, "spatial_position", H5P_DEFAULT);
  if (spatial_group < 0) {
    fprintf(stderr, "cndbplugin: Cannot open spatial_position group\n");
    free(cndb->genomic_pos);
    if (cndb->loop_pairs) free(cndb->loop_pairs);
    H5Gclose(cndb->replica_group);
    free(cndb->filepath);
    free(cndb);
    H5Fclose(file_id);
    return NULL;
  }

  H5Gget_info(spatial_group, &group_info);
  cndb->nframes = group_info.nlinks;
  printf("cndbplugin: Found %d frames\n", cndb->nframes);

  /* Allocate frame name and number arrays */
  cndb->frame_names = (char **)malloc(cndb->nframes * sizeof(char *));
  cndb->frame_numbers = (int *)malloc(cndb->nframes * sizeof(int));

  /* Read frame names and convert to numbers for sorting */
  for (i = 0; i < cndb->nframes; i++) {
    char fname[64];
    H5Lget_name_by_idx(spatial_group, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                       fname, sizeof(fname), H5P_DEFAULT);
    cndb->frame_names[i] = strdup(fname);
    cndb->frame_numbers[i] = atoi(fname);
  }

  /* Sort frames numerically */
  /* Use simple bubble sort with parallel swap of names */
  for (i = 0; i < cndb->nframes - 1; i++) {
    for (int j = 0; j < cndb->nframes - i - 1; j++) {
      if (cndb->frame_numbers[j] > cndb->frame_numbers[j+1]) {
        /* Swap numbers */
        int temp_num = cndb->frame_numbers[j];
        cndb->frame_numbers[j] = cndb->frame_numbers[j+1];
        cndb->frame_numbers[j+1] = temp_num;

        /* Swap names */
        char *temp_name = cndb->frame_names[j];
        cndb->frame_names[j] = cndb->frame_names[j+1];
        cndb->frame_names[j+1] = temp_name;
      }
    }
  }

  H5Gclose(spatial_group);

  return cndb;
}


/* Read structure information */
static int read_cndb_structure(void *mydata, int *optflags, molfile_atom_t *atoms) {
  cndb_handle *cndb = (cndb_handle *)mydata;
  int i;

  /* Set optional flags */
  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS |
              MOLFILE_BFACTOR | MOLFILE_OCCUPANCY;

  /* Fill in atom data */
  for (i = 0; i < cndb->natoms; i++) {
    molfile_atom_t *atom = &atoms[i];

    /* Atom properties - use bead types if available */
    if (cndb->bead_types && cndb->bead_types[i]) {
      /* Use bead type from dataset for both name and type */
      strncpy(atom->name, cndb->bead_types[i], 15);
      atom->name[15] = '\0';
      strncpy(atom->type, cndb->bead_types[i], 15);
      atom->type[15] = '\0';
    } else {
      /* Fallback to generic bead name */
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

    /* Use genomic position in beta/occupancy fields */
    atom->bfactor = (float)cndb->genomic_pos[i*2];     // Start position
    atom->occupancy = (float)cndb->genomic_pos[i*2+1]; // End position

    atom->mass = 1.0;
    atom->radius = 5.0;           // Larger radius for coarse-grained bead
    atom->atomicnumber = 0;       // No atomic number (not an atom)
  }

  return MOLFILE_SUCCESS;
}


/* Read bond information (chromatin loops) */
static int read_cndb_bonds(void *mydata, int *nbonds, int **from, int **to,
                           float **bondorder, int **bondtype,
                           int *nbondtypes, char ***bondtypename) {
  cndb_handle *cndb = (cndb_handle *)mydata;
  int i;

  if (cndb->nloops == 0) {
    *nbonds = 0;
    return MOLFILE_SUCCESS;
  }

  /* First pass: count valid bonds */
  int valid_bonds = 0;
  for (i = 0; i < cndb->nloops; i++) {
    int idx1 = cndb->loop_pairs[i*2];
    int idx2 = cndb->loop_pairs[i*2+1];

    /* Check if indices are within bounds */
    if (idx1 >= 0 && idx1 < cndb->natoms && idx2 >= 0 && idx2 < cndb->natoms) {
      valid_bonds++;
    } else {
      if (i < 5) {  /* Only print first few warnings */
        fprintf(stderr, "cndbplugin: Skipping loop %d with out-of-bounds indices: %d, %d (natoms=%d)\n",
                i, idx1, idx2, cndb->natoms);
      }
    }
  }

  if (valid_bonds == 0) {
    printf("cndbplugin: No valid loop constraints found (all out of bounds)\n");
    *nbonds = 0;
    return MOLFILE_SUCCESS;
  }

  printf("cndbplugin: Using %d valid loop constraints (out of %d total)\n", valid_bonds, cndb->nloops);
  *nbonds = valid_bonds;

  /* Allocate bond arrays */
  *from = (int *)malloc(valid_bonds * sizeof(int));
  *to = (int *)malloc(valid_bonds * sizeof(int));
  *bondorder = NULL;  // Not used
  *bondtype = NULL;   // Not used

  /* Fill in valid loop constraints as bonds */
  int bond_idx = 0;
  for (i = 0; i < cndb->nloops; i++) {
    int idx1 = cndb->loop_pairs[i*2];
    int idx2 = cndb->loop_pairs[i*2+1];

    /* Bounds check */
    if (idx1 >= 0 && idx1 < cndb->natoms && idx2 >= 0 && idx2 < cndb->natoms) {
      (*from)[bond_idx] = idx1;
      (*to)[bond_idx] = idx2;
      bond_idx++;
    }
  }

  *nbondtypes = 0;
  *bondtypename = NULL;

  return MOLFILE_SUCCESS;
}


/* Read next timestep */
static int read_cndb_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  cndb_handle *cndb = (cndb_handle *)mydata;
  char dataset_path[512];
  hid_t dset_id, space_id;
  hsize_t dims[2];
  float *coords;
  int i;

  if (cndb->current_frame >= cndb->nframes) {
    return MOLFILE_EOF;
  }

  /* Build path to spatial position dataset */
  snprintf(dataset_path, sizeof(dataset_path), "%s/spatial_position/%s",
           cndb->replica_name, cndb->frame_names[cndb->current_frame]);

  /* Open dataset */
  dset_id = H5Dopen2(cndb->file_id, dataset_path, H5P_DEFAULT);
  if (dset_id < 0) {
    fprintf(stderr, "cndbplugin: Cannot open dataset %s\n", dataset_path);
    return MOLFILE_ERROR;
  }

  /* Verify dimensions */
  space_id = H5Dget_space(dset_id);
  H5Sget_simple_extent_dims(space_id, dims, NULL);
  if (dims[0] != cndb->natoms || dims[1] != 3) {
    fprintf(stderr, "cndbplugin: Unexpected dimensions in frame %d\n", cndb->current_frame);
    H5Sclose(space_id);
    H5Dclose(dset_id);
    return MOLFILE_ERROR;
  }

  /* Allocate temporary buffer and read coordinates */
  coords = (float *)malloc(cndb->natoms * 3 * sizeof(float));
  H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords);

  /* Copy to timestep structure and convert nm to Angstroms (multiply by 10) */
  if (ts != NULL) {
    for (i = 0; i < cndb->natoms; i++) {
      ts->coords[i*3]     = coords[i*3]     * 10.0f;
      ts->coords[i*3 + 1] = coords[i*3 + 1] * 10.0f;
      ts->coords[i*3 + 2] = coords[i*3 + 2] * 10.0f;
    }
  }

  free(coords);
  H5Sclose(space_id);
  H5Dclose(dset_id);

  cndb->current_frame++;

  return MOLFILE_SUCCESS;
}


/* Close file and clean up */
static void close_cndb_read(void *mydata) {
  cndb_handle *cndb = (cndb_handle *)mydata;
  int i;

  if (cndb) {
    if (cndb->frame_names) {
      for (i = 0; i < cndb->nframes; i++) {
        free(cndb->frame_names[i]);
      }
      free(cndb->frame_names);
    }
    if (cndb->frame_numbers) free(cndb->frame_numbers);
    if (cndb->genomic_pos) free(cndb->genomic_pos);
    if (cndb->bead_types) {
      /* Free variable-length strings allocated by HDF5 */
      for (i = 0; i < cndb->natoms; i++) {
        if (cndb->bead_types[i]) {
          free(cndb->bead_types[i]);
        }
      }
      free(cndb->bead_types);
    }
    if (cndb->loop_pairs) free(cndb->loop_pairs);
    if (cndb->replica_group >= 0) H5Gclose(cndb->replica_group);
    if (cndb->file_id >= 0) H5Fclose(cndb->file_id);
    if (cndb->filepath) free(cndb->filepath);
    free(cndb);
  }
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
  plugin.minorv = 1;
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
