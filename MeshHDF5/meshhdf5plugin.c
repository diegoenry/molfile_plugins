/*
 * HDF5 Mesh (.h5mesh) graphics plugin for VMD
 *
 * Reads triangle meshes stored in HDF5 format and renders them as
 * VMD graphics objects (TRINORM + NORMS pairs), identical to how
 * VMD's built-in STL plugin renders meshes — but with proper per-face
 * normals read from the file for correct lighting.
 *
 * File structure expected:
 *   /vertices    [N, 3] float32  — vertex xyz coordinates
 *   /triangles   [M, 3] int32    — triangle vertex indices (0-based)
 *   /normals     [M, 3] float32  — per-face normals
 *
 *   Root attributes:
 *     format        = "h5mesh"  (string)
 *     version       = 1         (int)
 *     num_vertices  = N         (int)
 *     num_triangles = M        (int)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include "molfile_plugin.h"

typedef struct {
  hid_t file_id;
  int num_vertices;
  int num_triangles;
  molfile_graphics_t *graphics;  /* TRINORM/NORMS array, owned by plugin */
} h5mesh_t;


static void *open_h5mesh_read(const char *filepath, const char *filetype,
                               int *natoms) {
  hid_t file_id, attr_id, dset_id, space_id;
  hsize_t dims[2];
  h5mesh_t *h;

  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    fprintf(stderr, "h5meshplugin) Cannot open file %s\n", filepath);
    return NULL;
  }

  /* Verify format attribute */
  attr_id = H5Aopen(file_id, "format", H5P_DEFAULT);
  if (attr_id < 0) {
    fprintf(stderr, "h5meshplugin) No 'format' attribute in %s\n", filepath);
    H5Fclose(file_id);
    return NULL;
  }
  {
    char fmt[32];
    hid_t atype = H5Aget_type(attr_id);
    memset(fmt, 0, sizeof(fmt));
    H5Aread(attr_id, atype, fmt);
    H5Tclose(atype);
    H5Aclose(attr_id);
    if (strcmp(fmt, "h5mesh") != 0) {
      fprintf(stderr, "h5meshplugin) Not an h5mesh file (format='%s')\n", fmt);
      H5Fclose(file_id);
      return NULL;
    }
  }

  /* Read counts — try attributes first, fall back to dataset dims */
  int nv = 0, nt = 0;

  attr_id = H5Aopen(file_id, "num_vertices", H5P_DEFAULT);
  if (attr_id >= 0) {
    H5Aread(attr_id, H5T_NATIVE_INT, &nv);
    H5Aclose(attr_id);
  } else {
    dset_id = H5Dopen2(file_id, "vertices", H5P_DEFAULT);
    if (dset_id < 0) { H5Fclose(file_id); return NULL; }
    space_id = H5Dget_space(dset_id);
    H5Sget_simple_extent_dims(space_id, dims, NULL);
    nv = (int)dims[0];
    H5Sclose(space_id);
    H5Dclose(dset_id);
  }

  attr_id = H5Aopen(file_id, "num_triangles", H5P_DEFAULT);
  if (attr_id >= 0) {
    H5Aread(attr_id, H5T_NATIVE_INT, &nt);
    H5Aclose(attr_id);
  } else {
    dset_id = H5Dopen2(file_id, "triangles", H5P_DEFAULT);
    if (dset_id < 0) { H5Fclose(file_id); return NULL; }
    space_id = H5Dget_space(dset_id);
    H5Sget_simple_extent_dims(space_id, dims, NULL);
    nt = (int)dims[0];
    H5Sclose(space_id);
    H5Dclose(dset_id);
  }

  h = (h5mesh_t *)calloc(1, sizeof(h5mesh_t));
  if (!h) { H5Fclose(file_id); return NULL; }

  h->file_id = file_id;
  h->num_vertices = nv;
  h->num_triangles = nt;
  h->graphics = NULL;

  /* Graphics-only plugin: no atoms */
  *natoms = 0;

  printf("h5meshplugin) Opened %s: %d vertices, %d triangles\n",
         filepath, nv, nt);

  return h;
}


/*
 * Build the graphics array using MOLFILE_TRIANGLE (1 element per triangle).
 * This matches VMD's built-in STL plugin approach and is fast because VMD
 * batches TRIANGLE primitives efficiently. VMD computes flat-shading normals
 * automatically from the vertex winding order.
 */
static int read_h5mesh_rawgraphics(void *v, int *nelem,
                                    const molfile_graphics_t **data) {
  h5mesh_t *h = (h5mesh_t *)v;
  hid_t dset_id;
  int i;
  int nt = h->num_triangles;
  int nv = h->num_vertices;

  /* Read vertices [N, 3] */
  float *verts = (float *)malloc(nv * 3 * sizeof(float));
  if (!verts) return MOLFILE_ERROR;

  dset_id = H5Dopen2(h->file_id, "vertices", H5P_DEFAULT);
  if (dset_id < 0) { free(verts); return MOLFILE_ERROR; }
  H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, verts);
  H5Dclose(dset_id);

  /* Read triangles [M, 3] */
  int32_t *tris = (int32_t *)malloc(nt * 3 * sizeof(int32_t));
  if (!tris) { free(verts); return MOLFILE_ERROR; }

  dset_id = H5Dopen2(h->file_id, "triangles", H5P_DEFAULT);
  if (dset_id < 0) { free(verts); free(tris); return MOLFILE_ERROR; }
  H5Dread(dset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, tris);
  H5Dclose(dset_id);

  /*
   * Allocate graphics array: 1 COLOR element + nt TRIANGLE elements.
   * The leading COLOR element sets the mesh color.
   * Priority: H5MESH_COLOR env var > /color HDF5 attribute > default gray.
   *
   * From VMD Tcl console:
   *   set env(H5MESH_COLOR) "0.2,0.6,0.8"
   *   mol new file.h5mesh
   */
  int total_elems = 1 + nt;
  h->graphics = (molfile_graphics_t *)calloc(total_elems,
                                              sizeof(molfile_graphics_t));
  if (!h->graphics) {
    free(verts); free(tris);
    return MOLFILE_ERROR;
  }

  /* Element 0: COLOR */
  {
    float r = 0.7f, g = 0.7f, b = 0.7f;

    /* Try /color attribute in the file */
    hid_t color_attr = H5Aopen(h->file_id, "color", H5P_DEFAULT);
    if (color_attr >= 0) {
      float rgb[3];
      H5Aread(color_attr, H5T_NATIVE_FLOAT, rgb);
      H5Aclose(color_attr);
      r = rgb[0]; g = rgb[1]; b = rgb[2];
    }

    /* Environment variable overrides file attribute */
    const char *env_color = getenv("H5MESH_COLOR");
    if (env_color) {
      float er, eg, eb;
      if (sscanf(env_color, "%f,%f,%f", &er, &eg, &eb) == 3) {
        r = er; g = eg; b = eb;
      }
    }

    h->graphics[0].type = MOLFILE_COLOR;
    h->graphics[0].data[0] = r;
    h->graphics[0].data[1] = g;
    h->graphics[0].data[2] = b;
  }

  /* Elements 1..nt: TRIANGLE */
  for (i = 0; i < nt; i++) {
    int a = tris[i * 3];
    int b = tris[i * 3 + 1];
    int c = tris[i * 3 + 2];

    if (a < 0 || a >= nv || b < 0 || b >= nv || c < 0 || c >= nv)
      continue;

    molfile_graphics_t *g = &h->graphics[1 + i];
    g->type = MOLFILE_TRIANGLE;
    g->data[0] = verts[a * 3];
    g->data[1] = verts[a * 3 + 1];
    g->data[2] = verts[a * 3 + 2];
    g->data[3] = verts[b * 3];
    g->data[4] = verts[b * 3 + 1];
    g->data[5] = verts[b * 3 + 2];
    g->data[6] = verts[c * 3];
    g->data[7] = verts[c * 3 + 1];
    g->data[8] = verts[c * 3 + 2];
  }

  free(verts);
  free(tris);

  *nelem = total_elems;
  *data = h->graphics;

  printf("h5meshplugin) Created %d graphics elements (%d triangles)\n",
         total_elems, nt);

  return MOLFILE_SUCCESS;
}


static void close_h5mesh_read(void *v) {
  h5mesh_t *h = (h5mesh_t *)v;
  if (h) {
    if (h->graphics) free(h->graphics);
    if (h->file_id >= 0) H5Fclose(h->file_id);
    free(h);
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
  plugin.name = "h5mesh";
  plugin.prettyname = "HDF5 Mesh";
  plugin.author = "Diego Enry Barreto Gomes";
  plugin.majorv = 1;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "h5mesh";

  plugin.open_file_read = open_h5mesh_read;
  plugin.read_rawgraphics = read_h5mesh_rawgraphics;
  plugin.close_file_read = close_h5mesh_read;

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}
