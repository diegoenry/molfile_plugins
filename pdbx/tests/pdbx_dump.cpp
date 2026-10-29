/*
 * pdbx_dump.cpp -- standalone driver for the pdbx molfile plugin
 *
 * Compiles the plugin source directly and runs it the way VMD does:
 * open -> read_structure -> read_bonds -> read_next_timestep... -> close.
 * Buffers are sized exactly to the reported atom count so that a sanitizer
 * build (make asan) reports any out-of-bounds write. Bond indices are checked
 * to be within [1, natoms].
 *
 * Usage:
 *   pdbx_dump [-q] [-c] [-a N] [-w out.cif] file.cif
 *     -q         summary line only
 *     -c         machine-comparable dump of every atom and frame (no type field,
 *                which the writer normalizes to the element symbol)
 *     -a N       print the first N atoms (default 10, -1 = all)
 *     -w FILE    write the structure and all frames back out with the plugin
 *     -b N       start with an N-byte read buffer (tests refills at every token)
 *
 * Exit status: 0 ok, 1 open/read failure, 2 consistency error.
 */

#define PDBX_TEST_HOOKS
#include "../pdbxplugin.C"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static molfile_plugin_t *g_plugin = NULL;

static int register_cb(void *, vmdplugin_t *p) {
  g_plugin = (molfile_plugin_t *)p;
  return VMDPLUGIN_SUCCESS;
}

int main(int argc, char **argv) {
  bool quiet = false, compare = false;
  long show = 10;
  const char *out = NULL;
  const char *in = NULL;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-q")) quiet = true;
    else if (!strcmp(argv[i], "-c")) compare = quiet = true;
    else if (!strcmp(argv[i], "-a") && i + 1 < argc) show = atol(argv[++i]);
    else if (!strcmp(argv[i], "-w") && i + 1 < argc) out = argv[++i];
    else if (!strcmp(argv[i], "-b") && i + 1 < argc) g_read_buffer_bytes = (size_t)atol(argv[++i]);
    else in = argv[i];
  }
  if (!in) {
    fprintf(stderr, "usage: %s [-q] [-a N] [-w out.cif] file.cif\n", argv[0]);
    return 1;
  }

  VMDPLUGIN_init();
  VMDPLUGIN_register(NULL, register_cb);
  if (!g_plugin) return 1;

  int natoms = -2;
  void *h = g_plugin->open_file_read(in, "pdbx", &natoms);
  if (!h) { printf("RESULT open-failed\n"); return 1; }
  if (natoms < 0) { printf("RESULT bad-natoms %d\n", natoms); return 2; }

  /* exact-size allocations: any overrun is a sanitizer error */
  molfile_atom_t *atoms = (molfile_atom_t *)calloc(natoms > 0 ? (size_t)natoms : 1,
                                                   sizeof(molfile_atom_t));
  int optflags = 0;
  int rc = g_plugin->read_structure(h, &optflags, atoms);
  if (rc != MOLFILE_SUCCESS && !(rc == MOLFILE_NOSTRUCTUREDATA && natoms == 0)) {
    printf("RESULT structure-failed %d\n", rc);
    g_plugin->close_file_read(h);
    free(atoms);
    return 1;
  }

  int nbonds = 0, nbondtypes = 0, *from = NULL, *to = NULL, *btype = NULL;
  float *border = NULL;
  char **btypename = NULL;
  int status = 0;
  if (g_plugin->read_bonds(h, &nbonds, &from, &to, &border, &btype, &nbondtypes,
                           &btypename) != MOLFILE_SUCCESS) {
    printf("RESULT bonds-failed\n");
    status = 2;
  }
  for (int b = 0; b < nbonds; b++) {
    if (from[b] < 1 || from[b] > natoms || to[b] < 1 || to[b] > natoms ||
        from[b] == to[b]) {
      printf("RESULT bad-bond %d: %d-%d (natoms %d)\n", b, from[b], to[b], natoms);
      status = 2;
      break;
    }
  }

  std::vector<std::vector<float> > frames;
  for (;;) {
    molfile_timestep_t ts;
    memset(&ts, 0, sizeof(ts));
    float *coords = (float *)malloc(natoms > 0 ? 3 * (size_t)natoms * sizeof(float) : 4);
    ts.coords = coords;
    int r = g_plugin->read_next_timestep(h, natoms, &ts);
    if (r != MOLFILE_SUCCESS) { free(coords); break; }
    frames.push_back(std::vector<float>(coords, coords + 3 * (size_t)natoms));
    if (frames.size() == 1 && !quiet && ts.A > 0)
      printf("cell %.3f %.3f %.3f  %.2f %.2f %.2f\n", ts.A, ts.B, ts.C,
             ts.alpha, ts.beta, ts.gamma);
    free(coords);
  }

  int ngraphics = 0;
  const molfile_graphics_t *g = NULL;
  if (g_plugin->read_rawgraphics) g_plugin->read_rawgraphics(h, &ngraphics, &g);

  double sum = 0;
  for (size_t f = 0; f < frames.size(); f++)
    for (size_t k = 0; k < frames[f].size(); k++) sum += frames[f][k];

  int nZ0 = 0;
  for (int i = 0; i < natoms; i++) if (atoms[i].atomicnumber == 0) nZ0++;

  printf("RESULT ok natoms=%d frames=%lu bonds=%d orders=%s graphics=%d flags=0x%x "
         "elementX=%d coordsum=%.3f\n", natoms, (unsigned long)frames.size(), nbonds,
         border ? "yes" : "no", ngraphics, optflags, nZ0, sum);

  if (compare) {
    /* the writer always emits occupancy and B columns, and bonds are not
     * written, so those flags are not compared */
    printf("FLAGS 0x%x\n",
           optflags & ~(MOLFILE_BONDSSPECIAL | MOLFILE_OCCUPANCY | MOLFILE_BFACTOR));
    for (int i = 0; i < natoms; i++) {
      const molfile_atom_t &a = atoms[i];
      printf("A|%s|%s|%d|%s|%s|%s|%s|%.2f|%.2f|%.0f|%d\n", a.name, a.resname, a.resid,
             a.chain, a.segid, a.altloc, a.insertion, a.occupancy, a.bfactor, a.charge,
             a.atomicnumber);
    }
    for (size_t f = 0; f < frames.size(); f++)
      for (int i = 0; i < natoms; i++)
        printf("X|%lu|%d|%.3f|%.3f|%.3f\n", (unsigned long)f, i,
               frames[f][3 * (size_t)i] + 0.0f,          /* -0.000 -> 0.000 */
               frames[f][3 * (size_t)i + 1] + 0.0f, frames[f][3 * (size_t)i + 2] + 0.0f);
  }

  if (!quiet) {
    printf("flags:%s%s%s%s%s%s%s%s%s\n",
           (optflags & MOLFILE_INSERTION) ? " insertion" : "",
           (optflags & MOLFILE_OCCUPANCY) ? " occupancy" : "",
           (optflags & MOLFILE_BFACTOR) ? " bfactor" : "",
           (optflags & MOLFILE_MASS) ? " mass" : "",
           (optflags & MOLFILE_CHARGE) ? " charge" : "",
           (optflags & MOLFILE_RADIUS) ? " radius" : "",
           (optflags & MOLFILE_ALTLOC) ? " altloc" : "",
           (optflags & MOLFILE_ATOMICNUMBER) ? " atomicnumber" : "",
           (optflags & MOLFILE_BONDSSPECIAL) ? " bondsspecial" : "");
    long n = (show < 0 || show > natoms) ? natoms : show;
    if (n > 0)
      printf("%6s %-6s %-6s %-5s %6s %-2s %-7s %-2s %-2s %5s %6s %5s %3s %9s %9s %9s\n",
             "index", "name", "type", "resn", "resid", "ch", "segid", "al", "in",
             "occ", "beta", "chg", "Z", "x", "y", "z");
    for (long i = 0; i < n; i++) {
      const molfile_atom_t &a = atoms[i];
      const float *c = frames.empty() ? NULL : &frames[0][3 * (size_t)i];
      printf("%6ld %-6s %-6s %-5s %6d %-2s %-7s %-2s %-2s %5.2f %6.2f %5.1f %3d "
             "%9.3f %9.3f %9.3f\n", i, a.name, a.type, a.resname, a.resid, a.chain,
             a.segid, a.altloc, a.insertion, a.occupancy, a.bfactor, a.charge,
             a.atomicnumber, c ? c[0] : 0.0f, c ? c[1] : 0.0f, c ? c[2] : 0.0f);
    }
    int nb = nbonds < 10 ? nbonds : 10;
    for (int b = 0; b < nb; b++)
      printf("bond %d-%d %s%s/%d -- %s%s/%d order %.1f\n", from[b], to[b],
             atoms[from[b] - 1].name, atoms[from[b] - 1].altloc, atoms[from[b] - 1].resid,
             atoms[to[b] - 1].name, atoms[to[b] - 1].altloc, atoms[to[b] - 1].resid,
             border ? border[b] : 1.0f);
  }

  g_plugin->close_file_read(h);

  if (out && natoms > 0) {
    void *w = g_plugin->open_file_write(out, "pdbx", natoms);
    if (!w) { printf("RESULT write-open-failed\n"); status = 2; }
    else {
      if (g_plugin->write_structure(w, optflags, atoms) != MOLFILE_SUCCESS) status = 2;
      for (size_t f = 0; f < frames.size(); f++) {
        molfile_timestep_t ts;
        memset(&ts, 0, sizeof(ts));
        ts.coords = &frames[f][0];
        if (g_plugin->write_timestep(w, &ts) != MOLFILE_SUCCESS) status = 2;
      }
      g_plugin->close_file_write(w);
    }
  }

  free(atoms);
  return status;
}
