/*
 * fuzz_pdbx.cpp -- libFuzzer target for the pdbx molfile plugin
 *
 * Each input is written to a temporary file and read the way VMD reads it
 * (open, structure, bonds, every frame, graphics, close), with buffers sized
 * exactly to the reported atom count. Anything that loads is then written
 * back out with the plugin's writer and read again; the round trip must keep
 * the atom count, the frames, and every field the writer stores losslessly.
 * Each input is also read a second time with a read buffer of a few bytes,
 * so that every token straddles buffer refills; both readings must agree
 * exactly. A violation aborts, so libFuzzer saves the input as a crash.
 *
 * Build and run:  make fuzz   (needs clang with libFuzzer, e.g. Homebrew llvm)
 */

#define PDBX_TEST_HOOKS
#include "../pdbxplugin.C"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

molfile_plugin_t *g_plugin = NULL;
std::string g_in, g_out;

int register_cb(void *, vmdplugin_t *p) {
  g_plugin = (molfile_plugin_t *)p;
  return VMDPLUGIN_SUCCESS;
}

void fail(const char *what, long i = -1) {
  fprintf(stderr, "fuzz_pdbx: %s", what);
  if (i >= 0) fprintf(stderr, " (atom/frame %ld)", i);
  fprintf(stderr, "\n");
  abort();
}

struct Loaded {
  int natoms;
  int optflags;
  std::vector<molfile_atom_t> atoms;
  std::vector<std::vector<float> > frames;
  std::vector<int> bonds;              /* from, to pairs */
  std::vector<float> orders;
  std::vector<molfile_graphics_t> graphics;
  float cell[6];
};

bool load(const char *path, Loaded &L) {
  memset(L.cell, 0, sizeof(L.cell));
  int natoms = -2;
  void *h = g_plugin->open_file_read(path, "pdbx", &natoms);
  if (!h) return false;
  if (natoms < 0) fail("negative atom count");
  L.natoms = natoms;
  L.optflags = 0;
  L.atoms.assign(natoms > 0 ? (size_t)natoms : 1, molfile_atom_t());
  L.atoms.shrink_to_fit();
  int rc = g_plugin->read_structure(h, &L.optflags, &L.atoms[0]);
  if (rc != MOLFILE_SUCCESS && !(rc == MOLFILE_NOSTRUCTUREDATA && natoms == 0))
    fail("read_structure failed after a successful open");
  if (natoms == 0) L.atoms.clear();

  int nbonds = 0, nbondtypes = 0, *from = NULL, *to = NULL, *btype = NULL;
  float *order = NULL;
  char **btypename = NULL;
  if (g_plugin->read_bonds(h, &nbonds, &from, &to, &order, &btype, &nbondtypes,
                           &btypename) != MOLFILE_SUCCESS)
    fail("read_bonds failed");
  if (nbonds < 0 || (nbonds > 0 && (!from || !to))) fail("bad bond arrays");
  for (int b = 0; b < nbonds; b++) {
    if (from[b] < 1 || from[b] > natoms || to[b] < 1 || to[b] > natoms || from[b] == to[b])
      fail("bond index out of range", b);
    if (order && !(order[b] > 0.0f && order[b] <= 4.0f)) fail("bad bond order", b);
    L.bonds.push_back(from[b]);
    L.bonds.push_back(to[b]);
    if (order) L.orders.push_back(order[b]);
  }

  for (;;) {
    molfile_timestep_t ts;
    memset(&ts, 0, sizeof(ts));
    std::vector<float> c(natoms > 0 ? 3 * (size_t)natoms : 1);
    ts.coords = &c[0];
    if (g_plugin->read_next_timestep(h, natoms, &ts) != MOLFILE_SUCCESS) break;
    if (natoms == 0) fail("frame for a structure without atoms");
    if (L.frames.empty()) {
      L.cell[0] = ts.A; L.cell[1] = ts.B; L.cell[2] = ts.C;
      L.cell[3] = ts.alpha; L.cell[4] = ts.beta; L.cell[5] = ts.gamma;
    }
    for (size_t k = 0; k < c.size(); k++)
      if (!std::isfinite(c[k])) fail("non-finite coordinate");
    L.frames.push_back(c);
    if (L.frames.size() > 1000000) fail("runaway frames");
  }

  int ng = 0;
  const molfile_graphics_t *g = NULL;
  if (g_plugin->read_rawgraphics(h, &ng, &g) != MOLFILE_SUCCESS) fail("rawgraphics failed");
  if (ng < 0 || (ng > 0 && !g)) fail("bad graphics array");
  for (int i = 0; i < ng; i++)
    if (!(g[i].size > 0.0f) || !std::isfinite(g[i].data[0])) fail("bad sphere", i);
  if (ng > 0) L.graphics.assign(g, g + ng);

  g_plugin->close_file_read(h);
  if (natoms > 0 && L.frames.empty()) fail("atoms without coordinates");
  return true;
}

/* a string field the writer stores without change */
bool clean(const char *s, size_t cap) {
  size_t n = strnlen(s, cap);
  if (n == 0 || s[0] == ' ' || s[n - 1] == ' ') return false;
  for (size_t i = 0; i < n; i++) {
    unsigned char c = (unsigned char)s[i];
    if (c < 32 || c >= 127) return false;   /* not representable in CIF 1.1 */
  }
  return true;
}

bool close_enough(double a, double b, double abs_tol) {
  return std::fabs(a - b) <= abs_tol + 1e-6 * std::fabs(a);
}

void roundtrip(const Loaded &A) {
  void *w = g_plugin->open_file_write(g_out.c_str(), "pdbx", A.natoms);
  if (!w) fail("open_file_write failed");
  if (g_plugin->write_structure(w, A.optflags, &A.atoms[0]) != MOLFILE_SUCCESS)
    fail("write_structure failed");
  for (size_t f = 0; f < A.frames.size(); f++) {
    molfile_timestep_t ts;
    memset(&ts, 0, sizeof(ts));
    ts.coords = const_cast<float *>(&A.frames[f][0]);
    if (g_plugin->write_timestep(w, &ts) != MOLFILE_SUCCESS) fail("write_timestep failed");
  }
  g_plugin->close_file_write(w);

  Loaded B;
  if (!load(g_out.c_str(), B)) fail("written file does not load");
  if (B.natoms != A.natoms) fail("atom count changed on round trip");
  if (B.frames.size() != A.frames.size()) fail("frame count changed on round trip");
  for (int i = 0; i < A.natoms; i++) {
    const molfile_atom_t &a = A.atoms[(size_t)i], &b = B.atoms[(size_t)i];
    if (a.resid != b.resid) fail("resid changed", i);
    if (a.atomicnumber != b.atomicnumber) fail("atomic number changed", i);
    if (clean(a.name, sizeof(a.name)) && strcmp(a.name, b.name) != 0) fail("name changed", i);
    if (clean(a.resname, sizeof(a.resname)) && strcmp(a.resname, b.resname) != 0)
      fail("resname changed", i);
    if ((A.optflags & MOLFILE_ALTLOC) && a.altloc[0] > ' ' && a.altloc[0] < 127 &&
        a.altloc[0] != b.altloc[0])
      fail("altloc changed", i);
    if ((A.optflags & MOLFILE_INSERTION) && a.insertion[0] > ' ' && a.insertion[0] < 127 &&
        a.insertion[0] != b.insertion[0])
      fail("insertion code changed", i);
    if ((A.optflags & MOLFILE_OCCUPANCY) && !close_enough(a.occupancy, b.occupancy, 0.0051))
      fail("occupancy changed", i);
    if ((A.optflags & MOLFILE_BFACTOR) && !close_enough(a.bfactor, b.bfactor, 0.0051))
      fail("B factor changed", i);
    if ((A.optflags & MOLFILE_CHARGE) && a.charge == std::floor(a.charge) &&
        std::fabs(a.charge) < 100.0f && a.charge != b.charge)
      fail("formal charge changed", i);
  }
  for (size_t f = 0; f < A.frames.size(); f++)
    for (size_t k = 0; k < A.frames[f].size(); k++)
      if (!close_enough(A.frames[f][k], B.frames[f][k], 0.00051))
        fail("coordinate changed", (long)f);
}

/* the two readings of one input must be identical */
void same(const Loaded &A, const Loaded &B) {
  if (A.natoms != B.natoms || A.optflags != B.optflags) fail("buffer size changed atoms/flags");
  if (A.natoms > 0 && memcmp(&A.atoms[0], &B.atoms[0], A.atoms.size() * sizeof(A.atoms[0])))
    fail("buffer size changed atom fields");
  if (A.frames != B.frames) fail("buffer size changed coordinates");
  if (A.bonds != B.bonds || A.orders != B.orders) fail("buffer size changed bonds");
  if (memcmp(A.cell, B.cell, sizeof(A.cell))) fail("buffer size changed the cell");
  if (A.graphics.size() != B.graphics.size() ||
      (!A.graphics.empty() && memcmp(&A.graphics[0], &B.graphics[0],
                                     A.graphics.size() * sizeof(A.graphics[0]))))
    fail("buffer size changed graphics");
}

}  // namespace

extern "C" int LLVMFuzzerInitialize(int *, char ***) {
  VMDPLUGIN_init();
  VMDPLUGIN_register(NULL, register_cb);
  const char *tmp = getenv("TMPDIR");
  std::string dir = (tmp && *tmp) ? tmp : "/tmp";
  if (dir[dir.size() - 1] != '/') dir += '/';
  char pid[32];
  snprintf(pid, sizeof(pid), "%ld", (long)getpid());
  g_in = dir + "pdbx_fuzz_in_" + pid + ".cif";
  g_out = dir + "pdbx_fuzz_out_" + pid + ".cif";
  if (!freopen("/dev/null", "w", stdout)) return 0;   /* plugin messages */
  return 0;
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  FILE *fp = fopen(g_in.c_str(), "wb");
  if (!fp) return 0;
  if (size) fwrite(data, 1, size, fp);
  fclose(fp);
  Loaded A, T;
  g_read_buffer_bytes = (size_t)1 << 16;
  bool ok = load(g_in.c_str(), A);
  g_read_buffer_bytes = 1 + (size ? data[size / 2] % 7 : 0);   /* 1..7 bytes */
  bool ok_tiny = load(g_in.c_str(), T);
  g_read_buffer_bytes = (size_t)1 << 16;
  if (ok != ok_tiny) fail("buffer size changed whether the file loads");
  if (ok) same(A, T);
  if (ok && A.natoms > 0) roundtrip(A);
  return 0;
}
