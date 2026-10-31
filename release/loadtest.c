/* loadtest: load a molfile plugin the way VMD does and print a fingerprint
 * of what it reads, so builds for different platforms can be compared.
 *
 *   loadtest plugin.so                     list the formats it registers
 *   loadtest plugin.so TYPE FILE [FRAMES]  read FILE (at most FRAMES frames)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include "molfile_plugin.h"

typedef int (*init_fn)(void);
typedef int (*register_fn)(void *, vmdplugin_register_cb);

static molfile_plugin_t *plugins[16];
static int nplugins = 0;

static int collect(void *v, vmdplugin_t *p) {
  (void)v;
  if (!strcmp(p->type, MOLFILE_PLUGIN_TYPE) && nplugins < 16)
    plugins[nplugins++] = (molfile_plugin_t *)p;
  return VMDPLUGIN_SUCCESS;
}

static unsigned long long fnv(unsigned long long h, const void *data, size_t n) {
  const unsigned char *b = (const unsigned char *)data;
  for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 1099511628211ULL; }
  return h;
}
static unsigned long long fnv_str(unsigned long long h, const char *s) {
  return fnv(h, s, strlen(s) + 1);
}
/* floats are hashed as text, so -0.0 and rounding noise in the last bit
   do not count as differences */
static unsigned long long fnv_num(unsigned long long h, double x) {
  char buf[64];
  snprintf(buf, sizeof buf, "%.4f", x);
  if (!strcmp(buf, "-0.0000")) strcpy(buf, "0.0000");
  return fnv_str(h, buf);
}

int main(int argc, char **argv) {
  if (argc != 2 && argc != 4 && argc != 5) {
    fprintf(stderr, "usage: loadtest plugin.so [TYPE FILE [FRAMES]]\n");
    return 2;
  }
#ifdef _WIN32
  HMODULE lib = LoadLibraryA(argv[1]);
  if (!lib) { printf("LOAD FAILED (error %lu)\n", GetLastError()); return 1; }
  init_fn init = (init_fn)GetProcAddress(lib, "vmdplugin_init");
  register_fn reg = (register_fn)GetProcAddress(lib, "vmdplugin_register");
  init_fn fini = (init_fn)GetProcAddress(lib, "vmdplugin_fini");
#else
  void *lib = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if (!lib) { printf("LOAD FAILED: %s\n", dlerror()); return 1; }
  init_fn init = (init_fn)dlsym(lib, "vmdplugin_init");
  register_fn reg = (register_fn)dlsym(lib, "vmdplugin_register");
  init_fn fini = (init_fn)dlsym(lib, "vmdplugin_fini");
#endif
  if (!init || !reg || !fini) { printf("MISSING ENTRY POINTS\n"); return 1; }
  if (init() != VMDPLUGIN_SUCCESS) { printf("INIT FAILED\n"); return 1; }
  reg(NULL, collect);
  if (argc == 2) {
    for (int i = 0; i < nplugins; i++)
      printf("plugin %s %d.%d abi=%d ext=%s\n", plugins[i]->name, plugins[i]->majorv,
             plugins[i]->minorv, plugins[i]->abiversion,
             plugins[i]->filename_extension ? plugins[i]->filename_extension : "");
    fini();
    return nplugins > 0 ? 0 : 1;
  }

  const char *type = argv[2], *file = argv[3];
  long maxframes = argc == 5 ? atol(argv[4]) : 1000000L;
  molfile_plugin_t *p = NULL;
  for (int i = 0; i < nplugins; i++)
    if (!strcmp(plugins[i]->name, type)) p = plugins[i];
  if (!p) { printf("NO PLUGIN NAMED %s\n", type); return 1; }

  int natoms = MOLFILE_NUMATOMS_UNKNOWN;
  void *h = p->open_file_read(file, type, &natoms);
  if (!h) { printf("OPEN FAILED\n"); return 1; }
  printf("%s %d.%d natoms=%d", p->name, p->majorv, p->minorv, natoms);

  if (p->read_structure && natoms > 0) {
    molfile_atom_t *a = (molfile_atom_t *)calloc((size_t)natoms, sizeof *a);
    int flags = 0;
    if (p->read_structure(h, &flags, a) != MOLFILE_SUCCESS) { printf(" STRUCTURE FAILED\n"); return 1; }
    unsigned long long hs = 1469598103934665603ULL;
    for (int i = 0; i < natoms; i++) {
      hs = fnv_str(hs, a[i].name); hs = fnv_str(hs, a[i].type); hs = fnv_str(hs, a[i].resname);
      hs = fnv_str(hs, a[i].segid); hs = fnv_str(hs, a[i].chain); hs = fnv_str(hs, a[i].altloc);
      hs = fnv_str(hs, a[i].insertion); hs = fnv_num(hs, a[i].resid);
      if (flags & MOLFILE_OCCUPANCY) hs = fnv_num(hs, a[i].occupancy);
      if (flags & MOLFILE_BFACTOR) hs = fnv_num(hs, a[i].bfactor);
      if (flags & MOLFILE_MASS) hs = fnv_num(hs, a[i].mass);
      if (flags & MOLFILE_CHARGE) hs = fnv_num(hs, a[i].charge);
      if (flags & MOLFILE_RADIUS) hs = fnv_num(hs, a[i].radius);
      if (flags & MOLFILE_ATOMICNUMBER) hs = fnv_num(hs, a[i].atomicnumber);
    }
    printf(" flags=0x%x atoms=%016llx", flags, hs);
    free(a);
  }

  if (p->read_bonds && natoms > 0) {
    int nb = 0, nbt = 0, *from = NULL, *to = NULL, *bt = NULL;
    float *bo = NULL;
    char **btn = NULL;
    if (p->read_bonds(h, &nb, &from, &to, &bo, &bt, &nbt, &btn) != MOLFILE_SUCCESS) {
      printf(" BONDS FAILED\n");
      return 1;
    }
    unsigned long long hb = 1469598103934665603ULL;
    int bad = 0;
    for (int i = 0; i < nb; i++) {
      if (from[i] < 1 || from[i] > natoms || to[i] < 1 || to[i] > natoms) bad++;
      hb = fnv_num(hb, from[i]); hb = fnv_num(hb, to[i]);
      if (bo) hb = fnv_num(hb, bo[i]);
    }
    printf(" bonds=%d%s bad=%d hash=%016llx", nb, bo ? "+orders" : "", bad, hb);
  }

  if (p->read_rawgraphics) {
    int n = 0;
    const molfile_graphics_t *g = NULL;
    if (p->read_rawgraphics(h, &n, &g) != MOLFILE_SUCCESS) { printf(" GRAPHICS FAILED\n"); return 1; }
    unsigned long long hg = 1469598103934665603ULL;
    for (int i = 0; i < n; i++) {
      hg = fnv_num(hg, g[i].type); hg = fnv_num(hg, g[i].style); hg = fnv_num(hg, g[i].size);
      for (int k = 0; k < 9; k++) hg = fnv_num(hg, g[i].data[k]);
    }
    printf(" graphics=%d hash=%016llx", n, hg);
  }
  printf("\n");

  if (p->read_next_timestep && natoms > 0) {
    molfile_timestep_t ts;
    memset(&ts, 0, sizeof ts);
    ts.coords = (float *)malloc(3 * sizeof(float) * (size_t)natoms);
    long nf = 0;
    unsigned long long hc = 1469598103934665603ULL;
    while (nf < maxframes && p->read_next_timestep(h, natoms, &ts) == MOLFILE_SUCCESS) {
      double sum = 0;
      for (long i = 0; i < 3L * natoms; i++) { sum += ts.coords[i]; hc = fnv_num(hc, ts.coords[i]); }
      if (nf == 0)
        printf("frame0 coordsum=%.3f box=%.3f,%.3f,%.3f,%.2f,%.2f,%.2f time=%.3f\n", sum, ts.A, ts.B,
               ts.C, ts.alpha, ts.beta, ts.gamma, ts.physical_time);
      nf++;
    }
    printf("frames=%ld coords=%016llx\n", nf, hc);
    free(ts.coords);
  }
  p->close_file_read(h);
  fini();
  return 0;
}
