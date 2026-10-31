# Release builds

These scripts build the five plugins for every VMD platform on one Mac, test them, and package one archive per platform:

| Archive | VMD platform | Runs on |
|---|---|---|
| `molfile_plugins-TAG-MACOSXARM64.tar.gz` | `MACOSXARM64` | macOS 11 or later, Apple silicon |
| `molfile_plugins-TAG-MACOSXX86_64.tar.gz` | `MACOSXX86_64` | macOS 10.13 or later, Intel |
| `molfile_plugins-TAG-LINUXAMD64.tar.gz` | `LINUXAMD64` | Linux x86_64, glibc 2.17 or later |
| `molfile_plugins-TAG-LINUXARM64.tar.gz` | `LINUXARM64` | Linux ARM64, glibc 2.17 or later |
| `molfile_plugins-TAG-WIN64.zip` | `WIN64` | 64-bit VMD on Windows 10 or later |

Each archive contains:
- `molfile/` with the five plugins;
- a `README.txt` that tells users how to load them;
- the plugin license;
- the licenses of the libraries built into the plugins.

HDF5 and zlib are built from source and linked into the plugins, so users install nothing. Each plugin exports only VMD's three entry points, so the copies of HDF5 inside different plugins don't interfere.

## Making a release

```sh
brew install cmake ninja zig mingw-w64    # once; also needs Xcode's command line tools
release/make_release.sh v1.0              # TAG names the archives; default: the commit hash
```

The script works through these steps:
1. Finds VMD's plugin headers in `/Applications`. If it can't, set `VMD_PLUGIN_INCLUDE` to the folder with `molfile_plugin.h`.
2. Downloads HDF5 and zlib, and checks their SHA-256. The versions and checksums are in `config.sh`.
3. Exports the committed state of `HEAD`. Uncommitted changes are not included.
4. Builds all five platforms at once:
   - **macOS:** Apple clang;
   - **Linux:** zig, targeting glibc 2.17;
   - **Windows:** MinGW-w64.

   The first run also builds HDF5 for each platform, which takes a few minutes. Later runs reuse those builds (`REBUILD_DEPS=1` forces a rebuild).
5. Checks every plugin:
   - it exports only `vmdplugin_init`, `vmdplugin_register` and `vmdplugin_fini`;
   - it needs only system libraries;
   - on Linux, it needs no glibc newer than 2.17.
6. Runs the tests this Mac can run (see below).
7. Writes the archives and `molfile_plugins-TAG-SHA256SUMS.txt` to `dist/`.

Other settings:
- `ARCHES="MACOSXARM64 WIN64"` builds only some platforms.
- `W` sets the work folder (default `~/.cache/molfile-release`, about 250 MB).
- `DIST` sets the output folder.

## Testing

`loadtest.c` loads a plugin the way VMD does. It then prints a fingerprint of what the plugin read: atom count, a hash of all atom fields, the bonds, every frame's coordinates and box, and the graphics. `run_tests.sh` runs it on:
- every pdbx test file, plain, BinaryCIF, `.gz` and the deliberately damaged ones;
- `example.sdf` and `equil.h5md`;
- wrong and missing files for each HDF5 plugin.

All platforms must print the same lines. The repository has no CNDB or mesh file, so pass a folder with some (`*.cndb`, `*.h5md`, `*.h5mesh`):

```sh
EXTRA_TESTFILES=~/my-test-files release/make_release.sh v1.0
```

What gets tested where:
- **macOS:** `make_release.sh` runs the tests on the Mac builds; on Apple silicon, the Intel build runs under Rosetta.
- **Linux:** run `release/test_linux.sh`. It runs the same tests in a container for each Linux build, with no network and the work folder read-only, then compares the results with the Mac's. It needs Docker, Docker Desktop or colima; the top of the script has the settings, including how to turn on x86_64 emulation on Apple silicon.
- **In VMD:** `vmd_check.tcl` loads the plugins into a real VMD and opens a file with each:

  ```sh
  V="/path/to/VMD.app/Contents/vmd"   # the folder with vmd_MACOSXARM64
  VMDDIR="$V" "$V/vmd_MACOSXARM64" -dispdev text -e release/vmd_check.tcl \
      -args ~/.cache/molfile-release/out/MACOSXARM64/molfile ~/.cache/molfile-release
  ```

- **Windows:** the builds can't run on a Mac. `make_release.sh` checks their exports and imported DLLs. For a quick check on Windows, copy `~/.cache/molfile-release/loadtest/loadtest-WIN64.exe` next to the extracted `molfile` folder, and run `loadtest-WIN64.exe molfile\pdbxplugin.so`. It must list the `pdbx` and `mmcif` formats; `loadtest-WIN64.exe molfile\pdbxplugin.so pdbx file.cif` reads a file.

## macOS signing

VMD 1.9.4 from the VMD website is notarized with Apple's hardened runtime, so macOS lets it load only plugins signed with the VMD developers' certificate. It refuses these plugins, and any other locally built plugin, with "different Team IDs". The archives' README tells users how to re-sign their copy of VMD. Signing the plugins with the VMD team's Developer ID would remove that step. VMD builds without the hardened runtime load the plugins as they are.

## Files

| File | Purpose |
|---|---|
| `make_release.sh` | the whole release: sources, builds, checks, tests, archives |
| `build.sh` | builds one platform (`build.sh LINUXAMD64`) |
| `package.sh` | packages one platform |
| `run_tests.sh`, `loadtest.c` | the cross-platform tests |
| `test_linux.sh` | runs the tests for the Linux builds in containers |
| `vmd_check.tcl` | loads the plugins into a real VMD |
| `config.sh` | library versions, checksums, minimum OS versions |
| `exports_*` | the symbols each platform's plugins may export |
