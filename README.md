# VMD molfile plugins

Temporary repo for new VMD molfile plugins.
Once stable they will be part of VMD 2.x official distribution.

| Folder | Reads (and writes) | VMD type | File extensions | Needs |
|---|---|---|---|---|
| [`pdbx`](pdbx/README.md) | mmCIF/PDBx and BinaryCIF, plain or `.gz`; also writes mmCIF | `pdbx` (also `mmcif`) | `.cif` `.mmcif` `.bcif` | C++11; zlib (optional, for `.gz`) |
| `sdf_reader_writer` | SDF/MOL V2000; also writes | `sdf` | `.sdf` `.mol` | C |
| [`CNDB_hdf5_reader`](CNDB_hdf5_reader/README.md) | OpenMiChroM CNDB | `cndb` | `.cndb` | C, HDF5 |
| [`gromacs_h5md`](gromacs_h5md/README.md) | GROMACS H5MD trajectories | `h5md_gromacs` | `.h5` `.h5md` | C, HDF5 |
| `MeshHDF5` | HDF5 meshes | `h5mesh` | `.h5mesh` | C, HDF5 |

The repository also has two tools that are not plugins:
- `MeshOptimize` builds `mesh_optimize`, a command-line STL smoother and decimator. It needs C++17; HDF5 and OpenMP are optional.
- [`gromacs_h5md_webBrowser`](gromacs_h5md_webBrowser/README.md) is a Python web viewer for H5MD files, with nothing to compile.

`mmcif/` only points to `pdbx`, which replaced the old mmCIF reader.

## Building

### What every platform needs

- **VMD's plugin headers:** `molfile_plugin.h` and `vmdplugin.h`, in `plugins/include` inside the VMD folder. They are the same on every platform. If an installation lacks them, take them from VMD's plugin source or from another VMD installation.
- **A plugin that matches VMD:** to find VMD's architecture and folder, type this in VMD's console:

  ```tcl
  vmdinfo arch          ;# MACOSXARM64, MACOSXX86_64, LINUXAMD64, LINUXARM64, WIN64, ...
  puts $env(VMDDIR)     ;# the VMD folder
  ```

- **HDF5:** only for `CNDB_hdf5_reader`, `gromacs_h5md` and `MeshHDF5`.

Each plugin builds one file, for example `pdbxplugin.so`. VMD uses that name on every platform, Windows included.

Release archives for all platforms, with HDF5 built in, are made on a Mac by `release/make_release.sh`; see [release/README.md](release/README.md).

### macOS

1. Install the compiler: `xcode-select --install`.
2. For the HDF5 plugins, install HDF5: `brew install hdf5` ([Homebrew](https://brew.sh)). zlib comes with macOS.
3. Build each plugin in its folder:

   ```sh
   cd pdbx              && make && cd ..
   cd sdf_reader_writer && make && cd ..
   cd CNDB_hdf5_reader  && make && cd ..
   cd MeshHDF5          && make && cd ..
   ```

4. The `gromacs_h5md` Makefile cannot take a header path with spaces, and VMD's app name has them. Give it a link without spaces:

   ```sh
   ln -s "/Applications/VMD 1.9.4a57-arm64-Rev12.app/Contents/vmd/plugins/include" ~/vmd-plugin-include
   cd gromacs_h5md && make -f Makefile.h5md_gromacs VMD_PLUGIN_INCLUDE=$HOME/vmd-plugin-include
   ```

The Makefiles look for VMD in `/Applications`, `~/Applications` and `~/software`. If a build stops with `'molfile_plugin.h' file not found`, give the header folder yourself, quoted:

```sh
make VMD_PLUGIN_INCLUDE="/Applications/VMD 1.9.4a57-arm64-Rev12.app/Contents/vmd/plugins/include"
```

Build for the architecture of your VMD. On an Apple silicon Mac the default is arm64, for `MACOSXARM64`. The Intel build of VMD (`MACOSXX86_64`) needs x86_64 plugins, for example `make CC="cc -arch x86_64"` (pdbx: `make CXX="c++ -arch x86_64"`), and x86_64 HDF5 libraries for the HDF5 plugins. `file pdbxplugin.so` shows what you built.

`MeshOptimize`: `make` builds it single-threaded with Apple's compiler. For OpenMP, `brew install gcc`; the Makefile then uses Homebrew's `g++` by itself.

### Linux

1. Install a compiler, zlib and HDF5:

   ```sh
   sudo apt install build-essential pkg-config zlib1g-dev libhdf5-dev        # Debian, Ubuntu
   sudo dnf install gcc gcc-c++ make pkgconf-pkg-config zlib-devel hdf5-devel # Fedora, RHEL
   ```

2. Build each plugin in its folder, with the same commands as in step 3 for macOS. `gromacs_h5md` needs no link on Linux: `make -f Makefile.h5md_gromacs`.

The Makefiles look for VMD in `/usr/local/lib/vmd` and `~/vmd`, and find HDF5 with `pkg-config`. If they don't find them, give the folders yourself:

```sh
make VMD_PLUGIN_INCLUDE=/path/to/vmd/plugins/include      # VMD's headers (every plugin)
make HDF5_PREFIX=/path/to/hdf5                            # HDF5 with include/ and lib/
```

### Windows

The Makefiles do not run on Windows, so call the compiler directly. Three things matter on Windows:
- **Export the entry points:** define `VMDPLUGIN_EXPORTS`, so that `vmdplugin.h` exports VMD's entry points from the DLL.
- **Name the DLL `.so`:** VMD loads `*.so` files from its plugin folder. `plugin dlopen` takes any name.
- **Build 64-bit:** for the 64-bit VMD (`WIN64`).

Set the header folder first (in `cmd`; the path is an example):

```bat
set VMDINC=C:\Program Files\VMD\plugins\include
```

**With Visual Studio** (in the "x64 Native Tools Command Prompt"; VMD itself is built with MSVC):

```bat
cl /nologo /O2 /MT /EHsc /TP /DVMDPLUGIN_EXPORTS /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_DEPRECATE ^
   /I"%VMDINC%" /LD pdbx\pdbxplugin.C /Fe:pdbxplugin.so
cl /nologo /O2 /MT /DVMDPLUGIN_EXPORTS /D_CRT_SECURE_NO_WARNINGS ^
   /I"%VMDINC%" /LD sdf_reader_writer\sdfplugin.c /Fe:sdfplugin.so
```

- **`/TP` (pdbx):** `cl` would compile `pdbxplugin.C` as C otherwise.
- **`/MT`:** builds the C runtime into the plugin, so no Visual C++ runtime DLL is needed.
- **`.cif.gz` support:** install zlib with [vcpkg](https://vcpkg.io) (`vcpkg install zlib:x64-windows-static`). Then add `/D_USE_ZLIB /I"%VCPKG_ROOT%\installed\x64-windows-static\include"` and, at the end, `"%VCPKG_ROOT%\installed\x64-windows-static\lib\zlib.lib"`.
- **HDF5 plugins:** install HDF5 (from the HDF Group, or `vcpkg install hdf5:x64-windows`). Add `/I<hdf5>\include` and `<hdf5>\lib\hdf5.lib`. With an HDF5 DLL, older HDF5 versions also need `/DH5_BUILT_AS_DYNAMIC_LIB`, and the HDF5 DLLs must be next to `vmd.exe` or on `PATH`.

**With MSYS2 / MinGW-w64** (in an MSYS2 "UCRT64" shell; set `VMDINC` with forward slashes):

```sh
pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-zlib mingw-w64-ucrt-x86_64-hdf5
g++ -std=c++11 -O2 -DVMDPLUGIN_EXPORTS -D_USE_ZLIB -I"$VMDINC" -shared -static \
    -o pdbxplugin.so pdbx/pdbxplugin.C -lz
gcc -O2 -DVMDPLUGIN_EXPORTS -I"$VMDINC" -shared -static -o sdfplugin.so sdf_reader_writer/sdfplugin.c
gcc -O2 -DVMDPLUGIN_EXPORTS -I"$VMDINC" -shared -static-libgcc \
    -o cndbplugin.so CNDB_hdf5_reader/cndbplugin.c -lhdf5        # likewise gromacs_h5md, MeshHDF5
g++ -std=c++17 -O2 -static -o mesh_optimize.exe MeshOptimize/mesh_optimize.cpp
```

- **`-static`:** puts the GCC runtime and zlib inside the plugin, so it needs no MinGW DLLs.
- **HDF5 plugins:** they use MSYS2's HDF5 DLLs, so add `C:\msys64\ucrt64\bin` to `PATH`, or copy the DLLs it lists next to `vmd.exe`. To list them: `ldd cndbplugin.so`.

### Tested

- **macOS** (Apple silicon, VMD 1.9.4 headers):
  - every plugin builds with the commands above;
  - every plugin loads through `dlopen` as VMD loads it, and registers its format;
  - `pdbx` was also run inside VMD 2.0.
- **Linux** (Ubuntu 24.04, ARM64, GCC 13, HDF5 1.10): the same, and `pdbx`'s `make test` passes.
- **Windows:** not run in a Windows VMD yet.
  - pdbx (without zlib) and sdf were cross-compiled with MinGW-w64. The DLLs export VMD's entry points and need only standard Windows system DLLs.
  - `mesh_optimize` also cross-compiles.
  - The HDF5 plugins were only compile-checked for Windows.
  - The MSVC commands follow VMD's own Windows build flags but were not run.

## Loading the plugins in VMD

- **One session:** `plugin dlopen /path/to/pdbxplugin.so` in VMD's console.
- **Every session:** put the same line in VMD's startup file:
  - **macOS and Linux:** `~/.vmdrc`.
  - **Windows:** `vmd.rc`. VMD reads it from the folder it starts in, then from `%HOME%` if set, then from the VMD folder.
- **Installed:** copy the `.so` file into VMD's plugin folder; you may need admin rights. That folder is:

  ```tcl
  file join $env(VMDDIR) plugins [vmdinfo arch] molfile
  ```

  On macOS and Linux, `make install` does this (`make -f Makefile.h5md_gromacs install` for gromacs_h5md).
- **Plugin path (macOS and Linux):** set `VMDPLUGINPATH=dir1:dir2` before starting VMD. Each folder holds a `molfile/` subfolder with the `.so` files. This does not work on Windows, because VMD splits the variable at `:`, which breaks drive letters.

To see what VMD loaded, type `plugin list "mol file reader"` in VMD's console. `plugin info "mol file reader" pdbx p; parray p` shows the version in use. Each plugin's README describes its options.

## Tests

`pdbx` has a full test suite (`cd pdbx && make test`, macOS and Linux; see [its README](pdbx/README.md)). The other plugins have VMD test scripts in their folders, for example `sdf_reader_writer/test_sdf_plugin.tcl`.

## License

Everything in this folder is under the University of Illinois Open Source License
Copyright 2026 Theoretical and Computational Biophysics Group,
All rights reserved.

Developed by:           Theoretical and Computational Biophysics Group
                        University of Illinois at Urbana-Champaign
                        http://www.ks.uiuc.edu/
