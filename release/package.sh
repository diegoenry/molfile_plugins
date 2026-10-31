#!/bin/sh
# Package one built architecture as DIST/molfile_plugins-TAG-ARCH.tar.gz
# (.zip for Windows), with a README, the plugin license and the licenses of
# the libraries built into the plugins.
#
#   release/package.sh ARCH TAG DIST
set -eu
R=$(cd "$(dirname "$0")" && pwd)
. "$R/config.sh"
ARCH=$1; TAG=$2; DIST=$3
NAME=molfile_plugins-$TAG-$ARCH
STAGE=$W/stage/$NAME
rm -rf "$STAGE"; mkdir -p "$STAGE/molfile" "$DIST"
cp "$W/out/$ARCH/molfile/"*.so "$STAGE/molfile/"
cp "$W/src/LICENSE" "$STAGE/LICENSE.txt"

case $ARCH in
  MACOSXARM64)  REQ="macOS $MACOS_ARM64_MIN or later, Apple silicon (VMD's MACOSXARM64 build)" ;;
  MACOSXX86_64) REQ="macOS $MACOS_X86_64_MIN or later, Intel (VMD's MACOSXX86_64 build)" ;;
  LINUXAMD64)   REQ="Linux x86_64 with glibc $GLIBC_MIN or later (VMD's LINUXAMD64 build)" ;;
  LINUXARM64)   REQ="Linux ARM64 with glibc $GLIBC_MIN or later (VMD's LINUXARM64 build)" ;;
  WIN64)        REQ="64-bit VMD for Windows (WIN64), Windows 10 or later" ;;
esac
case $ARCH in
  MACOSX*) BUNDLED="HDF5 $HDF5_VER is built into the HDF5 plugins (zlib comes with macOS)" ;;
  *)       BUNDLED="HDF5 $HDF5_VER and zlib $ZLIB_VER are built into the plugins" ;;
esac
COMMIT=$(cat "$W/src/.release-commit")
LINE="=================================================================================="

{
  printf 'HDF5 %s - built into cndbplugin, h5mdplugin_gromacs and meshhdf5plugin\n%s\n\n' "$HDF5_VER" "$LINE"
  cat "$W/deps/$ARCH/hdf5-LICENSE.txt"
  if [ -f "$W/deps/$ARCH/zlib-LICENSE.txt" ]; then
    printf '\n\nzlib %s - built into pdbxplugin and the HDF5 plugins\n%s\n\n' "$ZLIB_VER" "$LINE"
    cat "$W/deps/$ARCH/zlib-LICENSE.txt"
  fi
  if [ "$ARCH" = WIN64 ]; then
    printf '\n\nmingw-w64 winpthreads - built into every plugin except sdfplugin\n%s\n\n' "$LINE"
    cat "$W/dl/winpthreads-COPYING.txt"
    printf '\n\nmingw-w64 runtime - DLL startup code in every plugin\n%s\n\n' "$LINE"
    cat "$W/dl/mingw-w64-COPYING.txt"
  fi
} > "$STAGE/THIRD_PARTY_LICENSES.txt"

cat > "$STAGE/README.txt" <<EOF
VMD molfile plugins - $ARCH
Built from https://github.com/diegoenry/molfile_plugins, commit $COMMIT

Plugins (in molfile/)
  File                    VMD type          Version  Reads
  pdbxplugin.so           pdbx (and mmcif)  2.1      mmCIF/PDBx and BinaryCIF: .cif .mmcif .bcif,
                                                     also gzip-compressed; also writes mmCIF
  sdfplugin.so            sdf               1.1      SDF/MOL V2000: .sdf .mol; also writes them
  cndbplugin.so           cndb              0.2      OpenMiChroM CNDB: .cndb
  h5mdplugin_gromacs.so   h5md_gromacs      0.1      GROMACS H5MD: .h5md
  meshhdf5plugin.so       h5mesh            1.0      HDF5 triangle meshes: .h5mesh

Requirements
  $REQ.
  $BUNDLED;
  nothing else needs to be installed.

Using the plugins
  One session: in VMD's console,
      plugin dlopen /path/to/molfile/pdbxplugin.so
  Every session: put those lines in VMD's startup file
      (~/.vmdrc on macOS and Linux, vmd.rc on Windows).
  Installed for everyone: copy molfile/*.so into the folder this prints in VMD:
      file join \$env(VMDDIR) plugins [vmdinfo arch] molfile
EOF
case $ARCH in
  WIN64) ;;
  *) echo "  Or start VMD with VMDPLUGINPATH set to this folder." >> "$STAGE/README.txt" ;;
esac
cat >> "$STAGE/README.txt" <<'EOF'

  Check what VMD loaded: plugin list "mol file reader"
  VMD guesses the type from the last extension only, so name the type for
  compressed files:  mol new 1abc.cif.gz type pdbx
  A newer version of a plugin replaces an older one with the same name.
EOF
case $ARCH in
  MACOSX*)
    if [ "$ARCH" = MACOSXARM64 ]; then exe=vmd_MACOSXARM64; else exe=vmd_MACOSXX86_64; fi
    cat >> "$STAGE/README.txt" <<EOF

macOS
  The plugins are not notarized. If VMD refuses to load one after you
  downloaded this archive with a browser, remove the quarantine flag:
      xattr -dr com.apple.quarantine /path/to/molfile

  VMD 1.9.4 from the VMD website is notarized with Apple's hardened runtime,
  which only lets it load plugins signed by the VMD developers. "plugin dlopen"
  then fails, and VMD prints "mapping process and mapped file (non-platform)
  have different Team IDs". To use these plugins with that VMD, re-sign its
  program and then the application on your Mac. This removes the hardened
  runtime from your copy only; use sudo if you can't write to /Applications:
      codesign --force --sign - "/Applications/VMD 1.9.4.app/Contents/vmd/$exe"
      codesign --force --deep --sign - "/Applications/VMD 1.9.4.app"
  Use your VMD's application name. VMD builds without the hardened runtime
  load the plugins as they are.
EOF
    ;;
esac
cat >> "$STAGE/README.txt" <<'EOF'

Licenses
  The plugins: University of Illinois Open Source License (LICENSE.txt).
  Libraries built into them: THIRD_PARTY_LICENSES.txt.
EOF

cd "$W/stage"
case $ARCH in
  WIN64) rm -f "$DIST/$NAME.zip"; zip -qr -X "$DIST/$NAME.zip" "$NAME" ;;
  *)     COPYFILE_DISABLE=1 tar --no-mac-metadata -czf "$DIST/$NAME.tar.gz" "$NAME" ;;  # no ._ files
esac
