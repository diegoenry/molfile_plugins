#!/bin/sh
# Build the release plugins for one VMD architecture.
#
#   release/build.sh ARCH     ARCH: MACOSXARM64 MACOSXX86_64 LINUXAMD64 LINUXARM64 WIN64
#
# Normally run by make_release.sh, which prepares the work folder $W first:
# $W/src (the repository at the commit being released), $W/inc (VMD's plugin
# headers) and $W/dl (the HDF5 and zlib sources).
#
# Runs on macOS: Apple clang for macOS, zig for Linux (targeting an old glibc,
# so the plugins run on old and new distributions), MinGW-w64 for Windows.
# HDF5 and zlib are built from source as static libraries and linked into the
# plugins, so users need nothing installed. Only VMD's entry points
# (vmdplugin_init/register/fini) are exported, so the copies of HDF5 inside
# different plugins can't clash. Output: $W/out/ARCH/molfile/*.so
set -eu
R=$(cd "$(dirname "$0")" && pwd)
. "$R/config.sh"

ARCH=${1:?usage: build.sh ARCH}
JOBS=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)

B=$W/bld/$ARCH          # build trees
D=$W/deps/$ARCH         # static HDF5 and zlib
O=$W/out/$ARCH/molfile  # plugins
# HDF5 and zlib are built once per architecture; REBUILD_DEPS=1 rebuilds them
REBUILD_DEPS=${REBUILD_DEPS:-0}
[ -f "$D/.done-hdf5-$HDF5_VER-zlib-$ZLIB_VER" ] || REBUILD_DEPS=1
if [ "$REBUILD_DEPS" = 1 ]; then rm -rf "$B" "$D"; fi
rm -rf "$O"
mkdir -p "$B/bin" "$D" "$O"

case $ARCH in
  MACOSXARM64)  KIND=mac; MARCH=arm64;  MINOS=$MACOS_ARM64_MIN ;;
  MACOSXX86_64) KIND=mac; MARCH=x86_64; MINOS=$MACOS_X86_64_MIN ;;
  LINUXAMD64)   KIND=linux; ZT=x86_64-linux-gnu.$GLIBC_MIN;  CPU=x86_64 ;;
  LINUXARM64)   KIND=linux; ZT=aarch64-linux-gnu.$GLIBC_MIN; CPU=aarch64 ;;
  WIN64)        KIND=win ;;
  *) echo "unknown ARCH $ARCH (use one of: $ALL_ARCHES)" >&2; exit 2 ;;
esac

# --- toolchain ---------------------------------------------------------------
TC=$B/toolchain.cmake
case $KIND in
  mac)
    export MACOSX_DEPLOYMENT_TARGET=$MINOS
    CC="clang -arch $MARCH -mmacosx-version-min=$MINOS"
    CXX="clang++ -arch $MARCH -mmacosx-version-min=$MINOS"
    CMAKE_TC="-DCMAKE_OSX_ARCHITECTURES=$MARCH -DCMAKE_OSX_DEPLOYMENT_TARGET=$MINOS"
    # -x: no local symbols in the output
    PLUGIN_LD="-bundle -Wl,-x -Wl,-exported_symbols_list,$R/exports_macos.txt"
    PLUGIN_DEF=""
    ZLIB_LIB="-lz"                       # part of macOS
    SYS_LIBS="-lm"
    ;;
  linux)
    for t in cc c++ ar ranlib; do
      case $t in cc|c++) cmd="zig $t -target $ZT" ;; *) cmd="zig $t" ;; esac
      printf '#!/bin/sh\nexec %s "$@"\n' "$cmd" > "$B/bin/$t"; chmod +x "$B/bin/$t"
    done
    CC="$B/bin/cc"; CXX="$B/bin/c++"
    cat > "$TC" <<EOF
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR $CPU)
set(CMAKE_C_COMPILER $B/bin/cc)
set(CMAKE_CXX_COMPILER $B/bin/c++)
set(CMAKE_AR $B/bin/ar CACHE FILEPATH "")
set(CMAKE_RANLIB $B/bin/ranlib CACHE FILEPATH "")
set(CMAKE_FIND_ROOT_PATH $D)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
EOF
    CMAKE_TC="-DCMAKE_TOOLCHAIN_FILE=$TC"
    PLUGIN_LD="-shared -s -Wl,--version-script=$R/exports_linux.map"
    PLUGIN_DEF=""
    ZLIB_LIB="$D/lib/libz.a"
    SYS_LIBS="-lm -ldl -lpthread"
    ;;
  win)
    CC="x86_64-w64-mingw32-gcc"; CXX="x86_64-w64-mingw32-g++"
    cat > "$TC" <<EOF
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR x86_64)
set(CMAKE_C_COMPILER $(command -v x86_64-w64-mingw32-gcc))
set(CMAKE_CXX_COMPILER $(command -v x86_64-w64-mingw32-g++))
set(CMAKE_RC_COMPILER $(command -v x86_64-w64-mingw32-windres))
set(CMAKE_AR $(command -v x86_64-w64-mingw32-ar) CACHE FILEPATH "")
set(CMAKE_RANLIB $(command -v x86_64-w64-mingw32-ranlib) CACHE FILEPATH "")
set(CMAKE_FIND_ROOT_PATH $D)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
EOF
    CMAKE_TC="-DCMAKE_TOOLCHAIN_FILE=$TC"
    # -static: the GCC, C++ and winpthreads runtimes go inside the DLL.
    # The .def file exports only the entry points; without it MinGW also
    # exports the C++ runtime and zlib from pdbxplugin.
    PLUGIN_LD="-shared -static -s $R/exports_win.def"
    PLUGIN_DEF="-DVMDPLUGIN_EXPORTS"
    ZLIB_LIB="$D/lib/libz.a"
    SYS_LIBS=""
    ;;
esac

# --- zlib (Linux, Windows; macOS has its own) --------------------------------
if [ "$KIND" != mac ] && [ "$REBUILD_DEPS" = 1 ]; then
  echo "== zlib $ZLIB_VER"
  mkdir -p "$B/zlib" "$D/include" "$D/lib"
  tar -xzf "$W/dl/zlib-$ZLIB_VER.tar.gz" -C "$B/zlib"
  zs=$B/zlib/zlib-$ZLIB_VER
  zdef=""; [ "$KIND" = linux ] && zdef="-DHAVE_UNISTD_H -D_LARGEFILE64_SOURCE=1"
  objs=""
  for f in adler32 compress crc32 deflate gzclose gzlib gzread gzwrite infback \
           inffast inflate inftrees trees uncompr zutil; do
    $CC -O2 -fPIC $zdef -c "$zs/$f.c" -o "$B/zlib/$f.o"
    objs="$objs $B/zlib/$f.o"
  done
  rm -f "$D/lib/libz.a"
  # shellcheck disable=SC2086
  case $KIND in
    linux) zig ar rcs "$D/lib/libz.a" $objs ;;
    win)   x86_64-w64-mingw32-ar rcs "$D/lib/libz.a" $objs ;;
  esac
  cp "$zs/zlib.h" "$zs/zconf.h" "$D/include/"
  cp "$zs/LICENSE" "$D/zlib-LICENSE.txt"
fi

# --- HDF5 --------------------------------------------------------------------
if [ "$REBUILD_DEPS" = 1 ]; then
  echo "== HDF5 $HDF5_VER (a few minutes)"
  mkdir -p "$B/hdf5"
  tar -xzf "$W/dl/hdf5-$HDF5_VER.tar.gz" -C "$B/hdf5"
  hs=$B/hdf5/hdf5-$HDF5_VER
  zopt=""
  [ "$KIND" != mac ] && zopt="-DZLIB_INCLUDE_DIR=$D/include -DZLIB_LIBRARY=$D/lib/libz.a"
  # The C library only, with the deflate filter (GROMACS H5MD files use it);
  # nothing is downloaded during the build.
  # shellcheck disable=SC2086
  cmake -S "$hs" -B "$B/hdf5/build" -G Ninja $CMAKE_TC \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$D" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DONLY_SHARED_LIBS=OFF \
    -DBUILD_TESTING=OFF -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_UTILS=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_HL_LIB=OFF -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_BUILD_FORTRAN=OFF -DHDF5_BUILD_JAVA=OFF -DHDF5_BUILD_DOC=OFF \
    -DHDF5_ENABLE_PARALLEL=OFF -DHDF5_ENABLE_THREADSAFE=OFF \
    -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DHDF5_ENABLE_SZIP_SUPPORT=OFF \
    -DHDF5_ENABLE_PLUGIN_SUPPORT=OFF -DHDF5_ALLOW_EXTERNAL_SUPPORT=NO \
    $zopt > "$B/hdf5/configure.log" 2>&1 || { tail -40 "$B/hdf5/configure.log"; exit 1; }
  cmake --build "$B/hdf5/build" -j "$JOBS" > "$B/hdf5/build.log" 2>&1 \
    || { tail -40 "$B/hdf5/build.log"; exit 1; }
  cmake --install "$B/hdf5/build" > "$B/hdf5/install.log" 2>&1
  grep -q 'define H5_HAVE_FILTER_DEFLATE 1' "$D/include/H5pubconf.h" \
    || { echo "HDF5 was built without the deflate filter" >&2; exit 1; }
  cp "$hs/LICENSE" "$D/hdf5-LICENSE.txt" 2>/dev/null || cp "$hs/COPYING" "$D/hdf5-LICENSE.txt"
  rm -rf "$B/hdf5"
  touch "$D/.done-hdf5-$HDF5_VER-zlib-$ZLIB_VER"
fi
HDF5_LIB=$(ls "$D"/lib/libhdf5.a "$D"/lib/libhdf5-static.a 2>/dev/null | head -1)
[ -n "$HDF5_LIB" ] || { echo "no static HDF5 library in $D/lib" >&2; exit 1; }

# --- plugins -----------------------------------------------------------------
S=$W/src
CFLAGS="-O2 -fPIC -I$W/inc $PLUGIN_DEF"
H5="-I$D/include"
H5LIBS="$HDF5_LIB $ZLIB_LIB $SYS_LIBS"
[ "$KIND" = win ] && H5LIBS="$H5LIBS -lshlwapi"
build() {  # build OUTPUT COMPILER SOURCE EXTRA...
  out=$1; comp=$2; src=$3; shift 3
  echo "== $ARCH $out"
  # shellcheck disable=SC2086
  $comp $CFLAGS -o "$O/$out" "$src" $PLUGIN_LD "$@"
}
# shellcheck disable=SC2086
build pdbxplugin.so "$CXX -std=c++11" "$S/pdbx/pdbxplugin.C" -D_USE_ZLIB -I"$D/include" $ZLIB_LIB $SYS_LIBS
# shellcheck disable=SC2086
build sdfplugin.so "$CC" "$S/sdf_reader_writer/sdfplugin.c" $SYS_LIBS
# shellcheck disable=SC2086
build cndbplugin.so "$CC" "$S/CNDB_hdf5_reader/cndbplugin.c" $H5 $H5LIBS
# shellcheck disable=SC2086
build h5mdplugin_gromacs.so "$CC" "$S/gromacs_h5md/h5mdplugin_gromacs.c" $H5 $H5LIBS
# shellcheck disable=SC2086
build meshhdf5plugin.so "$CC" "$S/MeshHDF5/meshhdf5plugin.c" $H5 $H5LIBS
