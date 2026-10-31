#!/bin/sh
# Build, check, test and package the molfile plugins for every VMD platform.
#
#   release/make_release.sh [TAG]
#
# TAG names the archives (default: the short commit hash). The committed state
# of HEAD is built; uncommitted changes are not included. Runs on macOS; see
# release/README.md for the tools it needs.
#
# Environment (all optional):
#   VMD_PLUGIN_INCLUDE  folder with VMD's molfile_plugin.h and vmdplugin.h
#                       (default: the first VMD found in /Applications)
#   EXTRA_TESTFILES     folder with more test files: *.cndb, *.h5md, *.h5mesh
#   ARCHES              architectures to build (default: all five)
#   DIST                where the archives go (default: dist/ in the repository)
#   W                   work folder (default: ~/.cache/molfile-release)
set -eu
export LC_ALL=C
R=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$R/.." && pwd)
. "$R/config.sh"
ARCHES=${ARCHES:-$ALL_ARCHES}
DIST=${DIST:-$REPO/dist}
COMMIT=$(git -C "$REPO" rev-parse --short HEAD)
TAG=${1:-$COMMIT}

die() { echo "make_release: $*" >&2; exit 1; }
has() { case " $ARCHES " in *" $1 "*) return 0 ;; esac; return 1; }
[ "$(uname -s)" = Darwin ] || die "this script runs on macOS"
case $W in *" "*) die "the work folder must not contain spaces: $W" ;; esac
for a in $ARCHES; do
  case " $ALL_ARCHES " in *" $a "*) ;; *) die "unknown architecture $a (use: $ALL_ARCHES)" ;; esac
done

# --- tools ---------------------------------------------------------------------
missing=""
need() { command -v "$1" >/dev/null 2>&1 || missing="$missing $1"; }
for t in clang clang++ cmake ninja curl shasum git zip gzip; do need $t; done
{ has LINUXAMD64 || has LINUXARM64; } && need zig
has WIN64 && { need x86_64-w64-mingw32-gcc; need x86_64-w64-mingw32-g++; need x86_64-w64-mingw32-objdump; }
[ -z "$missing" ] || die "missing tools:$missing
  (Xcode command line tools: xcode-select --install; the rest: brew install cmake ninja zig mingw-w64)"

# --- VMD's plugin headers ------------------------------------------------------
if [ -z "${VMD_PLUGIN_INCLUDE:-}" ]; then
  for d in /Applications/VMD*.app/Contents/vmd*/plugins/include \
           /Applications/VMD*.app/Contents/vmd*/lib/plugins/include \
           "$HOME"/Applications/VMD*.app/Contents/vmd*/plugins/include; do
    if [ -f "$d/molfile_plugin.h" ]; then VMD_PLUGIN_INCLUDE=$d; break; fi
  done
fi
[ -f "${VMD_PLUGIN_INCLUDE:-}/molfile_plugin.h" ] \
  || die "VMD's plugin headers not found; set VMD_PLUGIN_INCLUDE to the folder with molfile_plugin.h"
mkdir -p "$W/inc" "$W/dl"
cp "$VMD_PLUGIN_INCLUDE/molfile_plugin.h" "$VMD_PLUGIN_INCLUDE/vmdplugin.h" "$W/inc/"
echo "== headers: $VMD_PLUGIN_INCLUDE"

# --- sources, checked against their SHA-256 ----------------------------------------
fetch() {  # fetch FILE SHA256 URL...
  f=$W/dl/$1; sum=$2; shift 2
  if [ -f "$f" ] && echo "$sum  $f" | shasum -a 256 -c >/dev/null 2>&1; then return 0; fi
  for u in "$@"; do
    echo "== downloading $u"
    if curl -fL --retry 3 -o "$f.part" "$u" && echo "$sum  $f.part" | shasum -a 256 -c >/dev/null 2>&1; then
      mv "$f.part" "$f"; return 0
    fi
    echo "   download failed or checksum does not match"
  done
  rm -f "$f.part"
  die "could not get $1"
}
# shellcheck disable=SC2086
fetch "hdf5-$HDF5_VER.tar.gz" "$HDF5_SHA256" $HDF5_URLS
# shellcheck disable=SC2086
fetch "zlib-$ZLIB_VER.tar.gz" "$ZLIB_SHA256" $ZLIB_URLS
if has WIN64; then
  # shellcheck disable=SC2086
  fetch winpthreads-COPYING.txt "$WINPTHREADS_LICENSE_SHA256" $WINPTHREADS_LICENSE_URLS
  mc=${MINGW_COPYING:-$(brew --prefix mingw-w64 2>/dev/null || true)/COPYING}
  [ -f "$mc" ] || die "MinGW-w64's COPYING not found; set MINGW_COPYING to it"
  cp "$mc" "$W/dl/mingw-w64-COPYING.txt"
fi

# --- the source being released -------------------------------------------------
echo "== source: commit $COMMIT, archives named molfile_plugins-$TAG-*"
git -C "$REPO" diff --quiet HEAD -- 2>/dev/null || echo "   note: uncommitted changes are not included"
rm -rf "$W/src"; mkdir -p "$W/src"
git -C "$REPO" archive HEAD | tar -x -C "$W/src"
echo "$COMMIT" > "$W/src/.release-commit"

# --- test files ----------------------------------------------------------------
T=$W/testfiles
rm -rf "$T"; mkdir -p "$T"
gzip -9c "$W/src/pdbx/tests/data/chains.cif" > "$T/chains.cif.gz"
gzip -9c "$W/src/pdbx/tests/data/models.cif" > "$T/models.cif.gz"
gzip -9c "$W/src/pdbx/tests/bcif/chains.bcif" > "$T/chains.bcif.gz"
if [ -n "${EXTRA_TESTFILES:-}" ]; then
  [ -d "$EXTRA_TESTFILES" ] || die "EXTRA_TESTFILES is not a folder: $EXTRA_TESTFILES"
  for f in "$EXTRA_TESTFILES"/*.cndb "$EXTRA_TESTFILES"/*.h5md "$EXTRA_TESTFILES"/*.h5mesh; do
    [ -e "$f" ] && ln -s "$(readlink -f "$f")" "$T/"
  done
fi
for ext in cndb h5mesh; do
  ls "$T"/*.$ext >/dev/null 2>&1 \
    || echo "   note: no .$ext test file (EXTRA_TESTFILES); that plugin only gets the failure tests"
done

# --- build, in parallel ------------------------------------------------------------
echo "== building: $ARCHES (logs: $W/build-ARCH.log)"
for a in $ARCHES; do
  ( sh "$R/build.sh" "$a" > "$W/build-$a.log" 2>&1; echo $? > "$W/build-$a.status" ) &
done
wait
for a in $ARCHES; do
  [ "$(cat "$W/build-$a.status")" = 0 ] || { tail -30 "$W/build-$a.log"; die "build failed for $a"; }
done

# --- what each plugin exports and needs --------------------------------------------
echo "== checking exports and dependencies"
want="vmdplugin_fini vmdplugin_init vmdplugin_register"
for a in $ARCHES; do
  for so in "$W/out/$a/molfile/"*.so; do
    n=$(basename "$so")
    case $a in
      MACOSX*) ex=$(nm -gU "$so" | awk '{sub(/^_/, "", $3); print $3}')
               deps=$(otool -L "$so" | tail -n +2 | awk '{print $1}' | grep -v '^/usr/lib/' || true) ;;
      LINUX*)  ex=$(nm -D --defined-only "$so" | awk '{print $NF}')
               deps=$(objdump -p "$so" | awk '/NEEDED/{print $2}' \
                      | grep -vE '^(libc|libm|libdl|libpthread|ld-linux-[a-z0-9_-]+)\.so\.[0-9]+$' || true)
               g=$(objdump -T "$so" | grep -oE 'GLIBC_[0-9]+\.[0-9]+' | sed 's/GLIBC_//' \
                   | awk -F. -v m="$GLIBC_MIN" 'BEGIN{split(m,v,".")} $1>v[1] || ($1==v[1] && $2>v[2])' | head -1)
               [ -z "$g" ] || die "$a/$n needs glibc $g (newer than $GLIBC_MIN)" ;;
      WIN64)   ex=$(x86_64-w64-mingw32-objdump -p "$so" \
                    | awk '/\[Ordinal\/Name Pointer\] Table/{t=1; next} t && /\[/{print $NF} t && /^$/{t=0}')
               deps=$(x86_64-w64-mingw32-objdump -p "$so" | awk '/DLL Name:/{print $3}' \
                      | grep -viE '^(api-ms-win-crt-[a-z0-9-]+|kernel32|shlwapi)\.dll$' || true) ;;
    esac
    [ "$(echo $ex | tr ' ' '\n' | sort | tr '\n' ' ')" = "$want " ] \
      || die "$a/$n exports more than VMD's entry points: $(echo $ex | head -c 300)"
    [ -z "$deps" ] || die "$a/$n needs libraries users may not have: $deps"
  done
done
echo "   each plugin exports only $want and needs only system libraries"

# --- tests -------------------------------------------------------------------------
LT=$W/loadtest; mkdir -p "$LT"
has MACOSXARM64  && clang -O2 -arch arm64  -I"$W/inc" -o "$LT/loadtest-MACOSXARM64"  "$R/loadtest.c"
has MACOSXX86_64 && clang -O2 -arch x86_64 -I"$W/inc" -o "$LT/loadtest-MACOSXX86_64" "$R/loadtest.c"
has LINUXAMD64 && zig cc -target "x86_64-linux-gnu.$GLIBC_MIN"  -O2 -I"$W/inc" -o "$LT/loadtest-LINUXAMD64" "$R/loadtest.c" -ldl
has LINUXARM64 && zig cc -target "aarch64-linux-gnu.$GLIBC_MIN" -O2 -I"$W/inc" -o "$LT/loadtest-LINUXARM64" "$R/loadtest.c" -ldl
has WIN64 && x86_64-w64-mingw32-gcc -O2 -I"$W/inc" -o "$LT/loadtest-WIN64.exe" "$R/loadtest.c"

rosetta=no
[ "$(uname -m)" = arm64 ] && arch -x86_64 /usr/bin/true 2>/dev/null && rosetta=yes
ref=""
for a in MACOSXARM64 MACOSXX86_64; do
  has $a || continue
  case $a:$(uname -m):$rosetta in
    MACOSXARM64:arm64:*)  run="" ;;
    MACOSXX86_64:x86_64:*) run="" ;;
    MACOSXX86_64:arm64:yes) run="arch -x86_64" ;;
    *) echo "   $a: can't run on this Mac, not tested"; continue ;;
  esac
  # shellcheck disable=SC2086
  sh "$R/run_tests.sh" "$W/out/$a/molfile" "$LT/loadtest-$a" $run > "$W/results-$a.txt"
  n=$(wc -l < "$W/results-$a.txt" | tr -d ' ')
  ! grep -q ' exit=' "$W/results-$a.txt" || { grep ' exit=' "$W/results-$a.txt"; die "$a: a test crashed"; }
  if [ -z "$ref" ]; then
    ref=$a; echo "   $a: $n test cases, no crashes"
  elif cmp -s "$W/results-$ref.txt" "$W/results-$a.txt"; then
    echo "   $a: same results as $ref ($n test cases)"
  else
    diff "$W/results-$ref.txt" "$W/results-$a.txt" | head -20
    die "$a reads the test files differently from $ref"
  fi
done

# --- package -------------------------------------------------------------------------
mkdir -p "$DIST"
rm -f "$DIST/molfile_plugins-$TAG-"*
for a in $ARCHES; do sh "$R/package.sh" "$a" "$TAG" "$DIST"; done
( cd "$DIST" && shasum -a 256 "molfile_plugins-$TAG-"*.tar.gz "molfile_plugins-$TAG-"*.zip 2>/dev/null \
    > "molfile_plugins-$TAG-SHA256SUMS.txt" )
echo "== release files in $DIST:"
( cd "$DIST" && ls -l "molfile_plugins-$TAG-"* | awk '{printf "   %9d  %s\n", $5, $NF}' )
if has LINUXAMD64 || has LINUXARM64; then
  echo "Next: release/test_linux.sh runs the same tests on the Linux builds, in containers."
fi
if has WIN64; then
  echo "The Windows build can't run here; release/README.md has a quick check to run on Windows."
fi
