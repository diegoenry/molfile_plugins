#!/bin/sh
# Run the Linux builds' tests in containers and compare them with the macOS
# results. Run after make_release.sh.
#
#   release/test_linux.sh
#
# Environment (all optional):
#   DOCKER  the docker command (default: docker). With colima but no docker
#           CLI: DOCKER="limactl shell colima sudo docker" LIMA_HOME=~/.colima/_lima
#   IMAGE   Linux image to test in (default: ubuntu:24.04); any glibc
#           distribution works, and it doesn't need HDF5 or a compiler. It is
#           never downloaded here: docker pull --platform linux/amd64 IMAGE
#   ARCHES  which Linux builds to test (default: LINUXARM64 LINUXAMD64)
#
# The containers run without network and see the work folder read-only. On
# Apple silicon, x86_64 containers need emulation, which you can turn on until
# the next restart of the Docker VM with:
#   $DOCKER run --privileged --rm tonistiigi/binfmt --install amd64
set -u
export LC_ALL=C
R=$(cd "$(dirname "$0")" && pwd)
. "$R/config.sh"
DOCKER=${DOCKER:-docker}
IMAGE=${IMAGE:-ubuntu:24.04}
ARCHES=${ARCHES:-LINUXARM64 LINUXAMD64}

ref=""
for a in MACOSXARM64 MACOSXX86_64; do
  [ -s "$W/results-$a.txt" ] && { ref=$W/results-$a.txt; break; }
done
[ -n "$ref" ] || { echo "no macOS results to compare with; run make_release.sh first" >&2; exit 1; }

# The work folder, this folder, and the folders of linked test files
mounts="-v $W:$W:ro -v $R:$R:ro"
for f in "$W"/testfiles/*; do
  [ -L "$f" ] || continue
  d=$(dirname "$(readlink "$f")")
  case "$mounts" in *"-v $d:"*) ;; *) mounts="$mounts -v $d:$d:ro" ;; esac
done

status=0
for a in $ARCHES; do
  case $a in
    LINUXARM64) plat=linux/arm64 ;;
    LINUXAMD64) plat=linux/amd64 ;;
    *) echo "$a: not a Linux build" >&2; status=1; continue ;;
  esac
  [ -d "$W/out/$a/molfile" ] || { echo "$a: not built"; status=1; continue; }
  # shellcheck disable=SC2086
  if ! $DOCKER run --rm --pull never --network none --platform $plat --entrypoint true "$IMAGE" >/dev/null 2>&1; then
    echo "$a: can't run a $plat container of $IMAGE. Get the image with"
    echo "      $DOCKER pull --platform $plat $IMAGE"
    echo "    and, for x86_64 on Apple silicon, turn on emulation (see the top of this script)."
    status=1; continue
  fi
  # shellcheck disable=SC2086
  os=$($DOCKER run --rm --network none --platform $plat --entrypoint cat "$IMAGE" /etc/os-release \
       | sed -n 's/^PRETTY_NAME="*\([^"]*\)"*/\1/p')
  # shellcheck disable=SC2086
  $DOCKER run --rm --network none --platform $plat $mounts -e W="$W" --entrypoint sh "$IMAGE" \
    "$R/run_tests.sh" "$W/out/$a/molfile" "$W/loadtest/loadtest-$a" > "$W/results-$a.txt" 2>&1
  n=$(wc -l < "$W/results-$a.txt" | tr -d ' ')
  if grep -q ' exit=' "$W/results-$a.txt"; then
    echo "$a ($os): a test crashed:"; grep ' exit=' "$W/results-$a.txt"; status=1
  elif cmp -s "$ref" "$W/results-$a.txt"; then
    echo "$a ($os): same results as $(basename "$ref" .txt | sed 's/results-//') ($n test cases)"
  else
    echo "$a ($os): results differ from $(basename "$ref" .txt | sed 's/results-//'):"
    diff "$ref" "$W/results-$a.txt" | head -20; status=1
  fi
done
exit $status
