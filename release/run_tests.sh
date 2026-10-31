#!/bin/sh
# Read every test file with every plugin of one build, and print one
# fingerprint line per file. Builds for different platforms must print the
# same lines; make_release.sh and test_linux.sh compare them.
#
#   run_tests.sh PLUGIN_DIR LOADTEST [RUNNER...]
#
# LOADTEST is loadtest.c built for the platform; RUNNER runs it, for example
# "arch -x86_64" for Intel builds on Apple silicon. Test files: the pdbx test
# files, sdf_reader_writer/example.sdf and gromacs_h5md/equil.h5md from the
# repository, plus $W/testfiles (made by make_release.sh: .gz copies, and any
# .cndb, .h5md and .h5mesh files from EXTRA_TESTFILES). A crash shows as exit=N.
set -u
export LC_ALL=C   # same file order on every system
W=${W:-$HOME/.cache/molfile-release}
P=$1; LT=$2; shift 2
T=$W/testfiles
S=$W/src

t() {  # t LABEL LOADTEST-ARGS...
  label=$1; shift
  out=$("$@" 2>/dev/null </dev/null)
  rc=$?
  printf '%s: %s' "$label" "$(printf '%s' "$out" | tr '\n' ' ')"
  [ $rc -le 1 ] || printf ' exit=%s' "$rc"   # 1 = clean failure
  echo
}
each() {  # each PATTERN...: the files that exist
  for f in "$@"; do [ -e "$f" ] && printf '%s\n' "$f"; done
}

for so in pdbxplugin sdfplugin cndbplugin h5mdplugin_gromacs meshhdf5plugin; do
  t "register $so" "$@" "$LT" "$P/$so.so"
done
each "$S"/pdbx/tests/data/*.cif "$S"/pdbx/tests/bcif/*.bcif "$T"/*.cif.gz "$T"/*.bcif.gz |
while read -r f; do
  t "pdbx $(basename "$f")" "$@" "$LT" "$P/pdbxplugin.so" pdbx "$f"
done
t "mmcif alias" "$@" "$LT" "$P/pdbxplugin.so" mmcif "$S/pdbx/tests/data/chains.cif"
t "sdf example.sdf" "$@" "$LT" "$P/sdfplugin.so" sdf "$S/sdf_reader_writer/example.sdf"
each "$T"/*.cndb | while read -r f; do
  t "cndb $(basename "$f") (5 frames)" "$@" "$LT" "$P/cndbplugin.so" cndb "$f" 5
done
each "$S/gromacs_h5md/equil.h5md" "$T"/*.h5md | while read -r f; do
  t "h5md $(basename "$f")" "$@" "$LT" "$P/h5mdplugin_gromacs.so" h5md_gromacs "$f"
done
each "$T"/*.h5mesh | while read -r f; do
  t "h5mesh $(basename "$f")" "$@" "$LT" "$P/meshhdf5plugin.so" h5mesh "$f"
done
for pt in cndbplugin:cndb h5mdplugin_gromacs:h5md_gromacs meshhdf5plugin:h5mesh; do
  so=${pt%%:*}; ty=${pt#*:}
  t "wrong file $ty" "$@" "$LT" "$P/$so.so" "$ty" "$S/sdf_reader_writer/example.sdf"
  t "missing file $ty" "$@" "$LT" "$P/$so.so" "$ty" "$T/does-not-exist"
done
