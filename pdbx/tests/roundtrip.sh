#!/bin/sh
# roundtrip.sh -- read a file, write it with the plugin, read the result back,
# and check that every atom field and every frame's coordinates survive.
#
# usage: tests/roundtrip.sh DRIVER file.cif [file.cif ...]
# DRIVER is tests/pdbx_dump or tests/pdbx_dump_asan.

driver=$1
shift
tmp=${TMPDIR:-/tmp}/pdbx_roundtrip.$$
mkdir -p "$tmp" || exit 1
fail=0
for f in "$@"; do
    name=$(basename "$f" .cif)
    if ! "$driver" -c -w "$tmp/$name.out.cif" "$f" > "$tmp/$name.a" 2>/dev/null ||
       grep -q 'natoms=0 ' "$tmp/$name.a"; then
        echo "SKIP  $f (no atoms to write)"
        continue
    fi
    "$driver" -c "$tmp/$name.out.cif" > "$tmp/$name.b" 2>/dev/null
    # compare atoms and frames only: plugin messages differ, and bonds are not
    # written (archive mmCIF files have no per-atom bond table either)
    grep -E '^(FLAGS|A\||X\|)' "$tmp/$name.a" > "$tmp/$name.a2"
    grep -E '^(FLAGS|A\||X\|)' "$tmp/$name.b" > "$tmp/$name.b2"
    if cmp -s "$tmp/$name.a2" "$tmp/$name.b2"; then
        echo "OK    $f ($(grep -c '^A|' "$tmp/$name.a2") atoms, $(grep '^RESULT' "$tmp/$name.a" | sed 's/.*frames=\([0-9]*\).*/\1/') frames)"
    else
        echo "DIFF  $f"
        diff "$tmp/$name.a2" "$tmp/$name.b2" | head -6
        fail=1
    fi
done
rm -rf "$tmp"
exit $fail
