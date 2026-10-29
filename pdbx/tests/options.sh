#!/bin/sh
# options.sh -- chain selection ("file.cif?chain=A,B", MMCIF_CHAIN), .cif.gz
# reading and writing, and the plugin names the library registers.
#
# usage: tests/options.sh DRIVER ZLIB      (ZLIB: 1 if built with zlib)

driver=$1
zlib=${2:-1}
# the plugin's environment variables must not leak in from the caller
unset MMCIF_CHAIN VMDPDBXUSELABELIDS
tmp=${TMPDIR:-/tmp}/pdbx_options.$$
mkdir -p "$tmp" || exit 1
fail=0
data=tests/data

# run with a time limit (macOS has no timeout(1) by default)
limited() { perl -e 'alarm 60; exec @ARGV or exit 127' "$@"; }

ok()  { echo "OK    $1"; }
bad() { echo "FAIL  $1"; [ -f "$tmp/out" ] && sed 's/^/      /' "$tmp/out" | head -8; fail=1; }

# check NAME 'RESULT prefix' 'message pattern or -' DRIVER-ARGS...
check() {
    name=$1; want=$2; msg=$3; shift 3
    limited "$driver" -q "$@" > "$tmp/out" 2>&1
    if ! grep -q "^RESULT $want" "$tmp/out"; then bad "$name (expected 'RESULT $want')"; return; fi
    if [ "$msg" != "-" ] && ! grep -q -- "$msg" "$tmp/out"; then
        bad "$name (expected '$msg')"; return
    fi
    ok "$name"
}

# the atoms and frames of a reading, without the plugin's own messages
dump() { limited "$driver" -c "$@" 2>/dev/null | grep -v '^pdbxplugin)'; }

# same A B: both readings produced a summary line and are identical
same() { grep -q '^RESULT' "$1" && grep -q '^RESULT' "$2" && cmp -s "$1" "$2"; }

# ---- registration --------------------------------------------------------

"$driver" -l > "$tmp/out" 2>&1
if grep -q '^PLUGIN pdbx 2\.0 cif,mmcif mol file reader reader=yes writer=yes' "$tmp/out" &&
   grep -q '^PLUGIN mmcif 2\.0 mmcif mol file reader reader=yes writer=yes' "$tmp/out"; then
    ok "registers pdbx and the mmcif alias (2.0)"
else
    bad "plugin registration"
fi

# ---- chain selection -----------------------------------------------------

f=$data/chains.cif
check "no selection"            "ok natoms=10 frames=1 bonds=3" -  "$f"
check "?chain=A"                "ok natoms=5 frames=2 bonds=1" \
      "skipped 2 explicit bond record(s) to chains that were not loaded"  "$f?chain=A"
check "?chain=A,B"              "ok natoms=7 frames=1 bonds=2" -  "$f?chain=A,B"
check "?chain=A B (blank)"      "ok natoms=7 frames=1 bonds=2" -  "$f?chain=A B"
check "?chain=AA (long ID)"     "ok natoms=2 frames=2 bonds=1" -  "$f?chain=AA"
check "?chain=a (case matters)" "ok natoms=1 frames=2 bonds=0" -  "$f?chain=a"
check "?chain=Z (absent)"       "open-failed" \
      "chain(s) Z (selected by the file name) have no atoms with coordinates in the first model, whose chains are A, B, AA, a" \
      "$f?chain=Z"
check "?chain=A,Z"              "ok natoms=5 frames=2" "chain(s) Z have no atoms with coordinates in the first model" "$f?chain=A,Z"
check "?chains=A&chain=B"       "ok natoms=7" -  "$f?chains=A&chain=B"
check "unknown option"          "open-failed" "unknown option 'colour=red'" "$f?colour=red"
check "empty selection"         "open-failed" "no chain IDs after"          "$f?chain="
check "tiny buffer + selection" "ok natoms=5 frames=2 bonds=1" -  -b 1 "$f?chain=A"

# (sh keeps an assignment made before a function call, so unset each time)
export MMCIF_CHAIN=B
check "MMCIF_CHAIN=B" "ok natoms=2 frames=1" "selected by MMCIF_CHAIN" "$f"
check "file name wins over MMCIF_CHAIN" "ok natoms=5" "selected by the file name" "$f?chain=A"
export MMCIF_CHAIN=' , '
check "blank MMCIF_CHAIN is ignored" "ok natoms=10" - "$f"
unset MMCIF_CHAIN
export VMDPDBXUSELABELIDS=1
check "label numbering selects label_asym_id" "ok natoms=1 frames=2 bonds=0" \
      "1 of 3 explicit bond record(s) did not match the loaded atoms" "$f?chain=C"
unset VMDPDBXUSELABELIDS

cp "$f" "$tmp/real?chain=A"
check "a file whose name has '?' is read as named" "ok natoms=10" - "$tmp/real?chain=A"

check "chemical component ignores the selection" "ok natoms=4 frames=1 bonds=3 orders=yes" \
      "does not apply to a chemical component file" "$data/chem_comp.cif?chain=A"
check "IHM spheres by asym_id"  "ok natoms=0 frames=0 bonds=0 orders=no graphics=2" - \
      "$data/ihm_spheres.cif?chain=A"
check "IHM spheres, other chain" "open-failed" "chain(s) B (selected by the file name) have no atoms" \
      "$data/ihm_spheres.cif?chain=B"

# chain IDs longer than segname: selected in full, or by what segname shows;
# a chain whose rows have no coordinates is not "loaded"
cat > "$tmp/long.cif" <<'EOF'
data_long
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM N N  GLY A          1 1 1 1
ATOM N N  GLY B          1 ? ? ?
ATOM N N  GLY LONGCHAIN9 1 3 3 3
ATOM N N  GLY LONGCHAIN8 1 4 4 4
EOF
check "long ID in full"          "ok natoms=1" - "$tmp/long.cif?chain=LONGCHAIN9"
check "long ID as segname shows" "ok natoms=2" - "$tmp/long.cif?chain=LONGCHA"
check "chain without coordinates" "ok natoms=1" \
      "chain(s) B have no atoms with coordinates in the file" "$tmp/long.cif?chain=A,B"

# only the atom table that is loaded counts: chain B of the IHM starting
# model is not loaded when _atom_site is
cat > "$tmp/ihm_both.cif" <<'EOF'
data_ihm_both
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM N N  GLY A 1 1.0 2.0 3.0
ATOM C CA GLY A 1 1.5 2.0 3.0
#
loop_
_ihm_starting_model_coord.starting_model_id
_ihm_starting_model_coord.asym_id
_ihm_starting_model_coord.seq_id
_ihm_starting_model_coord.comp_id
_ihm_starting_model_coord.atom_id
_ihm_starting_model_coord.type_symbol
_ihm_starting_model_coord.Cartn_x
_ihm_starting_model_coord.Cartn_y
_ihm_starting_model_coord.Cartn_z
1 B 1 GLY N N 5.0 5.0 5.0
EOF
check "chains of a table not loaded are not counted" "ok natoms=2" \
      "chain(s) B have no atoms with coordinates" "$tmp/ihm_both.cif?chain=A,B"

# an unmatched bond record inside the selection is reported as unmatched
{ cat "$data/chains.cif"; printf 'loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n'
  printf '_struct_conn.ptnr1_auth_asym_id\n_struct_conn.ptnr1_auth_seq_id\n'
  printf '_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_auth_asym_id\n'
  printf '_struct_conn.ptnr2_auth_seq_id\n_struct_conn.ptnr2_label_atom_id\n'
  printf 'typo1 covale A 1 SG A 2 XX\n'; } > "$tmp/typo.cif"
check "unmatched bond record with every chain selected" "ok natoms=10" \
      "1 of 4 explicit bond record(s) did not match the loaded atoms" "$tmp/typo.cif?chain=A,B,AA,a"
if grep -q "not loaded" "$tmp/out"; then bad "no chain was left out, yet bonds 'not loaded'"; fi

# a '?' in a directory name is part of the path
mkdir -p "$tmp/a?b"
cp "$data/chains.cif" "$tmp/a?b/x.cif"
check "'?' in a directory, with a selection" "ok natoms=5" - "$tmp/a?b/x.cif?chain=A"
check "'?' in a directory, missing file" "open-failed" "Unable to open file" "$tmp/a?b/none.cif"
mkdir -p "$tmp/run?x=1"
cp "$data/chains.cif" "$tmp/run?x=1/x.cif"
check "a directory named like an option" "ok natoms=10" - "$tmp/run?x=1/x.cif"
# a chain ID may contain '/'
sed 's/ B  1$/ A\/B 1/; s/ B  2$/ A\/B 2/' "$data/chains.cif" > "$tmp/slash.cif"
check "chain ID with '/'" "ok natoms=2" - "$tmp/slash.cif?chain=A/B"

# ---- gzip ----------------------------------------------------------------

if [ "$zlib" != "1" ]; then
    gzip -c $data/models.cif > "$tmp/models.cif.gz"
    check "gzip without zlib" "open-failed" "this build has no zlib" "$tmp/models.cif.gz"
    limited "$driver" -q -w "$tmp/w.cif.gz" $data/models.cif > "$tmp/out" 2>&1
    if grep -q "cannot write .gz files" "$tmp/out"; then ok "writing .gz without zlib"
    else bad "writing .gz without zlib"; fi
    rm -rf "$tmp"
    exit $fail
fi

# every fixture: the .gz reads exactly like the plain file, also with tiny
# buffers (damaged loops are re-read, which rewinds the compressed stream)
n=0; nbad=0
for f in $data/*.cif; do
    gzip -c "$f" > "$tmp/t.cif.gz"
    dump "$f" > "$tmp/plain"
    for b in 65536 1 3; do
        dump -b $b "$tmp/t.cif.gz" > "$tmp/gz"
        if ! same "$tmp/plain" "$tmp/gz"; then
            echo "FAIL  $f.gz differs from $f (buffer $b)"; nbad=$((nbad + 1)); fail=1
        fi
    done
    n=$((n + 1))
done
[ $nbad -eq 0 ] && ok "all $n fixtures read the same from .gz (buffers 64K, 1, 3)"

# without rewinding (as when the re-read budget is spent) a damaged loop keeps
# its strict reading: different, but it loads
f=$data/damaged_missing_value.cif
gzip -c "$f" > "$tmp/d.cif.gz"
dump "$f" > "$tmp/plain"
dump -z "$tmp/d.cif.gz" > "$tmp/gz"
if grep -q '^RESULT ok natoms=6' "$tmp/plain" && grep -q '^RESULT ok natoms=3' "$tmp/gz"; then
    ok "no rewinding: damaged loop keeps the strict reading"
else
    bad "no rewinding"
fi

# A re-read needs two rewinds (there, and back if it is not better). Whatever
# the budget, a .gz must give either the plain file's reading or the strict
# one, never a third (a budget that allowed only the first rewind kept the
# worse re-read).
nbad=0
for f in $data/*.cif; do
    gzip -c "$f" > "$tmp/t.cif.gz"
    dump "$f" > "$tmp/plain"
    dump -z -b 1 "$tmp/t.cif.gz" > "$tmp/strict"
    size=$(wc -c < "$f")
    e=0
    while [ $e -le $((3 * size)) ]; do
        dump -b 1 -Z $e "$tmp/t.cif.gz" > "$tmp/gz"
        if ! same "$tmp/plain" "$tmp/gz" && ! same "$tmp/strict" "$tmp/gz"; then
            echo "FAIL  $f.gz with a budget of $e bytes"; nbad=$((nbad + 1)); fail=1
        fi
        e=$((e + size / 16 + 1))
    done
done
[ $nbad -eq 0 ] && ok "any rewind budget gives the plain or the strict reading (all fixtures)"

f=$data/models.cif
gzip -c "$f" > "$tmp/m.cif.gz"
size=$(wc -c < "$tmp/m.cif.gz")

head -c $((size / 2)) "$tmp/m.cif.gz" > "$tmp/trunc.cif.gz"
check "truncated .gz" "open-failed" "compressed data ends early" "$tmp/trunc.cif.gz"

perl -e 'local $/; $_ = <STDIN>; substr($_, 20, 1) = chr(ord(substr($_, 20, 1)) ^ 0x55);
         print' < "$tmp/m.cif.gz" > "$tmp/bad.cif.gz"
check "damaged .gz" "open-failed" "compressed data is damaged" "$tmp/bad.cif.gz"

# two gzip members (concatenated files), then trailing bytes
head -n 15 "$f" | gzip -c > "$tmp/cat.cif.gz"
tail -n +16 "$f" | gzip -c >> "$tmp/cat.cif.gz"
dump "$f" > "$tmp/plain"
dump "$tmp/cat.cif.gz" > "$tmp/gz"
if same "$tmp/plain" "$tmp/gz"; then ok "concatenated gzip members"; else bad "concatenated gzip members"; fi
cp "$tmp/m.cif.gz" "$tmp/pad.cif.gz"; head -c 1000 /dev/zero >> "$tmp/pad.cif.gz"
dump "$tmp/pad.cif.gz" > "$tmp/gz"
if same "$tmp/plain" "$tmp/gz"; then ok "zero padding after the gzip data is ignored"
else bad "zero padding after the gzip data"; fi
cp "$tmp/m.cif.gz" "$tmp/trail.cif.gz"; printf 'trailing bytes\n' >> "$tmp/trail.cif.gz"
check "other bytes after the gzip data" "open-failed" "unexpected bytes after a gzip member" \
      "$tmp/trail.cif.gz"
# a damaged header of the second member must not load only the first one
tail -n +16 "$f" | gzip -c > "$tmp/m2.gz"
perl -e 'local $/; $_ = <STDIN>; substr($_, 1, 1) = chr(ord(substr($_, 1, 1)) ^ 1); print' \
    < "$tmp/m2.gz" > "$tmp/m2bad.gz"
head -n 15 "$f" | gzip -c > "$tmp/catbad.cif.gz"; cat "$tmp/m2bad.gz" >> "$tmp/catbad.cif.gz"
check "damaged second gzip member" "open-failed" "compressed data is damaged" "$tmp/catbad.cif.gz"

gzip -c "$tmp/m.cif.gz" > "$tmp/twice.cif.gz"
check "gzip inside gzip" "open-failed" "decompressed data is gzip-compressed" "$tmp/twice.cif.gz"
printf '\203\244data\300\250encoder' | gzip -c > "$tmp/bcif.gz"
check "gzip-compressed BinaryCIF" "open-failed" "decompressed data is BinaryCIF" "$tmp/bcif.gz"
: | gzip -c > "$tmp/empty.cif.gz"
check "empty .gz" "open-failed" "No atom coordinates found" "$tmp/empty.cif.gz"

if command -v mkfifo > /dev/null; then
    mkfifo "$tmp/fifo.cif.gz"
    gzip -c $data/damaged_missing_value.cif > "$tmp/fifo.cif.gz" &
    writer=$!
    check ".gz from a named pipe" "ok" - "$tmp/fifo.cif.gz"
    kill $writer 2>/dev/null; wait $writer 2>/dev/null    # if nothing read the pipe
fi

# writing .gz: gzip magic, and it reads back like the plain output
f=$data/chains.cif
limited "$driver" -q -w "$tmp/w.cif" "$f" > /dev/null 2>&1
limited "$driver" -q -w "$tmp/w.cif.gz" "$f" > /dev/null 2>&1
magic=$(od -An -tx1 -N2 "$tmp/w.cif.gz" | tr -d ' \n')
dump "$tmp/w.cif" > "$tmp/plain"
dump "$tmp/w.cif.gz" > "$tmp/gz"
if [ "$magic" = "1f8b" ] && grep -q '^RESULT ok natoms=10' "$tmp/gz" && same "$tmp/plain" "$tmp/gz" &&
   gzip -dc "$tmp/w.cif.gz" | cmp -s - "$tmp/w.cif"; then
    ok "writes .cif.gz (gzip -d gives the same text as the plain output)"
else
    bad "writing .cif.gz"
fi

rm -rf "$tmp"
exit $fail
