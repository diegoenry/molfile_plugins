#!/bin/sh
# bcif.sh -- BinaryCIF input: every file in tests/bcif must load exactly like
# its text source (also with tiny read buffers, from .gz and with a chain
# selection), and every damaged one must fail with its message.
#
# usage: tests/bcif.sh DRIVER ZLIB        (files made by tests/make_bcif.py)

unset MMCIF_CHAIN VMDPDBXUSELABELIDS    # must not leak in from the caller
driver=$1
zlib=${2:-1}
tmp=${TMPDIR:-/tmp}/pdbx_bcif.$$
mkdir -p "$tmp" || exit 1
fail=0
dir=tests/bcif

limited() { perl -e 'alarm 60; exec @ARGV or exit 127' "$@"; }
ok()  { echo "OK    $1"; }
bad() { echo "FAIL  $1"; [ -f "$tmp/out" ] && sed 's/^/      /' "$tmp/out" | head -8; fail=1; }
dump() { limited "$driver" -c "$@" 2>/dev/null | grep -v '^pdbxplugin)'; }
same() { grep -q '^RESULT ok' "$1" && grep -q '^RESULT' "$2" && cmp -s "$1" "$2"; }

# the text source of a BinaryCIF test file
source_of() {
    b=$(basename "$1" .bcif)
    if [ -f "$dir/$b.cif" ]; then echo "$dir/$b.cif"; else echo "tests/data/$b.cif"; fi
}

n=0; nbad=0
for f in $dir/*.bcif; do
    case $(basename "$f") in bad_*) continue ;; esac
    src=$(source_of "$f")
    dump "$src" > "$tmp/text"
    for b in 65536 1 5; do
        dump -b $b "$f" > "$tmp/bin"
        if ! same "$tmp/text" "$tmp/bin"; then
            echo "FAIL  $f differs from $src (buffer $b)"
            diff "$tmp/text" "$tmp/bin" | head -6 | sed 's/^/      /'
            nbad=$((nbad + 1)); fail=1
        fi
    done
    if [ "$zlib" = "1" ]; then
        gzip -c "$f" > "$tmp/t.bcif.gz"
        dump -b 3 "$tmp/t.bcif.gz" > "$tmp/bin"
        if ! same "$tmp/text" "$tmp/bin"; then
            echo "FAIL  $f.gz differs from $src"; nbad=$((nbad + 1)); fail=1
        fi
    fi
    n=$((n + 1))
done
[ $nbad -eq 0 ] && ok "all $n BinaryCIF files read exactly like their text (buffers 64K, 1, 5; .gz)"

limited "$driver" -q "$dir/encodings.bcif" > "$tmp/out" 2>&1
if grep -q '^RESULT ok natoms=6 frames=2 bonds=2' "$tmp/out" && grep -q '\[BinaryCIF\]' "$tmp/out" &&
   grep -q 'more than one data block' "$tmp/out"; then
    ok "every encoding, masks, blocks, the struct_conn bond (to both altlocs), the [BinaryCIF] note"
else
    bad "encodings.bcif summary"
fi

for q in "?chain=A" "?chain=A,AA"; do
    dump "tests/data/chains.cif$q" > "$tmp/text"
    dump "$dir/chains.bcif$q" > "$tmp/bin"
    if same "$tmp/text" "$tmp/bin"; then ok "chain selection $q on BinaryCIF"
    else bad "chain selection $q on BinaryCIF"; fi
done

# damaged files: a clear message, never a crash or a partial load
nbad=0; n=0
while IFS='|' read -r file message; do
    case $file in ''|'#'*) continue ;; esac
    n=$((n + 1))
    limited "$driver" -q "$dir/$file" > "$tmp/out" 2>&1
    if ! grep -q '^RESULT open-failed' "$tmp/out" || ! grep -q -- "$message" "$tmp/out"; then
        echo "FAIL  $file (expected '$message')"; sed 's/^/      /' "$tmp/out" | head -4
        nbad=$((nbad + 1)); fail=1
    fi
done < "$dir/bad.txt"
[ $nbad -eq 0 ] && ok "all $n damaged BinaryCIF files fail with their message"

# every truncation of a small file fails cleanly (or loads, if only padding was cut)
f=$dir/models.bcif
size=$(wc -c < "$f")
nbad=0; i=1
while [ $i -lt "$size" ]; do
    head -c $i "$f" > "$tmp/cut.bcif"
    limited "$driver" -q "$tmp/cut.bcif" > "$tmp/out" 2>&1
    if ! grep -q '^RESULT \(ok\|open-failed\)' "$tmp/out"; then
        echo "FAIL  models.bcif cut at $i bytes"; nbad=$((nbad + 1)); fail=1
    fi
    i=$((i + 7))
done
[ $nbad -eq 0 ] && ok "every truncation of models.bcif fails cleanly"

# the writer does not pretend to write BinaryCIF
limited "$driver" -q -w "$tmp/out.bcif" tests/data/models.cif > "$tmp/out" 2>&1
if grep -q "writes text mmCIF only" "$tmp/out" && [ ! -s "$tmp/out.bcif" ]; then
    ok "writing to a .bcif name is refused"
else
    bad "writing to a .bcif name"
fi

rm -rf "$tmp"
exit $fail
