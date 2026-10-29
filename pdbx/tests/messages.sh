#!/bin/sh
# messages.sh -- error paths: inputs that must fail with a clear message, or
# load without hanging, and console output that must stay clean.
#
# usage: tests/messages.sh DRIVER

unset MMCIF_CHAIN VMDPDBXUSELABELIDS    # must not leak in from the caller
driver=$1
tmp=${TMPDIR:-/tmp}/pdbx_messages.$$
mkdir -p "$tmp" || exit 1
fail=0

# run with a time limit (macOS has no timeout(1) by default)
limited() { perl -e 'alarm 30; exec @ARGV or exit 127' "$@"; }

expect() {          # expect NAME PATTERN FILE
    limited "$driver" -q "$3" > "$tmp/out" 2>&1
    if grep -q "$2" "$tmp/out"; then
        echo "OK    $1"
    else
        echo "FAIL  $1 (expected '$2')"; sed 's/^/      /' "$tmp/out" | head -5; fail=1
    fi
}

# a compression format that is not read is named, not reported as "no
# coordinates" (gzip is read; see options.sh)
printf 'BZh91AY&SY\001\002\003\004' > "$tmp/models.cif"
expect "bzip2 input" "is bzip2-compressed" "$tmp/models.cif"

# a read error is an error, not an empty file
expect "directory as input" "read error" "$tmp"

# non-seekable input: the malformed-loop re-read is not possible, the strict
# reading is kept, and the load still finishes
if command -v mkfifo > /dev/null; then
    mkfifo "$tmp/fifo.cif"
    cat tests/data/nonstandard_elements.cif > "$tmp/fifo.cif" &
    writer=$!
    expect "named pipe" "RESULT ok" "$tmp/fifo.cif"
    kill $writer 2>/dev/null; wait $writer 2>/dev/null    # if nothing read the pipe
fi

# control characters from the file never reach the console
esc=$(printf '\033')
sed "s/ Cg / ${esc}[31mCg /" tests/data/nonstandard_elements.cif > "$tmp/escape.cif"
limited "$driver" -q "$tmp/escape.cif" > "$tmp/out" 2>&1
if grep -q "$esc" "$tmp/out"; then
    echo "FAIL  escape sequences are printed"; fail=1
elif grep -q '?\[31mCg' "$tmp/out"; then
    echo "OK    escape sequences are filtered"
else
    echo "FAIL  escape test: unexpected output"; sed 's/^/      /' "$tmp/out" | head -5; fail=1
fi

rm -rf "$tmp"
exit $fail
