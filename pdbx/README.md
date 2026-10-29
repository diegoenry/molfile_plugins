# pdbx — PDBx/mmCIF reader and writer for VMD (v2.1)

A from-scratch rewrite of VMD's `pdbx` molfile plugin, and the only mmCIF
plugin in this repository: the separate `mmcif` reader that used to live in
`../mmcif` is merged into it. It reads text mmCIF and BinaryCIF, plain or
gzip-compressed. It registers as `pdbx` (files ending in `.cif`, `.mmcif` or
`.bcif`) and also as `mmcif`, so scripts that say `type mmcif` keep
working. Both names have a higher version than the plugins they replace
(built-in pdbx 0.14, mmcif 1.x), so VMD uses this one once it is loaded, even
if an old `mmcifplugin.so` is still installed.

## Why

The original plugin aborts the whole load when **any** atom has an element
symbol that is not in the periodic table (coarse-grained beads, dummy atoms,
`X`, ...):

```
pdbxplugin) encountered 3 bad element indices!
ERROR) molecule_structure: severe error indicated by plugin aborting loading of molecule 0
```

It is also line-based rather than a CIF parser. It finds the atom table by
searching for the text `_atom_site.`, and in every PDB entry that has been
revised, that text first appears in the revision history. So it also fails on
ordinary entries such as 4hhb, 1hzh and 6vxx.

## What it does now

| | old pdbx 0.14 | new pdbx 2.1 |
|---|---|---|
| Unknown element symbol | whole load aborted | atom loads as element `X`; symbol kept in `type` |
| Revised PDB entries (4hhb, 1hzh, 6vxx, ...) | fail | load |
| NMR / multi-model files | all models merged into one frame | one frame per model |
| `resid`, `chain` | `label_seq_id`, `label_asym_id` (waters and ligands get resid 0) | `auth_seq_id`, `auth_asym_id`, same as the PDB-format file |
| `altloc` | model number | `label_alt_id` |
| `segname` | empty | full chain ID (see below) |
| Explicit bonds | `_pdbx_validate_rmsd_bond/angle` | those, plus `_struct_conn` (disulfide, covalent, metal coordination) |
| Unit cell | not read | `_cell` |
| Rows split over lines, quoted values with spaces, text fields | misparsed | handled |
| Rows with missing or extra values | everything after them shifts | each damaged row is matched to its columns by value type, or skipped if that is ambiguous, with a warning |
| Chemical component files (e.g. `ATP.cif`) | fail | load with all bonds and bond orders |
| Writer | malformed output for empty fields or names with spaces | valid, quoted, locale-independent mmCIF, one model per frame |
| `.cif.gz` (as the PDB distributes files) | not read | read and written directly (when built with zlib) |
| Loading only some chains | no | `file.cif?chain=A,B` or `MMCIF_CHAIN` |
| BinaryCIF (`.bcif`, as served by RCSB, PDBe and Mol*) | not read | read, also `.bcif.gz`, exactly like the text file |

For 1crn, 4hhb, 1d3z, 1ejg, 1hzh, 1bna, 2n9k and 6vxx, every atom field VMD
exposes matches what VMD's own PDB reader gives for the `.pdb` version of the
same entry: name, resname, resid, chain, insertion, altloc, element,
occupancy, beta, and coordinates in every frame.

### For users of the former `mmcif` plugin

Everything it did is here: unknown elements, the unit cell, `_struct_conn`
bonds, NMR models as frames, and the chain selection (`?chain=A,B`,
`MMCIF_CHAIN`). Some results differ, because they now match VMD's PDB reader:

| | mmcif 1.x | this plugin |
|---|---|---|
| Alternate locations | only the first (`A`) was kept | all are kept, as in a PDB file (select with `altloc "" A`) |
| `segname` | `label_asym_id` | the full author chain ID (`AA`, `A-2`, ...) |
| `type` | the atom name | the element symbol from the file (keeps bead types such as `Cg`) |
| Mass of an unknown element | 0 | guessed by VMD from the atom name |
| A model with different atoms | ended the trajectory | is skipped; later matching models still load |
| Rows with missing or extra values | the rest of the table shifted | repaired or skipped (see Damaged rows) |
| Bonds to symmetry copies | added between the wrong atoms | skipped |
| Several data blocks, IHM files, ligand files | not read | read |
| BinaryCIF | no | yes |
| Writing | no | yes |

## Build and load

```sh
make                 # builds pdbxplugin.so (auto-detects the VMD headers and zlib)
make config          # shows the detected paths and whether zlib is used
make VMD_PLUGIN_INCLUDE=/path/to/vmd/plugins/include   # if detection fails
make ZLIB=0          # without zlib: .cif.gz files are then refused with a message
```

Load it for one session (or put the line in `~/.vmdrc`):

```tcl
plugin dlopen /path/to/pdbxplugin.so
mol new myfile.cif                    ;# picked by extension, or: type pdbx
mol new 1abc.bcif                     ;# BinaryCIF, also picked by extension
mol new 1abc.cif.gz type pdbx         ;# VMD only looks at the last extension (.gz)
mol new "myfile.cif?chain=A,B" type pdbx
animate write pdbx out.cif.gz         ;# a name ending in .gz is written compressed
```

To make it permanent, `make install` copies it into VMD's molfile plugin
directory. VMD prefers the higher version number, so it shadows the built-in
plugin without removing it. To check which one is active:
`plugin info "mol file reader" pdbx p; parray p` shows version 2.1.

For a static VMD build, copy `pdbxplugin.C` over
`plugins/molfile_plugin/src/pdbxplugin.C`. It is a single C++11 file with no
dependencies beyond the molfile headers. For `.cif.gz` support, compile it
with `-D_USE_ZLIB` and link `-lz`, as VMD's LAMMPS and VTF plugins do.

## Loading only some chains

Put the chain IDs after a `?` at the end of the file name, or set
`MMCIF_CHAIN` (the file name wins when both are given):

```tcl
mol new "4v6x.cif?chain=A,B" type pdbx    ;# or "?chain=A B"
set env(MMCIF_CHAIN) "AA,A-2"; mol new big.cif
```

- IDs are chain IDs as in `segname` (`auth_asym_id`, or `label_asym_id`
  with `VMDPDBXUSELABELIDS`). They are case-sensitive: `a` and `A` are
  different chains in large assemblies. A chain ID longer than `segname`
  holds (7 characters) is selected by its full ID or by the 7 characters
  VMD shows, which selects every chain that shares them.
- The result is exactly the selected atoms of the full load. The first model
  still defines the structure; a later model becomes a frame if it matches on
  the selected chains (so a model that differs only in another chain is a
  frame).
- Bonds to chains that were not loaded are skipped, with a note; other bond
  records that match no loaded atom are reported as unmatched, as without a
  selection. IHM spheres are selected by their `asym_id`. Chemical component
  files have no chains; the selection is ignored with a note.
- A selected chain without atoms (not in the file, or without usable
  coordinates) is reported; if no atom is selected, the load fails and lists
  the chains that have atoms.
- A file whose name really contains `?` is opened as named, and a `?` in a
  directory name is part of the path. Any other option after `?` is an error
  rather than being ignored. IDs are separated by commas or blanks, and `&`
  starts another option, so an ID containing `&` or `?` can be selected only
  with `MMCIF_CHAIN`.
- `MMCIF_CHAIN` applies to every load until it is unset; each filtered load
  says so in the console.
- The whole file is still read (the chain is in every row), but only the
  selected atoms are kept, so memory drops: 3j3q with two of its 1,356
  chains needs about 50 MB instead of 590 MB.

## Field mapping (`_atom_site`)

| VMD | mmCIF | notes |
|---|---|---|
| name | `label_atom_id` | fallback `auth_atom_id` |
| type | `type_symbol` | fallback: the atom name |
| element | `type_symbol` | case-insensitive; `Fe2+`, `D`, `T` accepted; unknown → `X` |
| resname | `label_comp_id` | fallback `auth_comp_id` |
| resid | `auth_seq_id` | fallback `label_seq_id` |
| chain | first character of `auth_asym_id` | the molfile ABI stores one character |
| segname | full `auth_asym_id` | up to 7 characters |
| altloc | `label_alt_id` | |
| insertion | `pdbx_PDB_ins_code` | |
| occupancy / beta | `occupancy` / `B_iso_or_equiv` | |
| charge | `pdbx_formal_charge` | `1`, `-1`, `2+`, `1-` accepted |
| frames | `pdbx_PDB_model_num` | fallback `ihm_model_id` |

Notes:

- **Unknown elements.** The plugin reports them once
  (`Note: 3 atom(s) have element symbols that are not in VMD's periodic table: Cg (2), Dm (1)`).
  It deliberately leaves mass and radius unset, so VMD computes them per atom:
  from the element when it is known, and from the atom name otherwise. This is
  the same fallback VMD uses for any file without element data. The original
  symbol is still selectable, e.g. `type Cg`.
- **Long chain IDs.** Large assemblies use IDs such as `AA` or `A-2`. VMD's
  `chain` field holds only the first character, so the full ID is also stored
  in `segname`. Select those chains with `segname AA`.
- **Models.** The first model defines the structure. A later model becomes a
  frame only if it has the same atoms in the same order. Otherwise it is
  skipped with a warning (common in integrative models with varying
  composition). Models are compared as they are read, so a model that does
  not match uses no memory for its coordinates.
- **Damaged rows.** A row with missing values (or extra ones, such as an
  unquoted name with a blank) no longer shifts the rest of the loop. Each
  damaged row is read from its own line, and the plugin tries the likely
  places for the gap: the end of the row (a cut-off row), one position in
  the middle, or one run of extra values. Each placement is checked against
  the column types (coordinates are numbers, an insertion code is one
  character, ...) and against the complete rows already read. The row is
  used only if one placement fits best; otherwise it is skipped with a
  warning, because a missing atom is better than one at the wrong place.
- **Bonds.** Explicit bonds are added on top of VMD's distance-based bond
  search (`MOLFILE_BONDSSPECIAL`). Links to symmetry-related copies are skipped.
  Alternate locations bond within the same conformer. A partner given only as
  "atom O of chain C" (no residue number) is used only if exactly one residue
  matches, and rows without a connection type are ignored.
- **Label numbering.** Set `VMDPDBXUSELABELIDS=1` before loading to get
  `label_asym_id` / `label_seq_id` instead. In this mode, waters and ligands
  (which have no `label_seq_id`) fall back to `auth_seq_id`.
- **Compressed files.** gzip data is recognized by its content, not the
  file name, and decompressed while it is read, so memory use is the same as
  for the plain file (3j3q.cif.gz loads in 2.2 s, about 0.1 s more). Concatenated
  `.gz` files work, and zero padding after the data is ignored. A truncated
  or damaged file, or other bytes after the data (usually a damaged second
  part), fails with a message saying so rather than loading part of the
  file. Re-reading a damaged table means decompressing again from the
  start, so that work is capped (a few times the file size). The table is
  re-read only if the plugin can also go back to its normal reading; past
  the cap, as for a pipe, it keeps the normal reading.
- **BinaryCIF.** Recognized by its content, so a `.bcif.gz` or a `.bcif`
  with another name works too. The columns the plugin uses are decoded and
  turned back into the text a CIF file holds (a coordinate stored as
  17.047 reads as `17.047`), so everything above applies unchanged and a
  `.bcif` loads exactly like the `.cif` it was made from. The file is read
  whole, so loading needs more memory than text (3j3q in VMD: 1.24 GB
  instead of 1.04 GB) but is as fast (3.2 s; the file is 29 MB instead of
  241 MB). Damaged files fail with a message naming the column and the
  problem; every size and offset is checked, and a file cannot claim far
  more values than its size allows. All seven encodings of the format are
  read (ByteArray, FixedPoint, IntervalQuantization, RunLength, Delta,
  IntegerPacking, StringArray).
- **IHM.** `_ihm_starting_model_coord` is read when there is no `_atom_site`,
  and `_ihm_sphere_obj_site` spheres of the first model become VMD graphics,
  as in the original plugin. For coarse-grained spheres as real atoms, use
  `../ihmcif`.

## Limitations

- Only the first data block that contains coordinates is read.
- Archive mmCIF files do not list bonds inside a ligand. Stretched ones, such
  as heme Fe–N, can be missed by VMD's distance search; a PDB file gets them
  from CONECT records.
- The writer does not write bonds, entity information or `_struct_conn`.
  It writes `segname` as the chain ID when the segname starts with the chain
  letter (this keeps long chain IDs read by this plugin), so a CHARMM
  segment `PROA` on chain `P` is written as chain `PROA`. For an atom
  without an element, the writer puts its type in `type_symbol` only if that
  type reads back as unknown (a bead type such as `Cg`); a type such as `CA`
  is written as `?`, so it never comes back as calcium. Output is plain
  CIF 1.1: characters it does not allow (control characters, non-ASCII
  bytes) are written as `_`, and a value that neither quote style can hold is
  written as a text field.
- gzip is the only compression read (without zlib, not even that). bzip2,
  xz, zstd and zip are recognized and rejected with a message that says so.
- BinaryCIF is read, not written: a name ending in `.bcif` is refused by the
  writer. A BinaryCIF file larger than 1 GiB (decompressed) is refused; the
  largest PDB entries are far below that.
- VMD picks a plugin by the last extension only, so `.cif.gz` needs
  `type pdbx`. The plugin does not claim `.gz` itself, since that would
  route every compressed format to it.
- A single value larger than 256 MiB stops the load with an error. Only an
  unterminated text field can produce one; real values are a few MB at most.
- At most 1,000,000 models per file are considered for frames; rows of
  further models are ignored with a warning.

## Tests

```sh
make test     # 30 fixtures in tests/data: summary lines vs tests/expected.txt (also
              # with 1-7 byte read buffers), write/read round trip, error paths,
              # chain selection, .gz reading and writing, BinaryCIF (tests/bcif)
make asan     # the same driver built with AddressSanitizer + UBSan
make fuzz     # coverage-guided fuzzing (libFuzzer); FUZZ_TIME=300 seconds, FUZZ_JOBS=4
make ZLIB=0 test   # the build without zlib
```

- `tests/pdbx_dump.cpp` is a standalone driver that runs the plugin the way
  VMD does. It allocates buffers to the exact atom count so a sanitizer build
  catches any overrun, and it checks that bond indices are in range.
  `-a N` prints atoms, `-w out.cif` writes the file back out, and `-b N`
  starts with an N-byte read buffer so every token straddles a refill.
  `-z` forbids rewinding compressed input, and `-l` lists the registered
  plugin names. File names may carry options (`file.cif?chain=A`).
- `tests/roundtrip.sh` reads a file, writes it, reads it again, and compares
  every atom and frame.
- `tests/messages.sh` checks the error paths: compressed input, a read error
  (a directory), a named pipe (no rewinding possible), and that control
  characters from the file never reach the console.
- `tests/options.sh` checks the chain selection (file name, `MMCIF_CHAIN`,
  label numbering, long IDs, missing chains, bad options, bond messages),
  that every fixture reads the same from a `.gz` (also with tiny buffers),
  that any rewind budget gives either the plain or the normal reading,
  truncated, damaged, concatenated, padded and nested `.gz` files, `.gz`
  from a pipe, writing `.cif.gz`, and that both plugin names are registered.
  The test scripts unset `MMCIF_CHAIN` and `VMDPDBXUSELABELIDS` first.
- `tests/bcif.sh` checks BinaryCIF: every file in `tests/bcif` must load
  exactly like its text source (also with tiny buffers, from `.gz`, and with
  a chain selection), each damaged `bad_*.bcif` must fail with the message
  listed in `tests/bcif/bad.txt`, every truncation of a file must fail
  cleanly, and the writer must refuse a `.bcif` name.
- `tests/make_bcif.py` made the files in `tests/bcif`, which are committed,
  so the tests need no Python. Run it only to change them:
  `uv run --with mmcif --with msgpack python tests/make_bcif.py`. It converts
  text fixtures with RCSB's py-mmcif (the encoder behind
  `models.rcsb.org/<id>.bcif`), writes `encodings.bcif` and its text twin
  `encodings.cif` using every encoding (also those py-mmcif never uses),
  checks those two with py-mmcif's own decoder, and writes the damaged
  files.
- `tests/fuzz_pdbx.cpp` is a libFuzzer target. Each input is read twice, with
  the normal buffer and with one of a few bytes, and both readings must agree
  exactly. Every input that loads is also written and read back, and the
  round trip must keep the atom count, frames, coordinates and every field
  the writer stores losslessly. Each input is also gzip-compressed (sometimes
  in two members) and must read the same from the `.gz`, and the chain of
  one atom is selected with `?chain=`: that load must hold exactly that
  chain's atoms of the full load. The BinaryCIF files in `tests/bcif` are
  among the seeds, so the fuzzer also mutates BinaryCIF. It needs a clang
  that ships libFuzzer (Homebrew `llvm`; Apple's clang does not). Corpus and
  crash files go to `$TMPDIR/pdbx_fuzz`.
- `tests/compare_with_pdb.tcl` runs in VMD and compares `ID.cif` (this plugin)
  against `ID.pdb` (VMD's PDB reader), field by field.
- `tests/vmd_checks.tcl` runs end-to-end checks inside VMD: unknown elements,
  frames, writing, bonds, IHM graphics, label numbering, the `mmcif` name,
  chain selection, reading and writing `.cif.gz`, and BinaryCIF.

```sh
vmd -dispdev text -e tests/vmd_checks.tcl -args $PWD/pdbxplugin.so $PWD/tests/data
```

### Testing with real PDB entries

None of the test scripts download anything. To check the plugin against real
entries, download them yourself into any folder:

```sh
mkdir -p ~/pdbx_test && cd ~/pdbx_test
for id in 1crn 4hhb 1d3z 1ejg 1hzh 1bna 2n9k 6vxx; do
    curl -fL -o $id.pdb    https://files.rcsb.org/download/$id.pdb
    curl -fL -o $id.cif.gz https://files.rcsb.org/download/$id.cif.gz && gunzip -f $id.cif.gz
done
```

Then compare each `.cif` (this plugin) with its `.pdb` (VMD's PDB reader):

```sh
vmd -dispdev text -e tests/compare_with_pdb.tcl -args \
    $PWD/pdbxplugin.so ~/pdbx_test 1crn 4hhb 1d3z 1ejg 1hzh 1bna 2n9k 6vxx
```

Each entry prints `CMP <id> MATCH`, or `DIFF` followed by the fields that
differ.

The driver also loads any single file, which is useful for entries that are
too large to have a `.pdb` version (e.g. 4v6x, 3j3q):

```sh
make tests/pdbx_dump
./tests/pdbx_dump -q ~/pdbx_test/4v6x.cif
```

BinaryCIF versions of the same entries come from RCSB's model server. Each
must load exactly like the `.cif` (run in bash or zsh):

```sh
cd ~/pdbx_test
for id in 1crn 4hhb 1d3z 1ejg 1hzh 1bna 2n9k 6vxx; do
    curl -fL -o $id.bcif https://models.rcsb.org/$id.bcif
done
cd -                      # back to the pdbx directory
for id in 1crn 4hhb 1d3z 1ejg 1hzh 1bna 2n9k 6vxx; do
    if cmp -s <(./tests/pdbx_dump -c ~/pdbx_test/$id.cif  | grep -v '^pdbxplugin)') \
              <(./tests/pdbx_dump -c ~/pdbx_test/$id.bcif | grep -v '^pdbxplugin)')
    then echo "$id same"; else echo "$id DIFF"; fi
done
```

The server sometimes cuts off a large download (curl then exits with an
error, and the plugin reports the file as truncated). The model server's
query form, `https://models.rcsb.org/v1/<id>/full?encoding=bcif`, is written
by Mol*'s encoder instead; it holds the same atoms in a different order and
leaves out some tables (such as the validation bonds).

Chemical component files come from a different URL and can be passed to
`vmd_checks.tcl` as extra arguments. They print their bond orders:

```sh
curl -fL -o ~/pdbx_test/ATP.cif https://files.rcsb.org/ligands/download/ATP.cif
vmd -dispdev text -e tests/vmd_checks.tcl -args \
    $PWD/pdbxplugin.so $PWD/tests/data ~/pdbx_test/ATP.cif
```

Verified on 11 PDB entries (up to the 2.4-million-atom 3j3q), the 30 fixtures,
a 3,300-file mutation corpus, and about 2.6 million coverage-guided fuzzing
runs under AddressSanitizer/UBSan, the last 218,000 of them on the final code
with BinaryCIF seeds and the gzip and chain-selection checks, with no
sanitizer reports, crashes or hangs. gemmi reads the writer's output, plain
and `.gz`, exactly as this plugin does, including hostile names. Loading 3j3q
in VMD takes about the same time and memory as the original plugin (2.9 s,
1.0 GB), and reading 3J3Q.cif.gz as the PDB distributes it takes the same as
the plain file.

BinaryCIF was checked against both encoders in use: for the 8 entries above,
RCSB's `.bcif` files (py-mmcif) load identically to the `.cif`, field by
field and frame by frame, and the Mol* model server's files hold the same
atoms and coordinates. So does 3j3q converted with py-mmcif. 1,500 randomly
corrupted BinaryCIF files under AddressSanitizer all loaded or failed
cleanly.

Author: Diego Enry Barreto Gomes. University of Illinois Open Source License
(see the repository README).
