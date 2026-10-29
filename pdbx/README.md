# pdbx — PDBx/mmCIF reader and writer for VMD (v1.0)

A from-scratch rewrite of VMD's `pdbx` molfile plugin. It registers under the
same name (`pdbx`) with a higher version than the built-in plugin (0.14), so VMD
uses it automatically for `.cif` files once it is loaded.

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

| | old pdbx 0.14 | new pdbx 1.0 |
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

For 1crn, 4hhb, 1d3z, 1ejg, 1hzh, 1bna, 2n9k and 6vxx, every atom field VMD
exposes matches what VMD's own PDB reader gives for the `.pdb` version of the
same entry: name, resname, resid, chain, insertion, altloc, element,
occupancy, beta, and coordinates in every frame.

## Build and load

```sh
make                 # builds pdbxplugin.so (auto-detects the VMD headers)
make config          # shows the detected paths
make VMD_PLUGIN_INCLUDE=/path/to/vmd/plugins/include   # if detection fails
```

Load it for one session (or put the line in `~/.vmdrc`):

```tcl
plugin dlopen /path/to/pdbxplugin.so
mol new myfile.cif            ;# picked by extension, or: type pdbx
```

To make it permanent, `make install` copies it into VMD's molfile plugin
directory. VMD prefers the higher version number, so it shadows the built-in
plugin without removing it. To check which one is active:
`plugin info "mol file reader" pdbx p; parray p` shows version 1.0.

For a static VMD build, copy `pdbxplugin.C` over
`plugins/molfile_plugin/src/pdbxplugin.C`. It is a single C++11 file with no
dependencies beyond the molfile headers.

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
- No gzip support (decompress `.cif.gz` first). Compressed files (gzip,
  bzip2, xz, zstd, zip) and BinaryCIF are recognized and rejected with a
  message that says so.
- A single value larger than 256 MiB stops the load with an error. Only an
  unterminated text field can produce one; real values are a few MB at most.
- At most 1,000,000 models per file are considered for frames; rows of
  further models are ignored with a warning.

## Tests

```sh
make test     # 28 fixtures in tests/data: summary lines vs tests/expected.txt (also
              # with 1-7 byte read buffers), write/read round trip, error paths
make asan     # the same driver built with AddressSanitizer + UBSan
make fuzz     # coverage-guided fuzzing (libFuzzer); FUZZ_TIME=300 seconds, FUZZ_JOBS=4
```

- `tests/pdbx_dump.cpp` is a standalone driver that runs the plugin the way
  VMD does. It allocates buffers to the exact atom count so a sanitizer build
  catches any overrun, and it checks that bond indices are in range.
  `-a N` prints atoms, `-w out.cif` writes the file back out, and `-b N`
  starts with an N-byte read buffer so every token straddles a refill.
- `tests/roundtrip.sh` reads a file, writes it, reads it again, and compares
  every atom and frame.
- `tests/messages.sh` checks the error paths: compressed input, a read error
  (a directory), a named pipe (no rewinding possible), and that control
  characters from the file never reach the console.
- `tests/fuzz_pdbx.cpp` is a libFuzzer target. Each input is read twice, with
  the normal buffer and with one of a few bytes, and both readings must agree
  exactly. Every input that loads is also written and read back, and the
  round trip must keep the atom count, frames, coordinates and every field
  the writer stores losslessly. It needs a clang
  that ships libFuzzer (Homebrew `llvm`; Apple's clang does not). Corpus and
  crash files go to `$TMPDIR/pdbx_fuzz`.
- `tests/compare_with_pdb.tcl` runs in VMD and compares `ID.cif` (this plugin)
  against `ID.pdb` (VMD's PDB reader), field by field.
- `tests/vmd_checks.tcl` runs end-to-end checks inside VMD: unknown elements,
  frames, writing, bonds, IHM graphics, label numbering.

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

Chemical component files come from a different URL and can be passed to
`vmd_checks.tcl` as extra arguments. They print their bond orders:

```sh
curl -fL -o ~/pdbx_test/ATP.cif https://files.rcsb.org/ligands/download/ATP.cif
vmd -dispdev text -e tests/vmd_checks.tcl -args \
    $PWD/pdbxplugin.so $PWD/tests/data ~/pdbx_test/ATP.cif
```

Verified on 11 PDB entries (up to the 2.4-million-atom 3j3q), the 28 fixtures,
a 3,300-file mutation corpus, and about 1.8 million coverage-guided fuzzing
runs, all under AddressSanitizer/UBSan. There were no sanitizer reports,
crashes or hangs. gemmi reads the writer's output exactly as this plugin
does, including hostile names. Loading 3j3q in VMD takes about the same time
and memory as the original plugin (2.9 s, 1.0 GB).

Author: Diego Enry Barreto Gomes. University of Illinois Open Source License
(see the repository README).
