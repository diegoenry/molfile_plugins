# mmcif — merged into ../pdbx

The `mmcif` reader that lived here is now part of the single mmCIF plugin in
[`../pdbx`](../pdbx/README.md). That plugin registers under both names,
`pdbx` and `mmcif`, so scripts that say `mol new file.cif type mmcif` keep
working, and its higher version (2.0) makes VMD prefer it over an old
`mmcifplugin.so` that is still installed.

Everything this reader did is there, including chain selection with
`file.cif?chain=A,B` or `MMCIF_CHAIN`. Some results differ, because they now
match VMD's PDB reader (all alternate locations are kept, `segname` is the
author chain ID, `type` is the element symbol). The pdbx README has the
details.

The old sources (`mmcifplugin.cpp`, `cifparse.cpp`) are in the git history,
before the commit that removed them.
