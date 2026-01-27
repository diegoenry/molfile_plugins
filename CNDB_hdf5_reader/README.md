# OpenMiChroM CNDB Plugin for VMD

A VMD molfile plugin for reading OpenMiChroM chromatin dynamics simulation files in HDF5 format (.cndb).

## Features

- **Structure Reading**: Loads coarse-grained chromatin bead structures
- **Trajectory Support**: Reads multi-frame trajectories with thousands of timesteps
- **Genomic Coordinates**: Stores genomic positions (start/end) in beta/occupancy fields
- **Loop Constraints**: Imports chromatin loop constraints as bonds for visualization
- **HDF5 Backend**: Efficient reading of large trajectory files
- **Unit Conversion**: Automatically converts coordinates from nanometers to angstroms for VMD

## File Format

CNDB files contain chromatin dynamics data organized as:

```
/Header                         - Metadata attributes
/replica{N}_chr{M}/             - Replica groups (one per chromosome/replica)
  ├── genomic_position          - [natoms, 2] genomic coordinates
  ├── loops                     - [nloops, 2] loop constraint pairs
  └── spatial_position/         - Trajectory frames
      ├── 1                     - [natoms, 3] coordinates for frame 1
      ├── 2                     - [natoms, 3] coordinates for frame 2
      └── ...
```

## Installation

### Requirements

- VMD 1.9.3 or later
- HDF5 library (via Homebrew on macOS)
- C compiler (clang/gcc)

### Build

```bash
make
```

### Install to VMD

```bash
make install
```

This copies the plugin to your VMD plugins directory.

## Usage

### Loading in VMD

```tcl
# Load structure only
mol new your_file.cndb type cndb

# Load structure and first 100 frames
mol new your_file.cndb type cndb first 0 last 99 waitfor all

# Load all frames (may take time for large files)
mol new your_file.cndb type cndb waitfor all
```

### Command Line

```bash
vmd -cndb your_file.cndb
```

### Accessing Genomic Data

Genomic coordinates are stored in atom properties:

```tcl
set sel [atomselect top "all"]
$sel get beta       # Genomic start positions
$sel get occupancy  # Genomic end positions
```

### Selecting by Bead Type

```tcl
# Select active chromatin (A compartments)
set active [atomselect top "type A1 A2"]

# Select inactive chromatin (B compartments)
set inactive [atomselect top "type B1 B3"]

# Select specific compartment type
set a1_beads [atomselect top "type A1"]

# Count beads by type
set a1_count [atomselect top "type A1" num]
```

### Visualizing Chromatin Loops

Loop constraints are loaded as bonds:

```tcl
# Show loops as tubes
mol representation Bonds
mol color ColorID 3
mol addrep top
```

### Coloring by Compartment Type

```tcl
# Color by atom type (compartment)
mol color Type
mol representation VDW 5.0
mol addrep top

# Or use custom coloring
color Type A1 red
color Type A2 orange
color Type B1 blue
color Type B3 cyan
color Type NA gray
```

## Atom Properties

Each chromatin bead has the following properties:

- **name**: Bead type from CNDB file (e.g., `A1`, `A2`, `B1`, `B3`, `NA`)
- **type**: Same as name - chromatin compartment type
- **resname**: `CHR` (chromatin)
- **resid**: Bead index (1-based)
- **chain**: `A`
- **segid**: Replica name (truncated to 7 chars)
- **beta**: Genomic start position (bp)
- **occupancy**: Genomic end position (bp)
- **mass**: 1.0
- **radius**: 5.0 (for visualization)
- **atomicnumber**: 0 (not a real atom)

### Chromatin Compartment Types

The plugin reads bead types from the `types` dataset:

- **A1**, **A2**: Active chromatin compartments (euchromatin)
- **B1**, **B3**: Inactive chromatin compartments (heterochromatin)
- **NA**: Not applicable / constitutive heterochromatin

## Example Visualization Scripts

### Basic Beads and Loops

```tcl
# Load CNDB file
mol new Bcell_short.cndb type cndb first 0 last 99

# Color by genomic position
mol color Beta
mol representation VDW 5.0
mol material Opaque
mol addrep top

# Show chromatin loops
mol representation Bonds 2.0
mol color ColorID 3
mol selection "all"
mol addrep top

# Animate
animate goto 0
play
```

### Smooth Tube Representation

```tcl
# Load CNDB file (single frame for tube)
mol new Bcell_short.cndb type cndb first 0 last 0

# Load tube drawing functions
source chromatin_tube.tcl

# Draw smooth tube through all beads (defaults: r=10, seg=12, res=6)
draw_chromatin_tube top

# Draw tube colored by compartment type (rainbow effect)
draw_chromatin_tube top 10.0 12 6 blue Opaque 1

# Draw separate tubes for each compartment type
draw_chromatin_tube_by_type top

# Custom high-resolution tube
draw_chromatin_tube top 12.0 16 8 red Glossy
```

## Testing

Run the test script:

```bash
vmd -dispdev text -e test_cndb_full.tcl
```

## Benchmarking

The tube visualization includes comprehensive performance benchmarking:

```bash
# Run benchmark suite (tests multiple parameter combinations)
vmd -dispdev text -e benchmark_tube.tcl
```

Benchmarking is automatic when using `draw_chromatin_tube` - it reports:
- Time breakdown by operation (spline, RMF, geometry, graphics)
- Percentage distribution
- Triangle count
- Total execution time

See `BENCHMARKING.md` for detailed optimization strategies.

## File Structure

- `cndbplugin.c` - Main plugin source code
- `cndbplugin.so` - Compiled plugin
- `Makefile` - Build configuration
- `chromatin_tube.tcl` - Smooth tube visualization (B-spline based with benchmarking)
- `test_cndb_plugin.tcl` - Basic test script
- `test_cndb_full.tcl` - Full trajectory test
- `test_types.tcl` - Bead type distribution test
- `test_tube_colored.tcl` - Colored tube test
- `demo_all_tubes.tcl` - Demo of all tube visualization modes
- `benchmark_tube.tcl` - Performance benchmark suite
- `visualize_chromatin.tcl` - Interactive visualization
- `visualize_compartments.tcl` - Compartment-specific visualization
- `example_analysis.tcl` - Analysis examples
- `README.md` - This file
- `SUMMARY.md` - Technical implementation details
- `CHANGELOG.md` - Version history
- `BENCHMARKING.md` - Performance optimization guide

## Technical Details

### Performance

The plugin:
- Uses HDF5 for efficient data access
- Sorts frame datasets numerically (not lexicographically)
- Validates loop constraint indices against structure size
- Memory-efficient reading (doesn't load all frames at once)
- Converts coordinates from nanometers (CNDB format) to angstroms (VMD standard)

### Limitations

- Currently loads the first replica/chromosome found in the file
- Loop constraints with out-of-bounds indices are skipped
- Segment ID limited to 7 characters (VMD limitation)

### Multi-Replica Files

For files with multiple replicas (e.g., `replica1_chr10`, `replica2_chr10`), the plugin currently loads only the first replica found. To load other replicas, you can modify the source code or create separate files.

## Troubleshooting

### "Cannot open file"
- Ensure HDF5 library is installed: `brew install hdf5`
- Check that the file is a valid HDF5 file: `h5dump -n yourfile.cndb`

### "No replica groups found"
- Verify the file contains replica groups: `h5ls yourfile.cndb`
- Group names should match pattern: `replica*_chr*`

### Crash on load
- Check that loop indices are within bounds
- Verify HDF5 dataset dimensions match expected format

## Development

### Modifying the Plugin

After editing `cndbplugin.c`:

```bash
make clean
make
make install
```

### Debug Mode

Add debug output by modifying printf statements in the source code.

## Citation

If you use this plugin in your research, please cite:

- VMD: Humphrey, W., Dalke, A. and Schulten, K., "VMD - Visual Molecular Dynamics", J. Molec. Graphics, 1996, 14.1, 33-38.
- OpenMiChroM: [Add appropriate citation]

## License

This plugin is provided as-is for research use.

## Author

Diego Enry Barreto Gomes 

## Version

0.1 - Initial release
