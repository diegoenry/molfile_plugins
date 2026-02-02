# GROMACS H5MD Browser

A web-based file viewer and analyzer for H5MD (Hierarchical Data Format for Molecular Dynamics) files produced by GROMACS. Browse, explore, and visualize molecular dynamics trajectory data through a modern interface with interactive 3D structure viewing.

https://github.com/user-attachments/assets/57d6aa2f-08ac-4b7d-8ad0-dae802ebbe2a

## Features

- **File Browser**: Navigate directories and select H5MD files with real-time file listing
- **3D Visualization**: Interactive molecular structure viewing using Mol* (Molstar) with trajectory animation
- **Molecule Filtering**: Smart classification of solvent, ions, and proteins for selective loading
- **Multiple Data Views**:
  - Summary: H5MD metadata, topology overview, trajectory statistics
  - 3D View: Interactive visualization with playback controls
  - Topology: Molecule types, atom counts, residues, masses, charges
  - Particles: Box dimensions, position/velocity data, boundary conditions
  - Connectivity: Bond information with atom index preview
  - Frames: Frame-by-frame navigation with position bounds and timing
  - Structure: HDF5 file hierarchy explorer with expandable tree view
- **PDB Export**: Export single frames or full trajectories to PDB format
- **Theming**: 9 visual themes including Dark, Light, VMD, PyMOL, and Chimera color schemes

## Requirements

- Python 3.6+
- Modern web browser with WebGL support

## Installation

```bash
# Clone or navigate to the repository
cd gromacs_h5md_webBrowser

# Option 1: Use the launcher script (recommended)
./run_server.sh /path/to/h5md/files

# Option 2: Manual setup
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python app.py /path/to/h5md/files
```

## Usage

```bash
# Basic usage
./run_server.sh /path/to/h5md/files

# With custom port
./run_server.sh /path/to/h5md/files --port 8080

# Direct Python execution with all options
python app.py /path/to/h5md/files --port 5000 --host 0.0.0.0 --debug
```

### Command-line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `directory` | Directory containing H5MD files | Current directory |
| `--port` | Port to run server on | 5000 |
| `--host` | Host to bind to | 127.0.0.1 |
| `--debug` | Enable Flask debug mode | Off |

### Getting Started

1. Start the server with the command above
2. Open your browser to `http://localhost:5000`
3. Use the folder icon or click the directory path to browse for H5MD files
4. Select a file from the sidebar to load it
5. Explore data through the available tabs

### 3D Viewer Controls

- **Molecule Selection**: Choose which molecule types to load (solvent filtering available)
- **Trajectory Mode**: View single frames or animate full trajectories
- **Frame Step**: Control which frames to include when loading trajectories
- **Max Atoms**: Set limits to prevent browser overload with large systems
- **Standard Mol* Controls**: Rotate, zoom, and pan the structure

## Project Structure

```
gromacs_h5md_webBrowser/
├── app.py              # Flask backend server
├── requirements.txt    # Python dependencies
├── run_server.sh       # Launcher script
├── templates/
│   └── index.html      # Main HTML template
└── static/
    ├── app.js          # JavaScript frontend
    └── style.css       # CSS styling
```

## Dependencies

**Python packages** (installed automatically):
- `flask>=2.0.0` - Web framework
- `h5py>=3.0.0` - HDF5 file reading
- `numpy>=1.20.0` - Numerical computing

**External resources** (loaded from CDN):
- Mol* (Molstar) - 3D molecular visualization

## API Endpoints

| Endpoint | Description |
|----------|-------------|
| `/api/files` | List H5MD files in current directory |
| `/api/directory` | Get/set working directory |
| `/api/browse` | Browse filesystem directories |
| `/api/file/{filename}/summary` | File summary and metadata |
| `/api/file/{filename}/topology` | Molecular topology data |
| `/api/file/{filename}/particles` | Particle system information |
| `/api/file/{filename}/connectivity` | Bond/connectivity data |
| `/api/file/{filename}/structure` | HDF5 file hierarchy |
| `/api/file/{filename}/frame/{idx}` | Frame-specific data |
| `/api/file/{filename}/pdb/{idx}` | Export frame as PDB |
| `/api/file/{filename}/trajectory-pdb` | Export trajectory as PDB |
| `/api/file/{filename}/viewer-info` | 3D viewer initialization data |

## License

Part of the molfile_plugins project.
