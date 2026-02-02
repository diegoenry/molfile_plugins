// GROMACS H5MD Browser - Web Application JavaScript

let currentFile = null;
let currentFrame = 0;
let maxFrame = 0;

// Mol* viewer state
let molstarViewer = null;
let viewerFrame = 0;
let viewerMaxFrame = 0;
let viewerPlaying = false;
let viewerPlayInterval = null;
let viewerInitialized = false;
let moleculeTypes = [];
let selectedMolTypes = [];
let trajectoryLoaded = false;
let loadMode = 'single';  // 'single' or 'trajectory'

// Theme management
function initTheme() {
    const savedTheme = localStorage.getItem('h5md-theme') || 'dark';
    document.documentElement.setAttribute('data-theme', savedTheme);
    const themeSelect = document.getElementById('theme-select');
    if (themeSelect) {
        themeSelect.value = savedTheme;
    }
}

function setupThemeSwitcher() {
    const themeSelect = document.getElementById('theme-select');
    if (themeSelect) {
        themeSelect.addEventListener('change', (e) => {
            const theme = e.target.value;
            document.documentElement.setAttribute('data-theme', theme);
            localStorage.setItem('h5md-theme', theme);
        });
    }
}

// Folder browser state
let currentBrowsePath = '';

// Initialize application
document.addEventListener('DOMContentLoaded', () => {
    initTheme();
    setupThemeSwitcher();
    loadCurrentDirectory();
    loadFileList();
    setupTabs();
    setupFrameControls();
    setupViewerControls();
    setupFolderBrowser();
});

// Load and display current directory
async function loadCurrentDirectory() {
    try {
        const response = await fetch('/api/directory');
        const data = await response.json();
        const dirPath = document.getElementById('directory-path');
        if (dirPath) {
            dirPath.textContent = data.path;
            dirPath.title = data.path;
        }
    } catch (error) {
        console.error('Failed to load directory:', error);
    }
}

// Setup folder browser
function setupFolderBrowser() {
    const browseBtn = document.getElementById('browse-folder-btn');
    const currentDir = document.getElementById('current-directory');
    const modal = document.getElementById('folder-modal');
    const closeBtn = document.getElementById('close-folder-modal');
    const upBtn = document.getElementById('folder-up-btn');
    const goBtn = document.getElementById('folder-go-btn');
    const pathInput = document.getElementById('folder-path-input');
    const selectBtn = document.getElementById('folder-select-btn');

    // Open modal
    const openModal = async () => {
        modal.style.display = 'flex';
        // Start browsing from current directory
        const response = await fetch('/api/directory');
        const data = await response.json();
        currentBrowsePath = data.path;
        browseTo(currentBrowsePath);
    };

    if (browseBtn) browseBtn.addEventListener('click', openModal);
    if (currentDir) currentDir.addEventListener('click', openModal);

    // Close modal
    if (closeBtn) {
        closeBtn.addEventListener('click', () => {
            modal.style.display = 'none';
        });
    }

    // Close on backdrop click
    if (modal) {
        modal.addEventListener('click', (e) => {
            if (e.target === modal) {
                modal.style.display = 'none';
            }
        });
    }

    // Navigate up
    if (upBtn) {
        upBtn.addEventListener('click', async () => {
            const response = await fetch(`/api/browse?path=${encodeURIComponent(currentBrowsePath)}`);
            const data = await response.json();
            if (data.parent) {
                browseTo(data.parent);
            }
        });
    }

    // Go to path
    if (goBtn && pathInput) {
        const goToPath = () => {
            const path = pathInput.value.trim();
            if (path) {
                browseTo(path);
            }
        };
        goBtn.addEventListener('click', goToPath);
        pathInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') goToPath();
        });
    }

    // Select folder
    if (selectBtn) {
        selectBtn.addEventListener('click', async () => {
            await selectFolder(currentBrowsePath);
            modal.style.display = 'none';
        });
    }
}

// Browse to a directory
async function browseTo(path) {
    const folderList = document.getElementById('folder-list');
    const pathInput = document.getElementById('folder-path-input');
    const h5mdCount = document.getElementById('folder-h5md-count');

    folderList.innerHTML = '<p class="loading">Loading...</p>';

    try {
        const response = await fetch(`/api/browse?path=${encodeURIComponent(path)}`);
        const data = await response.json();

        if (data.error) {
            folderList.innerHTML = `<p class="placeholder">${data.error}</p>`;
            return;
        }

        currentBrowsePath = data.current;
        pathInput.value = data.current;

        // Update h5md count
        if (data.h5md_count > 0) {
            h5mdCount.innerHTML = `<strong>${data.h5md_count}</strong> H5MD file(s) in this folder`;
        } else {
            h5mdCount.textContent = 'No H5MD files in this folder';
        }

        // Render folder list
        if (data.directories.length === 0) {
            folderList.innerHTML = '<p class="placeholder">No subfolders</p>';
            return;
        }

        folderList.innerHTML = data.directories.map(dir => `
            <div class="folder-item" data-path="${dir.path}">
                <span class="folder-icon">📁</span>
                <span class="folder-name">${dir.name}</span>
                ${dir.h5md_count > 0
                    ? `<span class="folder-badge">${dir.h5md_count} H5MD</span>`
                    : '<span class="folder-badge empty">-</span>'
                }
            </div>
        `).join('');

        // Add click handlers for folders
        folderList.querySelectorAll('.folder-item').forEach(item => {
            item.addEventListener('dblclick', () => {
                browseTo(item.dataset.path);
            });
            item.addEventListener('click', () => {
                // Single click selects, updates path input
                folderList.querySelectorAll('.folder-item').forEach(i => i.classList.remove('selected'));
                item.classList.add('selected');
            });
        });

    } catch (error) {
        folderList.innerHTML = `<p class="placeholder">Error: ${error.message}</p>`;
    }
}

// Select a folder and update the app
async function selectFolder(path) {
    try {
        const response = await fetch('/api/directory', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ path: path })
        });

        const data = await response.json();

        if (data.error) {
            alert(`Error: ${data.error}`);
            return;
        }

        // Update UI
        const dirPath = document.getElementById('directory-path');
        if (dirPath) {
            dirPath.textContent = data.path;
            dirPath.title = data.path;
        }

        // Reload file list
        currentFile = null;
        document.getElementById('current-file').textContent = 'No file selected';
        await loadFileList();

    } catch (error) {
        alert(`Error: ${error.message}`);
    }
}

// Load file list from server
async function loadFileList() {
    const fileList = document.getElementById('file-list');

    try {
        const response = await fetch('/api/files');
        const files = await response.json();

        if (files.length === 0) {
            fileList.innerHTML = `
                <div class="no-files-message">
                    <p>No H5MD files found</p>
                    <p class="hint">Click 📁 or the folder path above to browse for a different folder</p>
                </div>
            `;
            return;
        }

        fileList.innerHTML = files.map(file => `
            <div class="file-item" data-file="${file.name}">
                <span class="name" title="${file.name}">${file.name}</span>
                <span class="size">${file.size}</span>
            </div>
        `).join('');

        // Add click handlers
        fileList.querySelectorAll('.file-item').forEach(item => {
            item.addEventListener('click', () => selectFile(item.dataset.file));
        });

        // Select first file
        if (files.length > 0) {
            selectFile(files[0].name);
        }
    } catch (error) {
        fileList.innerHTML = `<p class="placeholder">Error loading files: ${error.message}</p>`;
    }
}

// Select a file and load its data
async function selectFile(filename) {
    currentFile = filename;
    currentFrame = 0;
    viewerFrame = 0;

    // Stop playback if running
    if (viewerPlaying) {
        toggleViewerPlay();
    }

    // Update UI
    document.getElementById('current-file').textContent = filename;
    document.querySelectorAll('.file-item').forEach(item => {
        item.classList.toggle('selected', item.dataset.file === filename);
    });

    // Load all tabs
    await Promise.all([
        loadSummary(filename),
        loadTopology(filename),
        loadParticles(filename),
        loadConnectivity(filename),
        loadStructure(filename),
    ]);

    // Load first frame
    loadFrame(filename, 0);

    // Load viewer if tab is active
    if (document.getElementById('viewer3d').classList.contains('active')) {
        loadViewer(filename);
    }
}

// Tab handling
function setupTabs() {
    document.querySelectorAll('.tab').forEach(tab => {
        tab.addEventListener('click', () => {
            const tabId = tab.dataset.tab;

            // Update active tab
            document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
            tab.classList.add('active');

            // Update active pane
            document.querySelectorAll('.tab-pane').forEach(p => p.classList.remove('active'));
            document.getElementById(tabId).classList.add('active');

            // Initialize viewer when 3D View tab is selected
            if (tabId === 'viewer3d' && currentFile) {
                loadViewer(currentFile);
            }
        });
    });
}

// Load summary data
async function loadSummary(filename) {
    const content = document.getElementById('summary-content');
    content.innerHTML = '<p class="loading">Loading summary...</p>';

    try {
        const response = await fetch(`/api/file/${filename}/summary`);
        const data = await response.json();

        let html = '<div class="card-grid">';

        // H5MD Info
        html += `
            <div class="card">
                <h3>H5MD File Information</h3>
                <div class="info-item">
                    <span class="info-label">File:</span>
                    <span class="info-value">${data.filename}</span>
                </div>
                ${data.h5md_info.version ? `
                <div class="info-item">
                    <span class="info-label">H5MD Version:</span>
                    <span class="info-value">${data.h5md_info.version.join('.')}</span>
                </div>` : ''}
                ${data.h5md_info.author ? `
                <div class="info-item">
                    <span class="info-label">Author:</span>
                    <span class="info-value">${data.h5md_info.author}</span>
                </div>` : ''}
                ${data.h5md_info.creator ? `
                <div class="info-item">
                    <span class="info-label">Creator:</span>
                    <span class="info-value">${data.h5md_info.creator.name} ${data.h5md_info.creator.version || ''}</span>
                </div>` : ''}
                ${data.h5md_info.modules ? `
                <div class="info-item">
                    <span class="info-label">Modules:</span>
                    <span class="info-value">${data.h5md_info.modules.join(', ')}</span>
                </div>` : ''}
            </div>
        `;

        // Topology Summary
        if (!data.topology.error) {
            html += `
                <div class="card">
                    <h3>System Topology</h3>
                    ${data.topology.system_name ? `
                    <div class="info-item">
                        <span class="info-label">System:</span>
                        <span class="info-value">${truncate(data.topology.system_name, 50)}</span>
                    </div>` : ''}
                    <div class="info-item">
                        <span class="info-label">Total Atoms:</span>
                        <span class="info-value highlight">${formatNumber(data.topology.total_atoms)}</span>
                    </div>
                    <div class="info-item">
                        <span class="info-label">Molecule Types:</span>
                        <span class="info-value">${data.topology.molecule_block_names?.length || 0}</span>
                    </div>
                    ${data.topology.molecule_block_names?.map((name, i) => `
                    <div class="info-item">
                        <span class="info-label"></span>
                        <span class="info-value">• ${name}: ${formatNumber(data.topology.molecule_block_counts[i])} molecules</span>
                    </div>`).join('') || ''}
                </div>
            `;
        }

        // Trajectory Summary
        if (!data.trajectory.error) {
            const traj = data.trajectory;
            html += `
                <div class="card">
                    <h3>Trajectory Data</h3>
                    ${traj.position ? `
                    <div class="info-item">
                        <span class="info-label">Frames:</span>
                        <span class="info-value highlight">${formatNumber(traj.position.nframes)}</span>
                    </div>
                    <div class="info-item">
                        <span class="info-label">Atoms:</span>
                        <span class="info-value">${formatNumber(traj.position.natoms)}</span>
                    </div>
                    <div class="info-item">
                        <span class="info-label">Position Unit:</span>
                        <span class="info-value">${traj.position.unit || 'nm'}</span>
                    </div>` : ''}
                    ${traj.time ? `
                    <div class="info-item">
                        <span class="info-label">Time Range:</span>
                        <span class="info-value">${traj.time.start.toFixed(2)} - ${traj.time.end.toFixed(2)} ${traj.time.unit || 'ps'}</span>
                    </div>` : ''}
                    <div class="info-item">
                        <span class="info-label">Velocities:</span>
                        <span class="info-value">${traj.velocity ? 'Yes (' + (traj.velocity.unit || 'nm/ps') + ')' : 'No'}</span>
                    </div>
                    <div class="info-item">
                        <span class="info-label">Box:</span>
                        <span class="info-value">${traj.box ? 'Yes (' + (traj.box.unit || 'nm') + ')' : 'No'}</span>
                    </div>
                </div>
            `;

            // Update max frame for frame controls
            maxFrame = (traj.position?.nframes || 1) - 1;
            updateFrameLabel();
        }

        // Connectivity Summary
        html += `
            <div class="card">
                <h3>Connectivity</h3>
                <div class="info-item">
                    <span class="info-label">Total Bonds:</span>
                    <span class="info-value highlight">${formatNumber(data.connectivity.nbonds || 0)}</span>
                </div>
            </div>
        `;

        html += '</div>';
        content.innerHTML = html;

    } catch (error) {
        content.innerHTML = `<p class="placeholder">Error loading summary: ${error.message}</p>`;
    }
}

// Load topology data
async function loadTopology(filename) {
    const content = document.getElementById('topology-content');
    content.innerHTML = '<p class="loading">Loading topology...</p>';

    try {
        const response = await fetch(`/api/file/${filename}/topology`);
        const data = await response.json();

        if (data.error) {
            content.innerHTML = `<p class="placeholder">${data.error}</p>`;
            return;
        }

        let html = `
            <div class="card">
                <h3>Molecule Types</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Type</th>
                            <th>Count</th>
                            <th>Atoms/Mol</th>
                            <th>Total Atoms</th>
                            <th>Residues</th>
                            <th>Mass (u)</th>
                            <th>Net Charge</th>
                        </tr>
                    </thead>
                    <tbody>
        `;

        data.molecule_types?.forEach((mol, i) => {
            const count = data.molecule_block_counts[i];
            const totalAtoms = (mol.particle_count || 0) * count;
            html += `
                <tr>
                    <td>${mol.name}</td>
                    <td>${formatNumber(count)}</td>
                    <td>${formatNumber(mol.particle_count || 0)}</td>
                    <td>${formatNumber(totalAtoms)}</td>
                    <td>${mol.residue_count || '-'}</td>
                    <td>${mol.total_mass?.toFixed(1) || '-'}</td>
                    <td>${mol.net_charge?.toFixed(2) || '-'}</td>
                </tr>
            `;
        });

        html += `
                    </tbody>
                </table>
            </div>
        `;

        content.innerHTML = html;

    } catch (error) {
        content.innerHTML = `<p class="placeholder">Error loading topology: ${error.message}</p>`;
    }
}

// Load particles data
async function loadParticles(filename) {
    const content = document.getElementById('particles-content');
    content.innerHTML = '<p class="loading">Loading particles info...</p>';

    try {
        const response = await fetch(`/api/file/${filename}/particles`);
        const data = await response.json();

        if (data.error) {
            content.innerHTML = `<p class="placeholder">${data.error}</p>`;
            return;
        }

        let html = '';

        data.systems?.forEach(sysName => {
            const sys = data[sysName];
            if (!sys) return;

            html += `<div class="card"><h3>Particle System: ${sysName}</h3>`;

            // Box info
            if (sys.box) {
                html += `<h4 style="color: var(--success); margin: 1rem 0 0.5rem;">Box</h4>`;
                if (sys.box.dimension) {
                    html += `<div class="info-item"><span class="info-label">Dimension:</span><span class="info-value">${sys.box.dimension}D</span></div>`;
                }
                if (sys.box.boundary) {
                    html += `<div class="info-item"><span class="info-label">Boundary:</span><span class="info-value">${sys.box.boundary.join(', ')}</span></div>`;
                }
                if (sys.box.edges) {
                    html += `<div class="info-item"><span class="info-label">Edges Dataset:</span><span class="info-value">[${sys.box.edges.shape.join(' x ')}] (${sys.box.edges.dtype})</span></div>`;
                    if (sys.box.edges.unit) {
                        html += `<div class="info-item"><span class="info-label">Unit:</span><span class="info-value">${sys.box.edges.unit}</span></div>`;
                    }
                    if (sys.box.edges.first_frame) {
                        const box = sys.box.edges.first_frame;
                        html += `<div class="info-item"><span class="info-label">First Frame:</span><span class="info-value">[${box[0][0].toFixed(3)}, ${box[1][1].toFixed(3)}, ${box[2][2].toFixed(3)}]</span></div>`;
                    }
                }
            }

            // Position info
            if (sys.position?.value) {
                html += `<h4 style="color: var(--success); margin: 1rem 0 0.5rem;">Position</h4>`;
                html += `<div class="info-item"><span class="info-label">Dataset:</span><span class="info-value">[${sys.position.value.shape.join(' x ')}] (${sys.position.value.dtype})</span></div>`;
                html += `<div class="info-item"><span class="info-label">Frames:</span><span class="info-value">${formatNumber(sys.position.value.nframes)}</span></div>`;
                html += `<div class="info-item"><span class="info-label">Atoms:</span><span class="info-value">${formatNumber(sys.position.value.natoms)}</span></div>`;
                if (sys.position.value.unit) {
                    html += `<div class="info-item"><span class="info-label">Unit:</span><span class="info-value">${sys.position.value.unit}</span></div>`;
                }
            }

            // Velocity info
            if (sys.velocity?.value) {
                html += `<h4 style="color: var(--success); margin: 1rem 0 0.5rem;">Velocity</h4>`;
                html += `<div class="info-item"><span class="info-label">Dataset:</span><span class="info-value">[${sys.velocity.value.shape.join(' x ')}] (${sys.velocity.value.dtype})</span></div>`;
                if (sys.velocity.value.unit) {
                    html += `<div class="info-item"><span class="info-label">Unit:</span><span class="info-value">${sys.velocity.value.unit}</span></div>`;
                }
            }

            html += '</div>';
        });

        content.innerHTML = html || '<p class="placeholder">No particle data available</p>';

    } catch (error) {
        content.innerHTML = `<p class="placeholder">Error loading particles: ${error.message}</p>`;
    }
}

// Load connectivity data
async function loadConnectivity(filename) {
    const content = document.getElementById('connectivity-content');
    content.innerHTML = '<p class="loading">Loading connectivity...</p>';

    try {
        const response = await fetch(`/api/file/${filename}/connectivity`);
        const data = await response.json();

        let html = '<div class="card"><h3>Bond Information</h3>';

        const info = data.info;
        html += `
            <div class="info-item"><span class="info-label">Total Bonds:</span><span class="info-value highlight">${formatNumber(info.nbonds || 0)}</span></div>
        `;

        if (info.bonds_shape) {
            html += `<div class="info-item"><span class="info-label">Dataset Shape:</span><span class="info-value">[${info.bonds_shape.join(' x ')}]</span></div>`;
        }
        if (info.min_atom_index !== undefined) {
            html += `<div class="info-item"><span class="info-label">Atom Index Range:</span><span class="info-value">${info.min_atom_index} - ${info.max_atom_index}</span></div>`;
        }
        if (info.bonded_atoms_count) {
            html += `<div class="info-item"><span class="info-label">Bonded Atoms:</span><span class="info-value">${formatNumber(info.bonded_atoms_count)}</span></div>`;
        }

        html += '</div>';

        // Bond preview table
        if (data.bonds_preview?.length > 0) {
            html += `
                <div class="card">
                    <h3>Bond Preview (first ${data.bonds_preview.length} bonds)</h3>
                    <table>
                        <thead>
                            <tr>
                                <th>Index</th>
                                <th>Atom 1</th>
                                <th>Atom 2</th>
                                <th>Name 1</th>
                                <th>Name 2</th>
                            </tr>
                        </thead>
                        <tbody>
            `;

            data.bonds_preview.forEach(bond => {
                html += `
                    <tr>
                        <td>${bond.index}</td>
                        <td>${bond.atom1_idx}</td>
                        <td>${bond.atom2_idx}</td>
                        <td>${bond.atom1_name}</td>
                        <td>${bond.atom2_name}</td>
                    </tr>
                `;
            });

            html += '</tbody></table></div>';
        }

        content.innerHTML = html;

    } catch (error) {
        content.innerHTML = `<p class="placeholder">Error loading connectivity: ${error.message}</p>`;
    }
}

// Load structure data
async function loadStructure(filename) {
    const content = document.getElementById('structure-content');
    content.innerHTML = '<p class="loading">Loading structure...</p>';

    try {
        const response = await fetch(`/api/file/${filename}/structure`);

        if (!response.ok) {
            throw new Error(`HTTP error: ${response.status}`);
        }

        const structure = await response.json();
        console.log('Structure loaded, paths:', Object.keys(structure).length);

        // Build a nested tree structure from flat paths
        const tree = buildTree(structure);

        // Render the tree with expand/collapse controls
        let html = `
            <div class="card">
                <h3>HDF5 File Structure</h3>
                <div class="tree-controls">
                    <button id="expand-all-btn" class="tree-ctrl-btn">Expand All</button>
                    <button id="collapse-all-btn" class="tree-ctrl-btn">Collapse All</button>
                </div>
                <div class="structure-tree">
                    ${renderTreeNode(filename, tree, true)}
                </div>
            </div>
        `;
        content.innerHTML = html;

        // Setup tree interaction handlers
        setupTreeHandlers();

    } catch (error) {
        console.error('Structure loading error:', error);
        content.innerHTML = `<p class="placeholder">Error loading structure: ${error.message}</p>`;
    }
}

// Build nested tree from flat path structure
function buildTree(structure) {
    const tree = { children: {}, attrs: {}, type: 'group' };

    // Get root attributes
    if (structure['/'] && structure['/'].attrs) {
        tree.attrs = structure['/'].attrs;
    }

    Object.keys(structure).forEach(path => {
        if (path === '/') return;

        const parts = path.split('/').filter(p => p);
        let current = tree;

        parts.forEach((part, index) => {
            if (!current.children[part]) {
                current.children[part] = { children: {}, attrs: {}, type: 'group' };
            }
            if (index === parts.length - 1) {
                // This is the actual node - copy its properties
                const info = structure[path];
                current.children[part].type = info.type;
                current.children[part].attrs = info.attrs || {};
                if (info.type === 'dataset') {
                    current.children[part].shape = info.shape;
                    current.children[part].dtype = info.dtype;
                }
            }
            current = current.children[part];
        });
    });

    return tree;
}

// Render a tree node recursively
function renderTreeNode(name, node, isRoot = false, expanded = true) {
    const hasChildren = Object.keys(node.children).length > 0;
    const hasAttrs = Object.keys(node.attrs).length > 0;
    const isDataset = node.type === 'dataset';
    const canExpand = hasChildren || hasAttrs;

    // Determine icon and classes
    let icon = '';
    let nodeClass = '';
    let expandState = expanded ? 'expanded' : 'collapsed';

    if (isDataset) {
        icon = '<span class="tree-icon dataset-icon">◆</span>';
        nodeClass = 'structure-dataset';
    } else if (hasChildren) {
        icon = `<span class="tree-icon folder-icon">${expanded ? '▼' : '▶'}</span>`;
        nodeClass = 'structure-group';
    } else {
        icon = '<span class="tree-icon empty-icon">○</span>';
        nodeClass = 'structure-group';
    }

    // Build the node label
    let label = `<span class="${nodeClass}">${name}${isDataset ? '' : '/'}</span>`;
    if (isDataset && node.shape) {
        const shape = node.shape.join(' × ');
        label += ` <span class="structure-shape">[${shape}]</span> <span class="structure-dtype">(${node.dtype})</span>`;
    }

    // Start building HTML
    let html = `<div class="tree-node ${isRoot ? 'root' : ''} ${expandState}" data-expanded="${expanded}">`;
    html += `<div class="tree-node-header${canExpand ? ' expandable' : ''}">`;
    html += icon + label;
    html += '</div>';

    // Children container
    if (canExpand) {
        html += `<div class="tree-node-children" style="${expanded ? '' : 'display: none;'}">`;

        // Render attributes first
        if (hasAttrs) {
            html += '<div class="tree-attrs">';
            Object.entries(node.attrs).forEach(([attrName, attrVal]) => {
                let valStr = formatAttrValue(attrVal);
                html += `<div class="tree-attr">`;
                html += `<span class="tree-icon attr-icon">@</span>`;
                html += `<span class="attr-name">${attrName}</span>: `;
                html += `<span class="attr-value">${valStr}</span>`;
                html += `</div>`;
            });
            html += '</div>';
        }

        // Render child nodes (groups first, then datasets, alphabetically)
        const childNames = Object.keys(node.children).sort((a, b) => {
            const aIsGroup = node.children[a].type === 'group';
            const bIsGroup = node.children[b].type === 'group';
            if (aIsGroup !== bIsGroup) return bIsGroup ? 1 : -1;
            return a.localeCompare(b);
        });

        childNames.forEach(childName => {
            // Start with first two levels expanded
            const childExpanded = isRoot;
            html += renderTreeNode(childName, node.children[childName], false, childExpanded);
        });

        html += '</div>';
    }

    html += '</div>';
    return html;
}

// Format attribute value for display
function formatAttrValue(value) {
    if (Array.isArray(value)) {
        if (value.length > 5) {
            return `[${value.slice(0, 3).join(', ')}, ... (${value.length} items)]`;
        }
        return `[${value.join(', ')}]`;
    }
    let str = String(value);
    if (str.length > 50) {
        return str.substring(0, 47) + '...';
    }
    return str;
}

// Setup tree interaction handlers
function setupTreeHandlers() {
    // Toggle individual nodes
    document.querySelectorAll('.tree-node-header.expandable').forEach(header => {
        header.addEventListener('click', (e) => {
            const node = header.closest('.tree-node');
            const children = node.querySelector('.tree-node-children');
            const icon = header.querySelector('.folder-icon');
            const isExpanded = node.dataset.expanded === 'true';

            if (children) {
                if (isExpanded) {
                    children.style.display = 'none';
                    node.dataset.expanded = 'false';
                    node.classList.remove('expanded');
                    node.classList.add('collapsed');
                    if (icon) icon.textContent = '▶';
                } else {
                    children.style.display = 'block';
                    node.dataset.expanded = 'true';
                    node.classList.remove('collapsed');
                    node.classList.add('expanded');
                    if (icon) icon.textContent = '▼';
                }
            }
        });
    });

    // Expand all button
    const expandAllBtn = document.getElementById('expand-all-btn');
    if (expandAllBtn) {
        expandAllBtn.addEventListener('click', () => {
            document.querySelectorAll('.tree-node').forEach(node => {
                const children = node.querySelector('.tree-node-children');
                const icon = node.querySelector('.folder-icon');
                if (children) {
                    children.style.display = 'block';
                    node.dataset.expanded = 'true';
                    node.classList.remove('collapsed');
                    node.classList.add('expanded');
                    if (icon) icon.textContent = '▼';
                }
            });
        });
    }

    // Collapse all button
    const collapseAllBtn = document.getElementById('collapse-all-btn');
    if (collapseAllBtn) {
        collapseAllBtn.addEventListener('click', () => {
            document.querySelectorAll('.tree-node:not(.root)').forEach(node => {
                const children = node.querySelector('.tree-node-children');
                const icon = node.querySelector('.folder-icon');
                if (children) {
                    children.style.display = 'none';
                    node.dataset.expanded = 'false';
                    node.classList.remove('expanded');
                    node.classList.add('collapsed');
                    if (icon) icon.textContent = '▶';
                }
            });
        });
    }
}

// Frame controls
function setupFrameControls() {
    document.getElementById('frame-first').addEventListener('click', () => {
        currentFrame = 0;
        loadFrame(currentFile, currentFrame);
    });

    document.getElementById('frame-prev').addEventListener('click', () => {
        if (currentFrame > 0) {
            currentFrame--;
            loadFrame(currentFile, currentFrame);
        }
    });

    document.getElementById('frame-next').addEventListener('click', () => {
        if (currentFrame < maxFrame) {
            currentFrame++;
            loadFrame(currentFile, currentFrame);
        }
    });

    document.getElementById('frame-last').addEventListener('click', () => {
        currentFrame = maxFrame;
        loadFrame(currentFile, currentFrame);
    });
}

function updateFrameLabel() {
    document.getElementById('frame-label').textContent = `Frame: ${currentFrame} / ${maxFrame}`;
}

async function loadFrame(filename, frameIdx) {
    if (!filename) return;

    const content = document.getElementById('frames-content');
    updateFrameLabel();

    try {
        const response = await fetch(`/api/file/${filename}/frame/${frameIdx}`);
        const data = await response.json();

        if (data.error) {
            content.innerHTML = `<p class="placeholder">${data.error}</p>`;
            return;
        }

        let html = '<div class="card-grid">';

        html += `<div class="card"><h3>Frame ${frameIdx}</h3>`;

        if (data.time !== undefined) {
            html += `<div class="info-item"><span class="info-label">Time:</span><span class="info-value">${data.time.toFixed(4)} ps</span></div>`;
        }

        if (data.position) {
            html += `
                <h4 style="color: var(--success); margin: 1rem 0 0.5rem;">Position Bounds (nm)</h4>
                <div class="info-item"><span class="info-label">X:</span><span class="info-value">${data.position.min[0].toFixed(3)} to ${data.position.max[0].toFixed(3)}</span></div>
                <div class="info-item"><span class="info-label">Y:</span><span class="info-value">${data.position.min[1].toFixed(3)} to ${data.position.max[1].toFixed(3)}</span></div>
                <div class="info-item"><span class="info-label">Z:</span><span class="info-value">${data.position.min[2].toFixed(3)} to ${data.position.max[2].toFixed(3)}</span></div>
                <div class="info-item"><span class="info-label">Center:</span><span class="info-value">(${data.position.center[0].toFixed(3)}, ${data.position.center[1].toFixed(3)}, ${data.position.center[2].toFixed(3)})</span></div>
            `;
        }

        if (data.box) {
            html += `
                <h4 style="color: var(--success); margin: 1rem 0 0.5rem;">Box Dimensions (nm)</h4>
                <div class="info-item"><span class="info-label">X:</span><span class="info-value">${data.box[0][0].toFixed(3)}</span></div>
                <div class="info-item"><span class="info-label">Y:</span><span class="info-value">${data.box[1][1].toFixed(3)}</span></div>
                <div class="info-item"><span class="info-label">Z:</span><span class="info-value">${data.box[2][2].toFixed(3)}</span></div>
            `;
        }

        html += '</div></div>';
        content.innerHTML = html;

    } catch (error) {
        content.innerHTML = `<p class="placeholder">Error loading frame: ${error.message}</p>`;
    }
}

// Utility functions
function formatNumber(num) {
    return num?.toLocaleString() || '0';
}

function truncate(str, maxLen) {
    if (!str) return '';
    return str.length > maxLen ? str.substring(0, maxLen - 3) + '...' : str;
}

// ==================== Mol* Viewer Functions ====================

function setupViewerControls() {
    document.getElementById('viewer-first').addEventListener('click', () => {
        viewerFrame = 0;
        loadViewerFrame();
    });

    document.getElementById('viewer-prev').addEventListener('click', () => {
        if (viewerFrame > 0) {
            viewerFrame--;
            loadViewerFrame();
        }
    });

    document.getElementById('viewer-next').addEventListener('click', () => {
        if (viewerFrame < viewerMaxFrame) {
            viewerFrame++;
            loadViewerFrame();
        }
    });

    document.getElementById('viewer-last').addEventListener('click', () => {
        viewerFrame = viewerMaxFrame;
        loadViewerFrame();
    });

    document.getElementById('viewer-play').addEventListener('click', toggleViewerPlay);

    document.getElementById('viewer-reset').addEventListener('click', () => {
        if (molstarViewer) {
            molstarViewer.plugin.managers.camera.reset();
        }
    });

    document.getElementById('viewer-max-atoms').addEventListener('change', () => {
        if (currentFile && selectedMolTypes.length > 0) {
            loadViewerFrame();
        }
    });

    // Molecule selection controls
    document.getElementById('mol-select-all').addEventListener('click', () => {
        document.querySelectorAll('.mol-type-checkbox').forEach(cb => cb.checked = true);
        updateSelectedAtomCount();
    });

    document.getElementById('mol-select-none').addEventListener('click', () => {
        document.querySelectorAll('.mol-type-checkbox').forEach(cb => cb.checked = false);
        updateSelectedAtomCount();
    });

    document.getElementById('mol-select-protein').addEventListener('click', () => {
        document.querySelectorAll('.mol-type-checkbox').forEach(cb => {
            const molType = moleculeTypes.find(m => m.name === cb.dataset.molType);
            cb.checked = molType && !molType.is_solvent && !molType.is_ion;
        });
        updateSelectedAtomCount();
    });

    document.getElementById('mol-load-btn').addEventListener('click', loadSelectedMolecules);

    document.getElementById('viewer-change-selection').addEventListener('click', showMoleculeSelector);

    // Trajectory mode selection
    document.querySelectorAll('input[name="load-mode"]').forEach(radio => {
        radio.addEventListener('change', (e) => {
            loadMode = e.target.value;
            const detailPanel = document.getElementById('traj-options-detail');
            detailPanel.style.display = loadMode === 'trajectory' ? 'flex' : 'none';
            updateTrajectoryFrameCount();
        });
    });

    document.getElementById('traj-step').addEventListener('change', updateTrajectoryFrameCount);
}

function updateTrajectoryFrameCount() {
    const step = parseInt(document.getElementById('traj-step').value);
    const frameCount = Math.ceil((viewerMaxFrame + 1) / step);
    document.getElementById('traj-frame-count').textContent = `(${frameCount} frames will be loaded)`;
}

function showMoleculeSelector() {
    document.getElementById('mol-selector-panel').style.display = 'flex';
    document.getElementById('viewer-panel').style.display = 'none';

    // Stop playback if running
    if (viewerPlaying) {
        toggleViewerPlay();
    }

    // Reset trajectory loaded flag
    trajectoryLoaded = false;

    // Update trajectory frame count
    updateTrajectoryFrameCount();
}

function showViewerPanel() {
    document.getElementById('mol-selector-panel').style.display = 'none';
    document.getElementById('viewer-panel').style.display = 'flex';
}

function updateSelectedAtomCount() {
    let totalAtoms = 0;
    document.querySelectorAll('.mol-type-checkbox:checked').forEach(cb => {
        const molType = moleculeTypes.find(m => m.name === cb.dataset.molType);
        if (molType) {
            totalAtoms += molType.total_atoms;
        }
    });
    document.getElementById('mol-selected-atoms').textContent = `Selected atoms: ${formatNumber(totalAtoms)}`;
}

function loadSelectedMolecules() {
    // Get selected molecule types
    selectedMolTypes = [];
    document.querySelectorAll('.mol-type-checkbox:checked').forEach(cb => {
        selectedMolTypes.push(cb.dataset.molType);
    });

    if (selectedMolTypes.length === 0) {
        alert('Please select at least one molecule type');
        return;
    }

    showViewerPanel();

    if (loadMode === 'trajectory') {
        loadTrajectory();
    } else {
        trajectoryLoaded = false;
        loadViewerWithSelection();
    }
}

function toggleViewerPlay() {
    viewerPlaying = !viewerPlaying;
    const playBtn = document.getElementById('viewer-play');

    if (viewerPlaying) {
        playBtn.innerHTML = '&#10074;&#10074;'; // Pause symbol
        playBtn.classList.add('playing');

        if (trajectoryLoaded && molstarViewer) {
            // Use Mol*'s animation for trajectory
            try {
                const plugin = molstarViewer.plugin;
                // Try to start Mol*'s animation
                const animState = plugin.managers.animation;
                if (animState) {
                    plugin.managers.animation.play({ name: 'model-loop' }, { speed: 5 });
                }
            } catch (error) {
                console.log('Mol* animation not available, using manual animation');
                // Fall back to manual animation
                startManualAnimation();
            }
        } else {
            startManualAnimation();
        }
    } else {
        playBtn.innerHTML = '&#9654;'; // Play symbol
        playBtn.classList.remove('playing');

        if (trajectoryLoaded && molstarViewer) {
            try {
                molstarViewer.plugin.managers.animation.stop();
            } catch (error) {
                // Ignore
            }
        }

        if (viewerPlayInterval) {
            clearInterval(viewerPlayInterval);
            viewerPlayInterval = null;
        }
    }
}

function startManualAnimation() {
    viewerPlayInterval = setInterval(() => {
        viewerFrame++;
        if (viewerFrame > viewerMaxFrame) {
            viewerFrame = 0;
        }
        loadViewerFrame();
    }, 200); // 5 fps
}

function updateViewerFrameLabel() {
    document.getElementById('viewer-frame-label').textContent = `Frame: ${viewerFrame} / ${viewerMaxFrame}`;
}

async function initMolstarViewer() {
    if (viewerInitialized) return;

    const container = document.getElementById('molstar-viewer');
    if (!container) return;

    try {
        // Initialize Mol* Viewer
        molstarViewer = await molstar.Viewer.create(container, {
            layoutIsExpanded: false,
            layoutShowControls: true,
            layoutShowRemoteState: false,
            layoutShowSequence: false,
            layoutShowLog: false,
            layoutShowLeftPanel: false,
            viewportShowExpand: true,
            viewportShowSelectionMode: false,
            viewportShowAnimation: true,  // Enable animation controls for trajectories
            collapseLeftPanel: true,
            collapseRightPanel: true,
        });

        viewerInitialized = true;
        console.log('Mol* viewer initialized');
    } catch (error) {
        console.error('Failed to initialize Mol* viewer:', error);
        container.innerHTML = `<p class="placeholder" style="padding: 2rem;">Failed to initialize 3D viewer: ${error.message}</p>`;
    }
}

async function loadViewer(filename) {
    if (!filename) return;

    // Show selector panel first
    showMoleculeSelector();

    const listContainer = document.getElementById('mol-selector-list');
    listContainer.innerHTML = '<p class="loading">Loading molecule types...</p>';

    // Reset selected types
    selectedMolTypes = [];

    try {
        // Get viewer info including molecule types
        const response = await fetch(`/api/file/${filename}/viewer-info`);

        if (!response.ok) {
            throw new Error(`HTTP error: ${response.status}`);
        }

        const info = await response.json();

        if (info.error) {
            listContainer.innerHTML = `<p class="placeholder">Error: ${info.error}</p>`;
            return;
        }

        viewerMaxFrame = Math.max(0, (info.nframes || 1) - 1);
        viewerFrame = 0;
        moleculeTypes = info.molecule_types || [];

        console.log('Loaded molecule types:', moleculeTypes);

        // Populate molecule selector
        populateMoleculeSelector();

    } catch (error) {
        console.error('Failed to load viewer info:', error);
        listContainer.innerHTML = `<p class="placeholder">Failed to load molecule types: ${error.message}</p>`;
    }
}

function populateMoleculeSelector() {
    const listContainer = document.getElementById('mol-selector-list');

    if (moleculeTypes.length === 0) {
        listContainer.innerHTML = '<p class="placeholder">No molecule type information available</p>';
        return;
    }

    let html = '';
    moleculeTypes.forEach(mol => {
        const typeClass = mol.is_solvent ? 'solvent' : (mol.is_ion ? 'ion' : '');
        const typeLabel = mol.is_solvent ? ' (solvent)' : (mol.is_ion ? ' (ion)' : '');

        html += `
            <div class="mol-type-item">
                <input type="checkbox" class="mol-type-checkbox"
                       data-mol-type="${mol.name}"
                       ${mol.default_selected ? 'checked' : ''}>
                <div class="mol-type-info">
                    <div class="mol-type-name ${typeClass}">${mol.name}${typeLabel}</div>
                    <div class="mol-type-details">${formatNumber(mol.count)} molecules × ${mol.atoms_per_molecule} atoms</div>
                </div>
                <div class="mol-type-atoms">${formatNumber(mol.total_atoms)} atoms</div>
            </div>
        `;
    });

    listContainer.innerHTML = html;

    // Add change handlers
    listContainer.querySelectorAll('.mol-type-checkbox').forEach(cb => {
        cb.addEventListener('change', updateSelectedAtomCount);
    });

    // Update initial count
    updateSelectedAtomCount();

    // Update trajectory frame count
    updateTrajectoryFrameCount();
}

async function loadViewerWithSelection() {
    // Initialize viewer if needed
    await initMolstarViewer();
    if (!molstarViewer) return;

    // Show custom frame controls for single frame mode
    document.querySelector('.viewer-frame-controls').style.display = 'flex';

    updateViewerFrameLabel();

    // Load first frame with selection
    await loadViewerFrame();
}

async function loadTrajectory() {
    if (!currentFile) return;

    // Initialize viewer if needed
    await initMolstarViewer();
    if (!molstarViewer) return;

    const maxAtoms = document.getElementById('viewer-max-atoms').value;
    const step = document.getElementById('traj-step').value;

    // Update UI to show loading
    document.getElementById('viewer-atoms').textContent = 'Loading trajectory...';
    document.getElementById('viewer-time').textContent = '';

    // Hide custom frame controls for trajectory mode (use Mol*'s controls)
    document.querySelector('.viewer-frame-controls').style.display = 'none';

    try {
        // Clear previous structure
        await molstarViewer.plugin.clear();

        // Build URL with selected molecule types and trajectory params
        const molTypesParam = encodeURIComponent(selectedMolTypes.join(','));
        const trajUrl = `/api/file/${currentFile}/trajectory-pdb?max_atoms=${maxAtoms}&mol_types=${molTypesParam}&step=${step}`;

        console.log('Loading trajectory from:', trajUrl);

        await molstarViewer.loadStructureFromUrl(trajUrl, 'pdb', false, {
            representationParams: {
                theme: {
                    globalName: 'element-symbol',
                    carbonColor: { name: 'element-symbol' }
                }
            }
        });

        trajectoryLoaded = true;

        // Update atom count based on selection
        let selectedAtoms = 0;
        selectedMolTypes.forEach(name => {
            const mol = moleculeTypes.find(m => m.name === name);
            if (mol) selectedAtoms += mol.total_atoms;
        });
        document.getElementById('viewer-atoms').textContent = `Atoms: ${formatNumber(Math.min(selectedAtoms, parseInt(maxAtoms)))}`;

        const loadedFrames = Math.ceil((viewerMaxFrame + 1) / parseInt(step));
        document.getElementById('viewer-time').textContent = `Trajectory: ${loadedFrames} frames (use Mol* controls to animate)`;

    } catch (error) {
        console.error('Failed to load trajectory:', error);
        document.getElementById('viewer-atoms').textContent = 'Error loading trajectory';
    }
}

async function loadViewerFrame() {
    if (!currentFile || !molstarViewer) return;
    if (selectedMolTypes.length === 0) return;

    // If trajectory is loaded, use Mol*'s model navigation instead
    if (trajectoryLoaded) {
        try {
            // Navigate to the specific model in the trajectory
            const plugin = molstarViewer.plugin;
            const models = plugin.managers.structure.hierarchy.current.structures[0]?.models || [];
            if (models.length > viewerFrame) {
                // This would require more complex Mol* API usage
                // For now, just update the label
                updateViewerFrameLabel();
            }
        } catch (error) {
            console.error('Error navigating trajectory:', error);
        }
        return;
    }

    const maxAtoms = document.getElementById('viewer-max-atoms').value;
    updateViewerFrameLabel();

    try {
        // Clear previous structure
        await molstarViewer.plugin.clear();

        // Build URL with selected molecule types
        const molTypesParam = encodeURIComponent(selectedMolTypes.join(','));
        const pdbUrl = `/api/file/${currentFile}/pdb/${viewerFrame}?max_atoms=${maxAtoms}&mol_types=${molTypesParam}`;

        await molstarViewer.loadStructureFromUrl(pdbUrl, 'pdb', false, {
            representationParams: {
                theme: {
                    globalName: 'element-symbol',
                    carbonColor: { name: 'element-symbol' }
                }
            }
        });

        // Update atom count based on selection
        let selectedAtoms = 0;
        selectedMolTypes.forEach(name => {
            const mol = moleculeTypes.find(m => m.name === name);
            if (mol) selectedAtoms += mol.total_atoms;
        });
        document.getElementById('viewer-atoms').textContent = `Atoms: ${formatNumber(Math.min(selectedAtoms, parseInt(maxAtoms)))}`;

        // Get frame time info
        const frameResponse = await fetch(`/api/file/${currentFile}/frame/${viewerFrame}`);
        const frameData = await frameResponse.json();

        if (frameData.time !== undefined) {
            document.getElementById('viewer-time').textContent = `Time: ${frameData.time.toFixed(2)} ps`;
        }

    } catch (error) {
        console.error('Failed to load frame:', error);
    }
}
