#!/usr/bin/env python3
"""
GROMACS H5MD File Browser - Web Version

A Flask-based web application for browsing H5MD molecular dynamics files.

Usage:
    python app.py [directory] [--port PORT]
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple

import h5py
import numpy as np
from flask import Flask, render_template, jsonify, request, Response

app = Flask(__name__)

# Global configuration
H5MD_DIRECTORY = "."
CURRENT_FILE: Optional[str] = None
CURRENT_READER = None


class H5MDReader:
    """Reader for GROMACS H5MD files."""

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.file: Optional[h5py.File] = None
        self._open()

    def _open(self):
        """Open the H5MD file."""
        self.file = h5py.File(self.filepath, 'r')

    def close(self):
        """Close the file."""
        if self.file:
            self.file.close()
            self.file = None

    def get_h5md_info(self) -> Dict[str, Any]:
        """Get H5MD metadata."""
        info = {}
        if 'h5md' not in self.file:
            return {"error": "Not a valid H5MD file"}

        h5md = self.file['h5md']

        if 'version' in h5md.attrs:
            info['version'] = [int(x) for x in h5md.attrs['version']]

        if 'author' in h5md:
            author = h5md['author']
            if 'name' in author.attrs:
                name = author.attrs['name']
                info['author'] = name.decode('utf-8') if isinstance(name, bytes) else name

        if 'creator' in h5md:
            creator = h5md['creator']
            creator_info = {}
            if 'name' in creator.attrs:
                name = creator.attrs['name']
                creator_info['name'] = name.decode('utf-8') if isinstance(name, bytes) else name
            if 'version' in creator.attrs:
                version = creator.attrs['version']
                creator_info['version'] = version.decode('utf-8') if isinstance(version, bytes) else version
            info['creator'] = creator_info

        if 'modules' in h5md:
            info['modules'] = list(h5md['modules'].keys())

        return info

    def get_gromacs_topology(self) -> Dict[str, Any]:
        """Get GROMACS topology information."""
        topo = {}

        try:
            topo_group = self.file['h5md/modules/gromacs_topology']
        except KeyError:
            return {"error": "No GROMACS topology found"}

        if 'system_name' in topo_group.attrs:
            name = topo_group.attrs['system_name']
            topo['system_name'] = name.decode('utf-8') if isinstance(name, bytes) else name

        if 'version' in topo_group.attrs:
            topo['version'] = [int(x) for x in topo_group.attrs['version']]

        if 'molecule_block_names' in topo_group.attrs:
            names = topo_group.attrs['molecule_block_names']
            if isinstance(names[0], bytes):
                names = [n.decode('utf-8') for n in names]
            topo['molecule_block_names'] = list(names)

        if 'molecule_block_counts' in topo_group.attrs:
            topo['molecule_block_counts'] = [int(x) for x in topo_group.attrs['molecule_block_counts']]

        molecule_types = []
        if 'molecule_block_names' in topo:
            for name in topo['molecule_block_names']:
                if name in topo_group:
                    mol_info = self._get_molecule_type_info(topo_group[name], name)
                    molecule_types.append(mol_info)
        topo['molecule_types'] = molecule_types

        total_atoms = 0
        if 'molecule_block_names' in topo and 'molecule_block_counts' in topo:
            for mol_type, count in zip(topo['molecule_types'], topo['molecule_block_counts']):
                if 'particle_count' in mol_type:
                    total_atoms += mol_type['particle_count'] * count
        topo['total_atoms'] = total_atoms

        return topo

    def _get_molecule_type_info(self, group: h5py.Group, name: str) -> Dict[str, Any]:
        """Get information about a molecule type."""
        info = {'name': name}

        if 'particle_count' in group.attrs:
            info['particle_count'] = int(group.attrs['particle_count'])

        if 'residue_count' in group.attrs:
            info['residue_count'] = int(group.attrs['residue_count'])

        if 'particle_name_table' in group:
            names = group['particle_name_table'][:]
            if isinstance(names[0], bytes):
                names = [n.decode('utf-8') for n in names]
            info['atom_names'] = list(names)

        if 'residue_name_table' in group:
            names = group['residue_name_table'][:]
            if isinstance(names[0], bytes):
                names = [n.decode('utf-8') for n in names]
            info['residue_names'] = list(names)

        if 'mass' in group:
            masses = group['mass'][:]
            info['mass_range'] = [float(np.min(masses)), float(np.max(masses))]
            info['total_mass'] = float(np.sum(masses))

        if 'charge' in group:
            charges = group['charge'][:]
            info['charge_range'] = [float(np.min(charges)), float(np.max(charges))]
            info['net_charge'] = float(np.sum(charges))

        return info

    def get_trajectory_info(self) -> Dict[str, Any]:
        """Get trajectory information."""
        traj = {}

        try:
            particles = self.file['particles']
        except KeyError:
            return {"error": "No particles group found"}

        system_name = list(particles.keys())[0] if particles.keys() else None
        if not system_name:
            return {"error": "No particle system found"}

        traj['system_name'] = system_name
        system = particles[system_name]

        if 'position' in system:
            pos_group = system['position']
            if 'value' in pos_group:
                pos_data = pos_group['value']
                traj['position'] = {
                    'shape': list(pos_data.shape),
                    'dtype': str(pos_data.dtype),
                    'nframes': pos_data.shape[0],
                    'natoms': pos_data.shape[1],
                }
                if 'unit' in pos_data.attrs:
                    unit = pos_data.attrs['unit']
                    traj['position']['unit'] = unit.decode('utf-8') if isinstance(unit, bytes) else unit

            if 'time' in pos_group:
                time_data = pos_group['time'][:]
                traj['time'] = {
                    'start': float(time_data[0]),
                    'end': float(time_data[-1]),
                    'nframes': len(time_data),
                }
                if 'unit' in pos_group['time'].attrs:
                    unit = pos_group['time'].attrs['unit']
                    traj['time']['unit'] = unit.decode('utf-8') if isinstance(unit, bytes) else unit

        if 'velocity' in system:
            vel_group = system['velocity']
            if 'value' in vel_group:
                vel_data = vel_group['value']
                traj['velocity'] = {
                    'shape': list(vel_data.shape),
                    'dtype': str(vel_data.dtype),
                }
                if 'unit' in vel_data.attrs:
                    unit = vel_data.attrs['unit']
                    traj['velocity']['unit'] = unit.decode('utf-8') if isinstance(unit, bytes) else unit

        if 'box' in system:
            box_group = system['box']
            if 'edges' in box_group:
                edges_group = box_group['edges']
                if 'value' in edges_group:
                    edges_data = edges_group['value']
                    traj['box'] = {
                        'shape': list(edges_data.shape),
                        'dtype': str(edges_data.dtype),
                    }
                    first_box = edges_data[0]
                    last_box = edges_data[-1]
                    traj['box']['first_frame'] = first_box.tolist()
                    traj['box']['last_frame'] = last_box.tolist()
                    if 'unit' in edges_data.attrs:
                        unit = edges_data.attrs['unit']
                        traj['box']['unit'] = unit.decode('utf-8') if isinstance(unit, bytes) else unit

        return traj

    def get_particles_info(self) -> Dict[str, Any]:
        """Get detailed particles group information."""
        info = {}

        try:
            particles = self.file['particles']
        except KeyError:
            return {"error": "No particles group found"}

        info['systems'] = list(particles.keys())

        for sys_name in info['systems']:
            system = particles[sys_name]
            sys_info = {}

            if 'box' in system:
                box = system['box']
                box_info = {}

                if 'dimension' in box.attrs:
                    box_info['dimension'] = int(box.attrs['dimension'])
                if 'boundary' in box.attrs:
                    boundary = box.attrs['boundary']
                    if isinstance(boundary[0], bytes):
                        boundary = [b.decode('utf-8') for b in boundary]
                    box_info['boundary'] = list(boundary)

                if 'edges' in box:
                    edges = box['edges']
                    if 'value' in edges:
                        edge_data = edges['value']
                        box_info['edges'] = {
                            'shape': list(edge_data.shape),
                            'dtype': str(edge_data.dtype),
                        }
                        if 'unit' in edge_data.attrs:
                            unit = edge_data.attrs['unit']
                            box_info['edges']['unit'] = unit.decode('utf-8') if isinstance(unit, bytes) else unit
                        box_info['edges']['first_frame'] = edge_data[0].tolist()
                        box_info['edges']['last_frame'] = edge_data[-1].tolist()

                    if 'step' in edges:
                        step_data = edges['step'][:]
                        box_info['step'] = {
                            'first': int(step_data[0]),
                            'last': int(step_data[-1]),
                            'count': len(step_data),
                        }

                    if 'time' in edges:
                        time_data = edges['time'][:]
                        time_info = {
                            'first': float(time_data[0]),
                            'last': float(time_data[-1]),
                            'count': len(time_data),
                        }
                        if 'unit' in edges['time'].attrs:
                            unit = edges['time'].attrs['unit']
                            time_info['unit'] = unit.decode('utf-8') if isinstance(unit, bytes) else unit
                        box_info['time'] = time_info

                sys_info['box'] = box_info

            if 'position' in system:
                pos = system['position']
                pos_info = {}

                if 'value' in pos:
                    pos_data = pos['value']
                    pos_info['value'] = {
                        'shape': list(pos_data.shape),
                        'dtype': str(pos_data.dtype),
                        'nframes': pos_data.shape[0],
                        'natoms': pos_data.shape[1],
                    }
                    if 'unit' in pos_data.attrs:
                        unit = pos_data.attrs['unit']
                        pos_info['value']['unit'] = unit.decode('utf-8') if isinstance(unit, bytes) else unit

                sys_info['position'] = pos_info

            if 'velocity' in system:
                vel = system['velocity']
                vel_info = {}

                if 'value' in vel:
                    vel_data = vel['value']
                    vel_info['value'] = {
                        'shape': list(vel_data.shape),
                        'dtype': str(vel_data.dtype),
                    }
                    if 'unit' in vel_data.attrs:
                        unit = vel_data.attrs['unit']
                        vel_info['value']['unit'] = unit.decode('utf-8') if isinstance(unit, bytes) else unit

                sys_info['velocity'] = vel_info

            info[sys_name] = sys_info

        return info

    def get_connectivity_info(self) -> Dict[str, Any]:
        """Get bond connectivity information."""
        conn = {}

        try:
            connectivity = self.file['connectivity']
        except KeyError:
            return {"nbonds": 0, "note": "No connectivity data"}

        if 'bond_count' in connectivity.attrs:
            conn['nbonds'] = int(connectivity.attrs['bond_count'])

        if 'bonds' in connectivity:
            bonds = connectivity['bonds']
            conn['bonds_shape'] = list(bonds.shape)
            conn['bonds_dtype'] = str(bonds.dtype)

            nbonds = bonds.shape[0]
            preview_count = min(20, nbonds)
            conn['first_bonds'] = bonds[:preview_count].tolist()
            if nbonds > preview_count:
                conn['last_bonds'] = bonds[-5:].tolist()

            all_bonds = bonds[:]
            conn['min_atom_index'] = int(all_bonds.min())
            conn['max_atom_index'] = int(all_bonds.max())
            conn['bonded_atoms_count'] = len(np.unique(all_bonds.flatten()))

        return conn

    def get_bond_preview_with_names(self, num_bonds: int = 30) -> List[Dict[str, Any]]:
        """Get bond preview with atom names."""
        bonds_preview = []

        try:
            connectivity = self.file['connectivity']
            if 'bonds' not in connectivity:
                return bonds_preview

            bonds = connectivity['bonds']
            nbonds = min(num_bonds, bonds.shape[0])
            bond_data = bonds[:nbonds]

            atom_names = self._build_atom_name_map()

            for i, (atom1, atom2) in enumerate(bond_data):
                bond_info = {
                    'index': i,
                    'atom1_idx': int(atom1),
                    'atom2_idx': int(atom2),
                    'atom1_name': atom_names.get(int(atom1), '?'),
                    'atom2_name': atom_names.get(int(atom2), '?'),
                }
                bonds_preview.append(bond_info)

        except KeyError:
            pass

        return bonds_preview

    def _build_atom_name_map(self) -> Dict[int, str]:
        """Build a map of global atom index to atom name."""
        atom_names = {}

        try:
            topo_group = self.file['h5md/modules/gromacs_topology']
        except KeyError:
            return atom_names

        if 'molecule_block_names' not in topo_group.attrs:
            return atom_names

        mol_names = topo_group.attrs['molecule_block_names']
        if isinstance(mol_names[0], bytes):
            mol_names = [n.decode('utf-8') for n in mol_names]

        mol_counts = topo_group.attrs.get('molecule_block_counts', [])

        global_idx = 0
        for mol_name, mol_count in zip(mol_names, mol_counts):
            if mol_name not in topo_group:
                continue

            mol_group = topo_group[mol_name]
            particle_count = mol_group.attrs.get('particle_count', 0)

            if 'particle_name' in mol_group and 'particle_name_table' in mol_group:
                name_indices = mol_group['particle_name'][:]
                name_table = mol_group['particle_name_table'][:]
                if isinstance(name_table[0], bytes):
                    name_table = [n.decode('utf-8') for n in name_table]

                for _ in range(mol_count):
                    for local_idx in range(particle_count):
                        name_idx = name_indices[local_idx]
                        atom_names[global_idx] = name_table[name_idx]
                        global_idx += 1
            else:
                global_idx += particle_count * mol_count

        return atom_names

    def get_frame_data(self, frame_idx: int) -> Dict[str, Any]:
        """Get data for a specific frame."""
        data = {'frame': frame_idx}

        try:
            system = self.file['particles/system']
        except KeyError:
            return {"error": "No particle system found"}

        if 'position' in system and 'value' in system['position']:
            pos = system['position']['value'][frame_idx]
            data['position'] = {
                'min': pos.min(axis=0).tolist(),
                'max': pos.max(axis=0).tolist(),
                'center': pos.mean(axis=0).tolist(),
            }

        if 'box' in system and 'edges' in system['box'] and 'value' in system['box']['edges']:
            box = system['box']['edges']['value'][frame_idx]
            data['box'] = box.tolist()

        if 'position' in system and 'time' in system['position']:
            data['time'] = float(system['position']['time'][frame_idx])

        return data

    def get_hdf5_structure(self) -> Dict[str, Any]:
        """Get the full HDF5 file structure."""
        def visit_item(name, obj):
            path = "/" + name
            info = {"type": "group" if isinstance(obj, h5py.Group) else "dataset"}
            if isinstance(obj, h5py.Dataset):
                info["shape"] = list(obj.shape)
                info["dtype"] = str(obj.dtype)
            if obj.attrs:
                info["attrs"] = {k: self._convert_attr(v) for k, v in obj.attrs.items()}
            structure[path] = info

        structure = {"/": {"type": "group", "attrs": {}}}
        self.file.visititems(visit_item)
        return structure

    def _convert_attr(self, value):
        """Convert HDF5 attribute to Python-friendly format."""
        if isinstance(value, bytes):
            return value.decode('utf-8')
        if isinstance(value, np.ndarray):
            if value.dtype.kind == 'S':
                return [v.decode('utf-8') if isinstance(v, bytes) else v for v in value]
            return value.tolist()
        if isinstance(value, (np.integer, np.floating)):
            return value.item()
        return value

    def get_molecule_types_for_viewer(self) -> List[Dict[str, Any]]:
        """Get molecule type information for the viewer selection UI."""
        mol_types = []

        try:
            topo_group = self.file['h5md/modules/gromacs_topology']
        except KeyError:
            # Try alternate path
            try:
                topo_group = self.file['h5md/gromacs_topology']
            except KeyError:
                return mol_types

        if 'molecule_block_names' not in topo_group.attrs:
            return mol_types

        mol_names = topo_group.attrs['molecule_block_names']
        if len(mol_names) > 0 and isinstance(mol_names[0], bytes):
            mol_names = [n.decode('utf-8') for n in mol_names]
        mol_names = list(mol_names)

        mol_counts = list(topo_group.attrs.get('molecule_block_counts', []))

        for i, (name, count) in enumerate(zip(mol_names, mol_counts)):
            particle_count = 0
            if name in topo_group:
                particle_count = int(topo_group[name].attrs.get('particle_count', 0))

            total_atoms = particle_count * int(count)

            # Guess if this is water/solvent based on common names
            is_solvent = name.upper() in ['SOL', 'WAT', 'HOH', 'TIP3', 'TIP4', 'SPC', 'SPCE', 'TIP3P', 'TIP4P', 'WATER', 'T3C', 'T3P', 'T4P', 'T4E', 'T5P']
            is_ion = name.upper() in ['NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'NA+', 'CL-', 'K+', 'MG2+', 'CA2+', 'ION', 'POT', 'SOD', 'CLA']

            mol_types.append({
                'index': i,
                'name': name,
                'count': int(count),
                'atoms_per_molecule': particle_count,
                'total_atoms': total_atoms,
                'is_solvent': is_solvent,
                'is_ion': is_ion,
                'default_selected': not is_solvent,  # Select non-solvent by default
            })

        return mol_types

    def get_pdb_frame(self, frame_idx: int, max_atoms: int = 50000,
                      selected_mol_types: Optional[List[str]] = None) -> str:
        """Generate PDB format string for a frame.

        Args:
            frame_idx: Frame index to export
            max_atoms: Maximum atoms to include (for performance)
            selected_mol_types: List of molecule type names to include (None = all)

        Returns:
            PDB format string
        """
        lines = []
        lines.append(f"REMARK   Frame {frame_idx} from {Path(self.filepath).name}")

        try:
            system = self.file['particles/system']
        except KeyError:
            return "REMARK   No particle system found\nEND\n"

        if 'position' not in system or 'value' not in system['position']:
            return "REMARK   No position data found\nEND\n"

        # Get all coordinates (convert nm to Angstroms)
        pos_data = system['position']['value']
        total_atoms = pos_data.shape[1]
        all_coords = pos_data[frame_idx, :, :] * 10.0  # nm to Angstrom

        # Get atom info with molecule type filtering
        atom_info, atom_indices = self._build_atom_info_for_pdb_filtered(
            total_atoms, selected_mol_types, max_atoms
        )

        if not atom_indices:
            lines.append("REMARK   No atoms selected")
            lines.append("END")
            return "\n".join(lines)

        lines.append(f"REMARK   {len(atom_indices)} atoms selected")

        # Get box dimensions if available
        if 'box' in system and 'edges' in system['box'] and 'value' in system['box']['edges']:
            box = system['box']['edges']['value'][frame_idx] * 10.0  # nm to Angstrom
            a, b, c = box[0][0], box[1][1], box[2][2]
            lines.append(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}  90.00  90.00  90.00 P 1           1")

        # Build mapping from old indices to new PDB indices
        old_to_new = {old_idx: new_idx + 1 for new_idx, old_idx in enumerate(atom_indices)}

        # Generate ATOM records
        for pdb_idx, file_idx in enumerate(atom_indices):
            info = atom_info.get(file_idx, {})
            atom_name = info.get('name', 'X')[:4]
            res_name = info.get('resname', 'UNK')[:3]
            res_id = info.get('resid', 1) % 10000
            chain = info.get('chain', 'A')
            element = info.get('element', atom_name[0] if atom_name else 'X')[:2]

            x, y, z = all_coords[file_idx]

            atom_name_fmt = f"{atom_name:<4}" if len(atom_name) < 4 else atom_name[:4]
            lines.append(
                f"ATOM  {(pdb_idx + 1) % 100000:5d} {atom_name_fmt:4s} {res_name:3s} {chain:1s}{res_id:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
            )

        # Add connectivity if available and not too many atoms
        if len(atom_indices) <= 10000:
            bonds = self._get_bonds_for_pdb_filtered(atom_indices, old_to_new)
            for new_idx, bonded_new_indices in bonds.items():
                if bonded_new_indices:
                    conect_line = f"CONECT{new_idx:5d}"
                    for bonded in bonded_new_indices[:4]:
                        conect_line += f"{bonded:5d}"
                    lines.append(conect_line)

        lines.append("END")
        return "\n".join(lines)

    def _build_atom_info_for_pdb(self, max_atoms: int) -> Dict[int, Dict[str, Any]]:
        """Build atom info dictionary for PDB generation."""
        atom_info = {}

        try:
            topo_group = self.file['h5md/modules/gromacs_topology']
        except KeyError:
            return atom_info

        if 'molecule_block_names' not in topo_group.attrs:
            return atom_info

        mol_names = topo_group.attrs['molecule_block_names']
        if isinstance(mol_names[0], bytes):
            mol_names = [n.decode('utf-8') for n in mol_names]

        mol_counts = list(topo_group.attrs.get('molecule_block_counts', []))

        # Element mapping from common atom name prefixes
        element_map = {
            'C': 'C', 'N': 'N', 'O': 'O', 'H': 'H', 'S': 'S', 'P': 'P',
            'CA': 'C', 'CB': 'C', 'CG': 'C', 'CD': 'C', 'CE': 'C', 'CZ': 'C',
            'NA': 'NA', 'CL': 'CL', 'MG': 'MG', 'ZN': 'ZN', 'FE': 'FE',
            'OW': 'O', 'HW': 'H', 'MW': 'X',  # Water
        }

        global_idx = 0
        chain_idx = 0
        chains = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        for mol_name, mol_count in zip(mol_names, mol_counts):
            if mol_name not in topo_group:
                continue

            mol_group = topo_group[mol_name]
            particle_count = int(mol_group.attrs.get('particle_count', 0))

            # Get name tables
            name_table = []
            if 'particle_name_table' in mol_group:
                name_table = mol_group['particle_name_table'][:]
                if len(name_table) > 0 and isinstance(name_table[0], bytes):
                    name_table = [n.decode('utf-8') for n in name_table]

            name_indices = []
            if 'particle_name' in mol_group:
                name_indices = mol_group['particle_name'][:]

            resname_table = []
            if 'residue_name_table' in mol_group:
                resname_table = mol_group['residue_name_table'][:]
                if len(resname_table) > 0 and isinstance(resname_table[0], bytes):
                    resname_table = [n.decode('utf-8') for n in resname_table]

            resname_indices = []
            if 'residue_name' in mol_group:
                resname_indices = mol_group['residue_name'][:]

            residue_ids = []
            if 'residue_id' in mol_group:
                residue_ids = mol_group['residue_id'][:]

            chain = chains[chain_idx % len(chains)]

            for mol_copy in range(mol_count):
                for local_idx in range(particle_count):
                    if global_idx >= max_atoms:
                        return atom_info

                    info = {}

                    # Atom name
                    if len(name_indices) > local_idx and len(name_table) > 0:
                        name_idx = name_indices[local_idx]
                        info['name'] = name_table[name_idx] if name_idx < len(name_table) else 'X'
                    else:
                        info['name'] = 'X'

                    # Residue name
                    if len(resname_indices) > local_idx and len(resname_table) > 0:
                        resname_idx = resname_indices[local_idx]
                        info['resname'] = resname_table[resname_idx] if resname_idx < len(resname_table) else 'UNK'
                    else:
                        info['resname'] = mol_name[:3]

                    # Residue ID
                    if len(residue_ids) > local_idx:
                        info['resid'] = int(residue_ids[local_idx])
                    else:
                        info['resid'] = mol_copy + 1

                    info['chain'] = chain

                    # Element from atom name
                    atom_name = info['name'].upper()
                    if atom_name in element_map:
                        info['element'] = element_map[atom_name]
                    elif len(atom_name) > 0:
                        # First letter that's a valid element
                        info['element'] = atom_name[0]
                    else:
                        info['element'] = 'X'

                    atom_info[global_idx] = info
                    global_idx += 1

            chain_idx += 1

        return atom_info

    def _get_bonds_for_pdb(self, max_atoms: int) -> Dict[int, List[int]]:
        """Get bonds dictionary for PDB CONECT records."""
        bonds = {}

        try:
            connectivity = self.file['connectivity']
            if 'bonds' not in connectivity:
                return bonds

            bond_data = connectivity['bonds'][:]

            for atom1, atom2 in bond_data:
                atom1, atom2 = int(atom1), int(atom2)
                if atom1 < max_atoms and atom2 < max_atoms:
                    if atom1 not in bonds:
                        bonds[atom1] = []
                    if atom2 not in bonds:
                        bonds[atom2] = []
                    bonds[atom1].append(atom2)
                    bonds[atom2].append(atom1)

        except KeyError:
            pass

        return bonds

    def _build_atom_info_for_pdb_filtered(
        self, total_atoms: int, selected_mol_types: Optional[List[str]], max_atoms: int
    ) -> Tuple[Dict[int, Dict[str, Any]], List[int]]:
        """Build atom info dictionary with molecule type filtering.

        Returns:
            Tuple of (atom_info dict, list of selected atom indices)
        """
        atom_info = {}
        atom_indices = []

        try:
            topo_group = self.file['h5md/modules/gromacs_topology']
        except KeyError:
            # No topology, return all atoms up to max
            return {i: {} for i in range(min(total_atoms, max_atoms))}, list(range(min(total_atoms, max_atoms)))

        if 'molecule_block_names' not in topo_group.attrs:
            return {i: {} for i in range(min(total_atoms, max_atoms))}, list(range(min(total_atoms, max_atoms)))

        mol_names = topo_group.attrs['molecule_block_names']
        if isinstance(mol_names[0], bytes):
            mol_names = [n.decode('utf-8') for n in mol_names]

        mol_counts = list(topo_group.attrs.get('molecule_block_counts', []))

        # Element mapping
        element_map = {
            'C': 'C', 'N': 'N', 'O': 'O', 'H': 'H', 'S': 'S', 'P': 'P',
            'CA': 'C', 'CB': 'C', 'CG': 'C', 'CD': 'C', 'CE': 'C', 'CZ': 'C',
            'NA': 'NA', 'CL': 'CL', 'MG': 'MG', 'ZN': 'ZN', 'FE': 'FE',
            'OW': 'O', 'HW': 'H', 'MW': 'X',
        }

        global_idx = 0
        chain_idx = 0
        chains = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        for mol_name, mol_count in zip(mol_names, mol_counts):
            # Check if this molecule type is selected
            include_mol = (selected_mol_types is None) or (mol_name in selected_mol_types)

            if mol_name not in topo_group:
                continue

            mol_group = topo_group[mol_name]
            particle_count = int(mol_group.attrs.get('particle_count', 0))

            if not include_mol:
                # Skip this molecule type's atoms
                global_idx += particle_count * mol_count
                continue

            # Get name tables
            name_table = []
            if 'particle_name_table' in mol_group:
                name_table = mol_group['particle_name_table'][:]
                if len(name_table) > 0 and isinstance(name_table[0], bytes):
                    name_table = [n.decode('utf-8') for n in name_table]

            name_indices = []
            if 'particle_name' in mol_group:
                name_indices = mol_group['particle_name'][:]

            resname_table = []
            if 'residue_name_table' in mol_group:
                resname_table = mol_group['residue_name_table'][:]
                if len(resname_table) > 0 and isinstance(resname_table[0], bytes):
                    resname_table = [n.decode('utf-8') for n in resname_table]

            resname_indices = []
            if 'residue_name' in mol_group:
                resname_indices = mol_group['residue_name'][:]

            residue_ids = []
            if 'residue_id' in mol_group:
                residue_ids = mol_group['residue_id'][:]

            chain = chains[chain_idx % len(chains)]

            for mol_copy in range(mol_count):
                for local_idx in range(particle_count):
                    if len(atom_indices) >= max_atoms:
                        return atom_info, atom_indices

                    info = {}

                    # Atom name
                    if len(name_indices) > local_idx and len(name_table) > 0:
                        name_idx = name_indices[local_idx]
                        info['name'] = name_table[name_idx] if name_idx < len(name_table) else 'X'
                    else:
                        info['name'] = 'X'

                    # Residue name
                    if len(resname_indices) > local_idx and len(resname_table) > 0:
                        resname_idx = resname_indices[local_idx]
                        info['resname'] = resname_table[resname_idx] if resname_idx < len(resname_table) else 'UNK'
                    else:
                        info['resname'] = mol_name[:3]

                    # Residue ID
                    if len(residue_ids) > local_idx:
                        info['resid'] = int(residue_ids[local_idx])
                    else:
                        info['resid'] = mol_copy + 1

                    info['chain'] = chain
                    info['mol_type'] = mol_name

                    # Element from atom name
                    atom_upper = info['name'].upper()
                    if atom_upper in element_map:
                        info['element'] = element_map[atom_upper]
                    elif len(atom_upper) > 0:
                        info['element'] = atom_upper[0]
                    else:
                        info['element'] = 'X'

                    atom_info[global_idx] = info
                    atom_indices.append(global_idx)
                    global_idx += 1

            chain_idx += 1

        return atom_info, atom_indices

    def _get_bonds_for_pdb_filtered(
        self, atom_indices: List[int], old_to_new: Dict[int, int]
    ) -> Dict[int, List[int]]:
        """Get bonds dictionary for filtered atoms with remapped indices."""
        bonds = {}
        atom_set = set(atom_indices)

        try:
            connectivity = self.file['connectivity']
            if 'bonds' not in connectivity:
                return bonds

            bond_data = connectivity['bonds'][:]

            for atom1, atom2 in bond_data:
                atom1, atom2 = int(atom1), int(atom2)
                # Only include bonds where both atoms are selected
                if atom1 in atom_set and atom2 in atom_set:
                    new1 = old_to_new[atom1]
                    new2 = old_to_new[atom2]
                    if new1 not in bonds:
                        bonds[new1] = []
                    if new2 not in bonds:
                        bonds[new2] = []
                    bonds[new1].append(new2)
                    bonds[new2].append(new1)

        except KeyError:
            pass

        return bonds

    def get_trajectory_pdb(self, max_atoms: int = 50000,
                           selected_mol_types: Optional[List[str]] = None,
                           start_frame: int = 0,
                           end_frame: Optional[int] = None,
                           step: int = 1) -> str:
        """Generate multi-model PDB format string for trajectory.

        Args:
            max_atoms: Maximum atoms to include (for performance)
            selected_mol_types: List of molecule type names to include (None = all)
            start_frame: First frame to include
            end_frame: Last frame to include (None = all frames)
            step: Frame step (1 = every frame, 2 = every other frame, etc.)

        Returns:
            PDB format string with multiple MODEL records
        """
        lines = []
        lines.append(f"REMARK   Trajectory from {Path(self.filepath).name}")

        try:
            system = self.file['particles/system']
        except KeyError:
            return "REMARK   No particle system found\nEND\n"

        if 'position' not in system or 'value' not in system['position']:
            return "REMARK   No position data found\nEND\n"

        pos_data = system['position']['value']
        total_frames = pos_data.shape[0]
        total_atoms = pos_data.shape[1]

        # Determine frame range
        if end_frame is None or end_frame > total_frames:
            end_frame = total_frames
        frame_indices = list(range(start_frame, end_frame, step))

        lines.append(f"REMARK   {len(frame_indices)} frames, step={step}")

        # Get atom info with molecule type filtering (computed once)
        atom_info, atom_indices = self._build_atom_info_for_pdb_filtered(
            total_atoms, selected_mol_types, max_atoms
        )

        if not atom_indices:
            lines.append("REMARK   No atoms selected")
            lines.append("END")
            return "\n".join(lines)

        lines.append(f"REMARK   {len(atom_indices)} atoms per frame")

        # Build mapping from old indices to new PDB indices
        old_to_new = {old_idx: new_idx + 1 for new_idx, old_idx in enumerate(atom_indices)}

        # Get box dimensions if available (use first frame for CRYST1)
        if 'box' in system and 'edges' in system['box'] and 'value' in system['box']['edges']:
            box = system['box']['edges']['value'][start_frame] * 10.0  # nm to Angstrom
            a, b, c = float(box[0][0]), float(box[1][1]), float(box[2][2])
            lines.append(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}  90.00  90.00  90.00 P 1           1")

        # Get time data if available
        time_data = None
        if 'position' in system and 'time' in system['position']:
            time_data = system['position']['time'][:]

        # Generate MODEL records for each frame
        for model_num, frame_idx in enumerate(frame_indices, 1):
            lines.append(f"MODEL     {model_num:4d}")

            if time_data is not None:
                lines.append(f"REMARK   Time: {float(time_data[frame_idx]):.4f} ps")

            # Get coordinates for this frame (convert nm to Angstroms)
            coords = pos_data[frame_idx, :, :] * 10.0

            # Generate ATOM records
            for pdb_idx, file_idx in enumerate(atom_indices):
                info = atom_info.get(file_idx, {})
                atom_name = info.get('name', 'X')[:4]
                res_name = info.get('resname', 'UNK')[:3]
                res_id = info.get('resid', 1) % 10000
                chain = info.get('chain', 'A')
                element = info.get('element', atom_name[0] if atom_name else 'X')[:2]

                x, y, z = coords[file_idx]

                atom_name_fmt = f"{atom_name:<4}" if len(atom_name) < 4 else atom_name[:4]
                lines.append(
                    f"ATOM  {(pdb_idx + 1) % 100000:5d} {atom_name_fmt:4s} {res_name:3s} {chain:1s}{res_id:4d}    "
                    f"{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}  1.00  0.00          {element:>2s}"
                )

            lines.append("ENDMDL")

        # Add connectivity once at the end (applies to all models)
        if len(atom_indices) <= 10000:
            bonds = self._get_bonds_for_pdb_filtered(atom_indices, old_to_new)
            for new_idx, bonded_new_indices in bonds.items():
                if bonded_new_indices:
                    conect_line = f"CONECT{new_idx:5d}"
                    for bonded in bonded_new_indices[:4]:
                        conect_line += f"{bonded:5d}"
                    lines.append(conect_line)

        lines.append("END")
        return "\n".join(lines)

    def get_trajectory_info_for_viewer(self) -> Dict[str, Any]:
        """Get trajectory info needed for the 3D viewer."""
        info = {'nframes': 0, 'natoms': 0}

        try:
            system = self.file['particles/system']
            if 'position' in system and 'value' in system['position']:
                pos = system['position']['value']
                info['nframes'] = int(pos.shape[0])
                info['natoms'] = int(pos.shape[1])

            if 'box' in system and 'edges' in system['box'] and 'value' in system['box']['edges']:
                box = system['box']['edges']['value'][0] * 10.0  # nm to Angstrom
                info['box'] = [float(box[0][0]), float(box[1][1]), float(box[2][2])]
        except KeyError:
            pass

        return info


def get_reader(filepath: str) -> H5MDReader:
    """Get or create a reader for the given file."""
    global CURRENT_FILE, CURRENT_READER

    if CURRENT_FILE != filepath:
        if CURRENT_READER:
            CURRENT_READER.close()
        CURRENT_READER = H5MDReader(filepath)
        CURRENT_FILE = filepath

    return CURRENT_READER


@app.route('/')
def index():
    """Main page."""
    return render_template('index.html')


@app.route('/api/files')
def list_files():
    """List H5MD files in the directory."""
    files = []
    dir_path = Path(H5MD_DIRECTORY)
    for f in sorted(dir_path.glob("*.h5md")):
        size = f.stat().st_size
        if size > 1024*1024*1024:
            size_str = f"{size/(1024*1024*1024):.1f} GB"
        elif size > 1024*1024:
            size_str = f"{size/(1024*1024):.1f} MB"
        elif size > 1024:
            size_str = f"{size/1024:.1f} KB"
        else:
            size_str = f"{size} B"
        files.append({
            'name': f.name,
            'path': str(f),
            'size': size_str,
            'size_bytes': size
        })
    return jsonify(files)


@app.route('/api/directory')
def get_directory():
    """Get the current H5MD directory."""
    return jsonify({
        'path': H5MD_DIRECTORY,
        'name': Path(H5MD_DIRECTORY).name or H5MD_DIRECTORY
    })


@app.route('/api/directory', methods=['POST'])
def set_directory():
    """Set the H5MD directory."""
    global H5MD_DIRECTORY, CURRENT_FILE, CURRENT_READER

    data = request.get_json()
    new_path = data.get('path', '')

    if not new_path:
        return jsonify({"error": "No path provided"}), 400

    # Resolve and validate path
    new_path = os.path.abspath(os.path.expanduser(new_path))

    if not os.path.isdir(new_path):
        return jsonify({"error": f"Not a valid directory: {new_path}"}), 400

    # Close current file if open
    if CURRENT_READER:
        CURRENT_READER.close()
        CURRENT_READER = None
        CURRENT_FILE = None

    H5MD_DIRECTORY = new_path

    return jsonify({
        'path': H5MD_DIRECTORY,
        'name': Path(H5MD_DIRECTORY).name or H5MD_DIRECTORY
    })


@app.route('/api/browse')
def browse_directory():
    """Browse directories for folder selection."""
    path = request.args.get('path', '')

    if not path:
        # Default to home directory
        path = os.path.expanduser('~')
    else:
        path = os.path.abspath(os.path.expanduser(path))

    if not os.path.isdir(path):
        return jsonify({"error": "Invalid directory"}), 400

    # Get parent directory
    parent = str(Path(path).parent)
    if parent == path:
        parent = None  # At root

    # List subdirectories
    dirs = []
    try:
        for item in sorted(os.listdir(path)):
            item_path = os.path.join(path, item)
            if os.path.isdir(item_path) and not item.startswith('.'):
                # Count h5md files in this directory
                h5md_count = len(list(Path(item_path).glob("*.h5md")))
                dirs.append({
                    'name': item,
                    'path': item_path,
                    'h5md_count': h5md_count
                })
    except PermissionError:
        return jsonify({"error": "Permission denied"}), 403

    # Count h5md files in current directory
    h5md_count = len(list(Path(path).glob("*.h5md")))

    return jsonify({
        'current': path,
        'current_name': Path(path).name or path,
        'parent': parent,
        'directories': dirs,
        'h5md_count': h5md_count
    })


@app.route('/api/file/<path:filename>/summary')
def get_summary(filename):
    """Get file summary."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return jsonify({"error": "File not found"}), 404

    reader = get_reader(str(filepath))
    return jsonify({
        'filename': filename,
        'h5md_info': reader.get_h5md_info(),
        'topology': reader.get_gromacs_topology(),
        'trajectory': reader.get_trajectory_info(),
        'connectivity': reader.get_connectivity_info(),
    })


@app.route('/api/file/<path:filename>/topology')
def get_topology(filename):
    """Get topology details."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return jsonify({"error": "File not found"}), 404

    reader = get_reader(str(filepath))
    return jsonify(reader.get_gromacs_topology())


@app.route('/api/file/<path:filename>/particles')
def get_particles(filename):
    """Get particles information."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return jsonify({"error": "File not found"}), 404

    reader = get_reader(str(filepath))
    return jsonify(reader.get_particles_info())


@app.route('/api/file/<path:filename>/connectivity')
def get_connectivity(filename):
    """Get connectivity information."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return jsonify({"error": "File not found"}), 404

    reader = get_reader(str(filepath))
    return jsonify({
        'info': reader.get_connectivity_info(),
        'bonds_preview': reader.get_bond_preview_with_names(30),
    })


@app.route('/api/file/<path:filename>/structure')
def get_structure(filename):
    """Get HDF5 structure."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return jsonify({"error": "File not found"}), 404

    reader = get_reader(str(filepath))
    return jsonify(reader.get_hdf5_structure())


@app.route('/api/file/<path:filename>/frame/<int:frame_idx>')
def get_frame(filename, frame_idx):
    """Get frame data."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return jsonify({"error": "File not found"}), 404

    reader = get_reader(str(filepath))
    return jsonify(reader.get_frame_data(frame_idx))


@app.route('/api/file/<path:filename>/pdb/<int:frame_idx>')
def get_pdb(filename, frame_idx):
    """Get PDB format structure for a frame."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return Response("REMARK File not found\nEND\n", mimetype='chemical/x-pdb'), 404

    max_atoms = request.args.get('max_atoms', 50000, type=int)
    # Get selected molecule types from query parameter (comma-separated)
    mol_types_param = request.args.get('mol_types', None)
    selected_mol_types = None
    if mol_types_param:
        selected_mol_types = [m.strip() for m in mol_types_param.split(',') if m.strip()]

    reader = get_reader(str(filepath))
    pdb_data = reader.get_pdb_frame(frame_idx, max_atoms=max_atoms, selected_mol_types=selected_mol_types)
    return Response(pdb_data, mimetype='chemical/x-pdb')


@app.route('/api/file/<path:filename>/trajectory-pdb')
def get_trajectory_pdb(filename):
    """Get PDB format trajectory with multiple models (frames)."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return Response("REMARK File not found\nEND\n", mimetype='chemical/x-pdb'), 404

    max_atoms = request.args.get('max_atoms', 50000, type=int)
    start_frame = request.args.get('start', 0, type=int)
    end_frame = request.args.get('end', None, type=int)
    step = request.args.get('step', 1, type=int)

    # Get selected molecule types from query parameter (comma-separated)
    mol_types_param = request.args.get('mol_types', None)
    selected_mol_types = None
    if mol_types_param:
        selected_mol_types = [m.strip() for m in mol_types_param.split(',') if m.strip()]

    reader = get_reader(str(filepath))
    pdb_data = reader.get_trajectory_pdb(
        max_atoms=max_atoms,
        selected_mol_types=selected_mol_types,
        start_frame=start_frame,
        end_frame=end_frame,
        step=step
    )
    return Response(pdb_data, mimetype='chemical/x-pdb')


@app.route('/api/file/<path:filename>/viewer-info')
def get_viewer_info(filename):
    """Get info needed for 3D viewer initialization."""
    filepath = Path(H5MD_DIRECTORY) / filename
    if not filepath.exists():
        return jsonify({"error": "File not found"}), 404

    reader = get_reader(str(filepath))
    info = reader.get_trajectory_info_for_viewer()
    info['molecule_types'] = reader.get_molecule_types_for_viewer()
    return jsonify(info)


def main():
    global H5MD_DIRECTORY

    parser = argparse.ArgumentParser(description='H5MD File Browser Web Server')
    parser.add_argument('directory', nargs='?', default='.', help='Directory containing H5MD files')
    parser.add_argument('--port', type=int, default=5000, help='Port to run server on')
    parser.add_argument('--host', default='127.0.0.1', help='Host to bind to')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')

    args = parser.parse_args()

    H5MD_DIRECTORY = os.path.abspath(args.directory)

    if not os.path.isdir(H5MD_DIRECTORY):
        print(f"Error: {H5MD_DIRECTORY} is not a valid directory")
        sys.exit(1)

    print(f"H5MD Browser Web Server")
    print(f"  Directory: {H5MD_DIRECTORY}")
    print(f"  URL: http://{args.host}:{args.port}")
    print()

    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == '__main__':
    main()
