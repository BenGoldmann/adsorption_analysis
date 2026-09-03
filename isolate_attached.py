import sys
import re
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from read_parameters import read_parameters

def write_attached_polymers_lammpstrj():
    """
    Writes a LAMMPS dump file containing only the bottom surface atoms and the polymers attached to the bottom surface.
    """
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python script.py parameters.in")

    parameters_path = Path(sys.argv[1])
    if not parameters_path.is_file():
        raise SystemExit(f"ERROR: Parameters file not found: {sys.argv[1]}")

    params = read_parameters(parameters_path)

    working_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
    infile = working_dir / "dump_nosolvent.lammpstrj"

    if not infile.is_file():
        raise SystemExit(
            f"ERROR: Solvent-removed file does not exist: {infile}. "
            "Use remove_solvent.py to create it first."
        )

    for k in ("SURFACE_ID", "ADDITIVE_NUM", "ADDITIVE_ATOMS_NUM"):
        if k not in params:
            raise SystemExit(f"ERROR: {k} not found in parameters.in (expected line like: {k} = N)")

    try:
        surface_id = [int(i) for i in re.split(r"[,\s]+", params["SURFACE_ID"].strip()) if i]
    except ValueError:
        raise SystemExit("ERROR: SURFACE_ID must be integers.")

    try:
        additive_num = int(params["ADDITIVE_NUM"])
        additive_atoms_num = int(params["ADDITIVE_ATOMS_NUM"])
    except ValueError:
        raise SystemExit("ERROR: ADDITIVE_NUM and ADDITIVE_ATOMS_NUM must be integers.")

    print("Parameters read successfully")

    # Load fresh universe for this operation
    u = mda.Universe(infile, format="LAMMPSDUMP")
    
    # Jump to the specified frame
    u.trajectory[-1]  # Last frame
    
    # Get surface atoms
    seltext = "type " + " ".join(str(i) for i in surface_id)
    surface_atoms_group = u.select_atoms(seltext)
    
    # Get surface z-coordinates
    surface_z = surface_atoms_group.positions[:, 2].copy()
    surface_z.sort()
    
    # Get bottom surface z
    bottom_surface_z_actual = surface_z[len(surface_z) // 2 - 1]
    
    # Select bottom surface atoms (those below the midpoint)
    bottom_surface_atoms = surface_atoms_group.select_atoms(f"prop z < {bottom_surface_z_actual + 2.5}")
    
    # Identify polymers attached to the bottom surface
    attached_polymer_atoms = []
    for i in range(additive_num):
        start = len(surface_atoms_group) + i * additive_atoms_num
        stop = start + additive_atoms_num
        polymer_atoms = u.atoms[start:stop]
        
        # Check if any atom in the polymer is close to bottom surface
        is_attached = np.any(polymer_atoms.positions[:, 2] < bottom_surface_z_actual + 5.0)
        
        if is_attached:
            attached_polymer_atoms.extend(polymer_atoms)
    
    # Combine all atoms to write
    all_atoms = list(bottom_surface_atoms) + attached_polymer_atoms
    
    if not all_atoms:
        print("Warning: No atoms found to write to dump file.")
        return
    
    # Get box dimensions
    dimensions = u.dimensions
    output_file = working_dir / "attached_only_frame.lammpstrj"
    
    # Write LAMMPS dump file
    with open(output_file, 'w') as f:
        # Write header
        f.write("ITEM: TIMESTEP\n")
        f.write(f"{u.trajectory.frame}\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{len(all_atoms)}\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        f.write(f"0.0 {dimensions[0]}\n")
        f.write(f"0.0 {dimensions[1]}\n")
        f.write(f"0.0 {dimensions[2]}\n")
        f.write("ITEM: ATOMS id type x y z\n")
        
        # Write atom data
        for i, atom in enumerate(all_atoms, 1):
            atom_id = atom.index + 1  # LAMMPS uses 1-based indexing
            atom_type = atom.type
            x, y, z = atom.position
            f.write(f"{atom_id} {atom_type} {x:.6f} {y:.6f} {z:.6f}\n")
    
    print(f"Wrote {len(all_atoms)} atoms to {output_file}")
    print(f"  - {len(bottom_surface_atoms)} bottom surface atoms")
    print(f"  - {len(attached_polymer_atoms)} atoms from attached polymers")

def main():
    write_attached_polymers_lammpstrj()

if __name__ == "__main__":
    main()