#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
from multiprocessing import Pool
import re
import pickle
import sys
from pathlib import Path

from read_parameters import read_parameters


def prepare_trajectory():
    # Read params file and check all paramters put in okay.
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

    # Load universe once (only for geometry + counts)
    u = mda.Universe(infile, format="LAMMPSDUMP")

    n_frames = len(u.trajectory)

    print(f"Loaded universe {infile} with {n_frames} timesteps")

    # Surface atoms and surface z positions (frame 0)
    seltext = "type " + " ".join(str(i) for i in surface_id)
    surface_atoms = u.select_atoms(seltext)

    surface_z = surface_atoms.positions[:, 2].copy()
    surface_z.sort()

    if len(surface_z) < 2:
        raise SystemExit("ERROR: Not enough surface atoms found to define two surfaces.")

    bottom_surface_z = surface_z[len(surface_z) // 2 - 1]
    top_surface_z = surface_z[len(surface_z) // 2]

    print(f"Bottom surface at z={bottom_surface_z}")
    print(f"Top surface at z={top_surface_z}")

    # We need a stable indexing scheme for polymers. Your original approach assumed:
    # [surface atoms first] then polymers packed contiguously.
    n_surface_atoms = len(surface_atoms)

    print(f"Surface atoms: {n_surface_atoms}")

    # Prepare worker args: let each process open its own Universe
    worker_args = [
        (infile, n_surface_atoms, i, additive_atoms_num, bottom_surface_z, top_surface_z, n_frames)
        for i in range(additive_num)
    ]

    return worker_args, additive_num, n_frames, working_dir


def process_polymer(worker_arg):
    (
        infile, n_surface_atoms, polymer_i, additive_atoms_num,
        bottom_surface_z, top_surface_z, n_frames
    ) = worker_arg

    # Each worker reads its own Universe (avoids pickling / shared trajectory issues)
    u = mda.Universe(infile, format="LAMMPSDUMP")

    start = n_surface_atoms + polymer_i * additive_atoms_num
    stop = start + additive_atoms_num  # python stop is exclusive

    polymer = u.atoms[start:stop]

    attachment_by_frame = []
    surface_by_frame = []

    for ts in u.trajectory:
        z = polymer.positions[:, 2]

        # attachment decision: any atom within 5 A of either surface
        attached_to = None
        if np.any(z < bottom_surface_z + 5.0):
            attached_to = "bottom"
        elif np.any(z > top_surface_z - 5.0):
            attached_to = "top"

        if attached_to is None:
            attachment_by_frame.append([])  # not attached
            surface_by_frame.append(None)
            continue

        if attached_to == "bottom":
            d = z - bottom_surface_z
        else:
            d = top_surface_z - z

        # Build list of (atom_type, distance)
        # MDAnalysis usually stores types as strings; keep as-is
        frame_list = [(atom.type, float(dist)) for atom, dist in zip(polymer, d)]
        attachment_by_frame.append(frame_list)
        surface_by_frame.append(attached_to)

    # (Optional sanity check)
    if len(attachment_by_frame) != n_frames:
        raise RuntimeError("Frame count mismatch while processing polymer.")

    return attachment_by_frame, surface_by_frame


def main():
    worker_args, additive_num, n_frames, working_dir = prepare_trajectory()

    print("Start of analysis")

    with Pool() as pool:
        # results shape: [polymer][(frame atoms, frame surface)]
        results = pool.map(process_polymer, worker_args)

    # Transpose to desired structure: [frame][polymer][...]
    final = [[results[p][0][f] for p in range(additive_num)] for f in range(n_frames)]
    surfaces = [[results[p][1][f] for p in range(additive_num)] for f in range(n_frames)]

    print(f"Built attachment list with {len(final)} frames and {len(final[0]) if final else 0} polymers.")

    out_path = working_dir / "attachment_master.pkl"
    with out_path.open("wb") as handle:
        pickle.dump(final, handle, protocol=pickle.HIGHEST_PROTOCOL)

    surface_path = working_dir / "attachment_surface.pkl"
    with surface_path.open("wb") as handle:
        pickle.dump(surfaces, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Saved output to {out_path}")
    print(f"Saved surface labels to {surface_path}")


if __name__ == "__main__":
    main()