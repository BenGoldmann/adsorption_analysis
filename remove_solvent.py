#!/usr/bin/env python3
import re
import sys
from pathlib import Path

from read_parameters import read_parameters


def filter_lammpstrj_by_type(infile, outfile, keep_types):
    """
    Filters a LAMMPS trajectory file to keep only atoms of specified types.
    Args:
        infile: path to the input LAMMPS trajectory file
        outfile: path to the output filtered LAMMPS trajectory file
        keep_types: set of atom types to keep (integers)
    """
    with open(infile, "r") as fin, open(outfile, "w") as fout:
        # The variable fin changes each time it is called, so loop goes through the entire file until EOF.
        while True:
            timestep_hdr = fin.readline()
            if not timestep_hdr:
                break  # EOF

            if not timestep_hdr.startswith("ITEM: TIMESTEP"):
                raise SystemExit(f"ERROR: Expected 'ITEM: TIMESTEP', got: {timestep_hdr.strip()!r}")

            ts = fin.readline()
            if not ts:
                raise SystemExit("ERROR: Unexpected EOF after ITEM: TIMESTEP")

            num_hdr = fin.readline()
            if not num_hdr or not num_hdr.startswith("ITEM: NUMBER OF ATOMS"):
                raise SystemExit(f"ERROR: Expected 'ITEM: NUMBER OF ATOMS', got: {num_hdr.strip()!r}")

            natoms_line = fin.readline()
            if not natoms_line:
                raise SystemExit("ERROR: Unexpected EOF reading number of atoms")
            try:
                natoms = int(natoms_line.strip())
            except ValueError:
                raise SystemExit(f"ERROR: Bad NUMBER OF ATOMS line: {natoms_line!r}")

            box_hdr = fin.readline()
            if not box_hdr or not box_hdr.startswith("ITEM: BOX BOUNDS"):
                raise SystemExit(f"ERROR: Expected 'ITEM: BOX BOUNDS', got: {box_hdr.strip()!r}")

            box_lines = []
            for _ in range(3):
                b = fin.readline()
                if not b:
                    raise SystemExit("ERROR: Unexpected EOF in BOX BOUNDS")
                box_lines.append(b)

            atoms_hdr = fin.readline()
            if not atoms_hdr or not atoms_hdr.startswith("ITEM: ATOMS"):
                raise SystemExit(f"ERROR: Expected 'ITEM: ATOMS', got: {atoms_hdr.strip()!r}")

            cols = atoms_hdr.strip().split()[2:]
            if "type" not in cols:
                raise SystemExit(f"ERROR: Dump does not contain a 'type' column. Columns are: {cols}")
            type_idx = cols.index("type")

            kept_lines = []
            for _ in range(natoms):
                atom_line = fin.readline()
                if not atom_line:
                    raise SystemExit("ERROR: Unexpected EOF inside ATOMS section")

                parts = atom_line.split()
                if len(parts) < len(cols):
                    raise SystemExit(f"ERROR: Atom line has fewer columns than header.\nLine: {atom_line!r}")

                try:
                    atype = int(parts[type_idx])
                except ValueError:
                    raise SystemExit(f"ERROR: Non-integer type in line: {atom_line!r}")

                if atype in keep_types:
                    kept_lines.append(atom_line)

            # Write the whole frame with updated atom count
            fout.write(timestep_hdr)
            fout.write(ts)
            fout.write(num_hdr)
            fout.write(f"{len(kept_lines)}\n")
            fout.write(box_hdr)
            fout.writelines(box_lines)
            fout.write(atoms_hdr)
            fout.writelines(kept_lines)


def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python script.py parameters.in")

    params_file = Path(sys.argv[1])
    if not params_file.is_file():
        raise SystemExit(f"ERROR: Parameters file not found: {params_file}")

    params = read_parameters(params_file)
    print("Parameter file successfully read")

    if "LAMMPSDUMP_DIR" not in params:
        raise SystemExit("ERROR: LAMMPSDUMP_DIR not found in parameters.in")
    working_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
    infile = working_dir / "dump.lammpstrj"
    if not infile.is_file():
        raise SystemExit(f"ERROR: Input file does not exist: {infile}")

    for k in ("SURFACE_ID", "ADDITIVE_HEAD_ID", "ADDITIVE_TAIL_ID"):
        if k not in params:
            raise SystemExit(f"ERROR: {k} not found (expected: {k} = ... )")

    try:
        surface_id = [int(i) for i in re.split(r"[,\s]+", params["SURFACE_ID"].strip()) if i]
        additive_head_id = [int(i) for i in re.split(r"[,\s]+", params["ADDITIVE_HEAD_ID"].strip()) if i]
        additive_tail_id = [int(i) for i in re.split(r"[,\s]+", params["ADDITIVE_TAIL_ID"].strip()) if i]
    except ValueError:
        raise SystemExit("ERROR: SURFACE_ID, ADDITIVE_HEAD_ID, ADDITIVE_TAIL_ID must be lists of integers.")

    # Combine atom types to keep.
    keep_types = set(surface_id + additive_head_id + additive_tail_id)
    print(f"Atom types to keep: {sorted(keep_types)}")

    # Set outfile.
    outfile = infile.with_name(f"{infile.stem}_nosolvent{infile.suffix}")
    print(f"Target file: {outfile}")

    # Filter LAMMPS file.
    print("Filtering LAMMPS trajectory")
    filter_lammpstrj_by_type(infile, outfile, keep_types)

    if outfile.is_file() and outfile.stat().st_size > 0:
        print(f"Successfully generated {outfile}")
    else:
        raise SystemExit(f"ERROR: Output file not generated (or empty): {outfile}")


if __name__ == "__main__":
    main()
