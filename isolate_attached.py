import argparse
import math
import pickle
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from z_scripts.helper_scripts.read_parameters import read_parameters


def load_attachment_counts(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: Attachment time series not found: {path}")
    try:
        counts = np.loadtxt(path, comments="#", ndmin=2)
    except ValueError as exc:
        raise SystemExit(f"ERROR: Could not parse attachment counts in {path}: {exc}") from exc
    if counts.shape[1] < 4:
        raise SystemExit(
            f"ERROR: Expected step, total, top, and bottom columns in {path}"
        )
    return counts


def load_attachment_statistics(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: Attachment statistics not found: {path}")
    values = {}
    for line in path.read_text().splitlines():
        key, separator, value = line.partition("=")
        if not separator:
            continue
        try:
            values[key.strip()] = float(value.strip())
        except ValueError as exc:
            raise SystemExit(
                f"ERROR: Invalid value for {key.strip()} in {path}"
            ) from exc
    return values


def read_frame(fin):
    """Read one LAMMPS dump frame and return its header and atom lines."""
    first = fin.readline()
    if not first:
        return None
    timestep = fin.readline()
    number_header = fin.readline()
    number_line = fin.readline()
    if not timestep or not number_header or not number_line:
        raise SystemExit("ERROR: Incomplete LAMMPS dump frame header")
    try:
        atom_count = int(number_line.strip())
    except ValueError as exc:
        raise SystemExit(f"ERROR: Invalid atom count: {number_line.strip()!r}") from exc
    box_header = fin.readline()
    box_lines = [fin.readline() for _ in range(3)]
    atom_header = fin.readline()
    atom_lines = [fin.readline() for _ in range(atom_count)]
    if not box_header or any(not line for line in box_lines) or not atom_header:
        raise SystemExit("ERROR: Incomplete LAMMPS dump frame")
    if any(not line for line in atom_lines):
        raise SystemExit("ERROR: Unexpected end of file in LAMMPS atom data")
    if not atom_header.startswith("ITEM: ATOMS"):
        raise SystemExit(f"ERROR: Expected ATOMS header, found {atom_header.strip()!r}")
    columns = atom_header.strip().split()[2:]
    if "id" not in columns:
        raise SystemExit("ERROR: Dump ATOMS section must include an id column")
    return (first, timestep, number_header, box_header, box_lines, atom_header,
            atom_lines, columns, columns.index("id"))


def write_selected_frame(dump_path, output_path, frame_index, selected_atom_ids):
    selected_atom_ids = {int(atom_id) for atom_id in selected_atom_ids}
    with dump_path.open() as source, output_path.open("w") as target:
        current_index = 0
        while True:
            frame = read_frame(source)
            if frame is None:
                break
            (first, timestep, number_header, box_header, box_lines, atom_header,
             atom_lines, _, id_index) = frame
            if current_index == frame_index:
                kept_lines = []
                for atom_line in atom_lines:
                    fields = atom_line.split()
                    try:
                        atom_id = int(fields[id_index])
                    except (IndexError, ValueError) as exc:
                        raise SystemExit(f"ERROR: Invalid atom ID in row: {atom_line!r}") from exc
                    if atom_id in selected_atom_ids:
                        kept_lines.append(atom_line)

                target.writelines((first, timestep, number_header, f"{len(kept_lines)}\n"))
                target.write(box_header)
                target.writelines(box_lines)
                target.write(atom_header)
                target.writelines(kept_lines)
                return len(kept_lines)
            current_index += 1
    raise SystemExit(f"ERROR: Requested frame index {frame_index} is not in {dump_path}")


def selected_surface_atom_indices(universe, surface_ids, side_surface_ids, side):
    all_surface_atoms = universe.select_atoms(
        "type " + " ".join(map(str, surface_ids))
    )
    if len(all_surface_atoms) < 2:
        raise SystemExit("ERROR: Not enough surface atoms to separate top and bottom")
    z_positions = np.sort(all_surface_atoms.positions[:, 2])
    middle = len(z_positions) // 2
    split_z = (z_positions[middle - 1] + z_positions[middle]) / 2
    surface_type_group = universe.select_atoms(
        "type " + " ".join(map(str, side_surface_ids))
    )
    if side == "bottom":
        side_mask = surface_type_group.positions[:, 2] < split_z
    else:
        side_mask = surface_type_group.positions[:, 2] >= split_z
    return surface_type_group.indices[side_mask]


def collect_selected_atom_ids(universe, params, attachment_data, surface_data,
                              frame_index, side):
    surface_ids = [int(value) for value in params["SURFACE_ID"].split()]
    side_surface_ids = [
        int(value) for value in params[f"SURFACE_{side.upper()}_ID"].split()
    ]
    surface_indices = selected_surface_atom_indices(
        universe, surface_ids, side_surface_ids, side
    )
    surface_count = len(universe.select_atoms("type " + " ".join(map(str, surface_ids))))
    additive_num = int(params["ADDITIVE_NUM"])
    atoms_per_additive = int(params["ADDITIVE_ATOMS_NUM"])
    if frame_index >= len(attachment_data) or frame_index >= len(surface_data):
        raise SystemExit("ERROR: Attachment-analysis frame count is shorter than trajectory")

    selected_indices = set(map(int, surface_indices))
    frame = attachment_data[frame_index]
    labels = surface_data[frame_index]
    if len(frame) != additive_num or len(labels) != additive_num:
        raise SystemExit("ERROR: Attachment data polymer count does not match parameters")
    for polymer_index, (polymer_data, attached_side) in enumerate(zip(frame, labels)):
        if polymer_data and attached_side == side:
            start = surface_count + polymer_index * atoms_per_additive
            stop = min(start + atoms_per_additive, len(universe.atoms))
            selected_indices.update(range(start, stop))
    if not selected_indices:
        raise SystemExit(f"ERROR: No {side}-surface or attached additive atoms selected")
    atom_ids = universe.atoms.ids[list(sorted(selected_indices))]
    return atom_ids, len(selected_indices) - len(surface_indices)


def write_surface_snapshot(params, surface, side, frame_index, output_path):
    working_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
    dump_path = working_dir / "dump_nosolvent.lammpstrj"
    universe = mda.Universe(dump_path, format="LAMMPSDUMP")
    if frame_index >= len(universe.trajectory):
        raise SystemExit(f"ERROR: Frame index {frame_index} exceeds trajectory length")
    universe.trajectory[frame_index]

    for key in ("SURFACE_ID", "ADDITIVE_NUM", "ADDITIVE_ATOMS_NUM"):
        if key not in params:
            raise SystemExit(f"ERROR: {key} missing from isolate parameters")
    side_surface_key = f"SURFACE_{side.upper()}_ID"
    if side_surface_key not in params:
        raise SystemExit(f"ERROR: {side_surface_key} missing from isolate parameters")
    with (working_dir / "attachment_master.pkl").open("rb") as handle:
        attachment_data = pickle.load(handle)
    with (working_dir / "attachment_surface.pkl").open("rb") as handle:
        surface_data = pickle.load(handle)
    if len(attachment_data) != len(universe.trajectory) or len(surface_data) != len(attachment_data):
        raise SystemExit("ERROR: Attachment-analysis data does not match trajectory frame count")

    atom_ids, additive_atom_count = collect_selected_atom_ids(
        universe, params, attachment_data, surface_data, frame_index, side
    )
    output_path = working_dir / output_path
    atom_count = write_selected_frame(dump_path, output_path, frame_index, atom_ids)
    output_path.with_suffix(".side").write_text(f"{side}\n")
    print(
        f"Wrote {atom_count} atoms for {surface} ({side}) at frame {frame_index} "
        f"to {output_path} ({additive_atom_count} additive atoms)"
    )


def target_attachment_count(values, description):
    """Round the post-equilibrium mean attachment count up to a whole polymer."""
    if len(values) == 0:
        raise SystemExit(f"ERROR: No post-equilibrium frames available for {description}")
    return math.ceil(float(np.mean(values)))


def latest_frame_with_count(counts, start_frame, column, target_count, side):
    """Find the latest frame at/after equilibrium with exactly target_count attached."""
    if start_frame < 0 or start_frame >= len(counts):
        raise SystemExit(
            f"ERROR: Equilibrium frame {start_frame} is outside the "
            f"{len(counts)}-frame attachment series"
        )
    matches = np.flatnonzero(counts[start_frame:, column] == target_count)
    if not len(matches):
        raise SystemExit(
            f"ERROR: No {side}-surface frame at or after equilibrium frame "
            f"{start_frame} has exactly {target_count} attached polymer(s)"
        )
    return start_frame + int(matches[-1])


def select_representative_frames(params, counts, statistics):
    outputs = params.get("OUTPUT_SURFACES", "").split()
    if not outputs:
        raise SystemExit("ERROR: OUTPUT_SURFACES missing from isolate parameters")
    if len(counts) == 0:
        raise SystemExit("ERROR: Attachment time series contains no data rows")

    if len(outputs) == 2:
        required = ("TOP_EQUILIBRIUM", "BOTTOM_EQUILIBRIUM")
        missing = [key for key in required if key not in statistics]
        if missing:
            raise SystemExit(
                "ERROR: Asymmetric surfaces require " + ", ".join(missing)
                + " in attachment_statistics.txt"
            )
        top_start = int(statistics["TOP_EQUILIBRIUM"])
        bottom_start = int(statistics["BOTTOM_EQUILIBRIUM"])
        top_target = target_attachment_count(
            counts[top_start:, 2], "top BaTiO3 surface"
        )
        bottom_target = target_attachment_count(
            counts[bottom_start:, 3], "bottom BaTiO3 surface"
        )
        top_frame = latest_frame_with_count(
            counts, top_start, 2, top_target, "top"
        )
        bottom_frame = latest_frame_with_count(
            counts, bottom_start, 3, bottom_target, "bottom"
        )
        print(
            f"BaTiO3 target post-equilibrium attachments: "
            f"top={top_target} (latest matching frame {top_frame}), "
            f"bottom={bottom_target} (latest matching frame {bottom_frame})"
        )
        return (("tio2", "top", top_frame), ("bao", "bottom", bottom_frame))

    if len(outputs) != 1:
        raise SystemExit("ERROR: OUTPUT_SURFACES must contain one or two surface names")
    if counts.shape[1] < 4:
        raise SystemExit("ERROR: Need top and bottom attachment columns to select saturation")
    if "EQUILIBRIUM" not in statistics:
        raise SystemExit("ERROR: EQUILIBRIUM missing from attachment_statistics.txt")

    # For equivalent surfaces, pool top and bottom counts after equilibrium,
    # round the mean up to a whole polymer, then choose the latest face/frame
    # with exactly that many attached polymers.
    start_frame = int(statistics["EQUILIBRIUM"])
    if start_frame < 0 or start_frame >= len(counts):
        raise SystemExit(
            f"ERROR: EQUILIBRIUM frame {start_frame} is outside the "
            f"{len(counts)}-frame attachment series"
        )
    pooled_face_counts = np.concatenate((counts[start_frame:, 2], counts[start_frame:, 3]))
    target_count = target_attachment_count(pooled_face_counts, outputs[0])
    top_matches = np.flatnonzero(counts[start_frame:, 2] == target_count)
    bottom_matches = np.flatnonzero(counts[start_frame:, 3] == target_count)
    candidates = [
        (start_frame + int(index), "top") for index in top_matches
    ] + [
        (start_frame + int(index), "bottom") for index in bottom_matches
    ]
    if not candidates:
        raise SystemExit(
            f"ERROR: No top or bottom frame at or after equilibrium frame "
            f"{start_frame} has exactly {target_count} attached polymer(s)"
        )
    frame_index, side = max(candidates, key=lambda candidate: candidate[0])
    print(
        f"{outputs[0]} target post-equilibrium attachment count={target_count}; "
        f"selected latest matching {side} frame {frame_index}"
    )
    return ((outputs[0], side, frame_index),)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parameters", type=Path, help="Path to isolate_attached_parameters.in")
    args = parser.parse_args()
    params = read_parameters(args.parameters)
    required = (
        "LAMMPSDUMP_DIR", "ADDITIVE_NUM", "ADDITIVE_ATOMS_NUM", "SURFACE_ID",
        "OUTPUT_SURFACES",
    )
    missing = [key for key in required if key not in params]
    if missing:
        raise SystemExit(f"ERROR: Missing isolate parameters: {', '.join(missing)}")
    working_dir = Path(params.get("LAMMPSDUMP_DIR", "")).expanduser()
    if not working_dir.is_dir():
        raise SystemExit(f"ERROR: Working directory not found: {working_dir}")
    counts = load_attachment_counts(working_dir / "attachment_vs_time.txt")
    dump_path = working_dir / "dump_nosolvent.lammpstrj"
    if not dump_path.is_file():
        raise SystemExit(f"ERROR: Solvent-removed trajectory not found: {dump_path}")
    trajectory = mda.Universe(dump_path, format="LAMMPSDUMP")
    if len(counts) != len(trajectory.trajectory):
        raise SystemExit(
            f"ERROR: Attachment time series has {len(counts)} rows but trajectory has "
            f"{len(trajectory.trajectory)} frames"
        )
    statistics = load_attachment_statistics(working_dir / "attachment_statistics.txt")
    for surface, side, frame_index in select_representative_frames(
        params, counts, statistics
    ):
        write_surface_snapshot(
            params,
            surface,
            side,
            frame_index,
            f"{surface}_saturated.lammpstrj",
        )

if __name__ == "__main__":
    main()