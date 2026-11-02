#!/usr/bin/env python3
"""Calculate all-additive and attached-additive density profiles without plotting."""

import argparse
import ast
import pickle
import re
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from z_scripts.helper_scripts.read_parameters import read_parameters


BIN_WIDTH_ANGSTROM = 0.2
MAX_DISTANCE_ANGSTROM = 40.0
AMU_KG = 1.66053906660e-27


def parse_ids(value, key):
    try:
        identifiers = [
            int(item) for item in re.split(r"[,\s]+", value.strip()) if item
        ]
    except ValueError as exc:
        raise SystemExit(f"ERROR: {key} must contain integer atom types") from exc
    if not identifiers:
        raise SystemExit(f"ERROR: {key} cannot be empty")
    return identifiers


def read_element_masses(path):
    module = ast.parse(path.read_text(), filename=str(path))
    for statement in module.body:
        if isinstance(statement, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "ELEMENTS"
            for target in statement.targets
        ):
            elements = ast.literal_eval(statement.value)
            return {
                name: float(properties["mass"])
                for name, properties in elements.items()
            }
    raise ValueError(f"ELEMENTS dictionary not found in {path}")


def parse_statistics(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: Attachment statistics not found: {path}")
    values = {}
    for line in path.read_text().splitlines():
        key, separator, value = line.partition("=")
        if separator:
            try:
                values[key.strip()] = float(value.strip())
            except ValueError as exc:
                raise SystemExit(
                    f"ERROR: Invalid numeric value for {key.strip()} in {path}"
                ) from exc
    return values


def get_equilibrium_starts(stats, symmetric):
    if "EQUILIBRIUM" not in stats:
        raise SystemExit("ERROR: EQUILIBRIUM missing from attachment_statistics.txt")
    if symmetric:
        start = int(stats["EQUILIBRIUM"])
        return start, start

    needed = ("TOP_EQUILIBRIUM", "BOTTOM_EQUILIBRIUM")
    missing = [key for key in needed if key not in stats]
    if missing:
        raise SystemExit(
            "ERROR: Asymmetric surfaces require " + ", ".join(missing)
            + " in attachment_statistics.txt"
        )
    return int(stats["TOP_EQUILIBRIUM"]), int(stats["BOTTOM_EQUILIBRIUM"])


def find_surface_planes(universe, surface_ids):
    surface_atoms = universe.select_atoms("type " + " ".join(map(str, surface_ids)))
    z_positions = np.sort(surface_atoms.positions[:, 2])
    if len(z_positions) < 2:
        raise SystemExit("ERROR: Not enough surface atoms to identify both surfaces")
    middle = len(z_positions) // 2
    bottom = float(z_positions[middle - 1])
    top = float(z_positions[middle])
    if top <= bottom:
        raise SystemExit("ERROR: Could not distinguish top and bottom surface planes")
    return bottom, top, len(surface_atoms)


def add_density_frame(head_hist, tail_hist, atom_types, distances, head_ids,
                      tail_ids, type_masses, area_angstrom2):
    """Accumulate one frame's head/tail mass densities for a surface side."""
    atom_types = np.asarray(atom_types, dtype=int)
    distances = np.asarray(distances, dtype=float)
    bin_count = int(round(MAX_DISTANCE_ANGSTROM / BIN_WIDTH_ANGSTROM))
    volume_m3 = area_angstrom2 * BIN_WIDTH_ANGSTROM * 1e-30

    head_mask = np.isin(atom_types, list(head_ids)) & (distances >= 0)
    tail_mask = np.isin(atom_types, list(tail_ids)) & (distances >= 0)

    def histogram(mask):
        if not np.any(mask):
            return np.zeros(bin_count, dtype=float)
        masses = np.array([type_masses.get(atom_type, 0.0) for atom_type in atom_types[mask]])
        return np.histogram(
            distances[mask],
            bins=bin_count,
            range=(0, MAX_DISTANCE_ANGSTROM),
            weights=masses * AMU_KG / volume_m3,
        )[0]

    head_hist += histogram(head_mask)
    tail_hist += histogram(tail_mask)


def add_mass_counts(head_hist, tail_hist, atom_types, distances, head_ids,
                    tail_ids, type_masses):
    """Accumulate raw atomic masses (amu) into attached-profile bins."""
    atom_types = np.asarray(atom_types, dtype=int)
    distances = np.asarray(distances, dtype=float)
    bin_count = len(head_hist)
    head_mask = np.isin(atom_types, list(head_ids)) & (distances >= 0)
    tail_mask = np.isin(atom_types, list(tail_ids)) & (distances >= 0)

    for mask, histogram in ((head_mask, head_hist), (tail_mask, tail_hist)):
        if np.any(mask):
            masses = np.array(
                [type_masses.get(atom_type, 0.0) for atom_type in atom_types[mask]]
            )
            histogram += np.histogram(
                distances[mask],
                bins=bin_count,
                range=(0, MAX_DISTANCE_ANGSTROM),
                weights=masses,
            )[0]


def collect_all_density(trajectory_path, params, type_masses, top_start, bottom_start):
    head_ids = set(parse_ids(params["ADDITIVE_HEAD_ID"], "ADDITIVE_HEAD_ID"))
    tail_ids = set(parse_ids(params["ADDITIVE_TAIL_ID"], "ADDITIVE_TAIL_ID"))
    surface_ids = parse_ids(params["SURFACE_ID"], "SURFACE_ID")
    polymer_count = int(params["ADDITIVE_NUM"])
    atoms_per_polymer = int(params["ADDITIVE_ATOMS_NUM"])
    if polymer_count <= 0 or atoms_per_polymer <= 0:
        raise SystemExit("ERROR: ADDITIVE_NUM and ADDITIVE_ATOMS_NUM must be positive")

    universe = mda.Universe(trajectory_path, format="LAMMPSDUMP")
    bottom_plane, top_plane, surface_count = find_surface_planes(universe, surface_ids)
    midpoint = (bottom_plane + top_plane) / 2
    bin_count = int(round(MAX_DISTANCE_ANGSTROM / BIN_WIDTH_ANGSTROM))
    top_head = np.zeros(bin_count)
    top_tail = np.zeros(bin_count)
    bottom_head = np.zeros(bin_count)
    bottom_tail = np.zeros(bin_count)
    top_frame_count = 0
    bottom_frame_count = 0

    for ts in universe.trajectory:
        frame_index = ts.frame
        if frame_index < top_start and frame_index < bottom_start:
            continue
        area = float(ts.dimensions[0] * ts.dimensions[1])
        if area <= 0:
            raise SystemExit(f"ERROR: Invalid XY box area in frame {frame_index}")

        top_types, top_distances = [], []
        bottom_types, bottom_distances = [], []
        for polymer_index in range(polymer_count):
            start = surface_count + polymer_index * atoms_per_polymer
            stop = start + atoms_per_polymer
            if stop > len(universe.atoms):
                raise SystemExit(
                    "ERROR: Polymer atom range exceeds trajectory atom count; "
                    "check SURFACE_ID, ADDITIVE_NUM, and ADDITIVE_ATOMS_NUM"
                )
            polymer = universe.atoms[start:stop]
            z = polymer.positions[:, 2]
            if not len(z):
                continue
            if float(np.mean(z)) >= midpoint:
                if frame_index >= top_start:
                    top_types.extend(polymer.types)
                    top_distances.extend(np.abs(z - top_plane))
            elif frame_index >= bottom_start:
                bottom_types.extend(polymer.types)
                bottom_distances.extend(np.abs(z - bottom_plane))

        if frame_index >= top_start:
            add_density_frame(top_head, top_tail, top_types, top_distances,
                              head_ids, tail_ids, type_masses, area)
            top_frame_count += 1
        if frame_index >= bottom_start:
            add_density_frame(bottom_head, bottom_tail, bottom_types, bottom_distances,
                              head_ids, tail_ids, type_masses, area)
            bottom_frame_count += 1

    if not top_frame_count or not bottom_frame_count:
        raise SystemExit("ERROR: No frames available after one or both equilibrium points")
    top_head /= top_frame_count
    top_tail /= top_frame_count
    bottom_head /= bottom_frame_count
    bottom_tail /= bottom_frame_count
    return top_head, top_tail, bottom_head, bottom_tail


def collect_attached_density(output_dir, params, type_masses, top_start, bottom_start):
    master_path = output_dir / "attachment_master.pkl"
    surface_path = output_dir / "attachment_surface.pkl"
    for path in (master_path, surface_path):
        if not path.is_file():
            raise SystemExit(f"ERROR: Attachment analysis data not found: {path}")
    with master_path.open("rb") as handle:
        attachment_data = pickle.load(handle)
    with surface_path.open("rb") as handle:
        surface_data = pickle.load(handle)

    if not attachment_data or len(surface_data) != len(attachment_data):
        raise SystemExit("ERROR: Attachment data and surface labels are empty or mismatched")
    if any(len(fr) != len(sr) for fr, sr in zip(attachment_data, surface_data)):
        raise SystemExit("ERROR: Attachment surface labels do not match polymer data")

    head_ids = set(parse_ids(params["ADDITIVE_HEAD_ID"], "ADDITIVE_HEAD_ID"))
    tail_ids = set(parse_ids(params["ADDITIVE_TAIL_ID"], "ADDITIVE_TAIL_ID"))
    bin_count = int(round(MAX_DISTANCE_ANGSTROM / BIN_WIDTH_ANGSTROM))
    top_head = np.zeros(bin_count)
    top_tail = np.zeros(bin_count)
    bottom_head = np.zeros(bin_count)
    bottom_tail = np.zeros(bin_count)
    top_frames = max(0, len(attachment_data) - top_start)
    bottom_frames = max(0, len(attachment_data) - bottom_start)
    if not top_frames or not bottom_frames:
        raise SystemExit("ERROR: No attachment frames available after equilibrium")

    for frame_index, (frame, labels) in enumerate(zip(attachment_data, surface_data)):
        top_types, top_distances = [], []
        bottom_types, bottom_distances = [], []
        for polymer, side in zip(frame, labels):
            if not polymer:
                continue
            for atom_type, distance in polymer:
                if side == "top" and frame_index >= top_start:
                    top_types.append(atom_type)
                    top_distances.append(distance)
                elif side == "bottom" and frame_index >= bottom_start:
                    bottom_types.append(atom_type)
                    bottom_distances.append(distance)
        if frame_index >= top_start:
            add_mass_counts(top_head, top_tail, top_types, top_distances,
                            head_ids, tail_ids, type_masses)
        if frame_index >= bottom_start:
            add_mass_counts(bottom_head, bottom_tail, bottom_types, bottom_distances,
                            head_ids, tail_ids, type_masses)

    # Attachment records are already distances from the associated surface;
    # retain raw mass counts so the caller can apply surface-volume normalization.
    return top_head, top_tail, bottom_head, bottom_tail, top_frames, bottom_frames


def normalize_attached(top_head, top_tail, bottom_head, bottom_tail,
                       top_frames, bottom_frames, surface_area_angstrom2):
    """Convert accumulated atomic masses to average mass density profiles."""
    volume = surface_area_angstrom2 * BIN_WIDTH_ANGSTROM * 1e-30
    top_head = top_head * AMU_KG / volume / top_frames
    top_tail = top_tail * AMU_KG / volume / top_frames
    bottom_head = bottom_head * AMU_KG / volume / bottom_frames
    bottom_tail = bottom_tail * AMU_KG / volume / bottom_frames
    return top_head, top_tail, bottom_head, bottom_tail


def make_output_columns(top_head, top_tail, bottom_head, bottom_tail):
    top_total = top_head + top_tail
    bottom_total = bottom_head + bottom_tail
    both_head = (top_head + bottom_head) / 2
    both_tail = (top_tail + bottom_tail) / 2
    both_total = both_head + both_tail
    return np.column_stack((
        (np.arange(len(top_head)) + 0.5) * BIN_WIDTH_ANGSTROM / 10,
        both_total,
        both_head,
        both_tail,
        top_total,
        top_head,
        top_tail,
        bottom_total,
        bottom_head,
        bottom_tail,
    ))


def write_density(path, columns):
    np.savetxt(
        path,
        columns,
        delimiter="\t",
        header=(
            "distance_nm\tall_total_kg_m3\tall_head_kg_m3\tall_tail_kg_m3\t"
            "top_total_kg_m3\ttop_head_kg_m3\ttop_tail_kg_m3\t"
            "bottom_total_kg_m3\tbottom_head_kg_m3\tbottom_tail_kg_m3"
        ),
        fmt="%.6e",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parameters", type=Path, help="Path to density_calculator_parameters.in")
    args = parser.parse_args()
    params = read_parameters(args.parameters)
    required = (
        "LAMMPSDUMP_DIR", "ADDITIVE_NUM", "ADDITIVE_ATOMS_NUM",
        "ADDITIVE_HEAD_ID", "ADDITIVE_TAIL_ID", "SURFACE_ID", "SURFACE_SYMMETRIC",
    )
    missing = [key for key in required if key not in params]
    if missing:
        raise SystemExit(f"ERROR: Missing parameters: {', '.join(missing)}")

    output_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
    trajectory_path = output_dir / "dump_nosolvent.lammpstrj"
    if not trajectory_path.is_file():
        raise SystemExit(f"ERROR: Solvent-removed trajectory not found: {trajectory_path}")
    stats_path = output_dir / "attachment_statistics.txt"
    stats = parse_statistics(stats_path)
    symmetric = params["SURFACE_SYMMETRIC"].strip().lower() in {"1", "true", "yes"}
    top_start, bottom_start = get_equilibrium_starts(stats, symmetric)

    masses = read_element_masses(Path(__file__).with_name("elements.data"))
    type_to_element = {
        int(key[5:]): value
        for key, value in params.items()
        if key.startswith("ELEM_") and value in masses
    }
    head_types = set(parse_ids(params["ADDITIVE_HEAD_ID"], "ADDITIVE_HEAD_ID"))
    tail_types = set(parse_ids(params["ADDITIVE_TAIL_ID"], "ADDITIVE_TAIL_ID"))
    needed_types = head_types | tail_types
    missing_types = sorted(needed_types - type_to_element.keys())
    if missing_types:
        raise SystemExit(f"ERROR: Missing ELEM mappings for atom types: {missing_types}")
    type_masses = {
        atom_type: masses[type_to_element[atom_type]] for atom_type in needed_types
    }

    top_head, top_tail, bottom_head, bottom_tail = collect_all_density(
        trajectory_path, params, type_masses, top_start, bottom_start
    )

    all_profile = make_output_columns(top_head, top_tail, bottom_head, bottom_tail)
    write_density(output_dir / "total_density.txt", all_profile)

    # Attachment records carry distances from their own surfaces. Use the same
    # xy area measured from the first trajectory frame for physical density units.
    universe = mda.Universe(trajectory_path, format="LAMMPSDUMP")
    surface_area = float(universe.dimensions[0] * universe.dimensions[1])
    attached_raw = collect_attached_density(
        output_dir, params, type_masses, top_start, bottom_start
    )
    attached_top_head, attached_top_tail, attached_bottom_head, attached_bottom_tail, n_top, n_bottom = attached_raw
    attached_profiles = normalize_attached(
        attached_top_head, attached_top_tail,
        attached_bottom_head, attached_bottom_tail,
        n_top, n_bottom, surface_area,
    )
    attached_profile = make_output_columns(*attached_profiles)
    write_density(output_dir / "attached_density.txt", attached_profile)

    print(f"Wrote total density profile to {output_dir / 'total_density.txt'}")
    print(f"Wrote attached density profile to {output_dir / 'attached_density.txt'}")


if __name__ == "__main__":
    main()
