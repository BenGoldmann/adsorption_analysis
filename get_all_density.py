#!/usr/bin/env python3
"""Plot head and tail density profiles for all polymers after equilibrium."""

import argparse
import ast
import re
from pathlib import Path

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np

from read_parameters import read_parameters


def parse_ids(value):
	return [int(item) for item in re.split(r"[,\s]+", value.strip()) if item]


def read_element_masses(path):
	module = ast.parse(path.read_text(), filename=str(path))
	for statement in module.body:
		if isinstance(statement, ast.Assign):
			if any(
				isinstance(target, ast.Name) and target.id == "ELEMENTS"
				for target in statement.targets
			):
				elements = ast.literal_eval(statement.value)
				return {name: float(properties["mass"]) for name, properties in elements.items()}
	raise ValueError(f"ELEMENTS dictionary not found in {path}")


def read_equilibrium_start(working_dir):
	path = working_dir / "attachment_statistics.txt"
	if not path.is_file():
		raise FileNotFoundError(f"Attachment statistics not found: {path}")

	for line in path.read_text().splitlines():
		key, separator, value = line.partition("=")
		if separator and key.strip() == "EQUILIBRIUM":
			return int(float(value.strip()))
	raise ValueError(f"EQUILIBRIUM not found in {path}")


def surface_positions(universe, surface_ids):
	surface_atoms = universe.select_atoms("type " + " ".join(map(str, surface_ids)))
	positions = np.sort(surface_atoms.positions[:, 2])
	if len(positions) < 2:
		raise ValueError("Not enough surface atoms found to define two surfaces")
	midpoint = len(positions) // 2
	return positions[midpoint - 1], positions[midpoint]


def plot_density(frames, params, output_dir):
	if not frames:
		raise ValueError("Trajectory contains no frames")

	head_ids = set(parse_ids(params["ADDITIVE_HEAD_ID"]))
	tail_ids = set(parse_ids(params["ADDITIVE_TAIL_ID"]))
	surface_area = float(params["SURFACE_AREA"])

	masses = read_element_masses(Path(__file__).with_name("elements.data"))
	atom_type_to_element = {
		int(key[5:]): value
		for key, value in params.items()
		if key.startswith("ELEM_") and value in masses
	}

	bin_size = 0.2
	n_frames = len(frames)
	scaling_factor = surface_area * bin_size * n_frames / 10**30
	head_positions, head_weights = [], []
	tail_positions, tail_weights = [], []

	for frame in frames:
		for atom_type, distance in frame:
			atom_type = int(atom_type)
			element = atom_type_to_element.get(atom_type)
			if element is None:
				continue
			weight = masses[element] * 1.66054e-27 / scaling_factor
			if atom_type in head_ids:
				head_positions.append(distance)
				head_weights.append(weight)
			elif atom_type in tail_ids:
				tail_positions.append(distance)
				tail_weights.append(weight)

	bins = int(max(head_positions) / bin_size)
	head_density, edges = np.histogram(
		head_positions, bins=bins, weights=head_weights, range=(0, max(head_positions))
	)
	tail_density, _ = np.histogram(
		tail_positions, bins=bins, weights=tail_weights, range=(0, max(tail_positions))
	)
	centers = (edges[:-1] + edges[1:]) / 2 / 10

	figure, axis = plt.subplots(figsize=(5, 5))
	axis.plot(centers, head_density, label="Head", color="tab:orange")
	axis.plot(centers, tail_density, label="Tail", color="tab:blue")
	axis.set_xlabel("$r$ (nm)", fontsize=12)
	axis.set_ylabel(r"$\rho(r)$ (kg m$^{-3}$)", fontsize=12)
	axis.set_ylim(0, 500)
	axis.tick_params(axis="both", labelsize=12)
	axis.legend()
	figure.tight_layout()
	figure.savefig(output_dir / "all_density.jpg", bbox_inches="tight", dpi=1000)
	plt.close(figure)

	np.savetxt(
		output_dir / "all_density.txt",
		np.column_stack((centers, head_density, tail_density)),
		header="r (nm)  rho_head (kg/m^3)  rho_tail (kg/m^3)",
		fmt=("%.4f", "%.6e", "%.6e"),
	)


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("parameters", type=Path, help="Path to parameters.in")
	args = parser.parse_args()
	params = read_parameters(args.parameters)
	required = (
		"LAMMPSDUMP_DIR",
		"SURFACE_ID",
		"ADDITIVE_NUM",
		"ADDITIVE_ATOMS_NUM",
		"ADDITIVE_HEAD_ID",
		"ADDITIVE_TAIL_ID",
		"SURFACE_AREA",
	)
	missing = [key for key in required if key not in params]
	if missing:
		raise SystemExit(f"ERROR: Missing parameters: {', '.join(missing)}")

	output_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
	trajectory_path = output_dir / "dump_nosolvent.lammpstrj"
	if not trajectory_path.is_file():
		raise SystemExit(f"ERROR: Solvent-removed trajectory not found: {trajectory_path}")

	equilibrium_start = read_equilibrium_start(output_dir)
	universe = mda.Universe(trajectory_path, format="LAMMPSDUMP")
	bottom_surface, _ = surface_positions(universe, parse_ids(params["SURFACE_ID"]))
	surface_count = len(
		universe.select_atoms("type " + " ".join(map(str, parse_ids(params["SURFACE_ID"]))))
	)
	polymer_count = int(params["ADDITIVE_NUM"])
	atoms_per_polymer = int(params["ADDITIVE_ATOMS_NUM"])

	frames = []
	for ts in universe.trajectory:
		if ts.frame < equilibrium_start:
			continue
		frame = []
		for polymer_index in range(polymer_count):
			start = surface_count + polymer_index * atoms_per_polymer
			polymer = universe.atoms[start : start + atoms_per_polymer]
			z = polymer.positions[:, 2]
			distances = abs(z - bottom_surface)
			frame.extend((atom.type, float(distance)) for atom, distance in zip(polymer, distances))
		frames.append(frame)
	print(f"Collected {len(frames)} frames after equilibrium (frame {equilibrium_start})")

	if not frames:
		raise ValueError("No trajectory frames remain after equilibrium")

	plot_density(frames, params, output_dir)
	print(f"Saved all-polymer density profile to {output_dir}")


if __name__ == "__main__":
	main()
