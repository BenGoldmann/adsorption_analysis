#!/usr/bin/env python3
"""Sum attached densities from all attachment_* subdirectories."""

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np


def read_density(filename: Path):
	"""Read position, head density, and tail density from a data file."""
	data = []
	with filename.open() as file:
		for line in file:
			line = line.strip()
			if not line or line.startswith(("#", "@")):
				continue
			try:
				values = [float(value) for value in line.replace(",", " ").split()]
			except ValueError:
				continue
			if len(values) >= 3:
				data.append(values[:3])

	if not data:
		raise ValueError(f"No numeric three-column data found in {filename}")
	data = np.asarray(data, dtype=float)
	return data[:, 0], data[:, 1], data[:, 2]


def main():
	parser = argparse.ArgumentParser(
		description="Sum head and tail densities in attachment_* directories."
	)
	parser.add_argument(
		"directory",
		nargs="?",
		type=Path,
		default=Path.cwd(),
		help="directory to search (default: current directory)",
	)
	args = parser.parse_args()
	directory = args.directory.resolve()
	files = sorted(directory.glob("attachment_*/attached_density.txt"))
	if not files:
		parser.error(f"No attachment_*/attached_density.txt files found in {directory}")

	position, head_sum, tail_sum = read_density(files[0])
	for filename in files[1:]:
		other_position, other_head, other_tail = read_density(filename)
		if len(other_position) != len(position) or not np.allclose(other_position, position):
			raise ValueError(f"Incompatible position grid in {filename}")
		head_sum += other_head
		tail_sum += other_tail

	output_txt = directory / "attached_density.txt"
	np.savetxt(
		output_txt,
		np.column_stack((position, head_sum, tail_sum)),
		header="position head_density tail_density",
	)

	output_jpg = directory / "attached_density.jpg"
	print(f"Density plot is sum of these files:\n" + "\n".join(str(f) for f in files))

	figure, axis = plt.subplots(figsize=(5, 5))
	axis.plot(position, head_sum, label="Head", color="tab:orange")
	axis.plot(position, tail_sum, label="Tail", color="tab:blue")
	axis.set_xlabel("$r$ (nm)", fontsize=12)
	axis.set_ylabel(r"$\rho(r)$ (kg m$^{-3}$)", fontsize=12)
	axis.set_xlim((0, 3.5))
	axis.set_ylim((0, 500))
	axis.tick_params(axis='both', labelsize=12)
	figure.tight_layout()
	figure.savefig(output_jpg, dpi=300)
	plt.close(figure)


if __name__ == "__main__":
	main()
