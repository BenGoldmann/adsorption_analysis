#!/usr/bin/env python3
"""Plot attached-polymer density profiles for NiO and BaTiO3 surfaces."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_density(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
	"""Read radius, head-density, and tail-density columns from a density file."""
	if not path.is_file():
		raise FileNotFoundError(f"Density file not found: {path}")

	data = np.loadtxt(path, comments="#", ndmin=2)
	if data.shape[1] < 3:
		raise ValueError(f"Expected at least 3 columns in {path}")
	return data[:, 0], data[:, 1], data[:, 2]


def read_attachment_statistics(directory: Path) -> tuple[float, float]:
	"""Read coverage and attachment rate from a surface's statistics file."""
	for filename in ("attachment_statistics.txt", "attachment_stats.txt"):
		path = directory / filename
		if path.is_file():
			statistics = {}
			for line in path.read_text().splitlines():
				key, separator, value = line.partition("=")
				if separator:
					statistics[key.strip()] = float(value.strip())
			try:
				return statistics["COVERAGE"], statistics["ATTACHMENT_RATE"]
			except KeyError as error:
				raise ValueError(f"Missing {error.args[0]} in {path}") from error

	raise FileNotFoundError(
		f"No attachment_statistics.txt or attachment_stats.txt found in {directory}"
	)


def plot_density_pair(root: Path, polymer: str, output: Path) -> None:
	"""Create a two-row NiO/BaTiO3 density comparison for one polymer."""
	surfaces = ("NiO", "BaTiO3")
	polymer_directory = polymer if polymer.endswith("_dhta") else f"{polymer}_dhta"
	paths = {
		surface: root / f"{surface.lower()}_{polymer_directory}" / "attached_density.txt"
		for surface in surfaces
	}
	profiles = {surface: read_density(path) for surface, path in paths.items()}
	statistics = {
		surface: read_attachment_statistics(paths[surface].parent) for surface in surfaces
	}

	figure, axes = plt.subplots(
		nrows=1,
		ncols=2,
		sharex=True,
		sharey=True,
		figsize=(5, 5),
		constrained_layout=True,
	)

	for axis, surface in zip(axes, surfaces):
		radius, head, tail = profiles[surface]
		coverage, attachment_rate = statistics[surface]
		axis.plot(radius, head, label="Head", color="tab:orange")
		axis.plot(radius, tail, label="Tail", color="tab:blue")
		axis.set_title(surface, size=12)
		axis.grid(True, alpha=0.25)
		axis.tick_params(axis='both', labelsize=12)
		axis.set_xlim(0, 3.5)
		axis.set_ylim(0, 500)
		axis.text(
			0.95,
			0.975,
			f"Coverage: {coverage:.2f} pol./nm$^2$\n"
			f"Attachment rate: {attachment_rate:.0f}%",
			transform=axis.transAxes,
			verticalalignment="top",
			horizontalalignment="right",
			bbox=dict(boxstyle="round", facecolor="wheat", alpha=1),
		)

	axes[0].set_xlabel("$r$ (nm)", size=12)
	axes[1].set_xlabel("$r$ (nm)", size=12)
	axes[0].set_ylabel(r"$\rho(r)$ (kg m$^{-3}$)", size=12)
	figure.savefig(output, dpi=1000, bbox_inches="tight")
	plt.close(figure)


def main() -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument(
		"polymer",
		help="Polymer name, for example pibsa or dpibsa (optional _dhta suffix)",
	)
	parser.add_argument(
		"--root",
		type=Path,
		default=Path(__file__).resolve().parents[1],
		help="Project root (default: parent of z_scripts)",
	)
	parser.add_argument(
		"-o",
		"--output",
		type=Path,
		default=None,
		help="Output image path (default: attached_density_<polymer>.png in project root)",
	)
	args = parser.parse_args()
	output = args.output or args.root / f"attached_density_{args.polymer}.png"
	plot_density_pair(args.root, args.polymer, output)
	print(f"Saved plot to {output}")


if __name__ == "__main__":
	main()
