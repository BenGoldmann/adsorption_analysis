#!/usr/bin/env python3
"""Plot total and attached density profiles grouped by solvent."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path


SURFACE_COLORS = {
	"nio": "#1e1e1e",
	"bao": "#969696",
	"tio2": "#9696b4",
}
SOLVENT_STYLES = {
	"dhta": "-",
	"terpineol": "--",
}
SURFACE_LABELS = {
	"nio": "NiO",
	"bao": "BaO",
	"tio2": "TiO$_2$",
}
SOLVENT_LABELS = {
	"dhta": "DHTA",
	"terpineol": "terpineol",
}
SURFACE_ORDER = ("nio", "bao", "tio2")
HEAD_COLOR = "tab:orange"
TAIL_COLOR = "tab:blue"
FONT_SIZE = 12
X_LIMITS = (0.0, 3.5)  # nm
Y_LIMITS = (0.0, 500.0)  # kg/m^3


def find_density_files(search_root, density_kind):
	"""Return (surface, solvent, path) records for one density profile kind."""
	records = []
	suffix = f"_{density_kind}_density.txt"
	for path in search_root.rglob(f"*{suffix}"):
		surface = path.name.removesuffix(suffix).lower()
		if surface not in SURFACE_COLORS:
			continue

		solvent = None
		for parent in (path.parent, *path.parent.parents):
			if parent == search_root.parent:
				break
			directory_name = parent.name.lower()
			solvent = next(
				(
					candidate
					for candidate in SOLVENT_STYLES
					if directory_name.endswith(f"_{candidate}")
				),
				None,
			)
			if solvent is not None or parent == search_root:
				break
		if solvent is not None:
			records.append((surface, solvent, path))
	return sorted(records, key=lambda record: (record[1], record[0], str(record[2])))


def load_density(path):
	try:
		data = np.loadtxt(path, comments="#", ndmin=2)
	except ValueError as exc:
		raise SystemExit(f"ERROR: Could not parse density data in {path}: {exc}") from exc
	if data.shape[1] < 4:
		raise SystemExit(
			f"ERROR: Expected distance, total, head, and tail columns in {path}; "
			f"found {data.shape[1]} columns"
		)
	return data[:, :4]


def get_simulation_dir(path, search_root):
	"""Find the simulation ancestor that contains attachment_statistics.txt."""
	for parent in (path.parent, *path.parent.parents):
		if parent == search_root.parent:
			break
		if (parent / "attachment_statistics.txt").is_file():
			return parent
		if parent == search_root:
			break
	raise SystemExit(f"ERROR: Could not find attachment statistics for {path}")


def read_attachment_summary(simulation_dir, surface):
	"""Return post-equilibrium coverage and the applicable attachment rate."""
	stats_path = simulation_dir / "attachment_statistics.txt"
	if not stats_path.is_file():
		raise SystemExit(f"ERROR: Attachment statistics not found: {stats_path}")
	statistics = {}
	for line in stats_path.read_text().splitlines():
		key, separator, value = line.partition("=")
		if separator:
			try:
				statistics[key.strip()] = float(value.strip())
			except ValueError as exc:
				raise SystemExit(f"ERROR: Invalid value for {key.strip()} in {stats_path}") from exc

	is_batio3 = simulation_dir.name.lower().startswith("batio3_")
	if is_batio3:
		if surface == "tio2":
			coverage_key, rate_key = "TOP_COVERAGE", "TOP_ATTACHMENT_RATE"
		elif surface == "bao":
			coverage_key, rate_key = "BOTTOM_COVERAGE", "BOTTOM_ATTACHMENT_RATE"
		else:
			raise SystemExit(f"ERROR: Unexpected BaTiO3 surface prefix '{surface}'")
		missing = [key for key in (coverage_key, rate_key) if key not in statistics]
		if missing:
			raise SystemExit(f"ERROR: Missing {', '.join(missing)} in {stats_path}")
		coverage = statistics[coverage_key]
		# The opposite face shares the same dispersant pool; estimate total
		# attachment relative to that pool by doubling this face's rate.
		rate = min(100.0, 2.0 * statistics[rate_key])
	else:
		missing = [key for key in ("COVERAGE", "ATTACHMENT_RATE") if key not in statistics]
		if missing:
			raise SystemExit(f"ERROR: Missing {', '.join(missing)} in {stats_path}")
		coverage = statistics["COVERAGE"]
		rate = statistics["ATTACHMENT_RATE"]
	return coverage, rate


def plot_solvent_profiles(records, solvent, density_kind, search_root, output_path):
	solvent_records = [record for record in records if record[1] == solvent]
	if not solvent_records:
		return None
	solvent_records.sort(key=lambda record: SURFACE_ORDER.index(record[0]))

	subplot_count = len(solvent_records)
	figure, axes = plt.subplots(
		1,
		subplot_count,
		# Sized so each axes matches the attachment_vs_time axes height (~4.1 in),
		# keeping 12 pt text at the same size relative to the plot. The extra
		# width is room for the y-axis label on the leftmost subplot only.
		figsize=(2.5 * subplot_count, 5),
		squeeze=False,
	)

	for axis, (surface, _, path) in zip(axes[0], solvent_records):
		data = load_density(path)
		simulation_dir = get_simulation_dir(path, search_root)
		coverage, attachment_rate = read_attachment_summary(simulation_dir, surface)
		axis.plot(data[:, 0], data[:, 2], color=HEAD_COLOR, linestyle=SOLVENT_STYLES[solvent])
		axis.plot(data[:, 0], data[:, 3], color=TAIL_COLOR, linestyle=SOLVENT_STYLES[solvent])
		axis.set_box_aspect(2)  # Height is twice the width.
		axis.set_title(
			f"{SURFACE_LABELS[surface]} ({SOLVENT_LABELS[solvent]})", fontsize=FONT_SIZE
		)
		axis.text(
			0.5,
			0.93,
			f"Coverage: {coverage:.2f} pol./nm$^2$\n"
			f"Attachment rate: {attachment_rate:.0f}%",
			transform=axis.transAxes,
			ha="center",
			va="center",
			multialignment="center",
			fontsize=10,
			bbox={"boxstyle": "round", "facecolor": "wheat", "edgecolor": "black"},
		)
		axis.set_xlim(*X_LIMITS)
		axis.set_ylim(*Y_LIMITS)
		axis.set_xlabel("$r$ (nm)", fontsize=FONT_SIZE)
		if axis is axes[0][0]:
			axis.set_ylabel(r"$\rho(r)$ (kg/m$^3$)", fontsize=FONT_SIZE)
		else:
			axis.tick_params(axis="y", labelleft=False)
		axis.grid(True, color="0.85", linewidth=0.8)
		axis.set_axisbelow(True)
		axis.set_xticks([0, 1, 2, 3])
		axis.tick_params(axis="both", labelsize=FONT_SIZE)

		for spine in axis.spines.values():
			spine.set_color(SURFACE_COLORS[surface])
			spine.set_linewidth(1.5)

	figure.tight_layout()
	figure.savefig(output_path, bbox_inches="tight", dpi=300)
	plt.close(figure)
	print(f"Plotted {subplot_count} {solvent} {density_kind}-density profile(s) to {output_path}")
	return output_path


def plot_density(search_root):
	search_root = Path(search_root).resolve()
	output_paths = []
	found_any = False
	for density_kind in ("attached", "total"):
		records = find_density_files(search_root, density_kind)
		for solvent in SOLVENT_STYLES:
			output_path = search_root / f"{density_kind}_density_{solvent}.jpg"
			result = plot_solvent_profiles(
				records, solvent, density_kind, search_root, output_path
			)
			if result is not None:
				found_any = True
				output_paths.append(result)
	if not found_any:
		raise SystemExit(
			f"ERROR: No supported *_attached_density.txt or *_total_density.txt files found below {search_root}"
		)
	return output_paths


def main():
	systeminfo = get_systeminfo_path(Path.cwd())
	plot_density(systeminfo.parent)


if __name__ == "__main__":
	main()
