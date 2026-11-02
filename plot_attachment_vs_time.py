#!/usr/bin/env python3
"""Plot sorted attachment-versus-time data for all simulations under SYSTEMINFO."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path
from z_scripts.helper_scripts.read_parameters import read_parameters


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
TIMESTEP_NS = 2e-6  # 2 fs per simulation step.
FONT_SIZE = 12


def get_solvent(simulation_dir):
    name = simulation_dir.name.lower()
    for solvent in SOLVENT_STYLES:
        if name.endswith(f"_{solvent}"):
            return solvent
    return None


def get_simulation_dir(path, search_root):
    """Find the simulation ancestor that names both surface and solvent."""
    for parent in (path.parent, *path.parent.parents):
        if parent == search_root.parent:
            break
        name = parent.name.lower()
        has_surface = name.startswith(("nio_", "bao_", "tio2_", "batio3_"))
        has_solvent = any(name.endswith(f"_{solvent}") for solvent in SOLVENT_STYLES)
        if has_surface and has_solvent:
            return parent
        if parent == search_root:
            break
    raise SystemExit(f"ERROR: Could not identify simulation directory for {path}")


def surface_area_nm2(systeminfo, surface):
    area_key = "NIO_SURFACE_AREA" if surface == "nio" else "BAO_SURFACE_AREA"
    if area_key not in systeminfo:
        raise SystemExit(f"ERROR: {area_key} is required in SYSTEMINFO for {surface} coverage")
    try:
        area = float(systeminfo[area_key]) / 100.0  # SYSTEMINFO area is in Å².
    except ValueError as exc:
        raise SystemExit(f"ERROR: {area_key} must be a number") from exc
    if area <= 0:
        raise SystemExit(f"ERROR: {area_key} must be positive")
    return area


def load_attachment_data(path):
    try:
        data = np.loadtxt(path, comments="#", ndmin=2)
    except ValueError as exc:
        raise SystemExit(f"ERROR: Could not parse attachment data in {path}: {exc}") from exc
    if data.shape[1] < 2:
        raise SystemExit(f"ERROR: Expected at least two columns in {path}")
    return data[:, :2]


def find_attachment_files(search_root):
    """Find sorted attachment files below the directory containing SYSTEMINFO."""
    found = []
    for path in search_root.rglob("*_attachment_vs_time.txt"):
        prefix = path.name.removesuffix("_attachment_vs_time.txt").lower()
        if prefix not in SURFACE_COLORS:
            continue
        solvent = get_solvent(path.parent)
        if solvent is None:
            continue
        found.append((prefix, solvent, path))
    return sorted(
        found, key=lambda item: (item[1], SURFACE_ORDER.index(item[0]), str(item[2]))
    )


def plot_attachment_vs_time(search_root, output_path=None):
    search_root = Path(search_root).resolve()
    systeminfo = read_parameters(get_systeminfo_path(search_root))
    records = find_attachment_files(search_root)
    if not records:
        raise SystemExit(
            f"ERROR: No nio_/bao_/tio2_ attachment time-series files found below {search_root}"
        )

    figure, axis = plt.subplots(figsize=(5, 5))
    seen_labels = set()
    max_time_ns = 0.0
    for surface, solvent, path in records:
        data = load_attachment_data(path)
        time_ns = data[:, 0] * TIMESTEP_NS
        max_time_ns = max(max_time_ns, time_ns.max())
        simulation_dir = get_simulation_dir(path, search_root)
        # A symmetric slab has two equivalent surfaces; divide total attached
        # polymers across both faces. BaTiO3 files are already top or bottom only.
        is_batio3 = simulation_dir.name.lower().startswith("batio3_")
        divisor = 1 if is_batio3 else 2
        area_nm2 = surface_area_nm2(systeminfo, surface)
        coverage = data[:, 1] / (divisor * area_nm2)
        label = f"{SURFACE_LABELS[surface]} ({SOLVENT_LABELS[solvent]})"
        if label in seen_labels:
            label = "_nolegend_"
        else:
            seen_labels.add(label)
        axis.plot(
            time_ns,
            coverage,
            color=SURFACE_COLORS[surface],
            linestyle=SOLVENT_STYLES[solvent],
            label=label,
        )

    axis.set_xlabel("Simulation time (ns)", fontsize=FONT_SIZE)
    axis.set_ylabel("Surface coverage (polymers/nm$^2$)", fontsize=FONT_SIZE)
    axis.set_xlim(0, max_time_ns)
    axis.set_ylim(bottom=0)
    axis.legend(title="Surface (solvent)", fontsize=FONT_SIZE, title_fontsize=FONT_SIZE, loc="upper left")
    axis.tick_params(axis="both", labelsize=FONT_SIZE)
    axis.xaxis.get_offset_text().set_fontsize(FONT_SIZE)
    axis.yaxis.get_offset_text().set_fontsize(FONT_SIZE)
    figure.tight_layout()

    if output_path is None:
        output_path = search_root / "attachment_vs_time.jpg"
    output_path = Path(output_path)
    figure.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(figure)
    print(f"Plotted {len(records)} attachment time series to {output_path}")
    return output_path


def main():
    systeminfo = get_systeminfo_path(Path.cwd())
    plot_attachment_vs_time(systeminfo.parent)


if __name__ == "__main__":
    main()