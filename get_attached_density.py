#!/usr/bin/env python3
"""Calculate an attached-polymer density profile from attachment analysis output."""

import argparse
import ast
import pickle
import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from read_parameters import read_parameters


def read_element_masses(path):
    module = ast.parse(path.read_text(), filename=str(path))
    for statement in module.body:
        if isinstance(statement, ast.Assign):
            if any(isinstance(target, ast.Name) and target.id == "ELEMENTS" for target in statement.targets):
                elements = ast.literal_eval(statement.value)
                return {name: float(properties["mass"]) for name, properties in elements.items()}
    raise ValueError(f"ELEMENTS dictionary not found in {path}")


def parse_ids(value):
    return [int(item) for item in re.split(r"[,\s]+", value.strip()) if item]


def read_statistics(working_dir):
    path = working_dir / "attachment_statistics.txt"
    if not path.is_file():
        raise FileNotFoundError(f"Attachment statistics not found: {path}")
    values = {}
    for line in path.read_text().splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key.strip()] = float(value.strip())
    required = ("EQUILIBRIUM", "COVERAGE", "ATTACHMENT_RATE")
    missing = [key for key in required if key not in values]
    if missing:
        raise ValueError(f"Missing {', '.join(missing)} in {path}")
    return int(values["EQUILIBRIUM"]), values["COVERAGE"], values["ATTACHMENT_RATE"]


def plot_density(final, equilibrium_start, avg_coverage, avg_attachment_rate, params, output_dir):
    additive_head_id = parse_ids(params["ADDITIVE_HEAD_ID"])
    additive_tail_id = parse_ids(params["ADDITIVE_TAIL_ID"])
    surface_area = float(params["SURFACE_AREA"])
    element_masses = read_element_masses(Path(__file__).with_name("elements.data"))
    atom_type_to_element = {
        int(key[5:]): value for key, value in params.items()
        if key.startswith("ELEM_") and value in element_masses
    }

    if not final or equilibrium_start >= len(final):
        raise ValueError("EQUILIBRIUM is beyond the available attachment frames")

    hist_binsize = 0.2
    num_equilibrium_frames = len(final) - equilibrium_start
    scaling_factor = 2 * surface_area * hist_binsize * num_equilibrium_frames / 10**30
    head_positions, head_weights = [], []
    tail_positions, tail_weights = [], []

    for frame in final[equilibrium_start:]:
        for polymer in frame:
            for atom_type, z_position in polymer:
                atom_type = int(atom_type)
                element = atom_type_to_element.get(atom_type)
                if element is None:
                    continue
                weight = element_masses[element] * 1.66054e-27 / scaling_factor
                if atom_type in additive_head_id:
                    head_positions.append(z_position)
                    head_weights.append(weight)
                elif atom_type in additive_tail_id:
                    tail_positions.append(z_position)
                    tail_weights.append(weight)

    n_bins = int(40 / hist_binsize)
    head_density, edges = np.histogram(head_positions, bins=n_bins, weights=head_weights, range=(0, 40))
    tail_density, _ = np.histogram(tail_positions, bins=n_bins, weights=tail_weights, range=(0, 40))
    centers = (edges[:-1] + edges[1:]) / 2 / 10

    figure, axis = plt.subplots(figsize=(5, 5))
    axis.plot(centers, head_density, label="Head", color="tab:orange")
    axis.plot(centers, tail_density, label="Tail", color="tab:blue")
    axis.set_xlabel("$r$ (nm)", fontsize=12)
    axis.set_ylabel(r"$\rho(r)$ (kg m$^{-3}$)", fontsize=12)
    axis.set_xlim((0, 3.5))
    axis.set_ylim((0, 500))
    axis.tick_params(axis='both', labelsize=12)
    axis.text(
        0.95, 0.95,
        f"Coverage: {avg_coverage:.2f} pol./nm$^2$\nAttachment rate: {avg_attachment_rate:.0f}%",
        transform=axis.transAxes, ha="right", va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=1), fontsize=12
    )
    figure.tight_layout()
    figure.savefig(output_dir / "attached_density.jpg", bbox_inches="tight", dpi=1000)
    plt.close(figure)

    np.savetxt(
        output_dir / "attached_density.txt",
        np.column_stack((centers, head_density, tail_density)),
        header="r (nm)  rho_head (kg/m^3)  rho_tail (kg/m^3)",
        fmt=("%.4f", "%.6e", "%.6e"),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parameters", type=Path, help="Path to parameters.in")
    args = parser.parse_args()
    params = read_parameters(args.parameters)
    required = ("LAMMPSDUMP_DIR", "ADDITIVE_HEAD_ID", "ADDITIVE_TAIL_ID", "SURFACE_AREA")
    missing = [key for key in required if key not in params]
    if missing:
        raise SystemExit(f"ERROR: Missing parameters: {', '.join(missing)}")
    working_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
    if not working_dir.is_dir():
        raise SystemExit(f"ERROR: Working directory not found: {working_dir}")

    statistics = read_statistics(working_dir)
    with (working_dir / "attachment_master.pkl").open("rb") as handle:
        final = pickle.load(handle)
    plot_density(final, *statistics, params, working_dir)


if __name__ == "__main__":
    main()