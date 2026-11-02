#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np

from z_scripts.helper_scripts.read_parameters import read_parameters
from z_scripts.helper_scripts.open_pickle import pickle_open
from z_scripts.helper_scripts.read_dump_steps import read_dump_steps


def find_equilibrium(data, threshold=0.0025):
    """
    Find equilibrium by calculating running means from different starting points.
    Returns the mean and starting frame when the value stabilizes.
    
    Args:
        data: list of values to analyze
        threshold: relative change threshold to consider equilibrated (default 0.1%)
    
    Returns:
        equilibrium_value: the stable average value
        equilibrium_start: the frame index where equilibrium begins
    """
    if not data:
        raise ValueError("Cannot find equilibrium in empty data")
    if len(data) < 2:
        return np.mean(data), 0
    
    means = []
    
    # Calculate means starting from different points to the end
    for start in range(len(data)):
        mean = np.mean(data[start:])
        means.append(mean)
    
    # Find when means stabilize
    for i in range(1, len(means)):
        if means[i-1] != 0:
            relative_change = abs(means[i] - means[i-1]) / abs(means[i-1])
        else:
            relative_change = abs(means[i] - means[i-1])
        
        if relative_change < threshold:
            # Found equilibrium
            return means[i], i
    
    # If no equilibrium found, return the last mean and its starting point
    return means[-1], len(means) - 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "parameters",
        type=Path,
        help="Path to statistics_calculator_parameters.in",
    )
    args = parser.parse_args()
    params = read_parameters(args.parameters)
    for key in ("LAMMPSDUMP_DIR", "SURFACE_AREA"):
        if key not in params:
            raise SystemExit(f"ERROR: {key} not found in {args.parameters}")
    output_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
    if not output_dir.is_dir():
        raise SystemExit(f"ERROR: Working directory not found: {output_dir}")
    try:
        surface_area = float(params["SURFACE_AREA"])
    except ValueError as exc:
        raise SystemExit("ERROR: SURFACE_AREA must be a number.") from exc
    if surface_area <= 0:
        raise SystemExit("ERROR: SURFACE_AREA must be positive.")

    master_path = output_dir / "attachment_master.pkl"
    if not master_path.is_file():
        raise SystemExit(f"ERROR: Attachment analysis output not found: {master_path}")
    final = pickle_open(master_path)
    surface_path = output_dir / "attachment_surface.pkl"
    if not surface_path.is_file():
        raise SystemExit(
            f"ERROR: Surface labels not found: {surface_path}. "
            "Run attachment_analysis.py again to generate them."
        )
    surfaces = pickle_open(surface_path)

    if not final or not final[0]:
        raise SystemExit("ERROR: Attachment analysis output contains no frames or polymers")
    if len(surfaces) != len(final) or any(
        len(surface_frame) != len(attachment_frame)
        for surface_frame, attachment_frame in zip(surfaces, final)
    ):
        raise SystemExit("ERROR: Attachment surface labels do not match attachment data")

    dump_path = output_dir / "dump_nosolvent.lammpstrj"
    if not dump_path.is_file():
        raise SystemExit(f"ERROR: Trajectory file not found: {dump_path}")
    steps = read_dump_steps(dump_path)
    if len(steps) != len(final):
        raise SystemExit(
            f"ERROR: Found {len(steps)} trajectory timesteps but {len(final)} attachment frames"
        )

    n_polymers = len(final[0])
    total_attached = []
    top_attached = []
    bottom_attached = []
    for attachment_frame, surface_frame in zip(final, surfaces):
        bottom_count = sum(
            attachment != [] and surface == "bottom"
            for attachment, surface in zip(attachment_frame, surface_frame)
        )
        top_count = sum(
            attachment != [] and surface == "top"
            for attachment, surface in zip(attachment_frame, surface_frame)
        )
        bottom_attached.append(bottom_count)
        top_attached.append(top_count)
        total_attached.append(bottom_count + top_count)

    # Keep the historical 100-frame block averaging used for equilibrium statistics.
    block_size = 100
    coverage_total = []
    coverage_bottom = []
    coverage_top = []
    attachment_rate = []
    attachment_rate_bottom = []
    attachment_rate_top = []
    for start in range(0, len(total_attached), block_size):
        stop = start + block_size
        total_mean = np.mean(total_attached[start:stop])
        coverage_total.append(total_mean / ((2 * surface_area) / 100))
        bottom_mean = np.mean(bottom_attached[start:stop])
        top_mean = np.mean(top_attached[start:stop])
        coverage_bottom.append(bottom_mean / (surface_area / 100))
        coverage_top.append(top_mean / (surface_area / 100))
        attachment_rate.append(total_mean / n_polymers)
        attachment_rate_bottom.append(bottom_mean / n_polymers)
        attachment_rate_top.append(top_mean / n_polymers)

    avg_coverage, equilibrium_block = find_equilibrium(coverage_total)
    avg_top_coverage, top_equilibrium_block = find_equilibrium(coverage_top)
    avg_bottom_coverage, bottom_equilibrium_block = find_equilibrium(coverage_bottom)
    avg_top_coverage = float(avg_top_coverage)
    avg_bottom_coverage = float(avg_bottom_coverage)
    if min(avg_top_coverage, avg_bottom_coverage) == 0:
        top_bottom_ratio = float("inf")
    else:
        top_bottom_ratio = max(avg_top_coverage, avg_bottom_coverage) / min(
            avg_top_coverage, avg_bottom_coverage
        )
    avg_attachment_rate, _ = find_equilibrium(attachment_rate)
    avg_top_attachment_rate, _ = find_equilibrium(attachment_rate_top)
    avg_bottom_attachment_rate, _ = find_equilibrium(attachment_rate_bottom)

    statistics_path = output_dir / "attachment_statistics.txt"
    with statistics_path.open("w") as handle:
        handle.write(
            f"EQUILIBRIUM={equilibrium_block * block_size}\n"
            f"TOP_EQUILIBRIUM={top_equilibrium_block * block_size}\n"
            f"BOTTOM_EQUILIBRIUM={bottom_equilibrium_block * block_size}\n"
            f"COVERAGE={avg_coverage:.2f}\n"
            f"TOP_COVERAGE={avg_top_coverage:.2f}\n"
            f"BOTTOM_COVERAGE={avg_bottom_coverage:.2f}\n"
            f"TOP_BOTTOM_RATIO={top_bottom_ratio:.2f}\n"
            f"ATTACHMENT_RATE={avg_attachment_rate * 100:.0f}\n"
            f"TOP_ATTACHMENT_RATE={avg_top_attachment_rate * 100:.0f}\n"
            f"BOTTOM_ATTACHMENT_RATE={avg_bottom_attachment_rate * 100:.0f}\n"
        )

    time_series_path = output_dir / "attachment_vs_time.txt"
    with time_series_path.open("w") as handle:
        handle.write("#step\ttotal_attachment\ttop_attachment\tbottom_attachment\n")
        for step, total, top, bottom in zip(
            steps, total_attached, top_attached, bottom_attached
        ):
            handle.write(f"{step}\t{total}\t{top}\t{bottom}\n")

    print(f"Wrote attachment statistics to {statistics_path}")
    print(f"Wrote attachment time series to {time_series_path}")


if __name__ == "__main__":
    main()