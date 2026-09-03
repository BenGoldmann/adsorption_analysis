#!/usr/bin/env python3
import argparse
import pickle
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from read_parameters import read_parameters


def pickle_open(path):
    with path.open("rb") as handle:
        return pickle.load(handle)


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


def plot_attachment_vs_time(final, params, output_dir):
    if not final or not final[0]:
        raise ValueError("Attachment analysis output contains no frames or polymers")

    try:
        surface_area = float(params["SURFACE_AREA"])
    except ValueError:
        raise SystemExit("ERROR: SURFACE_AREA must be a number.")

    surface_coverage_total = []
    attachment_rate_total = []
    
    num_polymers = len(final[0])
    
    # Count attached polymers per frame
    attached_per_frame = []
    for frame in range(len(final)):
        attached_count = sum(1 for pol in range(num_polymers) if final[frame][pol] != [])
        attached_per_frame.append(attached_count)
    
    # Block average
    block_size = 100
    for i in range(0, len(attached_per_frame), block_size):
        block = attached_per_frame[i:i+block_size]
        block_mean = np.mean(block)
        attachment_rate_total.append(block_mean / num_polymers)
        surface_coverage_total.append(block_mean / ((2 * surface_area) / 100))  # Two surfaces, convert area to nm^2

    avg_coverage, coverage_eq_start = find_equilibrium(surface_coverage_total)
    avg_attachment_rate, rate_eq_start = find_equilibrium(attachment_rate_total)
    
    fig, ax = plt.subplots(figsize=(5, 5))

    block_times = (np.arange(len(surface_coverage_total)) + 0.5) * 10
    ax.plot(block_times, surface_coverage_total)
    ax.axhline(avg_coverage, alpha=0.5)
    ax.axvline(block_times[coverage_eq_start], alpha=0.3, linestyle='--', color='gray')
    ax.set_xlabel('Time (ns)', size=12)
    ax.set_ylabel('Surface coverage (polymer/nm$^2$)', size=12)
    ax.tick_params(axis='both', labelsize=12)
    
    # Add text annotations
    ax.text(0.95, 0.05, f'Coverage: {avg_coverage:.2f} pol./nm$^2$\nAttachment rate: {avg_attachment_rate*100:.0f}%',
            transform=ax.transAxes, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=1))

    plt.tight_layout()
    plt.savefig(output_dir / 'attachment_vs_time.jpg', bbox_inches='tight', dpi=1000)

    with open(output_dir / "attachment_statistics.txt", "w") as handle:
        handle.write(
            f'EQUILIBRIUM={coverage_eq_start}\n'
            f'COVERAGE={avg_coverage:.2f}\n'
            f'ATTACHMENT_RATE={avg_attachment_rate * 100:.0f}\n'
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parameters", type=Path, help="Path to parameters.in")
    args = parser.parse_args()
    params = read_parameters(args.parameters)
    output_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
    if not output_dir.is_dir():
        raise SystemExit(f"ERROR: Working directory not found: {output_dir}")

    master_path = output_dir / "attachment_master.pkl"
    if not master_path.is_file():
        raise SystemExit(f"ERROR: Attachment analysis output not found: {master_path}")
    final = pickle_open(master_path)

    plot_attachment_vs_time(final, params, output_dir)

if __name__ == "__main__":
    main()