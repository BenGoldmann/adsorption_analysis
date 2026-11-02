#!/usr/bin/env python3
"""Split density and attachment data into files grouped by surface material."""

from pathlib import Path

import numpy as np


INPUT_COLUMNS = {
    "attachment_vs_time.txt": 4,
    "total_density.txt": 10,
    "attached_density.txt": 10,
}


def load_table(path, expected_columns):
    if not path.is_file():
        raise SystemExit(f"ERROR: Required input file not found: {path}")
    try:
        data = np.loadtxt(path, comments="#", ndmin=2)
    except ValueError as exc:
        raise SystemExit(f"ERROR: Could not parse numeric data in {path}: {exc}") from exc
    if data.shape[1] < expected_columns:
        raise SystemExit(
            f"ERROR: Expected at least {expected_columns} columns in {path}, "
            f"found {data.shape[1]}"
        )
    return data


def write_selected_columns(output_dir, output_name, data, columns, header):
    np.savetxt(
        output_dir / output_name,
        data[:, columns],
        delimiter="\t",
        header=header,
        fmt="%.8g",
    )


def sort_surface_data(working_dir):
    working_dir = Path(working_dir).resolve()
    name = working_dir.name.lower()
    if name.startswith("batio3_"):
        surface_type = "batio3"
    elif name.startswith("bao_"):
        surface_type = "bao"
    elif name.startswith("tio2_"):
        surface_type = "tio2"
    elif name.startswith("nio_"):
        surface_type = "nio"
    else:
        raise SystemExit(
            f"ERROR: Unrecognized simulation directory '{working_dir.name}'. "
            "Expected a name starting with nio_, bao_, tio2_, or batio3_."
        )

    attachment = load_table(working_dir / "attachment_vs_time.txt", 4)
    total_density = load_table(working_dir / "total_density.txt", 10)
    attached_density = load_table(working_dir / "attached_density.txt", 10)

    if surface_type == "batio3":
        outputs = (
            ("tio2_attachment_vs_time.txt", attachment, (0, 2), "step\ttop_attachment"),
            ("bao_attachment_vs_time.txt", attachment, (0, 3), "step\tbottom_attachment"),
            (
                "tio2_total_density.txt",
                total_density,
                (0, 4, 5, 6),
                "distance_nm\ttotal_kg_m3\thead_kg_m3\ttail_kg_m3",
            ),
            (
                "bao_total_density.txt",
                total_density,
                (0, 7, 8, 9),
                "distance_nm\ttotal_kg_m3\thead_kg_m3\ttail_kg_m3",
            ),
            (
                "tio2_attached_density.txt",
                attached_density,
                (0, 4, 5, 6),
                "distance_nm\ttotal_kg_m3\thead_kg_m3\ttail_kg_m3",
            ),
            (
                "bao_attached_density.txt",
                attached_density,
                (0, 7, 8, 9),
                "distance_nm\ttotal_kg_m3\thead_kg_m3\ttail_kg_m3",
            ),
        )
    else:
        prefix = surface_type
        outputs = (
            (f"{prefix}_attachment_vs_time.txt", attachment, (0, 1), "step\ttotal_attachment"),
            (
                f"{prefix}_total_density.txt",
                total_density,
                (0, 1, 2, 3),
                "distance_nm\ttotal_kg_m3\thead_kg_m3\ttail_kg_m3",
            ),
            (
                f"{prefix}_attached_density.txt",
                attached_density,
                (0, 1, 2, 3),
                "distance_nm\ttotal_kg_m3\thead_kg_m3\ttail_kg_m3",
            ),
        )

    for output_name, data, columns, header in outputs:
        write_selected_columns(working_dir, output_name, data, columns, header)
    return [working_dir / output_name for output_name, *_ in outputs]


def main():
    output_paths = sort_surface_data(Path.cwd())
    for path in output_paths:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
