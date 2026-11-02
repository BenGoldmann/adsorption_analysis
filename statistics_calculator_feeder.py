from pathlib import Path

from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path
from z_scripts.helper_scripts.read_parameters import read_parameters


def build_statistics_parameters(systeminfo, working_dir, output_path):
    """Write the minimal parameters needed by statistics_calculator.py."""
    working_dir = Path(working_dir).resolve()
    name = working_dir.name.lower()
    if name.startswith(("bao_", "batio3_", "tio2_")):
        area_key = "BAO_SURFACE_AREA"
    elif name.startswith("nio_"):
        area_key = "NIO_SURFACE_AREA"
    else:
        raise SystemExit(
            f"ERROR: Cannot determine surface area for '{working_dir.name}'. "
            "Expected a directory starting with bao_, batio3_, tio2_, or nio_."
        )

    if area_key not in systeminfo:
        raise SystemExit(f"ERROR: {area_key} not found in SYSTEMINFO")

    try:
        surface_area = float(systeminfo[area_key])
    except ValueError as exc:
        raise SystemExit(f"ERROR: {area_key} must be a number") from exc
    if surface_area <= 0:
        raise SystemExit(f"ERROR: {area_key} must be positive")

    with Path(output_path).open("w") as handle:
        handle.write(f"LAMMPSDUMP_DIR = {working_dir}\n")
        handle.write(f"SURFACE_AREA = {surface_area:g}\n")


def main():
    working_dir = Path.cwd()
    systeminfo_path = get_systeminfo_path(working_dir)
    systeminfo = read_parameters(systeminfo_path)
    output_path = working_dir / "statistics_calculator_parameters.in"
    build_statistics_parameters(systeminfo, working_dir, output_path)
    print(f"Created statistics calculator parameters: {output_path}")


if __name__ == "__main__":
    main()
