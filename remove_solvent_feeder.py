from pathlib import Path

from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path
from z_scripts.helper_scripts.read_parameters import read_parameters


SURFACE_TYPE_KEYS = {
    "nio": "NIO_IDS",
    "bao": "BATIO3_IDS",
    "tio2": "BATIO3_IDS",
    "batio3": "BATIO3_IDS",
}
SOLVENT_TYPE_KEYS = {
    "dhta": "DHTA_IDS",
    "terpineol": "TERPINEOL_IDS",
}


def build_remove_solvent_parameters(systeminfo, system_name, working_dir, output_path):
    """Write the dump-file path and non-solvent LAMMPS type IDs."""
    working_dir = Path(working_dir).resolve()
    name = system_name.lower()
    surface_prefix = next(
        (prefix for prefix in SURFACE_TYPE_KEYS if name.startswith(f"{prefix}_")),
        None,
    )
    solvent = next(
        (key for key in SOLVENT_TYPE_KEYS if name.endswith(f"_{key}")),
        None,
    )
    if surface_prefix is None:
        raise SystemExit(
            f"ERROR: Cannot determine surface for '{system_name}'. Expected nio_, bao_, "
            "tio2_, or batio3_."
        )
    if solvent is None:
        raise SystemExit(
            f"ERROR: Cannot determine solvent from '{system_name}'. "
            "Expected a _dhta or _terpineol suffix."
        )

    surface_key = SURFACE_TYPE_KEYS[surface_prefix]
    solvent_key = SOLVENT_TYPE_KEYS[solvent]
    id_keys = ("ADDITIVE_HEAD_IDS", "ADDITIVE_TAIL_IDS", surface_key, solvent_key)
    missing = [key for key in id_keys if key not in systeminfo]
    if missing:
        raise SystemExit(f"ERROR: Missing SYSTEMINFO entries: {', '.join(missing)}")

    atom_types = sorted({label for key in id_keys for label in systeminfo[key].split()})
    solvent_types = set(systeminfo[solvent_key].split())
    keep_types = [
        index + 1 for index, label in enumerate(atom_types) if label not in solvent_types
    ]

    dump_file = working_dir / "dump.lammpstrj"
    if not dump_file.is_file():
        raise SystemExit(f"ERROR: LAMMPS dump file not found: {dump_file}")

    with Path(output_path).open("w") as handle:
        handle.write(f"DUMP_FILE = {dump_file}\n")
        handle.write("KEEP_TYPES = " + " ".join(map(str, keep_types)) + "\n")


def main():
    working_dir = Path.cwd()
    systeminfo_path = get_systeminfo_path(working_dir)
    systeminfo = read_parameters(systeminfo_path)
    output_path = working_dir / "remove_solvent_parameters.in"
    build_remove_solvent_parameters(systeminfo, working_dir.name, working_dir, output_path)
    print(f"Created solvent-removal parameters: {output_path}")


if __name__ == "__main__":
    main()
