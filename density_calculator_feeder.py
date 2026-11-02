from pathlib import Path

from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path
from z_scripts.helper_scripts.read_parameters import read_parameters


def build_density_parameters(systeminfo, system_name, working_dir, output_path):
    """Write density-calculator settings using SYSTEMINFO atom labels."""
    name = system_name.lower()
    if name.startswith("batio3_"):
        surface_key = "BATIO3_IDS"
        top_surface_labels = ("Ti", "O")
        bottom_surface_labels = ("Ba", "O")
        symmetric = False
    elif name.startswith("bao_"):
        surface_key = "BATIO3_IDS" if "BATIO3_IDS" in systeminfo else "BAO_IDS"
        top_surface_labels = bottom_surface_labels = ("Ba", "O")
        symmetric = True
    elif name.startswith("tio2_"):
        surface_key = "BATIO3_IDS" if "BATIO3_IDS" in systeminfo else "TIO2_IDS"
        top_surface_labels = bottom_surface_labels = ("Ti", "O")
        symmetric = True
    elif name.startswith("nio_"):
        surface_key = "NIO_IDS"
        top_surface_labels = bottom_surface_labels = tuple(systeminfo.get(surface_key, "").split())
        symmetric = True
    else:
        raise SystemExit(
            f"ERROR: Cannot determine atom types for '{system_name}'. "
            "Expected a directory starting with bao_, batio3_, tio2_, or nio_."
        )

    surface_labels = tuple(systeminfo.get(surface_key, "").split())
    id_keys = (
        "ADDITIVE_HEAD_IDS",
        "ADDITIVE_TAIL_IDS",
        "DHTA_IDS",
        surface_key,
    )
    missing = [key for key in id_keys if key not in systeminfo]
    for key in ("ADDITIVE_ATOMS_NUM", "BATIO3_ADDITIVE_NUM" if surface_key != "NIO_IDS" else "NIO_ADDITIVE_NUM"):
        if key not in systeminfo:
            missing.append(key)
    if missing:
        raise SystemExit(f"ERROR: Missing SYSTEMINFO entries: {', '.join(missing)}")

    for label in (*top_surface_labels, *bottom_surface_labels):
        if label not in surface_labels:
            raise SystemExit(
                f"ERROR: Surface type '{label}' is not listed in {surface_key} for {system_name}"
            )

    atom_type_names = sorted(
        {label for key in id_keys for label in systeminfo[key].split()}
    )
    atom_type_ids = {label: index + 1 for index, label in enumerate(atom_type_names)}

    def resolve_ids(key):
        labels = systeminfo[key].split()
        unresolved = [label for label in labels if label not in atom_type_ids]
        if unresolved:
            raise SystemExit(
                f"ERROR: {key} references unknown atom types: {', '.join(unresolved)}"
            )
        return [atom_type_ids[label] for label in labels]

    missing_elements = [label for label in atom_type_names if label not in systeminfo]
    if missing_elements:
        raise SystemExit(
            "ERROR: Missing element mappings in SYSTEMINFO for atom types: "
            + ", ".join(missing_elements)
        )

    working_dir = Path(working_dir).resolve()
    output_path = Path(output_path)
    additive_num_key = "NIO_ADDITIVE_NUM" if surface_key == "NIO_IDS" else "BATIO3_ADDITIVE_NUM"

    def mapped_surface_ids(labels):
        return sorted(atom_type_ids[label] for label in labels)

    top_ids = mapped_surface_ids(top_surface_labels)
    bottom_ids = mapped_surface_ids(bottom_surface_labels)
    all_surface_ids = sorted(mapped_surface_ids(surface_labels))
    with output_path.open("w") as handle:
        handle.write(f"LAMMPSDUMP_DIR = {working_dir}\n")
        handle.write(f"ADDITIVE_NUM = {systeminfo[additive_num_key]}\n")
        handle.write(f"ADDITIVE_ATOMS_NUM = {systeminfo['ADDITIVE_ATOMS_NUM']}\n")
        handle.write(
            "ADDITIVE_HEAD_ID = "
            + " ".join(map(str, resolve_ids("ADDITIVE_HEAD_IDS")))
            + "\n"
        )
        handle.write(
            "ADDITIVE_TAIL_ID = "
            + " ".join(map(str, resolve_ids("ADDITIVE_TAIL_IDS")))
            + "\n"
        )
        handle.write("SURFACE_ID = " + " ".join(map(str, all_surface_ids)) + "\n")
        handle.write("SURFACE_TOP_ID = " + " ".join(map(str, top_ids)) + "\n")
        handle.write("SURFACE_BOTTOM_ID = " + " ".join(map(str, bottom_ids)) + "\n")
        handle.write(f"SURFACE_SYMMETRIC = {int(symmetric)}\n")
        for label in atom_type_names:
            handle.write(f"ELEM_{atom_type_ids[label]} = {systeminfo[label]}\n")


def main():
    working_dir = Path.cwd()
    systeminfo_path = get_systeminfo_path(working_dir)
    systeminfo = read_parameters(systeminfo_path)
    output_path = working_dir / "density_calculator_parameters.in"
    build_density_parameters(
        systeminfo, working_dir.name, working_dir, output_path
    )
    print(f"Created density calculator parameters: {output_path}")


if __name__ == "__main__":
    main()
