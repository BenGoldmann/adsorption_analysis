from pathlib import Path

from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path
from z_scripts.helper_scripts.read_parameters import read_parameters


def build_isolate_parameters(systeminfo, system_name, working_dir, output_path):
	"""Write the settings needed to isolate the most-saturated surface frame."""
	name = system_name.lower()
	if name.endswith("_dhta"):
		additive_id_key = "DHTA_IDS"
	elif name.endswith("_terpineol"):
		additive_id_key = "TERPINEOL_IDS"
	else:
		raise SystemExit(
			f"ERROR: Cannot determine additive from '{system_name}'. "
			"Expected a _dhta or _terpineol suffix."
		)

	if name.startswith("batio3_"):
		surface_key = "BATIO3_IDS"
		top_labels = bottom_labels = ("Ba", "Ti", "O")
		output_surfaces = ("tio2", "bao")
	elif name.startswith("bao_"):
		surface_key = "BATIO3_IDS" if "BATIO3_IDS" in systeminfo else "BAO_IDS"
		top_labels = bottom_labels = ("Ba", "Ti", "O")
		output_surfaces = ("bao",)
	elif name.startswith("tio2_"):
		surface_key = "BATIO3_IDS" if "BATIO3_IDS" in systeminfo else "TIO2_IDS"
		top_labels = bottom_labels = ("Ba", "Ti", "O")
		output_surfaces = ("tio2",)
	elif name.startswith("nio_"):
		surface_key = "NIO_IDS"
		top_labels = bottom_labels = tuple(systeminfo.get(surface_key, "").split())
		output_surfaces = ("nio",)
	else:
		raise SystemExit(
			f"ERROR: Cannot determine surface for '{system_name}'. "
			"Expected a directory starting with nio_, bao_, tio2_, or batio3_."
		)

	keys = ("ADDITIVE_HEAD_IDS", "ADDITIVE_TAIL_IDS", additive_id_key, surface_key)
	missing = [key for key in keys if key not in systeminfo]
	additive_num_key = "NIO_ADDITIVE_NUM" if surface_key == "NIO_IDS" else "BATIO3_ADDITIVE_NUM"
	for key in ("ADDITIVE_ATOMS_NUM", additive_num_key):
		if key not in systeminfo:
			missing.append(key)
	if missing:
		raise SystemExit(f"ERROR: Missing SYSTEMINFO entries: {', '.join(missing)}")

	surface_labels = systeminfo[surface_key].split()
	for label in (*top_labels, *bottom_labels):
		if label not in surface_labels:
			raise SystemExit(f"ERROR: Surface atom type '{label}' is missing from {surface_key}")

	atom_type_names = sorted(
		{label for key in keys for label in systeminfo[key].split()}
	)
	atom_type_ids = {label: index + 1 for index, label in enumerate(atom_type_names)}
	top_ids = sorted(atom_type_ids[label] for label in top_labels)
	bottom_ids = sorted(atom_type_ids[label] for label in bottom_labels)
	surface_ids = sorted(atom_type_ids[label] for label in surface_labels)
	working_dir = Path(working_dir).resolve()

	with Path(output_path).open("w") as handle:
		handle.write(f"LAMMPSDUMP_DIR = {working_dir}\n")
		handle.write(f"ADDITIVE_NUM = {systeminfo[additive_num_key]}\n")
		handle.write(f"ADDITIVE_ATOMS_NUM = {systeminfo['ADDITIVE_ATOMS_NUM']}\n")
		handle.write(f"SURFACE_ID = {' '.join(map(str, surface_ids))}\n")
		handle.write(f"SURFACE_TOP_ID = {' '.join(map(str, top_ids))}\n")
		handle.write(f"SURFACE_BOTTOM_ID = {' '.join(map(str, bottom_ids))}\n")
		handle.write(f"OUTPUT_SURFACES = {' '.join(output_surfaces)}\n")


def main():
	working_dir = Path.cwd()
	systeminfo_path = get_systeminfo_path(working_dir)
	systeminfo = read_parameters(systeminfo_path)
	output_path = working_dir / "isolate_attached_parameters.in"
	build_isolate_parameters(systeminfo, working_dir.name, working_dir, output_path)
	print(f"Created isolate-attached parameters: {output_path}")


if __name__ == "__main__":
	main()
