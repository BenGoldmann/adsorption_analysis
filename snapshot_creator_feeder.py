from pathlib import Path

from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path
from z_scripts.helper_scripts.read_parameters import read_parameters


SURFACE_TYPE_KEYS = {
	"nio": ("NIO_IDS", "NIO_COLOR"),
	"bao": ("BAO_IDS", "BAO_COLOR"),
	"tio2": ("TIO2_IDS", "TIO2_COLOR"),
}


def get_surface_labels(systeminfo, surface_prefix):
	"""Resolve output surface labels, supporting shared BaTiO3 atom types."""
	if surface_prefix not in SURFACE_TYPE_KEYS:
		raise SystemExit(f"ERROR: Unsupported saturated-snapshot prefix: {surface_prefix}")
	ids_key, color_key = SURFACE_TYPE_KEYS[surface_prefix]

	if surface_prefix in ("bao", "tio2") and "BATIO3_IDS" in systeminfo:
		labels = systeminfo["BATIO3_IDS"].split()
		if "Ba" not in labels or "Ti" not in labels or "O" not in labels:
			raise SystemExit(
				f"ERROR: BATIO3_IDS must contain Ba, Ti, and O to render {surface_prefix}"
			)
		labels = ["Ba", "Ti", "O"]
	elif ids_key in systeminfo:
		labels = systeminfo[ids_key].split()
	else:
		raise SystemExit(f"ERROR: Missing {ids_key} in SYSTEMINFO")
	return labels, color_key


def get_type_sequence_surface_labels(systeminfo, surface_prefix):
	"""Return all slab labels used when assigning dump atom-type numbers."""
	if surface_prefix == "nio":
		key = "NIO_IDS"
	else:
		key = "BATIO3_IDS" if "BATIO3_IDS" in systeminfo else SURFACE_TYPE_KEYS[surface_prefix][0]
	if key not in systeminfo:
		raise SystemExit(f"ERROR: Missing {key} in SYSTEMINFO")
	return systeminfo[key].split()


def build_snapshot_parameters(systeminfo, surface_prefix, output_dir, dump_path, sidecar_path, output_path, additive_name=None):
	"""Write VMD settings for one already-isolated saturated surface dump."""
	output_dir = Path(output_dir).resolve()
	surface_labels, surface_color_key = get_surface_labels(systeminfo, surface_prefix)
	all_surface_labels = get_type_sequence_surface_labels(systeminfo, surface_prefix)
	solvent_key = "TERPINEOL_IDS" if output_dir.name.lower().endswith("_terpineol") else "DHTA_IDS"
	id_keys = ("ADDITIVE_HEAD_IDS", "ADDITIVE_TAIL_IDS", solvent_key)
	missing = [key for key in id_keys if key not in systeminfo]
	if missing:
		raise SystemExit(f"ERROR: Missing SYSTEMINFO entries: {', '.join(missing)}")

	labels = sorted(
		{label for key in id_keys for label in systeminfo[key].split()}
		| set(all_surface_labels)
	)
	atom_type_ids = {label: index + 1 for index, label in enumerate(labels)}
	missing_elements = [label for label in labels if label not in systeminfo]
	if missing_elements:
		raise SystemExit(
			"ERROR: Missing SYSTEMINFO element mapping(s): " + ", ".join(missing_elements)
		)

	additive_prefix = (additive_name or "").upper()
	head_color_key = "ADDITIVE_HEAD_COLOR"
	tail_color_key = "ADDITIVE_TAIL_COLOR"
	if head_color_key not in systeminfo:
		head_color_key = f"{additive_prefix}_HEAD_COLOR"
	if tail_color_key not in systeminfo:
		tail_color_key = f"{additive_prefix}_TAIL_COLOR"
	color_keys = (surface_color_key, head_color_key, tail_color_key)
	missing_colors = [key for key in color_keys if key not in systeminfo]
	if missing_colors:
		raise SystemExit(f"ERROR: Missing SYSTEMINFO color(s): {', '.join(missing_colors)}")

	try:
		colors = {key: int(systeminfo[key]) for key in color_keys}
	except ValueError as exc:
		raise SystemExit("ERROR: SYSTEMINFO VMD colors must be integers") from exc

	side = sidecar_path.read_text().strip().lower() if sidecar_path.is_file() else "bottom"
	if side not in {"top", "bottom"}:
		raise SystemExit(f"ERROR: Expected 'top' or 'bottom' in {sidecar_path}")

	values = [
		("DUMP_FILE", Path(dump_path).resolve()),
		("OUTPUT_DIR", Path(output_dir).resolve()),
		("SURFACE_ID", [atom_type_ids[label] for label in surface_labels]),
		("ADDITIVE_HEAD_ID", [atom_type_ids[label] for label in systeminfo["ADDITIVE_HEAD_IDS"].split()]),
		("ADDITIVE_TAIL_ID", [atom_type_ids[label] for label in systeminfo["ADDITIVE_TAIL_IDS"].split()]),
		("SURFACE_COLOR", colors[surface_color_key]),
		("ADDITIVE_HEAD_COLOR", colors[head_color_key]),
		("ADDITIVE_TAIL_COLOR", colors[tail_color_key]),
		("ROTATE_TOP", int(side == "top")),
	]
	with Path(output_path).open("w") as handle:
		for label in labels:
			values.append((f"ELEM_{atom_type_ids[label]}", systeminfo[label]))
		for key, value in values:
			value_text = " ".join(map(str, value)) if isinstance(value, list) else str(value)
			handle.write(f"{key} = {value_text}\n")


def main():
	working_dir = Path.cwd()
	systeminfo_path = get_systeminfo_path(working_dir)
	systeminfo = read_parameters(systeminfo_path)
	search_root = systeminfo_path.parent
	found = []
	for surface_prefix in SURFACE_TYPE_KEYS:
		for dump_path in search_root.rglob(f"{surface_prefix}_saturated.lammpstrj"):
			output_path = dump_path.parent / f"{surface_prefix}_snapshot_parameters.in"
			build_snapshot_parameters(
				systeminfo,
				surface_prefix,
				dump_path.parent,
				dump_path,
				dump_path.with_suffix(".side"),
				output_path,
				search_root.name,
			)
			found.append(output_path)
	if not found:
		raise SystemExit(f"ERROR: No *_saturated.lammpstrj files found in {working_dir}")
	for path in found:
		print(f"Created snapshot parameters: {path}")


if __name__ == "__main__":
	main()
