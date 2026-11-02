from pathlib import Path

from z_scripts.helper_scripts.read_parameters import read_parameters
from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path


def build_attachment_parameters(systeminfo, system_name, output_path):
	"""Map SYSTEMINFO entries into attachment_analysis's required parameters."""
	name = system_name.lower()
	if name.startswith("bao_"):
		surface_key = "BATIO3"
		path_prefix = "BAO"
	elif name.startswith("tio2_"):
		surface_key = "BATIO3"
		path_prefix = "TIO2"
	elif name.startswith("batio3_"):
			surface_key = "BATIO3"
			path_prefix = "BATIO3"
	elif name.startswith("nio_"):
		surface_key = "NIO"
		path_prefix = "NIO"
	else:
		raise SystemExit(
			f"ERROR: Cannot determine surface for directory '{system_name}'. "
			"Expected a name starting with bao_, tio2_, batio3_, or nio_."
		)

	required_keys = (
		f"{path_prefix}_LAMMPSDUMP_DIR",
		"ADDITIVE_HEAD_IDS",
		"ADDITIVE_TAIL_IDS",
		"DHTA_IDS",
		f"{surface_key}_IDS",
		"ADDITIVE_ATOMS_NUM",
		f"{surface_key}_ADDITIVE_NUM",
	)
	missing = [key for key in required_keys if key not in systeminfo]
	if missing:
		raise SystemExit(f"ERROR: Missing SYSTEMINFO entries: {', '.join(missing)}")

	# Atom type numbers follow alphabetical order of the type identifiers.
	relevant_id_keys = (
		"ADDITIVE_HEAD_IDS",
		"ADDITIVE_TAIL_IDS",
		"DHTA_IDS",
		f"{surface_key}_IDS",
	)
	atom_type_names = sorted(
		{
			name
			for key in relevant_id_keys
			for name in systeminfo[key].split()
		}
	)
	atom_type_ids = {name: index + 1 for index, name in enumerate(atom_type_names)}

	def resolve_ids(key):
		labels = systeminfo[key].split()
		unknown = [label for label in labels if label not in atom_type_ids]
		if unknown:
			raise SystemExit(
				f"ERROR: {key} references undefined atom type(s): {', '.join(unknown)}"
			)
		return [atom_type_ids[label] for label in labels]

	output_values = (
		("LAMMPSDUMP_DIR", systeminfo[f"{path_prefix}_LAMMPSDUMP_DIR"]),
		("ADDITIVE_HEAD_ID", resolve_ids("ADDITIVE_HEAD_IDS")),
		("ADDITIVE_TAIL_ID", resolve_ids("ADDITIVE_TAIL_IDS")),
		("SURFACE_ID", resolve_ids(f"{surface_key}_IDS")),
		("ADDITIVE_ATOMS_NUM", systeminfo["ADDITIVE_ATOMS_NUM"]),
		("ADDITIVE_NUM", systeminfo[f"{surface_key}_ADDITIVE_NUM"]),
	)

	with Path(output_path).open("w") as handle:
		for key, value in output_values:
			value_text = " ".join(map(str, value)) if isinstance(value, list) else value
			handle.write(f"{key} = {value_text}\n")


def main():
	working_dir = Path.cwd()
	systeminfo_path = get_systeminfo_path(working_dir)
	systeminfo = read_parameters(systeminfo_path)
	output_path = working_dir / "attachment_analysis_parameters.in"
	build_attachment_parameters(systeminfo, working_dir.name, output_path)
	print(f"Created attachment analysis parameters: {output_path}")


if __name__ == "__main__":
	main()