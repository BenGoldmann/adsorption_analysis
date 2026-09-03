#!/usr/bin/env python3
"""Render the attachment frame using atom data from elements.data."""

import ast
import os
import re
import subprocess
import sys
from pathlib import Path

from remove_solvent import read_parameters


def read_element_data(path):
    module = ast.parse(path.read_text(), filename=str(path))
    for statement in module.body:
        if isinstance(statement, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "ELEMENTS"
            for target in statement.targets
        ):
            return ast.literal_eval(statement.value)
    raise ValueError(f"ELEMENTS dictionary not found in {path}")


def parse_ids(value):
    return [int(item) for item in re.split(r"[,\s]+", value.strip()) if item]


def render_attachment_frame():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python script.py parameters.in")

    parameters_path = Path(sys.argv[1])
    if not parameters_path.is_file():
        raise SystemExit(f"ERROR: Parameters file not found: {parameters_path}")

    params = read_parameters(parameters_path)
    elements = read_element_data(Path(__file__).with_name("elements.data"))
    working_dir = Path(params["LAMMPSDUMP_DIR"]).expanduser()
    infile = working_dir / "attached_only_frame.lammpstrj"

    if not os.path.isfile(infile):
        raise SystemExit(f"ERROR: Attached-only frame does not exist: {infile}")

    required = (
        "SURFACE_ID", "ADDITIVE_HEAD_ID", "ADDITIVE_TAIL_ID",
        "SURFACE_COLOR", "ADDITIVE_HEAD_COLOR", "ADDITIVE_TAIL_COLOR",
    )
    missing = [key for key in required if key not in params]
    if missing:
        raise SystemExit(f"ERROR: Missing parameters: {', '.join(missing)}")

    try:
        surface_id = parse_ids(params["SURFACE_ID"])
        additive_head_id = parse_ids(params["ADDITIVE_HEAD_ID"])
        additive_tail_id = parse_ids(params["ADDITIVE_TAIL_ID"])
        surface_color = int(params["SURFACE_COLOR"])
        additive_head_color = int(params["ADDITIVE_HEAD_COLOR"])
        additive_tail_color = int(params["ADDITIVE_TAIL_COLOR"])
    except ValueError as error:
        raise SystemExit(f"ERROR: Invalid parameter value: {error}") from error

    atom_type_to_element = {}
    for key, value in params.items():
        if key.startswith("ELEM_"):
            try:
                atom_type = int(key[5:])
            except ValueError:
                continue
            element_symbol = value.strip()
            if element_symbol not in elements:
                raise SystemExit(f"ERROR: Element '{element_symbol}' not found in elements.data.")
            atom_type_to_element[atom_type] = element_symbol

    tcl_script = generate_vmd_tcl_script(
        infile,
        atom_type_to_element,
        surface_id,
        additive_head_id,
        additive_tail_id,
        surface_color,
        additive_head_color,
        additive_tail_color,
        str(working_dir),
        elements,
    )

    tcl_file = working_dir / "render_attached.tcl"
    with tcl_file.open("w") as handle:
        handle.write(tcl_script)

    print(f"Generated VMD TCL script: {tcl_file}")
    try:
        subprocess.run(["vmd", "-dispdev", "text", "-e", str(tcl_file)], check=True)
        print("VMD rendering completed successfully")
        for output in (working_dir / "side", working_dir / "top", tcl_file):
            if output.exists():
                output.unlink()
        print("Temporary files removed")
    except subprocess.CalledProcessError as error:
        print(f"Error running VMD: {error}")
    except FileNotFoundError:
        print("Error: VMD not found in PATH. Please ensure VMD is installed.")


def generate_vmd_tcl_script(
    infile,
    atom_type_to_element,
    surface_id,
    additive_head_id,
    additive_tail_id,
    surface_color,
    additive_head_color,
    additive_tail_color,
    working_dir,
    elements,
):
    radius_commands = []
    for atom_type, element in atom_type_to_element.items():
        radius = elements[element]["radius"]
        radius_commands.append(f"    set_atom_radius {atom_type} {radius}")

    surface_types = " ".join(str(atom_type) for atom_type in surface_id)
    head_types = " ".join(str(atom_type) for atom_type in additive_head_id)
    tail_types = " ".join(str(atom_type) for atom_type in additive_tail_id)
    output_image_side = working_dir + "/side"
    output_image_top = working_dir + "/top"

    return f'''#!/usr/bin/env tclsh

mol new {infile}
display projection orthographic
display depth cue off
color Display Background white
axes location off
scale by 1.728

mol modselect 0 0 type {surface_types}
mol modstyle 0 0 VDW 1.000000 12.000000
mol modmaterial 0 0 Diffuse
mol modcolor 0 0 ColorID {surface_color}
mol modstyle 0 0 VDW 1.000000 50.000000

mol representation VDW 1.000000 50.000000
mol selection type {head_types}
mol material Diffuse
mol addrep 0
mol modcolor 1 0 ColorID {additive_head_color}

mol representation VDW 1.000000 50.000000
mol selection type {tail_types}
mol material Diffuse
mol addrep 0
mol modcolor 2 0 ColorID {additive_tail_color}

color change rgb 6 0.800000 0.800000 0.800000

proc set_atom_radius {{type radius}} {{
    set sel [atomselect top "type $type"]
    $sel set radius $radius
    $sel delete
}}

{chr(10).join(radius_commands)}

rotate x by -90
after 500
render Tachyon {output_image_side} "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1920 1080 -o %s.tga
puts "Rendered image saved to {output_image_side}"

rotate x by 90
rotate z by 180
after 500
render Tachyon {output_image_top} "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1920 1080 -o %s.tga
puts "Rendered image saved to {output_image_top}"
quit
'''


if __name__ == "__main__":
    render_attachment_frame()
