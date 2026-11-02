#!/usr/bin/env python3
"""Render VMD snapshots for available saturated-surface dump files."""

import ast
import subprocess
from pathlib import Path

from z_scripts.helper_scripts.get_systeminfo_path import get_systeminfo_path
from z_scripts.helper_scripts.read_parameters import read_parameters
from z_scripts.adsorption_analysis.snapshot_creator_feeder import (
    SURFACE_TYPE_KEYS,
    build_snapshot_parameters,
)


def parse_ids(value):
    try:
        return [int(item) for item in value.replace(",", " ").split()]
    except ValueError as exc:
        raise SystemExit(f"ERROR: Invalid atom type ID list: {value!r}") from exc


def read_element_data(path):
    module = ast.parse(path.read_text(), filename=str(path))
    for statement in module.body:
        if isinstance(statement, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "ELEMENTS"
            for target in statement.targets
        ):
            return ast.literal_eval(statement.value)
    raise ValueError(f"ELEMENTS dictionary not found in {path}")


def render_available_snapshots():
    systeminfo_path = get_systeminfo_path(Path.cwd())
    search_root = systeminfo_path.parent
    systeminfo = read_parameters(systeminfo_path)
    elements = read_element_data(Path(__file__).with_name("elements.data"))
    dump_paths = sorted(
        path
        for surface in SURFACE_TYPE_KEYS
        for path in search_root.rglob(f"{surface}_saturated.lammpstrj")
    )
    if not dump_paths:
        raise SystemExit(f"ERROR: No saturated surface dumps found below {search_root}")

    for infile in dump_paths:
        surface = infile.name.removesuffix("_saturated.lammpstrj")
        working_dir = infile.parent
        parameters_path = working_dir / f"{surface}_snapshot_parameters.in"
        build_snapshot_parameters(
            systeminfo,
            surface,
            working_dir,
            infile,
            infile.with_suffix(".side"),
            parameters_path,
            search_root.name,
        )
        params = read_parameters(parameters_path)
        atom_type_to_element = {
            int(key[5:]): value
            for key, value in params.items()
            if key.startswith("ELEM_")
        }
        for element in atom_type_to_element.values():
            if element not in elements:
                raise SystemExit(f"ERROR: Element '{element}' not found in elements.data")

        output_base = working_dir / f"{surface}_saturated"
        tcl_script = generate_vmd_tcl_script(
            infile=infile,
            atom_type_to_element=atom_type_to_element,
            surface_id=parse_ids(params["SURFACE_ID"]),
            additive_head_id=parse_ids(params["ADDITIVE_HEAD_ID"]),
            additive_tail_id=parse_ids(params["ADDITIVE_TAIL_ID"]),
            surface_color=int(params["SURFACE_COLOR"]),
            additive_head_color=int(params["ADDITIVE_HEAD_COLOR"]),
            additive_tail_color=int(params["ADDITIVE_TAIL_COLOR"]),
            output_base=str(output_base),
            elements=elements,
            rotate_top=params.get("ROTATE_TOP", "0") == "1",
        )
        tcl_file = working_dir / f"{surface}_render.tcl"
        tcl_file.write_text(tcl_script)
        print(f"Rendering {infile} with VMD")
        subprocess.run(["vmd", "-dispdev", "text", "-e", str(tcl_file)], check=True)
        tcl_file.unlink(missing_ok=True)
        parameters_path.unlink(missing_ok=True)


def generate_vmd_tcl_script(
    infile,
    atom_type_to_element,
    surface_id,
    additive_head_id,
    additive_tail_id,
    surface_color,
    additive_head_color,
    additive_tail_color,
    output_base,
    elements,
    rotate_top=False,
):
    radius_commands = []
    for atom_type, element in atom_type_to_element.items():
        radius = elements[element]["radius"]
        radius_commands.append(f"    set_atom_radius {atom_type} {radius}")

    surface_types = " ".join(str(atom_type) for atom_type in surface_id)
    head_types = " ".join(str(atom_type) for atom_type in additive_head_id)
    tail_types = " ".join(str(atom_type) for atom_type in additive_tail_id)
    output_image_side = output_base + "_side"
    output_image_top = output_base + "_top"
    rotate_side = "rotate z by 180\n" if rotate_top else ""
    undo_rotate_side = "rotate z by -180\n" if rotate_top else ""

    return f'''#!/usr/bin/env tclsh

mol new {{{infile}}}
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
color change rgb 7 0.800000 0.800000 0.900000

proc set_atom_radius {{type radius}} {{
    set sel [atomselect top "type $type"]
    $sel set radius $radius
    $sel delete
}}

{chr(10).join(radius_commands)}

rotate x by -90
{rotate_side}after
translate by 0 -0.7 0
render Tachyon {{{output_image_side}}} "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1920 1080 -o %s.tga
puts "Rendered image saved to {output_image_side}"
translate by 0 0.7 0

rotate x by 90
rotate z by 180
after 500
render Tachyon {{{output_image_top}}} "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1920 1080 -o %s.tga
puts "Rendered image saved to {output_image_top}"
quit
'''


if __name__ == "__main__":
    render_available_snapshots()
