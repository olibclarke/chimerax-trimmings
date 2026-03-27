import csv
import os
from uuid import uuid4

import numpy as np


def _model_id_key(model):
    return tuple(model.id)


def _model_id_spec(model_id):
    if isinstance(model_id, int):
        model_id = (model_id,)
    return "#" + ".".join(str(part) for part in model_id)


def _model_ids_spec(model_ids):
    parts = []
    for model_id in model_ids:
        if isinstance(model_id, int):
            model_id = (model_id,)
        parts.append(".".join(str(part) for part in model_id))
    return "#" + ",".join(parts)


def _find_model_by_id(session, model_id):
    if isinstance(model_id, int):
        wanted = (model_id,)
    else:
        wanted = tuple(model_id)
    for model in session.models.list():
        if tuple(model.id) == wanted:
            return model
    return None


def _quote_path(path):
    return '"' + path.replace('"', '\\"') + '"'


def _quote_string(value):
    return '"' + str(value).replace('"', '\\"') + '"'


def _next_unused_model_id(session, reserved_ids, start=9200):
    used = {model.id[0] for model in session.models.list()}
    used.update(reserved_ids)
    model_id = start
    while model_id in used:
        model_id += 1
    reserved_ids.add(model_id)
    return model_id


def _polymer_residues_for_chain(chain):
    from chimerax.atomic import Residues, Residue

    residues = [res for res in chain.existing_residues if res.polymer_type != Residue.PT_NONE]
    return Residues(residues)


def _ss_kind(residue):
    if getattr(residue, "is_helix", False):
        return "helix"
    if getattr(residue, "is_strand", False):
        return "strand"
    return None


def _window_ranges(residues, window_size):
    windows = []
    start = 0
    count = len(residues)
    while start < count:
        end = min(start + window_size - 1, count - 1)
        ss_kind = _ss_kind(residues[end])
        if ss_kind is not None:
            while end + 1 < count and _ss_kind(residues[end + 1]) == ss_kind:
                end += 1
        windows.append((start, end))
        start = end + 1

    # Fold tiny terminal leftovers into the adjacent panel instead of exporting
    # a nearly empty N- or C-terminal segment on its own.
    min_terminal_size = max(1, window_size // 2)
    if len(windows) > 1:
        first_start, first_end = windows[0]
        if first_end - first_start + 1 < min_terminal_size:
            _, next_end = windows[1]
            windows[1] = (first_start, next_end)
            windows.pop(0)
    if len(windows) > 1:
        last_start, last_end = windows[-1]
        if last_end - last_start + 1 < min_terminal_size:
            prev_start, _ = windows[-2]
            windows[-2] = (prev_start, last_end)
            windows.pop()
    return windows


def _principal_axis_transform(coords):
    from chimerax.geometry import Place, translation

    center = coords.mean(axis=0)
    centered = coords - center
    covariance = np.dot(centered.T, centered)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    major = eigenvectors[:, order[0]]
    second = eigenvectors[:, order[1]]

    major /= np.linalg.norm(major)
    second -= major * np.dot(second, major)
    if np.linalg.norm(second) < 1.0e-6:
        second = np.array((1.0, 0.0, 0.0))
        if abs(np.dot(second, major)) > 0.9:
            second = np.array((0.0, 0.0, 1.0))
        second -= major * np.dot(second, major)
    second /= np.linalg.norm(second)
    third = np.cross(second, major)
    third /= np.linalg.norm(third)

    base_rotation = np.column_stack((second, major, third)).T
    aligned = centered @ base_rotation.T

    # After aligning the major axis vertically, choose the rotation around
    # that vertical axis that maximizes the atoms' horizontal extent.
    best_angle = 0.0
    best_span = None
    for angle in np.linspace(0.0, np.pi, 181):
        c = np.cos(angle)
        s = np.sin(angle)
        x = c * aligned[:, 0] + s * aligned[:, 2]
        span = float(x.max() - x.min())
        if best_span is None or span > best_span:
            best_span = span
            best_angle = angle

    c = np.cos(best_angle)
    s = np.sin(best_angle)
    y_rotation = np.array(
        (
            (c, 0.0, s),
            (0.0, 1.0, 0.0),
            (-s, 0.0, c),
        )
    )
    rotation_matrix = y_rotation @ base_rotation

    rotation_place = Place(
        (
            (rotation_matrix[0, 0], rotation_matrix[0, 1], rotation_matrix[0, 2], 0.0),
            (rotation_matrix[1, 0], rotation_matrix[1, 1], rotation_matrix[1, 2], 0.0),
            (rotation_matrix[2, 0], rotation_matrix[2, 1], rotation_matrix[2, 2], 0.0),
        )
    )
    return translation(center) * rotation_place * translation(-center)


def _bounds_size(bounds):
    xyz_min = np.array(bounds.xyz_min)
    xyz_max = np.array(bounds.xyz_max)
    return xyz_max - xyz_min


def _residue_id_text(residue):
    insertion = residue.insertion_code if residue.insertion_code else ""
    return f"{residue.number}{insertion}"


def _residue_key(residue):
    insertion = residue.insertion_code if residue.insertion_code else ""
    return (residue.chain_id, residue.number, insertion)


def _display_state(session):
    return {_model_id_key(model): model.display for model in session.models.list()}


def _restore_display_state(session, states):
    for model in session.models.list():
        key = _model_id_key(model)
        if key in states:
            model.display = states[key]


def _view_name():
    return f"map_salami_{uuid4().hex}"


def _atom_display_state(structure):
    return np.array(structure.atoms.displays, copy=True)


def _restore_atom_display_state(structure, displays):
    structure.atoms.displays = displays


def _format_sigfigs(value, sigfigs=3):
    if isinstance(value, (int, np.integer)):
        return str(value)
    if isinstance(value, (float, np.floating)):
        return format(float(value), f".{sigfigs}g")
    return value


def map_salami(
    session,
    model=None,
    zone=2.0,
    map=None,
    patches=None,
    segment_size=40,
    image_height=2000,
    output_dir=os.path.expanduser("~/Desktop"),
):
    '''
    Export one PNG per local model-map fit segment and a matching CSV metadata file.

    Parameters
    ----------
    model : AtomicStructure
        Atomic model to segment into residue windows.
    zone : float, optional
        Surface-zoning distance in Angstroms. Default 2.0.
    map : Volume
        Map whose displayed surface/mesh will be zoned around each residue window.
    patches : int, optional
        If given, pass `maxComponents` to `surface zone` to limit the number
        of disconnected density patches kept.
    segment_size : int, optional
        Target number of residues per segment before extending to avoid breaking a
        secondary-structure element, with tiny terminal leftovers merged into
        the neighboring segment. Default 40.
    image_height : int, optional
        Height of each saved PNG in pixels. Width is computed automatically to
        preserve the starting scene aspect ratio. If omitted, the starting
        scene aspect ratio is preserved and the default height is 2000 pixels.
    output_dir : str, optional
        Folder for output PNGs and map_salami_metadata.csv. Default ~/Desktop.
    '''
    from chimerax.atomic import concise_residue_spec
    from chimerax.core.commands import run
    from chimerax.core.errors import UserError
    from chimerax.core.objects import Objects
    from chimerax.map import Volume
    from chimerax.std_commands.view import NamedView, view as view_command

    structure = model
    map_model = map
    panel_pad = 0.01
    start_window_width, start_window_height = session.main_view.window_size
    start_aspect = start_window_width / max(start_window_height, 1)

    if segment_size < 1:
        raise UserError("segment_size must be at least 1.")
    if zone <= 0:
        raise UserError("zone must be greater than 0.")
    if patches is not None and patches < 1:
        raise UserError("patches must be at least 1.")
    if image_height is not None and image_height < 1:
        raise UserError("image_height must be at least 1.")
    if structure is None:
        raise UserError("A model must be supplied with the model keyword.")
    if map_model is None:
        raise UserError("A map model must be supplied with the map keyword.")
    if not os.path.isdir(output_dir):
        raise UserError(f'Output directory does not exist: "{output_dir}"')
    if not isinstance(map_model, Volume):
        raise UserError("The map argument must be a volume model.")

    image_height = int(image_height)
    image_width = max(1, int(round(image_height * start_aspect)))

    display_state = _display_state(session)
    atom_display_state = _atom_display_state(structure)
    metadata_rows = []
    metadata_fieldnames = [
        "image_file",
        "chain_id",
        "start_residue",
        "end_residue",
        "residue_count",
        "zone_offset_a",
        "surface_levels",
        "image_levels",
        "image_pixel_size_a",
        "panel_width_a",
        "panel_height_a",
    ]
    if patches is not None:
        metadata_fieldnames.append("patches")
    panel_index = 1
    saved_view = NamedView(session.main_view, session.main_view.center_of_rotation, session.models.list())
    starting_scene_name = "map_salami_scene_start"

    run(session, f"scenes save {_quote_string(starting_scene_name)}")

    try:
        for chain in structure.chains:
            residues = _polymer_residues_for_chain(chain)
            if len(residues) == 0:
                continue

            for start_index, end_index in _window_ranges(residues, segment_size):
                window_residues = residues[start_index : end_index + 1]
                residue_spec = concise_residue_spec(session, window_residues)
                panel_atoms = window_residues.atoms
                if len(panel_atoms) == 0:
                    continue

                if len(map_model.surfaces) == 0:
                    raise UserError(
                        "The supplied map has no displayed surfaces to zone. Please show it as a surface or mesh first."
                    )

                zone_cmd = f"surface zone {map_model.atomspec} near {residue_spec} distance {zone}"
                if patches is not None:
                    zone_cmd += f" maxComponents {int(patches)}"
                run(session, zone_cmd)

                for model in session.models.list():
                    if tuple(model.id) not in (tuple(structure.id), tuple(map_model.id)):
                        model.display = False
                structure.display = True
                map_model.display = True
                structure.atoms.displays = False
                panel_atoms.displays = True
                run(session, f"~rib {structure.atomspec}; ~surf {structure.atomspec}; show {residue_spec} atoms")
                for surface in map_model.surfaces:
                    surface.display = True
                session.update_loop.update_graphics_now()

                coords = panel_atoms.scene_coords
                if len(coords) == 0:
                    raise UserError(f"Residue window for panel {panel_index} has no displayed atoms.")

                transform = _principal_axis_transform(coords)
                original_structure_position = structure.scene_position
                structure.scene_position = transform * structure.scene_position
                original_map_position = map_model.scene_position
                map_model.scene_position = transform * map_model.scene_position

                try:
                    view_command(
                        session,
                        Objects(atoms=panel_atoms, models=[map_model]),
                        clip=True,
                        cofr=True,
                        orient=True,
                        pad=panel_pad,
                        need_undo=False,
                    )
                    session.update_loop.update_graphics_now()

                    atom_bounds = panel_atoms.scene_bounds
                    if atom_bounds is not None:
                        session.main_view.camera.view_all(
                            atom_bounds, window_size=(int(image_width), int(image_height)), pad=panel_pad
                        )
                        session.main_view.center_of_rotation = atom_bounds.center()
                        session.update_loop.update_graphics_now()

                    panel_center = panel_atoms.scene_coords.mean(axis=0)
                    image_pixel_size_a = float(session.main_view.pixel_size(panel_center))

                    atom_xyz_min = panel_atoms.scene_coords.min(axis=0)
                    atom_xyz_max = panel_atoms.scene_coords.max(axis=0)
                    size = atom_xyz_max - atom_xyz_min
                    if len(size) == 3:
                        width_a, height_a, depth_a = [float(v) for v in size]
                    else:
                        width_a = height_a = depth_a = 0.0

                    scene_name = f"map_salami_scene_{panel_index}"
                    run(session, f"scenes save {_quote_string(scene_name)}")

                    image_path = os.path.join(output_dir, f"{panel_index}.png")
                    run(
                        session,
                        f"save {_quote_path(image_path)} width {int(image_width)} height {int(image_height)}",
                    )

                    start_residue = window_residues[0]
                    end_residue = window_residues[-1]
                    row = {
                        "image_file": os.path.basename(image_path),
                        "chain_id": chain.chain_id,
                        "start_residue": _residue_id_text(start_residue),
                        "end_residue": _residue_id_text(end_residue),
                        "residue_count": len(window_residues),
                        "zone_offset_a": _format_sigfigs(zone),
                        "surface_levels": ";".join(str(surface.level) for surface in map_model.surfaces),
                        "image_levels": ";".join(
                            f"{level}:{brightness}" for level, brightness in map_model.image_levels
                        ),
                        "image_pixel_size_a": _format_sigfigs(image_pixel_size_a),
                        "panel_width_a": _format_sigfigs(width_a),
                        "panel_height_a": _format_sigfigs(height_a),
                    }
                    if patches is not None:
                        row["patches"] = int(patches)
                    metadata_rows.append(row)
                finally:
                    structure.scene_position = original_structure_position
                    map_model.scene_position = original_map_position
                    run(session, f"surface unzone {map_model.atomspec}")
                    _restore_atom_display_state(structure, atom_display_state)
                panel_index += 1

        if not metadata_rows:
            raise UserError("No polymer residue windows were found in the supplied structure.")

        csv_path = os.path.join(output_dir, "map_salami_metadata.csv")
        with open(csv_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=metadata_fieldnames)
            writer.writeheader()
            writer.writerows(metadata_rows)

        session.logger.info(
            f"Exported {len(metadata_rows)} fit panels to {output_dir} and wrote {csv_path}."
        )
    finally:
        try:
            run(session, f"scenes restore {_quote_string(starting_scene_name)}")
        except Exception:
            _restore_atom_display_state(structure, atom_display_state)
            _restore_display_state(session, display_state)
            saved_view.set_view(session.main_view, session.models.list())


def register_command(logger):
    from chimerax.atomic import AtomicStructureArg
    from chimerax.core.commands import CmdDesc, FloatArg, IntArg, SaveFolderNameArg, register
    from chimerax.map import MapArg

    desc = CmdDesc(
        keyword=[
            ("model", AtomicStructureArg),
            ("map", MapArg),
            ("zone", FloatArg),
            ("patches", IntArg),
            ("segment_size", IntArg),
            ("image_height", IntArg),
            ("output_dir", SaveFolderNameArg),
        ],
        required_arguments=("model", "map"),
        synopsis="Export segment-wise map-fit panels for one structure and one map",
    )
    register("map_salami", desc, map_salami, logger=logger)


register_command(session.logger)
