_REORIENT_STATE_VERSION = 2
_HYDROPHOBIC_RESIDUE_NAMES = "ALA,VAL,ILE,LEU,MET,PHE,TRP,PRO"


def _sasa_atom_areas(atoms, probe_radius):
    import numpy as np
    from chimerax.core.errors import UserError
    from chimerax.surface import spheres_surface_area

    areas = spheres_surface_area(atoms.scene_coords, atoms.radii + probe_radius)
    # ChimeraX reports failed sphere calculations as -1; tiny negative values
    # from roundoff can instead be treated as zero exposed area.
    if not np.all(np.isfinite(areas)) or np.any(areas < -1.0e-8):
        raise UserError("ChimeraX could not calculate SASA for one or more atoms.")
    return np.maximum(areas, 0)


def _protein_sidechain_sasa(structure, probe_radius, context=None):
    import numpy as np
    from chimerax.atomic import Residue

    residues = structure.residues
    residues = residues.filter(residues.polymer_types == Residue.PT_PROTEIN)
    atoms = structure.atoms
    if context is not None:
        atoms = atoms | context
    atoms = atoms.filter(~np.isin(atoms.residues.names, list(Residue.water_res_names)))
    if not len(atoms):
        return residues, np.zeros(len(residues)), np.zeros(len(residues))

    # Calculate burial in the full model first; restricting this calculation
    # to hydrophobics or the requested chain would expose buried interfaces.
    areas = _sasa_atom_areas(atoms, probe_radius)
    residue_indices = residues.indices(atoms.residues)
    sidechain = atoms.is_side_onlys & (residue_indices >= 0)
    totals = np.bincount(
        residue_indices[sidechain], weights=areas[sidechain], minlength=len(residues)
    )
    references = np.zeros(len(residues))
    for index, residue in enumerate(residues):
        residue_atoms = residue.atoms
        sidechain = residue_atoms.is_side_onlys
        if sidechain.any():
            # Remove other residues while retaining this residue's backbone,
            # conformation, hydrogens, and radii for a matched reference.
            references[index] = _sasa_atom_areas(residue_atoms, probe_radius)[sidechain].sum()
    return residues, totals, references


def exposed_hydrophobics(
    session, atoms, minimum_exposure=0.25, minimum_area=0.0,
    residue_names=_HYDROPHOBIC_RESIDUE_NAMES, probe_radius=1.4,
):
    """Select hydrophobic residues with relative sidechain exposure above a cutoff.

    The atom specification identifies candidate residues, not the burial context.
    Each containing model is calculated independently using all non-water atoms,
    including hidden chains and ligands. Native sideonly atoms exclude backbone
    atoms and CA. Divide model sidechain SASA by that residue's isolated sidechain
    SASA with its own backbone retained, using the same radii and probe radius.
    This is a conformation-specific reference, not a tabulated tripeptide maximum.
    minimum_exposure is a fraction (default 0.25); minimum_area is an optional
    additional cutoff in square Angstroms. Residues without accessible reference
    sidechain area have zero relative exposure. Store absolute, reference, and
    relative areas on all protein residues for reuse; recompute after edits.
    """
    from math import isfinite
    import numpy as np
    from chimerax.atomic import Residue, Residues
    from chimerax.core.errors import UserError
    from chimerax.core.objects import Objects
    from chimerax.std_commands.select import select

    if not isfinite(minimum_exposure) or not 0 <= minimum_exposure <= 1:
        raise UserError("minimum_exposure must be a finite fraction between 0 and 1.")
    if not isfinite(minimum_area) or minimum_area < 0:
        raise UserError("minimum_area must be finite and at least 0 (square Angstroms).")
    if not isfinite(probe_radius) or probe_radius <= 0:
        raise UserError("probe_radius must be finite and greater than 0 (Angstroms).")
    names = set(residue_names.replace(",", " ").upper().split())
    if not names:
        raise UserError("residue_names must contain at least one residue name.")
    candidates = atoms.unique_residues
    candidates = candidates.filter(candidates.polymer_types == Residue.PT_PROTEIN)
    if not len(candidates):
        raise UserError("The atom specification must match at least one protein residue.")

    results = [
        _protein_sidechain_sasa(structure, probe_radius)
        for structure in candidates.atoms.unique_structures
    ]
    for name in ("sidechain_sasa", "sidechain_sasa_reference", "relative_sidechain_sasa"):
        Residue.register_attr(session, name, "Trimmings SASA", attr_type=float)
    for residues, areas, references in results:
        relative = np.divide(areas, references, out=np.zeros_like(areas), where=references > 0)
        relative = np.clip(relative, 0, 1)
        for residue, area, reference, fraction in zip(residues, areas, references, relative):
            residue.sidechain_sasa = float(area)
            residue.sidechain_sasa_reference = float(reference)
            residue.relative_sidechain_sasa = float(fraction)

    selected = Residues([
        residue for residue in candidates
        if residue.name.upper() in names
        and residue.relative_sidechain_sasa > minimum_exposure
        and residue.sidechain_sasa > minimum_area
    ])
    selected_atoms = selected.atoms
    select(session, Objects(atoms=selected_atoms, bonds=selected_atoms.intra_bonds))
    message = (
        f"Selected {len(selected)} hydrophobic residues with relative sidechain exposure "
        f"> {minimum_exposure:.1%} and SASA > {minimum_area:g} A^2 "
        f"(probe radius {probe_radius:g} A)."
    )
    session.logger.info(message)
    session.logger.status(message)
    return selected


def _patch_sphere_radius(binder_radius, binder_mass):
    from math import isfinite, pi
    from chimerax.core.errors import UserError

    if (binder_radius is None) == (binder_mass is None):
        raise UserError("Supply exactly one of binder_radius (Angstroms) or binder_mass (kDa).")
    size = binder_radius if binder_radius is not None else binder_mass
    if not isfinite(size) or size <= 0:
        raise UserError("Binder size must be finite and greater than zero.")
    if binder_radius is not None:
        return binder_radius
    # Equivalent dry sphere: 1 Da = 1.66054 A^3 g/cm^3; protein density 1.35 g/cm^3.
    return (3 * binder_mass * 1000 * 1.66054 / (4 * pi * 1.35)) ** (1 / 3)


def _sphere_has_clear_approach(center, direction, xyz, radii, binder_radius):
    import numpy as np

    offsets = xyz - center
    along = np.maximum(offsets @ direction, 0)
    # Exact minimum distance to the outward half-ray, including its starting
    # point. This tests a straight entry from infinity as well as the final pose.
    distance_squared = np.maximum(np.einsum("ij,ij->i", offsets, offsets) - along * along, 0)
    return not np.any(distance_squared < (radii + binder_radius - 1.0e-4) ** 2)


def _hydrophobic_patch_poses(
    atoms, atom_areas, residue_indices, environment, radius,
    contact_distance, link_distance, samples, minimum_residues, approach,
):
    import numpy as np
    from scipy.spatial import cKDTree

    xyz, radii = atoms.scene_coords, atoms.radii
    env_xyz, env_radii = environment.scene_coords, environment.radii
    env_tree, patch_tree = cKDTree(env_xyz), cKDTree(xyz)
    collision_range = radius + float(env_radii.max())
    contact_range = radius + float(radii.max()) + contact_distance
    neighbors = {int(index): set() for index in residue_indices}
    for a, b in patch_tree.query_pairs(link_distance):
        ra, rb = int(residue_indices[a]), int(residue_indices[b])
        if ra != rb:
            neighbors[ra].add(rb)
            neighbors[rb].add(ra)

    # Deterministic, approximately uniform directions; no random orientation
    # search. Finite sampling can miss a valid placement, so this is a screen.
    index = np.arange(samples)
    z = 1 - 2 * (index + 0.5) / samples
    angle = index * (np.pi * (3 - np.sqrt(5)))
    ring = np.sqrt(1 - z * z)
    directions = np.column_stack((ring * np.cos(angle), ring * np.sin(angle), z))
    poses = {}
    for anchor in np.flatnonzero(atoms.element_numbers > 1):
        centers = xyz[anchor] + (radius + radii[anchor]) * directions
        nearby_sets = env_tree.query_ball_point(centers, collision_range)
        for center, direction, nearby in zip(centers, directions, nearby_sets):
            nearby = np.asarray(nearby, dtype=int)
            distances = np.linalg.norm(env_xyz[nearby] - center, axis=1)
            if np.any(distances < radius + env_radii[nearby] - 1.0e-4):
                continue
            if approach and not _sphere_has_clear_approach(
                center, direction, env_xyz, env_radii, radius
            ):
                continue
            close = np.asarray(patch_tree.query_ball_point(center, contact_range), dtype=int)
            gaps = np.linalg.norm(xyz[close] - center, axis=1) - radii[close] - radius
            close = close[(gaps >= -1.0e-4) & (gaps <= contact_distance)]
            touching = set(int(i) for i in residue_indices[close])
            component, pending = set(), [int(residue_indices[anchor])]
            while pending:
                residue_index = pending.pop()
                if residue_index in component or residue_index not in touching:
                    continue
                component.add(residue_index)
                pending.extend(neighbors[residue_index] & (touching - component))
            if len(component) < minimum_residues:
                continue
            members = tuple(sorted(component))
            score = float(atom_areas[close[np.isin(residue_indices[close], members)]].sum())
            if members not in poses or score > poses[members][0]:
                poses[members] = (score, center.copy(), direction.copy())
    return poses


def hydrophobic_patches(
    session, atoms, binder_radius=None, binder_mass=None, context=None,
    minimum_exposure=0.25, minimum_residues=3, contact_distance=3.0,
    link_distance=6.0, samples=64, max_patches=10, approach=True, highlight=True,
    residue_names=_HYDROPHOBIC_RESIDUE_NAMES, probe_radius=1.4, gui=True,
):
    """Screen exposed hydrophobic clusters for access by a globular binder sphere.

    Supply a radius in Angstroms or mass in kDa (equivalent sphere at 1.35 g/cm^3).
    All patch members must lie near one shared clash-free sphere placement.
    By default its outward ray must also admit a straight collision-free entry.
    Rank by nearby water-probe sidechain SASA; this is not predicted affinity or
    designability. Optional context adds steric blockers from separate models.
    In GUI sessions, open or refresh the patch browser unless gui is false.
    """
    from math import isfinite
    import numpy as np
    from chimerax.atomic import Residue, Residues, concise_residue_spec
    from chimerax.core.errors import UserError
    from chimerax.core.objects import Objects
    from chimerax.std_commands.select import select

    radius = _patch_sphere_radius(binder_radius, binder_mass)
    if not isfinite(minimum_exposure) or not 0 <= minimum_exposure <= 1:
        raise UserError("minimum_exposure must be a fraction between 0 and 1.")
    if not isfinite(probe_radius) or probe_radius <= 0:
        raise UserError("probe_radius must be finite and positive (Angstroms).")
    if any(not isfinite(value) or value <= 0 for value in (contact_distance, link_distance)):
        raise UserError("contact_distance and link_distance must be finite and positive.")
    if minimum_residues < 1 or samples < 12 or max_patches < 1:
        raise UserError("Require minimum_residues >= 1, samples >= 12, and max_patches >= 1.")
    names = set(residue_names.replace(",", " ").upper().split())
    if not names:
        raise UserError("residue_names must contain at least one residue name.")
    candidates = atoms.unique_residues
    candidates = candidates.filter(candidates.polymer_types == Residue.PT_PROTEIN)
    if not len(candidates):
        raise UserError("The atom specification must match protein residues.")

    records = []
    analyzed_residues = []
    for structure in candidates.atoms.unique_structures:
        session.logger.status(f"Screening hydrophobic patches in #{structure.id_string}...")
        environment = structure.atoms if context is None else structure.atoms | context
        environment = environment.filter(~np.isin(
            environment.residues.names, list(Residue.water_res_names)
        ))
        residues, areas, references = _protein_sidechain_sasa(structure, probe_radius, context=context)
        analyzed_residues.extend(residues)
        fractions = np.divide(areas, references, out=np.zeros_like(areas), where=references > 0)
        fractions = np.clip(fractions, 0, 1)
        eligible = residues.filter(
            residues.mask(candidates) & np.isin(residues.names, list(names))
            & (fractions > minimum_exposure)
        )
        if len(eligible) < minimum_residues:
            continue
        env_areas = _sasa_atom_areas(environment, probe_radius)
        keep = (environment.residues.mask(eligible) & environment.is_side_onlys & (env_areas > 1.0e-8))
        patch_atoms = environment.filter(keep)
        if not len(patch_atoms):
            continue
        poses = _hydrophobic_patch_poses(
            patch_atoms, env_areas[keep], residues.indices(patch_atoms.residues),
            environment, radius, contact_distance, link_distance, samples,
            minimum_residues, approach,
        )
        for members, (score, center, direction) in poses.items():
            records.append({
                "structure": structure,
                "residues": residues[np.asarray(members, dtype=int)],
                "nearby_hydrophobic_sasa": score, "center_scene": center,
                "approach_direction": direction, "binder_radius": radius,
                "center_local": structure.scene_position.inverse() * center,
                "direction_local": structure.scene_position.inverse().transform_vector(direction),
            })

    records.sort(key=lambda record: (
        -record["nearby_hydrophobic_sasa"],
        tuple((r.structure.id, r.chain_id, r.number, r.insertion_code) for r in record["residues"]),
    ))
    # Report distinct patches rather than many nearly identical placements.
    patches, used = [], set()
    for record in records:
        members = set(record["residues"])
        if members & used:
            continue
        patches.append(record)
        used.update(members)
        if len(patches) == max_patches:
            break

    Residue.register_attr(session, "hydrophobic_patch_rank", "Trimmings patches", attr_type=int)
    for residue in analyzed_residues:
        residue.hydrophobic_patch_rank = 0
    with session.undo.aggregate("highlight hydrophobic patches"):
        for rank, patch in enumerate(patches, start=1):
            for residue in patch["residues"]:
                residue.hydrophobic_patch_rank = rank
            if highlight:
                from colorsys import hsv_to_rgb
                from chimerax.core.colors import Color
                from chimerax.std_commands.color import color

                rgba = (*hsv_to_rgb((0.08 + (rank - 1) * 0.618034) % 1, 0.8, 0.95), 1)
                color(session, Objects(atoms=patch["residues"].atoms), Color(rgba), target="ac")
            spec = concise_residue_spec(session, patch["residues"])
            session.logger.info(
                f"Patch {rank}: {spec}; {len(patch['residues'])} hydrophobics; "
                f"nearby sidechain SASA {patch['nearby_hydrophobic_sasa']:.3g} A^2."
            )
        selected = Residues([residue for patch in patches for residue in patch["residues"]]).atoms
        select(session, Objects(atoms=selected, bonds=selected.intra_bonds))
    session._trimmings_hydrophobic_patches = patches
    session.logger.info(
        f"Found {len(patches)} distinct candidate patches for a {radius:.3g}-A-radius sphere "
        f"using {samples} sampled directions per anchor. Geometric screen only; "
        "not a prediction of affinity, specificity, or binder design success."
    )
    if gui and session.ui.is_gui:
        patch_browser(session)
    return patches


def _create_patch_browser(session):
    # Follow ResNotes' list-and-navigation layout, without requiring its bundle.
    from Qt.QtCore import Qt
    from Qt.QtWidgets import QCheckBox, QHBoxLayout, QLabel, QListWidget, QPushButton, QVBoxLayout
    from chimerax.core.models import REMOVE_MODELS
    from chimerax.core.tools import ToolInstance
    from chimerax.std_commands.view import NamedView
    from chimerax.ui import MainToolWindow

    class HydrophobicPatchBrowser(ToolInstance):
        SESSION_SAVE = False
        SESSION_ENDURING = False

        def __init__(self):
            super().__init__(session, "Hydrophobic Patches")
            self._closed = False
            self._rows = []
            view = session.main_view
            # Save camera/clipping only: restoring this view must not move models.
            self._starting_view = NamedView(view, view.center_of_rotation.copy(), [])
            self._starting_cofr = view.center_of_rotation.copy()
            self._starting_cofr_method = view.center_of_rotation_method
            self.tool_window = MainToolWindow(self)
            root = QVBoxLayout(self.tool_window.ui_area)
            self.summary = QLabel("No patches")
            root.addWidget(self.summary)
            self.patch_list = QListWidget()
            self.patch_list.setAlternatingRowColors(True)
            self.patch_list.setMinimumWidth(290)
            self.patch_list.currentRowChanged.connect(self._activate_current)
            root.addWidget(self.patch_list, 1)
            self.details = QLabel()
            self.details.setWordWrap(True)
            self.details.setTextFormat(Qt.TextFormat.PlainText)
            self.details.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            root.addWidget(self.details)
            navigation = QHBoxLayout()
            self.previous_button = QPushButton("< Previous")
            self.next_button = QPushButton("Next >")
            self.previous_button.clicked.connect(lambda _checked=False: self._step(-1))
            self.next_button.clicked.connect(lambda _checked=False: self._step(1))
            navigation.addWidget(self.previous_button)
            navigation.addWidget(self.next_button)
            root.addLayout(navigation)
            self.auto_center = QCheckBox("Auto-center on patch")
            self.auto_center.setChecked(True)
            self.from_binder = QCheckBox("View from binder approach")
            self.from_binder.setChecked(True)
            for checkbox in (self.auto_center, self.from_binder):
                checkbox.toggled.connect(self._activate_current)
                root.addWidget(checkbox)
            actions = QHBoxLayout()
            self.focus_button = QPushButton("Focus patch")
            self.focus_button.clicked.connect(lambda _checked=False: self.focus_current())
            restore = QPushButton("Restore view")
            restore.clicked.connect(lambda _checked=False: self.restore_view())
            refresh = QPushButton("Refresh list")
            refresh.clicked.connect(lambda _checked=False: self.refresh())
            for button in (self.focus_button, restore, refresh):
                actions.addWidget(button)
            root.addLayout(actions)
            self.status = QLabel("Geometric screen only; rescreen after structural edits.")
            self.status.setWordWrap(True)
            root.addWidget(self.status)
            self._removed_handler = session.triggers.add_handler(REMOVE_MODELS, self._models_removed)
            self.tool_window.manage("side")
            self.refresh()

        def _live(self, patch):
            model = patch.get("structure")
            if model is None:
                try:
                    if len(patch["residues"]):
                        model = patch["residues"][0].structure
                        patch["structure"] = model
                except (ReferenceError, RuntimeError):
                    return False
            return model is not None and not model.deleted and model.id is not None

        def refresh(self):
            source = getattr(session, "_trimmings_hydrophobic_patches", [])
            self._rows = [(rank, patch) for rank, patch in enumerate(source, 1) if self._live(patch)]
            self.patch_list.blockSignals(True)
            self.patch_list.clear()
            for rank, patch in self._rows:
                self.patch_list.addItem(
                    f"Patch {rank} | {len(patch['residues'])} residues | "
                    f"{patch['nearby_hydrophobic_sasa']:.3g} A^2"
                )
            self.patch_list.setCurrentRow(0 if self._rows else -1)
            self.patch_list.blockSignals(False)
            self.summary.setText(f"{len(self._rows)} candidate patches (ranked by nearby hydrophobic SASA)")
            for button in (self.previous_button, self.next_button, self.focus_button):
                button.setEnabled(bool(self._rows))
            self._activate_current()

        def _current(self):
            row = self.patch_list.currentRow()
            return self._rows[row] if 0 <= row < len(self._rows) else None

        def _step(self, direction):
            if self._rows:
                row = (self.patch_list.currentRow() + direction) % len(self._rows)
                if row == self.patch_list.currentRow():
                    self._activate_current()
                else:
                    self.patch_list.setCurrentRow(row)
                    self.patch_list.scrollToItem(self.patch_list.currentItem())

        def _activate_current(self, *_args):
            from chimerax.atomic import concise_residue_spec
            from chimerax.core.objects import Objects
            from chimerax.std_commands.select import select

            current = self._current()
            if current is None:
                self.details.setText("No live patch results. Run hydrophobicpatches to generate candidates.")
                return
            rank, patch = current
            if not self._live(patch):
                self.refresh()
                return
            atoms = patch["residues"].atoms
            if not len(atoms):
                self.details.setText("The patch atoms were deleted. Rescreen the target.")
                return
            select(session, Objects(atoms=atoms, bonds=atoms.intra_bonds))
            self.details.setText(
                f"Patch {rank} ({self.patch_list.currentRow() + 1} of {len(self._rows)})\n"
                f"{concise_residue_spec(session, patch['residues'])}\n"
                f"Nearby sidechain SASA: {patch['nearby_hydrophobic_sasa']:.3g} A^2\n"
                f"Binder sphere radius: {patch['binder_radius']:.3g} A"
            )
            self.status.setText("Geometric screen only; rescreen after structural edits.")
            if self.auto_center.isChecked():
                self.focus_current()

        def focus_current(self):
            from chimerax.geometry import orthonormal_frame

            current = self._current()
            if current is None:
                return
            _rank, patch = current
            if not self._live(patch):
                self.refresh()
                return
            model = patch["structure"]
            if not model.visible:
                self.status.setText(f"Model #{model.id_string} is hidden; show it to inspect this patch.")
                return
            bounds = patch["residues"].atoms.scene_bounds
            if bounds is None:
                return
            view, center = session.main_view, bounds.center()
            if self.from_binder.isChecked():
                direction = patch["approach_direction"]
                if "direction_local" in patch:
                    direction = model.scene_position.transform_vector(patch["direction_local"])
                # The outward approach points toward the camera, whose view is -z.
                _right, up, _back = view.camera.position.axes()
                view.camera.position = orthonormal_frame(direction, ydir=up, origin=center)
            view.view_all(bounds, pad=0.15)
            normal = view.camera.view_direction()
            depth = max(bounds.radius(), 1.0)
            view.clip_planes.set_clip_position("near", center - depth * normal, view)
            view.clip_planes.set_clip_position("far", center + depth * normal, view)
            view.center_of_rotation = center

        def restore_view(self):
            view = session.main_view
            self._starting_view.set_view(view, [])
            view.center_of_rotation = self._starting_cofr.copy()
            view.center_of_rotation_method = self._starting_cofr_method
            self.status.setText("Restored the camera and clipping from when the browser opened.")

        def _models_removed(self, *_args):
            self.refresh()

        def delete(self):
            if self._closed:
                return
            self._closed = True
            self._removed_handler.remove()
            if getattr(session, "_trimmings_patch_browser", None) is self:
                session._trimmings_patch_browser = None
            super().delete()

    return HydrophobicPatchBrowser()


def patch_browser(session):
    """Open or refresh the click-through browser for the latest patch search."""
    from chimerax.core.errors import UserError

    if not session.ui.is_gui:
        raise UserError("The patch browser requires the ChimeraX GUI.")
    tool = getattr(session, "_trimmings_patch_browser", None)
    if tool is None or tool._closed:
        tool = _create_patch_browser(session)
        session._trimmings_patch_browser = tool
    else:
        tool.refresh()
        tool.display(True)
    return tool


def _motion_units(session):
    from chimerax.atomic import Structure
    from chimerax.core.models import Surface
    from chimerax.map import Volume

    candidates = {
        model for model in session.models.list()
        if isinstance(model, (Structure, Surface, Volume))
        and getattr(model, "has_scene_bounds", True)
    }
    # A structure/map owns its surfaces. Move the owner once, never both.
    units = set()
    for model in candidates:
        parent = model.parent
        while parent is not None and parent not in candidates:
            parent = parent.parent
        if parent is None:
            units.add(model)
    return units


def _requested_motion_units(session, models):
    from chimerax.core.errors import UserError

    requested = set(models)
    units = set()
    for model in _motion_units(session):
        ancestor = model
        while ancestor is not None:
            if ancestor in requested:
                units.add(model)
                break
            ancestor = ancestor.parent
    if not units:
        raise UserError("Specify structures, maps, independent surfaces, or their groups.")
    return units


def _rotate_motion_models(mode, axis, angle):
    from chimerax.geometry import rotation

    models = mode.models()
    if not models:
        return
    scene_axis = mode.camera_position.transform_vector(axis)
    # View.rotate(..., models) substitutes the models' bounds center. Apply
    # the transform directly so cofr/cofron's actual pivot is respected.
    transform = rotation(
        scene_axis, angle * mode.speed, center=mode.view.center_of_rotation
    )
    mode.view.move(transform, models)


def _active_motion_modes(session, controller):
    from chimerax.mouse_modes.std_modes import RotateMouseMode, TranslateMouseMode

    class ActiveModelsMixin:
        icon_file = None
        _uses_view_pivot = True

        def models(self):
            return controller.targets() if controller.enabled else []

        def mouse_down(self, event):
            controller.finish_gesture()
            super().mouse_down(event)

        def _rotate(self, axis, angle):
            _rotate_motion_models(self, axis, angle)

        def _translate(self, shift):
            if self.models():
                super()._translate(shift)

        def touchpad_two_finger_trans(self, event):
            if controller.begin_gesture():
                super().touchpad_two_finger_trans(event)

        def touchpad_two_finger_twist(self, event):
            if controller.begin_gesture():
                super().touchpad_two_finger_twist(event)

        def touchpad_three_finger_trans(self, event):
            if not controller.begin_gesture():
                return
            if "control" in event.modifiers:
                _dx, dy = event.three_finger_trans
                trackpad = session.ui.mouse_modes.trackpad
                scale = trackpad.trackpad_speed * session.main_view.window_size[0]
                # _translate expects camera-space pixels and converts to scene units.
                self._translate((0, 0, scale * dy))
            else:
                super().touchpad_three_finger_trans(event)

    class ActiveRotate(ActiveModelsMixin, RotateMouseMode):
        name = "rotate active models"

    class ActiveRotateAndSelect(ActiveRotate):
        name = "rotate active models and select"
        click_to_select = True

    class ActiveTranslate(ActiveModelsMixin, TranslateMouseMode):
        name = "translate active models"

    return ActiveRotate(session), ActiveRotateAndSelect(session), ActiveTranslate(session)


class _ModelMotion:
    def __init__(self, session):
        from Qt.QtCore import QTimer

        self.session = session
        self.enabled = False
        self.active = set()
        self.saved_bindings = {}
        self.rotate, self.rotate_select, self.translate = _active_motion_modes(session, self)
        for mode in (self.rotate, self.rotate_select, self.translate):
            session.ui.mouse_modes.add_mode(mode)
        self.gesture_positions = {}
        # Native touch events have no mouse-up callback. Group a continuous
        # gesture into one undo action after 300 ms without movement updates.
        self.gesture_timer = QTimer()
        self.gesture_timer.setSingleShot(True)
        self.gesture_timer.setInterval(300)
        self.gesture_timer.timeout.connect(self.finish_gesture)

    def targets(self):
        self.active = {model for model in self.active if not model.deleted and model.id is not None}
        return sorted(self.active, key=lambda model: model.id)

    def begin_gesture(self):
        models = self.targets() if self.enabled else []
        if not models:
            return False
        if not self.gesture_positions:
            self.gesture_positions = {model: model.position for model in models}
        self.gesture_timer.start()
        return True

    def finish_gesture(self):
        self.gesture_timer.stop()
        positions, self.gesture_positions = self.gesture_positions, {}
        changed = [
            model for model, position in positions.items()
            if not model.deleted and model.position != position
        ]
        if changed:
            from chimerax.core.undo import UndoState

            undo = UndoState("move active models")
            for model in changed:
                undo.add(model, "position", positions[model], model.position)
            self.session.undo.register(undo)

    def _bind(self, kind, key, modifiers, mode):
        manager = self.session.ui.mouse_modes
        binding = (kind, key, tuple(sorted(modifiers)))
        if binding not in self.saved_bindings:
            lookup = manager.mode if kind == "mouse" else manager.trackpad_mode
            self.saved_bindings[binding] = lookup(key, modifiers, exact=True)
        if kind == "mouse":
            manager.bind_mouse_mode(mouse_button=key, mouse_modifiers=modifiers, mode=mode)
        else:
            manager.bind_mouse_mode(trackpad_action=key, trackpad_modifiers=modifiers, mode=mode)

    def enable(self):
        if self.enabled:
            return
        from itertools import combinations

        manager = self.session.ui.mouse_modes
        replacements = {
            "rotate": self.rotate,
            "rotate and select": self.rotate_select,
            "rotate selected models": self.rotate,
            "translate": self.translate,
            "translate selected models": self.translate,
        }
        for binding in list(manager.bindings):
            replacement = replacements.get(binding.mode.name)
            if replacement is not None:
                replacement.speed = binding.mode.speed
                self._bind("mouse", binding.button, binding.modifiers, replacement)

        # Trackpad bindings are independent of mouse-button bindings, including
        # modifier-specific overrides that could otherwise move frozen models.
        for count in range(5):
            for modifiers in combinations(("alt", "command", "control", "shift"), count):
                for action in ("two finger swipe", "twist", "three finger swipe"):
                    previous = manager.trackpad_mode(action, list(modifiers), exact=True)
                    replacement = replacements.get(getattr(previous, "name", None))
                    if replacement is not None:
                        self._bind("trackpad", action, list(modifiers), replacement)
        for action, modifiers, mode in (
            ("two finger swipe", [], self.rotate),
            ("twist", [], self.rotate),
            ("three finger swipe", [], self.translate),
            ("three finger swipe", ["control"], self.translate),
        ):
            self._bind("trackpad", action, modifiers, mode)
        self.enabled = True

    def disable(self):
        self.finish_gesture()
        for (kind, key, modifiers), previous in self.saved_bindings.items():
            manager = self.session.ui.mouse_modes
            if kind == "mouse":
                manager.bind_mouse_mode(mouse_button=key, mouse_modifiers=list(modifiers), mode=previous)
            else:
                manager.bind_mouse_mode(trackpad_action=key, trackpad_modifiers=list(modifiers), mode=previous)
        self.saved_bindings.clear()
        self.enabled = False


def _model_motion_controller(session):
    from chimerax.core.errors import UserError

    if not session.ui.is_gui:
        raise UserError("Model movement controls require the ChimeraX GUI.")
    controller = getattr(session, "_trimmings_model_motion", None)
    if controller is None:
        controller = _ModelMotion(session)
        session._trimmings_model_motion = controller
    elif not getattr(controller.rotate, "_uses_view_pivot", False):
        from types import MethodType

        # Update our existing mode instances on reload without losing active
        # models or the original bindings needed by modelmotion off.
        controller.finish_gesture()
        for mode in (controller.rotate, controller.rotate_select, controller.translate):
            mode._rotate = MethodType(_rotate_motion_models, mode)
            mode._uses_view_pivot = True
    return controller


def _watch_new_motion_models(session):
    from chimerax.core.models import ADD_MODELS

    if getattr(session, "_trimmings_motion_added_handler", None) is not None:
        return

    def activate_new_models(_trigger_name, added_models):
        controller = getattr(session, "_trimmings_model_motion", None)
        if controller is None or not controller.enabled:
            return
        units = _motion_units(session)
        new_units = units.intersection(added_models)
        if new_units:
            # Finish the old movement before changing its targets. Owned map/
            # structure surfaces are excluded, so frozen owners stay frozen.
            controller.finish_gesture()
            controller.active.intersection_update(units)
            controller.active.update(new_units)

    session._trimmings_motion_added_handler = session.triggers.add_handler(
        ADD_MODELS, activate_new_models
    )


def _change_model_motion(session, models, operation):
    units = _motion_units(session) if models is None else _requested_motion_units(session, models)
    controller = _model_motion_controller(session)
    controller.finish_gesture()
    if not controller.enabled:
        controller.active = _motion_units(session)
    if operation == "replace":
        controller.active = units
    elif operation == "freeze":
        controller.active.difference_update(units)
    else:
        controller.active.update(units)
    controller.enable()
    _watch_new_motion_models(session)
    model_motion_status(session)


def model_motion_active(session, models):
    """Enable movement for only the specified models, independently of selection."""
    _change_model_motion(session, models, "replace")


def model_motion_freeze(session, models):
    _change_model_motion(session, models, "freeze")


def model_motion_unfreeze(session, models):
    _change_model_motion(session, models, "unfreeze")


def model_motion_all(session):
    _change_model_motion(session, None, "replace")


def model_motion_off(session):
    controller = getattr(session, "_trimmings_model_motion", None)
    if controller is not None and controller.enabled:
        controller.disable()
    handler = getattr(session, "_trimmings_motion_added_handler", None)
    if handler is not None:
        handler.remove()
        session._trimmings_motion_added_handler = None
    session.logger.status("Model movement controls off; previous input bindings restored.")


def model_motion_status(session):
    controller = getattr(session, "_trimmings_model_motion", None)
    if controller is None or not controller.enabled:
        session.logger.info("Model movement controls are off.")
        return
    targets = controller.targets()
    label = ", ".join(_model_status_label(model) for model in targets) or "none (all frozen)"
    session.logger.info(f"Active for movement: {label}")
    session.logger.status(f"Active for movement: {label}")


def next_model(session):
    _step_atomic_model_visibility(session, 1)


def previous_model(session):
    _step_atomic_model_visibility(session, -1)


def next_item(session):
    _step_model_or_map_visibility(session, 1)


def previous_item(session):
    _step_model_or_map_visibility(session, -1)


def lower_map_threshold(session):
    _adjust_displayed_map_thresholds(session, -0.5)


def raise_map_threshold(session):
    _adjust_displayed_map_thresholds(session, 0.5)


def reorient(session):
    axes, center, bounds = _visible_inertia_axes_and_bounds(session)
    if axes is None or center is None or bounds is None:
        session.logger.info("No visible atoms, surfaces, or contoured maps are available.")
        return

    target_direction, axis_label, up_axis = _next_reorient_direction(session, axes, center)
    _set_camera_view_direction(session, center, target_direction, up_axis, bounds)
    session.logger.status(f"Reoriented along the {axis_label} inertia axis")


def _model_id_key(model):
    return tuple(model.id)


def _model_id_spec_from_keys(model_ids):
    return "#" + ",".join(".".join(str(part) for part in model_id) for model_id in model_ids)


def _model_spec(models):
    return "#" + ",".join(model.id_string for model in models)


def _model_status_label(model):
    from os.path import basename

    name = getattr(model, "name", None) or "(unnamed)"
    return f"#{model.id_string} {basename(name)}"


def _show_target_status(session, model):
    session.logger.status(f"Showing {_model_status_label(model)}")


def _visible_inertia_axes_and_bounds(session):
    from chimerax.core.models import Surface
    from chimerax.core.objects import all_objects
    from chimerax.map import Volume
    from chimerax.map.volume import VolumeSurface
    from chimerax.std_commands.measure_inertia import map_points_and_weights, moments_of_inertia

    displayed = all_objects(session).displayed()
    bounds = displayed.bounds()

    vw = []
    atoms = displayed.atoms
    if len(atoms) > 0:
        vw.append((atoms.scene_coords, atoms.elements.masses))

    # Use visible non-volume surfaces directly so surface-only scenes can be
    # reoriented, but keep volume surfaces separate because map inertia should
    # be measured from the contoured density rather than from the triangulation.
    from chimerax.surface import vertex_areas

    surfaces = [
        model
        for model in displayed.models
        if isinstance(model, Surface) and not isinstance(model, VolumeSurface)
    ]
    for surface in surfaces:
        if getattr(surface, "vertices", None) is None or getattr(surface, "triangles", None) is None:
            continue
        if len(surface.vertices) == 0 or len(surface.triangles) == 0:
            continue
        weights = vertex_areas(surface.vertices, surface.triangles)
        points = surface.scene_position.transform_points(surface.vertices)
        vw.append((points, weights))

    maps = [model for model in displayed.models if isinstance(model, Volume)]
    for map_model in maps:
        vw.append(map_points_and_weights(map_model, scene_coordinates=True))

    axes, _moments, center = moments_of_inertia(vw)
    return axes, center, bounds


def _preferred_signed_axis(axis, reference_direction):
    from chimerax.geometry import inner_product

    return axis if inner_product(axis, reference_direction) >= 0 else -axis


def _signed_up_axis(axis, reference_up):
    from chimerax.geometry import inner_product

    return axis if inner_product(axis, reference_up) >= 0 else -axis


def _reorient_signature(axes, center):
    rounded_center = tuple(round(float(value), 3) for value in center)
    rounded_axes = tuple(
        tuple(round(abs(float(value)), 3) for value in axis)
        for axis in axes
    )
    return rounded_center, rounded_axes


def _reorient_candidates(axes, current_direction):
    # moments_of_inertia() sorts axes by increasing moment. For an ellipsoid,
    # that means longest physical axis first and shortest physical axis last.
    longest_axis, middle_axis, shortest_axis = axes
    shortest_near = _preferred_signed_axis(shortest_axis, current_direction)
    longest_near = _preferred_signed_axis(longest_axis, current_direction)
    middle_near = _preferred_signed_axis(middle_axis, current_direction)

    return (
        (shortest_near, "shortest", longest_axis),
        (-shortest_near, "shortest (opposite)", longest_axis),
        (middle_near, "middle", longest_axis),
        (-middle_near, "middle (opposite)", longest_axis),
        (longest_near, "longest", middle_axis),
        (-longest_near, "longest (opposite)", middle_axis),
    )


def _next_reorient_direction(session, axes, center):
    from chimerax.geometry import inner_product

    current_direction = session.main_view.camera.view_direction()
    signature = _reorient_signature(axes, center)
    state = getattr(session, "_reorient_state", None)
    if (
        state is None
        or state.get("signature") != signature
        or "candidates" not in state
        or state.get("version") != _REORIENT_STATE_VERSION
    ):
        candidates = _reorient_candidates(axes, current_direction)
        state = {
            "version": _REORIENT_STATE_VERSION,
            "signature": signature,
            "candidates": candidates,
            "index": None,
        }
    else:
        candidates = state["candidates"]

    dot_tolerance = 0.985
    current_index = None
    for index, (direction, _label, _up_axis) in enumerate(candidates):
        if inner_product(current_direction, direction) >= dot_tolerance:
            current_index = index
            break

    if state.get("index") == current_index and current_index is not None:
        next_index = (current_index + 1) % len(candidates)
    else:
        next_index = 0

    session._reorient_state = {
        "version": _REORIENT_STATE_VERSION,
        "signature": signature,
        "candidates": candidates,
        "index": next_index,
    }
    return candidates[next_index]


def _set_camera_view_direction(session, center, target_direction, up_axis, bounds):
    from numpy.linalg import norm
    from chimerax.geometry import orthonormal_frame

    view = session.main_view
    camera = view.camera
    current_origin = camera.position.origin()
    _x_axis, current_up, _z_axis = camera.position.axes()
    current_distance = max(float(norm(center - current_origin)), 1.0)

    # ChimeraX cameras look along -z, so the camera frame z axis must point
    # opposite the desired view direction.
    camera_origin = center - current_distance * target_direction
    up_reference = _signed_up_axis(up_axis, current_up)
    camera.position = orthonormal_frame(
        -target_direction,
        ydir=up_reference,
        origin=camera_origin,
    )

    # Recenter on the inertia ellipsoid first, then refit using the new view.
    view.center_of_rotation = center
    view.view_all(bounds)


def _atomic_models_by_ids(session, model_ids):
    from chimerax.atomic import AtomicStructure

    wanted = {tuple(model_id) for model_id in model_ids}
    return [
        model
        for model in session.models.list(type=AtomicStructure)
        if _model_id_key(model) in wanted
    ]


def _visible_map_models(session):
    from chimerax.map import Volume

    maps = list(session.models.list(type=Volume))
    maps.sort(key=lambda m: tuple(m.id))
    return [m for m in maps if m.display], maps


def _displayed_surface_maps(session):
    visible_maps, _all_maps = _visible_map_models(session)
    return [map_model for map_model in visible_maps if map_model.surface_shown]


def _model_panel_skip_models(session):
    from chimerax.model_panel.tool import model_panel

    return model_panel(session, "Model Panel").skip_models


def _is_skipped_in_model_panel(session, model):
    skip_models = _model_panel_skip_models(session)
    return skip_models.setdefault(model, False)


def _toggle_model_set(session, model_type, state_key, empty_message):
    models = list(session.models.list(type=model_type))
    if not models:
        session.logger.info(empty_message)
        return

    visible = [model for model in models if model.display]
    hidden = [model for model in models if not model.display]

    if visible:
        setattr(session, state_key, [_model_id_key(model) for model in visible])
        for model in visible:
            model.display = False
        return

    saved_ids = getattr(session, state_key, None)
    if saved_ids:
        remaining = {_model_id_key(model): model for model in hidden}
        restored = False
        for model_id in saved_ids:
            model = remaining.get(tuple(model_id))
            if model is not None:
                model.display = True
                restored = True
        if restored:
            return

    for model in hidden:
        model.display = True


def toggle_models(session):
    from chimerax.atomic import AtomicStructure

    _toggle_model_set(
        session,
        AtomicStructure,
        "_togglemodels_hidden_ids",
        "No atomic models are open.",
    )


def toggle_maps(session):
    from chimerax.map import Volume

    _toggle_model_set(
        session,
        Volume,
        "_togglemaps_hidden_ids",
        "No volume maps are open.",
    )


def cycle_lighting(session):
    from chimerax.core.commands import run

    presets = ["simple", "soft", "full", "flat", "gentle"]
    current_index = getattr(session, "_cyclelighting_index", -1)
    next_index = (current_index + 1) % len(presets)
    run(session, f"lighting {presets[next_index]}")
    session._cyclelighting_index = next_index


def _visible_atomic_models(session):
    from chimerax.atomic import AtomicStructure

    return [m for m in session.models.list(type=AtomicStructure) if m.display]


def _atomic_models_in_id_order(session):
    from chimerax.atomic import AtomicStructure

    models = list(session.models.list(type=AtomicStructure))
    models.sort(key=lambda model: tuple(model.id))
    return models


def _steppable_models_in_id_order(session):
    from chimerax.atomic import AtomicStructure
    from chimerax.map import Volume

    models = list(session.models.list(type=(AtomicStructure, Volume)))
    models.sort(key=lambda model: tuple(model.id))
    return models


def _step_atomic_model_visibility(session, direction):
    all_models = _atomic_models_in_id_order(session)
    if not all_models:
        session.logger.info("No atomic models are open.")
        return

    models = [model for model in all_models if not _is_skipped_in_model_panel(session, model)]
    if not models:
        session.logger.info("No unskipped atomic models are available.")
        return

    visible_models = [model for model in models if model.display]
    if not visible_models:
        target = models[0 if direction > 0 else -1]
        target.display = True
        _show_target_status(session, target)
        return

    current = visible_models[-1] if direction > 0 else visible_models[0]
    current_index = models.index(current)
    next_index = (current_index + direction) % len(models)

    # Match the old one-at-a-time stepping behavior without relying on
    # ChimeraX's private Model Panel implementation.
    for model in visible_models:
        model.display = False
    target = models[next_index]
    target.display = True
    _show_target_status(session, target)


def _step_model_or_map_visibility(session, direction):
    all_models = _steppable_models_in_id_order(session)
    if not all_models:
        session.logger.info("No atomic models or volume maps are open.")
        return

    models = [model for model in all_models if not _is_skipped_in_model_panel(session, model)]
    if not models:
        session.logger.info("No unskipped atomic models or volume maps are available.")
        return

    visible_models = [model for model in models if model.display]
    if not visible_models:
        target = models[0 if direction > 0 else -1]
        target.display = True
        _show_target_status(session, target)
        return

    current = visible_models[-1] if direction > 0 else visible_models[0]
    current_index = models.index(current)
    next_index = (current_index + direction) % len(models)

    for model in visible_models:
        model.display = False
    target = models[next_index]
    target.display = True
    _show_target_status(session, target)


def _sign_preserving_level(level, magnitude_delta, min_magnitude):
    sign = 1 if level >= 0 else -1
    magnitude = max(abs(level) + magnitude_delta, min_magnitude)
    return sign * magnitude


def _adjust_displayed_map_thresholds(session, sigma_delta):
    maps = _displayed_surface_maps(session)
    if not maps:
        session.logger.info("No displayed surface maps are available.")
        return

    adjusted = []
    for map_model in maps:
        _mean, sd, _rms = map_model.mean_sd_rms()
        if sd == 0:
            continue

        magnitude_delta = sigma_delta * sd
        min_magnitude = max(abs(sd) * 1e-6, 1e-12)
        new_levels = tuple(
            _sign_preserving_level(surface.level, magnitude_delta, min_magnitude)
            for surface in map_model.surfaces
        )
        map_model.set_parameters(
            surface_levels=new_levels,
            threaded_surface_calculation=True,
        )
        adjusted.append(f"#{map_model.id_string}")

    if adjusted:
        direction = "Raised" if sigma_delta > 0 else "Lowered"
        session.logger.status(
            f"{direction} thresholds by {abs(sigma_delta):.1f} sigma for {' '.join(adjusted)}"
        )
    else:
        session.logger.info("Displayed surface maps have zero standard deviation; thresholds unchanged.")


def _next_unused_model_id(session, reserved_ids):
    used_ids = {model.id[0] for model in session.models.list()}
    used_ids.update(reserved_ids)
    model_id = 9000
    while model_id in used_ids:
        model_id += 1
    reserved_ids.add(model_id)
    return model_id


def _close_cycle_model_display_clipboards(session):
    from chimerax.core.commands import run

    state = getattr(session, "_cyclemodeldisplay_state", None)
    if state is None:
        return

    clipboard_ids = state.get("clipboard_ids", ())
    if clipboard_ids:
        run(session, f"close {_model_id_spec_from_keys(clipboard_ids)}")


def _create_cycle_model_display_state(session, models):
    from chimerax.core.commands import run

    reserved_ids = set()
    clipboard_ids = []
    for model in models:
        clipboard_model_id = _next_unused_model_id(session, reserved_ids)
        run(
            session,
            f"combine #{model.id_string} close false modelId {clipboard_model_id} name cyclemodeldisplay_clipboard",
        )
        run(session, f"hide #{clipboard_model_id} models")
        clipboard_ids.append((clipboard_model_id,))

    return {
        "model_ids": tuple(_model_id_key(model) for model in models),
        "clipboard_ids": tuple(clipboard_ids),
        "index": 0,
    }


def _apply_model_display_preset(session, models, preset_index):
    from chimerax.core.commands import run

    spec = _model_spec(models)
    commands = [
        f"graphics silhouettes false; ca_and_sidechains {spec}",
        f"graphics silhouettes false; ca_trace {spec}",
        (
            f"graphics silhouettes false; nucleotides {spec} atoms; ~rib {spec}; ~surf {spec}; "
            f"disp {spec}; ~disp @H*&{spec}; style {spec} stick; style ions sphere; style solvent ball; "
            f"size {spec} ballscale 0.2; size {spec} stickradius 0.07"
        ),
        f"graphics silhouettes false; rib {spec}; ~disp {spec}; ~surf {spec}",
        (
            f"graphics silhouettes false; rib {spec}; ~surf {spec}; ~disp {spec}; disp @CA&protein&{spec}; "
            f"disp @P&nucleic&{spec}; style {spec} stick; disp sidechain&{spec}; "
            f"disp ~backbone&nucleic&{spec}; size {spec} stickradius 0.1; "
            f"size {spec} pseudobondradius 0.1"
        ),
        f"graphics silhouettes false; rib {spec}; ~disp {spec}; ~surf {spec}; rainbow {spec} palette RdYlBu-5",
    ]
    run(session, commands[preset_index])


def cycle_model_display(session):
    from chimerax.core.commands import run

    state = getattr(session, "_cyclemodeldisplay_state", None)
    if state is not None:
        models = _atomic_models_by_ids(session, state["model_ids"])
    else:
        models = _visible_atomic_models(session)
    if not models:
        session.logger.info("No atomic models are currently shown.")
        return

    model_ids = tuple(_model_id_key(model) for model in models)
    if state is None or state["model_ids"] != model_ids:
        _close_cycle_model_display_clipboards(session)
        state = _create_cycle_model_display_state(session, models)

    next_index = (state["index"] + 1) % 7
    if next_index == 0:
        clipboard_spec = _model_id_spec_from_keys(state["clipboard_ids"])
        run(session, f"show {clipboard_spec} models")
        for source_id, target_id in zip(state["clipboard_ids"], state["model_ids"]):
            run(
                session,
                f"mcopy {_model_id_spec_from_keys((source_id,))} to {_model_id_spec_from_keys((target_id,))} settings csv",
            )
        run(session, f"hide {clipboard_spec} models")
    else:
        _apply_model_display_preset(session, models, next_index - 1)

    state["index"] = next_index
    session._cyclemodeldisplay_state = state


def _step_map_visibility(session, direction):
    visible_maps, all_maps = _visible_map_models(session)
    if not all_maps:
        session.logger.info("No volume maps are open.")
        return

    maps = [map_model for map_model in all_maps if not _is_skipped_in_model_panel(session, map_model)]
    if not maps:
        session.logger.info("No unskipped volume maps are available.")
        return

    visible_maps = [map_model for map_model in visible_maps if map_model in maps]

    if not visible_maps:
        target = maps[0 if direction > 0 else -1]
        target.display = True
        _show_target_status(session, target)
        return

    current = visible_maps[-1] if direction > 0 else visible_maps[0]
    current_index = maps.index(current)
    next_index = (current_index + direction) % len(maps)

    for map_model in visible_maps:
        map_model.display = False
    target = maps[next_index]
    target.display = True
    _show_target_status(session, target)


def next_map(session):
    _step_map_visibility(session, 1)


def previous_map(session):
    _step_map_visibility(session, -1)


# Shape alignment: full implementation lives here, with no sibling-script dependency.
from contextvars import ContextVar
from math import pi, sqrt

_SHAPE_BATCH_CACHE = ContextVar("shapealign_batch_cache", default=None)
_SHAPE_MAP_WORKERS = ContextVar("shapealign_map_workers", default=None)
_SHAPE_TILE_PROGRESS = ContextVar("shapealign_tile_progress", default=None)
_SHAPE_UI = ContextVar("shapealign_ui", default=None)
_SHAPE_INPUT_GUARD = ContextVar("shapealign_input_guard", default=None)
_SHAPE_CACHE_BYTES = 256 * 1024 ** 2  # Retained preprocessing only, not fitter scratch space.
_SHAPE_PARALLEL_MAP_WORK = 500_000  # Amortize native-worker startup on substantial refinements.
_SHAPE_AUTO_PARALLEL_MIN = 4
_MAX_GRID_VOXELS = 16_000_000  # 64 MB per float32 grid, before fitter scratch space.
_MAX_CORRELATION_VOXELS = 8_000_000  # Several moment arrays are needed at once.
_MAX_PARTIAL_FFT_VOXELS = 4_000_000
_PARTIAL_SIZE_RATIO = 0.75  # Atom-count proxy for unequal density mass/volume.
_CHAIN_BACKUP_CORRELATION = 0.4  # Weak partial fits may have unmatched parts on both sides.
_CHAIN_BACKUP_MIN_FRACTION = 0.2  # Ignore tiny chains and cap the extra searches below.
_MAX_CHAIN_BACKUP_STARTS = 4
_AUTO_FINE_LIMITS = (2.0, 8.0)  # Bound unnecessary atom-level detail and preserve helix-scale features.
_SHAPE_MAP_BEAM = 48  # Preserve competing sites/orientations between scales.
_SHAPE_SUBMAP_FITS = 4  # Finish several fits per submap before comparing submaps.


class _ShapeBatchCache:
    """Command-scoped LRU; retain reference work without keeping stale session data."""

    def __init__(self, limit=_SHAPE_CACHE_BYTES):
        from collections import OrderedDict

        self.limit = limit
        self.entries = OrderedDict()
        self.bytes = 0

    def __enter__(self):
        self.token = _SHAPE_BATCH_CACHE.set(self)
        return self

    def __exit__(self, *exc):
        _SHAPE_BATCH_CACHE.reset(self.token)
        self.entries.clear()
        self.bytes = 0

    def get(self, key):
        entry = self.entries.get(key)
        if entry is None:
            return None
        self.entries.move_to_end(key)
        return entry[0]

    def put(self, key, value):
        size = self._size(value)
        old = self.entries.pop(key, None)
        if old is not None:
            self.bytes -= old[1]
        if size <= self.limit:
            while self.entries and self.bytes + size > self.limit:
                _, (_, removed) = self.entries.popitem(last=False)
                self.bytes -= removed
            self.entries[key] = (value, size)
            self.bytes += size
        return value

    @classmethod
    def _size(cls, value):
        import numpy as np

        if isinstance(value, np.ndarray):
            return value.nbytes
        if isinstance(value, (tuple, list)):
            return sum(cls._size(item) for item in value)
        if isinstance(value, _ShapeDensityInput):
            return ((value.matrix.nbytes if value.is_map else value.xyz.nbytes) + value.landmarks.nbytes
                    + value.frame_landmarks.nbytes + (0 if value.weights is None else value.weights.nbytes))
        if callable(getattr(value, "matrix", None)):
            return value.matrix().nbytes
        return 0


class _ShapeUI:
    """Service the native GUI between numerical steps, never from a worker."""

    def __init__(self, session):
        self.session, self.next_update = session, 0.

    def __enter__(self):
        from chimerax.core.errors import UserError

        if getattr(self.session, '_shapealign_running', False):
            raise UserError('shapealign is already running; wait for it to finish.')
        self.session._shapealign_running = True
        self.token = _SHAPE_UI.set(self)
        return self

    def __exit__(self, *args):
        _SHAPE_UI.reset(self.token)
        self.session._shapealign_running = False

    def checkpoint(self):
        from time import perf_counter

        guard = _SHAPE_INPUT_GUARD.get()
        if not self.session.ui.is_gui or guard is None or perf_counter() < self.next_update:
            return
        self.next_update = perf_counter() + .1
        # A batch owns an open native Undo aggregate. Repaint it, but do not
        # accept unrelated commands into its pending transaction.
        self._process_events(allow_input=not self.session.undo._register_stack)
        guard.validate()
        self.next_update = perf_counter() + .1

    @staticmethod
    def _process_events(allow_input):
        from Qt.QtCore import QCoreApplication, QEventLoop

        flags = (QEventLoop.ProcessEventsFlag.AllEvents if allow_input else
                 QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents)
        QCoreApplication.processEvents(flags, 5)


def _shape_checkpoint():
    ui = _SHAPE_UI.get()
    if ui is not None:
        ui.checkpoint()


def _shape_map_status(logger, message):
    progress = _SHAPE_TILE_PROGRESS.get()
    if progress is not None:
        message = f"{progress[0]}/{progress[1]} tiles complete | {message}"
    logger.status(message, blank_after=0)
    _shape_checkpoint()


def _shape_map_fft_threads(voxels):
    workers = _SHAPE_MAP_WORKERS.get()
    if workers is None or workers.requested == 1 or voxels < 262_144:
        return 1
    cpus = _shape_available_cpus()
    return max(1, cpus - 1) if workers.requested == 0 else min(cpus, workers.requested)


def _shape_refine_grid_pose(points, weights, field, xyz_to_ijk, initial, denominator):
    """One unchanged native optimization, shared by serial and worker paths."""
    import numpy as np
    from chimerax.map_fit.fitmap import locate_maximum, volume_values

    delta, _ = locate_maximum(initial.transform_points(points), weights, field, xyz_to_ijk,
                              max_steps=240, ijk_step_size_min=.01, ijk_step_size_max=.5, metric='sum product')
    pose = delta * initial
    values = volume_values(points, xyz_to_ijk * pose, field)
    return float(np.dot(weights.astype(np.float64), values)) / denominator, pose


def _shape_map_worker_main(directory):
    """Native numerical worker: read-only memory maps in, score/transform out."""
    import json
    from pathlib import Path
    import sys
    import traceback
    import numpy as np
    from chimerax.geometry import Place

    def compute(job):
        root = Path(directory) / job['stage']
        # These mappings die before the result is sent, so the parent can
        # remove a completed stage's files on Windows as well as Unix.
        points = np.load(root / 'points.npy', mmap_mode='r', allow_pickle=False)
        weights = np.load(root / 'weights.npy', mmap_mode='r', allow_pickle=False)
        field = np.load(root / 'field.npy', mmap_mode='r', allow_pickle=False)
        score, pose = _shape_refine_grid_pose(points, weights, field,
                                             Place(matrix=np.array(job['xyz_to_ijk'])),
                                             Place(matrix=np.array(job['matrix'])), job['denominator'])
        return dict(index=job['index'], stage=job['stage'], score=score, matrix=pose.matrix.tolist())

    for line in sys.stdin:
        job = json.loads(line)
        try:
            result = compute(job)
        except Exception:
            result = dict(index=job['index'], stage=job['stage'], error=traceback.format_exc())
        print('SHAPE_MAP:' + json.dumps(result), flush=True)


class _ShapeMapWorkers:
    """Lazy command-scoped native processes; session objects never leave the GUI."""

    def __init__(self, session, requested):
        self.session, self.requested = session, requested
        self.processes, self.threads, self.diagnostics = [], [], []
        self.temporary = None
        self.disabled, self.used = requested == 1, False

    def __enter__(self):
        self.token = _SHAPE_MAP_WORKERS.set(self)
        return self

    def __exit__(self, *args):
        _SHAPE_MAP_WORKERS.reset(self.token)
        self.close()

    def close(self):
        import subprocess

        if hasattr(self, 'stop'):
            self.stop.set()
        for process in self.processes:
            if process.poll() is None:
                process.terminate()
        for process in self.processes:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        for thread in self.threads:
            thread.join(timeout=5)
        for process in self.processes:
            process.stdin.close()
            process.stdout.close()
        for diagnostic in self.diagnostics:
            diagnostic.close()
        self.processes.clear()
        self.threads.clear()
        self.diagnostics.clear()
        if self.temporary is not None:
            self.temporary.cleanup()
            self.temporary = None

    def start(self, count):
        import json
        import os
        from pathlib import Path
        import queue
        from shutil import copyfile
        import subprocess
        import sys
        import tempfile
        import threading
        from chimerax.core.errors import UserError

        try:
            from chimerax.core.python_utils import chimerax_python_executable
            executable = chimerax_python_executable()
        except ImportError:
            executable = ''
        if not Path(executable).is_file():
            executable = sys.executable if Path(sys.executable).name.lower().startswith('python') else ''
        script = Path(__file__).resolve()
        if not Path(executable).is_file() or not script.is_file():
            raise UserError('Cannot locate native Python/shapealign source for workers; use workers 1.')
        environment = os.environ.copy()
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                     'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            environment[name] = '1'
        environment.pop('PYTHONHOME', None)
        environment['PYTHONDONTWRITEBYTECODE'] = '1'
        self.stop = threading.Event()
        self.waiting, self.completed = queue.Queue(), queue.Queue()
        self.temporary = tempfile.TemporaryDirectory(prefix='shapealign_maps_')
        self.directory = Path(self.temporary.name)
        copyfile(script, self.directory / 'trimmings.py')
        (self.directory / 'path.json').write_text(json.dumps([os.path.abspath(p) for p in sys.path]))
        bootstrap = ('import json, pathlib, runpy, sys; root=pathlib.Path(sys.argv[1]); '
                     'sys.path[:]=json.loads((root/"path.json").read_text()); '
                     'runpy.run_path(str(root/"trimmings.py"))["_shape_map_worker_main"](root)')

        def exchange(process, diagnostic):
            try:
                while not self.stop.is_set():
                    try:
                        job = self.waiting.get(timeout=.1)
                    except queue.Empty:
                        continue
                    process.stdin.write(json.dumps(job) + '\n')
                    process.stdin.flush()
                    while True:
                        line = process.stdout.readline()
                        if not line:
                            raise RuntimeError('A native map-refinement worker exited unexpectedly.')
                        position = line.find('SHAPE_MAP:')
                        if position >= 0:
                            result = json.loads(line[position + len('SHAPE_MAP:'):])
                            if result['index'] != job['index'] or result['stage'] != job['stage']:
                                raise RuntimeError('A map-refinement worker returned an unexpected job.')
                            self.completed.put(result)
                            break
            except Exception as error:
                details = Path(diagnostic).read_text(errors='replace')[-4000:]
                self.completed.put(dict(error=str(error) + ('\n' + details if details else '')))

        for number in range(count):
            diagnostic = (self.directory / f'worker_{number}.log').open('w+')
            self.diagnostics.append(diagnostic)
            process = subprocess.Popen([executable, '-u', '-c', bootstrap, str(self.directory)],
                                       stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=diagnostic,
                                       text=True, encoding='utf-8', env=environment,
                                       **({'creationflags': subprocess.CREATE_NO_WINDOW} if sys.platform == 'win32' else {}))
            self.processes.append(process)
            thread = threading.Thread(target=exchange, args=(process, diagnostic.name), daemon=True)
            self.threads.append(thread)
            thread.start()

    def refine(self, points, weights, field, xyz_to_ijk, candidates, denominator):
        import queue
        import tempfile
        import numpy as np
        from pathlib import Path
        from time import perf_counter
        from chimerax.core.errors import UserError
        from chimerax.geometry import Place

        if self.disabled or len(candidates) < 2:
            return None
        if self.requested == 0 and len(points) * len(candidates) < _SHAPE_PARALLEL_MAP_WORK:
            return None
        count = _shape_worker_count(self.requested, len(candidates))
        if count < 2:
            return None
        if not self.processes:
            try:
                self.start(count)
            except (OSError, UserError) as error:
                self.close()
                if self.requested != 0:
                    raise
                self.disabled = True
                self.session.logger.warning(f'Parallel map refinement unavailable; continuing serially. {error}')
                return None
        self.used = True
        with tempfile.TemporaryDirectory(prefix='stage_', dir=self.directory) as temporary:
            directory = Path(temporary)
            for name, value in (('points', points), ('weights', weights), ('field', field)):
                np.save(directory / f'{name}.npy', value, allow_pickle=False)
            results = [None] * len(candidates)
            for index, initial in enumerate(candidates):
                self.waiting.put(dict(index=index, stage=directory.name, matrix=initial.matrix.tolist(),
                                      xyz_to_ijk=xyz_to_ijk.matrix.tolist(), denominator=denominator))
            completed, next_report = 0, 0.
            try:
                while completed < len(candidates):
                    if perf_counter() >= next_report:
                        _shape_map_status(self.session.logger, f'Map refinement: {completed}/{len(candidates)} starts '
                                           f'({len(self.processes)} workers)...')
                        next_report = perf_counter() + .1
                    try:
                        result = self.completed.get(timeout=.1)
                    except queue.Empty:
                        continue
                    if 'error' in result:
                        raise UserError('Native map refinement failed; no fit was applied.\n' + result['error'])
                    pose = Place(matrix=np.array(result['matrix'], dtype=float))
                    if not np.isfinite(pose.matrix).all() or not np.isfinite(result['score']):
                        raise UserError('A native map worker returned a nonfinite result; no fit was applied.')
                    results[result['index']] = result['score'], pose, result['index']
                    completed += 1
            except BaseException:
                # Stop mappings before stage cleanup, including cancellation.
                self.close()
                raise
        # Preserve original index order before stable score sorting, including ties.
        return sorted(results, key=lambda row: row[0], reverse=True)


class _ShapeMapInputGuard:
    """Reject stale fits if repaint-time callbacks change a numerical input."""

    def __init__(self, items):
        self.inputs, self.maps, self.changed = [], [], False
        for model, atoms in items:
            state = dict(model=model, scene=model.scene_position.matrix.copy(), atoms=atoms)
            if atoms is None:
                state['data'] = model.data
                self.maps.append(model.data)
            else:
                state.update(pointers=atoms.pointers.copy(), xyz=atoms.coords.copy(), bfactors=atoms.bfactors.copy())
            self.inputs.append(state)

    def _changed(self, reason):
        if reason in ('values changed', 'coordinates changed'):
            self.changed = True

    def __enter__(self):
        for data in self.maps:
            data.add_change_callback(self._changed)
        self.token = _SHAPE_INPUT_GUARD.set(self)
        return self

    def __exit__(self, *args):
        _SHAPE_INPUT_GUARD.reset(self.token)
        for data in self.maps:
            if self._changed in data.change_callbacks:
                data.remove_change_callback(self._changed)

    def validate(self):
        import numpy as np
        from chimerax.core.errors import UserError

        invalid = self.changed
        for state in self.inputs:
            model, atoms = state['model'], state['atoms']
            if model.deleted or not np.array_equal(model.scene_position.matrix, state['scene']):
                invalid = True
                break
            if atoms is None:
                invalid |= model.data is not state['data']
            else:
                if not np.isin(state['pointers'], model.atoms.pointers).all():
                    invalid = True
                    break
                invalid |= (not np.array_equal(atoms.pointers, state['pointers']) or
                            not np.array_equal(atoms.coords, state['xyz']) or
                            not np.array_equal(atoms.bfactors, state['bfactors']))
        if invalid:
            raise UserError('A fitting model or map changed during shapealign; no fit was applied.')


def _shape_array_key(array):
    """Exact snapshot identity, not a rounded geometric similarity test."""
    from hashlib import blake2b
    import numpy as np

    array = np.ascontiguousarray(array)
    return array.shape, array.dtype.str, blake2b(memoryview(array).cast("B"), digest_size=20).digest()


def _shape_cached(key, prepare):
    cache = _SHAPE_BATCH_CACHE.get()
    if cache is None:
        return prepare()
    value = cache.get(key)
    return cache.put(key, prepare()) if value is None else value


def _shape_reference_atom_grid(xyz, resolution, weights=None):
    def prepare():
        points, mass = _shape_points(xyz, resolution / 6, weights)
        return points, mass, _shape_grid(points, mass, resolution)

    if _SHAPE_BATCH_CACHE.get() is None:
        return prepare()
    return _shape_cached(("atom_grid", _shape_array_key(xyz), resolution,
                          None if weights is None else _shape_array_key(weights)), prepare)


def _shape_has_plddt(model):
    """Require explicit AlphaFold provenance, never a filename or numeric range."""
    if getattr(getattr(model, "database", None), "database", None) == "AlphaFold":
        return True
    metadata = getattr(model, "metadata", {})
    if not metadata.get("database_2") or not metadata.get("ma_qa_metric"):
        return False
    # AlphaFoldDB CIFs opened as local files lack the fetched-model marker.
    # Check both database provenance and a declared local pLDDT metric.
    from chimerax.mmcif import get_mmcif_tables_from_metadata

    database, metrics = get_mmcif_tables_from_metadata(model, ("database_2", "ma_qa_metric"), metadata=metadata)
    return (database is not None and metrics is not None
            and any(value.casefold() == "alphafolddb" for (value,) in
                    database.fields(("database_id",), allow_missing_fields=True))
            and any(mode.casefold() == "local" and kind.casefold() == "plddt" for mode, kind in
                    metrics.fields(("mode", "type"), allow_missing_fields=True)))


def _shape_confidence_weights(atoms, mode):
    """Only interpret B-factors as pLDDT with provenance or explicit consent."""
    import numpy as np
    from chimerax.core.errors import UserError

    if mode not in ("auto", "off", "plddt"):
        raise UserError("confidence must be auto, off, or plddt.")
    if mode == "off" or (mode == "auto" and (not len(atoms) or not _shape_has_plddt(atoms[0].structure))):
        return None
    values = atoms.bfactors.astype(np.float64)
    if not np.isfinite(values).all() or np.any((values < 0) | (values > 100)):
        raise UserError("pLDDT weighting requires finite B-factor confidence scores from 0 to 100; "
                        "use confidence off for experimental B-factors.")
    # Smoothstep: ignore <=50, half weight at 70, full weight from 90.
    t = np.clip((values - 50) / 40, 0, 1)
    return t * t * (3 - 2 * t)


def _shape_center(xyz, weights=None):
    import numpy as np

    return xyz.mean(axis=0) if weights is None else np.average(xyz, weights=weights, axis=0)


def _shape_effective_count(xyz, weights=None):
    """Uniformly scaling confidence must not turn a whole fit into a partial fit."""
    import numpy as np

    if weights is None:
        return len(xyz)
    scaled = weights / weights.max()
    return float(scaled.sum() ** 2 / np.dot(scaled, scaled))


def _shape_weighted_atoms(atoms, xyz, mode):
    """Drop zero-confidence atoms before allocating grids, not from the model."""
    import numpy as np
    from chimerax.core.errors import UserError

    weights = _shape_confidence_weights(atoms, mode)
    if weights is None:
        return atoms, xyz, None
    keep = weights > 0
    atoms, xyz, weights = atoms.filter(keep), xyz[keep], weights[keep]
    if len(xyz) < 4 or np.linalg.matrix_rank(xyz - _shape_center(xyz, weights)) < 2:
        raise UserError("Confidence weighting leaves fewer than four non-collinear fitting atoms; "
                        "expand the selection or use confidence off.")
    return atoms, xyz, weights


def _shape_log_confidence(session, model, weights, role):
    if weights is not None:
        session.logger.info(f"Shape alignment {role} #{model.id_string}: pLDDT weighting, "
                            f"{len(weights)} atoms above 50; effective atom weight {weights.sum():.3g}.")


def _shape_geometry(xyz, weights=None):
    """Six points with the same weighted centroid/covariance, for map-fit seeds."""
    import numpy as np

    if weights is None:
        return xyz
    center = _shape_center(xyz, weights)
    centered = xyz - center
    values, axes = np.linalg.eigh((centered.T * weights) @ centered / weights.sum())
    offsets = (axes * np.sqrt(3 * np.maximum(values, 0))).T
    return np.concatenate((offsets, -offsets)) + center


def _shape_reference_field(grid, objective):
    def prepare():
        import numpy as np

        if objective == "gradient":
            field, norm = _gradient_density_field(grid)
        else:
            field = grid.matrix()
            norm = float(np.einsum("ijk,ijk->", field, field, dtype=np.float64))
        # Retaining the grid prevents its identity being reused while a field
        # or FFT is cached. Counting shared arrays twice is conservative.
        return grid, field, norm

    return _shape_cached(("field", id(grid), objective), prepare)[1:]


def _shape_reference_fft(grid, field, fft_shape, objective, squared=False, workers=None):
    def prepare():
        from scipy.fft import rfftn

        return grid, rfftn(field * field if squared else field, fft_shape, workers=workers)

    return _shape_cached(("fft", id(grid), objective, tuple(fft_shape), squared), prepare)[1]


def _shape_spatial_size(xyz, weights=None):
    """Equivalent diameter of a uniform PCA ellipsoid, independent of atom count."""
    import numpy as np

    centered = xyz - _shape_center(xyz, weights)
    covariance = (centered.T @ centered / len(centered) if weights is None else
                  (centered.T * weights) @ centered / weights.sum())
    variances = np.linalg.eigvalsh(covariance)
    # A uniform solid ellipsoid has variance a^2/5 along semiaxis a. Its
    # equivalent-volume diameter uses all axes, not just the longest extent.
    diameters = 2 * np.sqrt(5 * np.maximum(variances, 0))
    # A planar atom cloud still has finite physical thickness at atomic scales.
    return float(np.exp(np.log(np.maximum(diameters, 2)).mean()))


def _shape_resolution_schedule(moving, reference, resolution=None, local=False,
                               moving_weights=None, reference_weights=None):
    """Derive both endpoints and geometric intermediate levels from fitted size."""
    import numpy as np

    moving_size = _shape_spatial_size(moving, moving_weights)
    reference_size = _shape_spatial_size(reference, reference_weights)
    size = min(moving_size, reference_size)
    # About four broad shape features across the object at the initial scale,
    # and twelve at the fine scale. These are size ratios, not molecule presets.
    fine = float(np.clip(size / 12, *_AUTO_FINE_LIMITS)) if resolution is None else float(resolution)
    coarse = max(fine, size / 4)
    if local:
        coarse = min(coarse, 2 * fine)  # Already aligned: avoid unnecessary broad blurring.
    steps = max(0, int(np.ceil(np.log2(coarse / fine))))
    levels = tuple(float(value) for value in np.geomspace(coarse, fine, steps + 1))
    return levels, {
        "resolution_mode": "automatic" if resolution is None else "explicit",
        "object_size": size, "moving_size_diameter": moving_size,
        "reference_size_diameter": reference_size,
    }


def _local_correlation_field(first, second, sigma):
    """Gaussian-weighted Pearson correlation and validity, both float32 grids.

    Consume the two density arrays as scratch space to limit peak memory.
    Undefined voxels contain zero in the field and zero in the validity grid.
    """
    import numpy as np
    from scipy.ndimage import gaussian_filter

    def smooth(array):
        gaussian_filter(array, sigma, output=array, mode="constant", truncate=3)
        return array

    # Independent positive scaling preserves Pearson correlation and makes the
    # absolute variance floor independent of Gaussian density normalization.
    for density in (first, second):
        peak = float(density.max())
        if peak > 0:
            density /= peak
    mean_first = smooth(first.copy())
    mean_second = smooth(second.copy())
    covariance = smooth(first * second)
    covariance -= mean_first * mean_second

    np.square(first, out=first)
    smooth(first)
    np.square(mean_first, out=mean_first)
    first -= mean_first
    mean_first *= 1e-5
    valid = (first > 1e-10) & (first > mean_first)

    np.square(second, out=second)
    smooth(second)
    np.square(mean_second, out=mean_second)
    second -= mean_second
    mean_second *= 1e-5
    valid &= (second > 1e-10) & (second > mean_second)

    np.multiply(first, second, out=first)
    np.maximum(first, 0, out=first)
    np.sqrt(first, out=first)
    np.divide(covariance, first, out=covariance, where=valid)
    covariance[~valid] = 0
    np.clip(covariance, -1, 1, out=covariance)
    first[:] = valid
    return covariance, first


def shape_correlations(moving, reference, resolution=8.0, window=12.0,
                       moving_weights=None, reference_weights=None):
    """Sample local density CC at moving atom centers in a shared coordinate frame.

    window is the FWHM of the Gaussian correlation neighbourhood, in Angstroms.
    Input arrays are untouched; output NaNs represent undefined correlations.
    """
    from time import perf_counter
    import numpy as np
    from chimerax.core.errors import UserError
    from chimerax.map.molmap import add_gaussians, bounding_grid
    from chimerax.map_data import ArrayGridData
    from chimerax.map_fit.fitmap import volume_values

    started = perf_counter()
    moving = np.asarray(moving, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    # Keep numerical coordinates small even for deposited models far from zero.
    origin = reference.mean(axis=0)
    moving, reference = moving - origin, reference - origin
    density_sigma = resolution / (pi * sqrt(2))
    window_sigma = window / sqrt(8 * np.log(2))
    step = resolution / 3
    # The full correlation kernel must fit around even the outermost atoms.
    pad = 4 * density_sigma + 3 * window_sigma + step
    bounds = np.concatenate((moving, reference))
    dimensions = np.ceil((np.ptp(bounds, axis=0) + 2 * pad) / step)
    if np.prod(dimensions) > _MAX_CORRELATION_VOXELS:
        raise UserError("The correlation grid is too large; increase correlation_resolution "
                        "or restrict the fitting regions.")
    first_grid = bounding_grid(bounds, step, pad)
    first = first_grid.matrix()
    second_grid = ArrayGridData(np.zeros_like(first), first_grid.origin, first_grid.step)
    for xyz, grid, mass in ((moving, first_grid, moving_weights), (reference, second_grid, reference_weights)):
        points, weights = _shape_points(xyz, resolution / 6, mass)
        add_gaussians(grid, points, weights, density_sigma, 4)
    field, valid = _local_correlation_field(first, second_grid.matrix(), window_sigma / step)
    points = moving.astype(np.float32)
    transform = first_grid.xyz_to_ijk_transform
    scores = volume_values(points, transform, field).copy()
    support = volume_values(points, transform, valid)
    # Do not interpolate missing values as if they were poor (zero) fits.
    scores[support < 0.999] = np.nan
    return scores, {
        "correlation_seconds": perf_counter() - started,
        "correlation_resolution": resolution, "correlation_window": window,
        "correlation_grid_voxels": int(field.size),
        "correlation_defined_atoms": int(np.isfinite(scores).sum()),
    }


class _ShapeCorrelationAttributes:
    """Bulk attribute undo that also restores the absence of an old attribute."""

    def __init__(self, items):
        self.items = items

    def restore(self, values):
        for item, value in zip(self.items, values):
            if item.deleted:
                continue
            if value is None:
                if hasattr(item, "shape_correlation"):
                    delattr(item, "shape_correlation")
            else:
                item.shape_correlation = value


def _assign_shape_correlations(session, model, atoms, scores, undo, color_fit):
    import numpy as np
    from chimerax.atomic import Atom, Residue
    from chimerax.core.colors import Colormap

    for cls in (Atom, Residue):
        cls.register_attr(session, "shape_correlation", "Shape alignment", attr_type=float)
    residues = atoms.unique_residues
    indices = residues.indices(atoms.residues)
    finite = np.isfinite(scores)
    sums = np.bincount(indices[finite], weights=scores[finite], minlength=len(residues))
    counts = np.bincount(indices[finite], minlength=len(residues))
    means = np.divide(sums, counts, out=np.full(len(residues), np.nan), where=counts > 0)
    atom_values = {atom: float(value) for atom, value in zip(atoms.filter(finite), scores[finite])}
    defined_residues = np.isfinite(means)
    residue_values = {r: float(value) for r, value in
                      zip(residues.filter(defined_residues), means[defined_residues])}
    # Replace the previous score set on this moving model, including stale
    # values outside a newly restricted ROI. Leave their colors untouched.
    for collection, current in ((model.atoms, atom_values), (model.residues, residue_values)):
        items = [item for item in collection if item in current or hasattr(item, "shape_correlation")]
        old = [getattr(item, "shape_correlation", None) for item in items]
        new = [current.get(item) for item in items]
        owner = _ShapeCorrelationAttributes(items)
        undo.add(owner, "restore", old, new, option="M")
        owner.restore(new)
    if not color_fit:
        return
    cmap = Colormap((0, 0.5, 1), ((0.85, 0.15, 0.1, 1), (1, 1, 1, 1), (0.1, 0.3, 0.9, 1)))

    def colors(values, old):
        good = np.isfinite(values)
        result = cmap.interpolated_rgba8(np.nan_to_num(values, nan=0).astype(np.float32))
        result[~good, :3] = 160  # Undefined is gray, not a low correlation.
        result[:, 3] = old[:, 3]
        return result

    for objects, attribute, values in ((atoms, "colors", scores), (residues, "ribbon_colors", means)):
        old = getattr(objects, attribute).copy()
        new = colors(values, old)
        undo.add(objects, attribute, old, new)
        setattr(objects, attribute, new)


def _shape_points(xyz, spacing, weights=None):
    """Collapse cells to mass-weighted centroids, retaining confidence mass."""
    import numpy as np

    cells = np.floor(xyz / spacing).astype(np.int64)
    cells -= cells.min(axis=0)
    keys = np.ravel_multi_index(cells.T, tuple(cells.max(axis=0) + 1))
    _, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
    if weights is not None:
        mass = np.bincount(inverse, weights=weights)
        keep = mass > 0
        centers = np.column_stack([
            np.bincount(inverse, weights=xyz[:, axis] * weights)[keep] / mass[keep]
            for axis in range(3)
        ])
        return centers.astype(np.float32), mass[keep].astype(np.float32)
    centers = np.column_stack([
        np.bincount(inverse, weights=xyz[:, axis]) / counts for axis in range(3)
    ])
    return centers.astype(np.float32), counts.astype(np.float32)


def _shape_grid(xyz, weights, resolution):
    import numpy as np
    from chimerax.core.errors import UserError
    from chimerax.map.molmap import add_gaussians, bounding_grid

    sigma = resolution / (pi * sqrt(2))  # Same convention as molmap.
    step, pad = resolution / 3, 4 * sigma
    dimensions = np.ceil((np.ptp(xyz, axis=0) + 2 * pad) / step)
    if np.prod(dimensions) > _MAX_GRID_VOXELS:
        raise UserError("The fitting grid is too large; increase resolution or restrict the fitting regions.")
    grid = bounding_grid(xyz, step, pad)
    add_gaussians(grid, xyz, weights, sigma, 4)
    return grid


def _shape_samples(grid, max_points):
    import numpy as np

    matrix = grid.matrix()
    # Include low-density outskirts, not only the dense molecular interior.
    indices = np.flatnonzero(matrix >= matrix.max() * 0.001)
    count = len(indices)
    if count > max_points:
        chosen = np.random.default_rng(0).choice(count, max_points, replace=False)
        indices = indices[np.sort(chosen)]
    k, j, i = np.unravel_index(indices, matrix.shape)
    points = np.column_stack((i, j, k)).astype(np.float32)
    grid.ijk_to_xyz_transform.transform_points(points, in_place=True)
    weights = matrix.ravel()[indices].copy()
    norm = float(np.einsum("ijk,ijk->", matrix, matrix, dtype=np.float64))
    return points, weights, count / len(indices), norm


def _gradient_density_field(grid):
    """Negative Laplacian and full-density gradient energy on a zero-padded grid.

    Integration by parts gives <grad m, grad r> = <m, -Laplacian r>.
    A scalar target field therefore lets the native fitter optimize the gradient
    objective directly, including its negative penalties around density peaks.
    """
    import numpy as np
    from scipy.ndimage import laplace
    from chimerax.core.errors import UserError

    density = grid.matrix()
    field = laplace(density, mode="constant", cval=0)
    field *= -1 / grid.step[0] ** 2
    energy = float(np.einsum("ijk,ijk->", density, field, dtype=np.float64))
    if not np.isfinite(energy) or energy <= 0:
        raise UserError("The simulated density has no finite gradient energy.")
    return field, energy


def _shape_starts(moving, reference, extra_starts=0, moving_weights=None, reference_weights=None):
    """Try every proper signed permutation of the principal-axis frames."""
    from itertools import permutations, product
    import numpy as np
    from chimerax.geometry import Place

    def frame(xyz, weights):
        center = xyz.mean(axis=0, dtype=np.float64) if weights is None else _shape_center(xyz, weights)
        centered = xyz - center
        _, axes = np.linalg.eigh(centered.T @ centered if weights is None else
                                 (centered.T * weights) @ centered)
        return center, axes

    cm, am = frame(moving, moving_weights)
    cr, ar = frame(reference, reference_weights)
    starts = []
    for order in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            rotation = (ar[:, order] * signs) @ am.T
            if np.linalg.det(rotation) > 0:
                starts.append(Place(matrix=np.column_stack((rotation, cr - rotation @ cm))))
    # A good existing orientation should not be discarded just because axes
    # become unstable in a nearly spherical or symmetric structure.
    starts.append(Place(origin=cr - cm))
    starts.append(Place())
    if extra_starts:
        # Reproducible broad SO(3) search for unstable principal axes or very
        # different compositions. These are rotations, never reflections.
        from chimerax.geometry import quaternion_rotation

        quaternions = np.random.default_rng(0).normal(size=(extra_starts, 4))
        quaternions /= np.linalg.norm(quaternions, axis=1)[:, None]
        for quaternion in quaternions:
            rotation = quaternion_rotation(quaternion).matrix[:, :3]
            starts.append(Place(matrix=np.column_stack((rotation, cr - rotation @ cm))))
    return starts


def _partial_shape_starts(moving, reference, resolution, extra_starts=0, status=None,
                          objective="density", max_starts=24,
                          moving_weights=None, reference_weights=None):
    """Search rotations and all coarse translations of a smaller density template.

    Density mode uses masked Pearson normalization. Gradient mode instead uses
    full-input gradient energies, with no local mask normalization or trimming.
    """
    import numpy as np
    from scipy.fft import rfftn, irfftn, next_fast_len
    from chimerax.core.errors import UserError
    from chimerax.geometry import Place, quaternion_rotation
    from chimerax.map.molmap import add_gaussians
    from chimerax.map_data import ArrayGridData

    step = resolution / 3
    sigma = resolution / (pi * sqrt(2))
    mx, mw = _shape_points(moving, step / 2, moving_weights)
    rx, rw, target_grid = _shape_reference_atom_grid(reference, resolution, reference_weights)
    target = target_grid.matrix()
    gradient = objective == "gradient"
    if gradient:
        target, target_energy = _shape_reference_field(target_grid, objective)
    radius = float(np.linalg.norm(moving, axis=1).max())
    side = int(np.ceil((2 * radius + 8 * sigma) / step))
    side += 1 - side % 2  # Odd side: the template center is an exact grid point.
    fft_shape = tuple(next_fast_len(size + side - 1) for size in target.shape)
    if np.prod(fft_shape) > _MAX_PARTIAL_FFT_VOXELS:
        raise UserError(f"The {'gradient' if gradient else 'partial'}-search FFT grid is too large; increase resolution "
                        "or restrict the fitting regions.")
    half = side // 2
    template = ArrayGridData(np.zeros((side,) * 3, np.float32),
                             (-step * half,) * 3, (step,) * 3)
    target_fft = _shape_reference_fft(target_grid, target, fft_shape, objective)
    target_squared_fft = None if gradient else _shape_reference_fft(
        target_grid, target, fft_shape, objective, squared=True)
    rotations = [pose.matrix[:, :3] for pose in _shape_starts(
        moving, reference, moving_weights=moving_weights, reference_weights=reference_weights)]
    def proper_axes(xyz, weights):
        _, axes = np.linalg.eigh(xyz.T @ xyz if weights is None else (xyz.T * weights) @ xyz)
        if np.linalg.det(axes) < 0:
            axes[:, 0] *= -1
        return axes

    moving_axes = proper_axes(moving, moving_weights)
    reference_axes = proper_axes(reference, reference_weights)
    quaternions = np.random.default_rng(0).normal(size=(96 + extra_starts, 4))
    quaternions /= np.linalg.norm(quaternions, axis=1)[:, None]
    rotations.extend(reference_axes @ quaternion_rotation(q).matrix[:, :3] @ moving_axes.T
                     for q in quaternions)
    # A subcomplex's principal axes need not resemble the whole target's axes.
    # Use template-sized target neighbourhoods as additional orientation hints;
    # every hint still scans all translations. This maintains angular capture
    # when a small object correctly gets a finer initial density resolution.
    anchors = []
    if len(reference) > len(moving):
        anchor = _shape_center(reference, reference_weights)
        nearest = np.full(len(rx), np.inf)
        for _ in range(min(8, len(rx))):
            anchors.append(anchor)
            squared = np.sum((reference - anchor) ** 2, axis=1)
            indices = np.argpartition(squared, len(moving) - 1)[:len(moving)]
            rotations.extend(pose.matrix[:, :3] for pose in
                             _shape_starts(moving, reference[indices], moving_weights=moving_weights,
                                           reference_weights=None if reference_weights is None else reference_weights[indices]))
            nearest = np.minimum(nearest, np.sum((rx - anchor) ** 2, axis=1))
            anchor = rx[np.argmax(nearest)]
    peaks = []
    evaluated = set()

    def evaluate(rotation):
        key = tuple(np.round(rotation, 6).ravel())
        if key in evaluated:
            return
        evaluated.add(key)
        density = template.matrix()
        density[:] = 0
        add_gaussians(template, mx @ rotation.T, mw, sigma, 4)
        if gradient:
            _, energy = _gradient_density_field(template)
            cross = irfftn(target_fft * rfftn(density[::-1, ::-1, ::-1], fft_shape), fft_shape)
            cross /= sqrt(energy * target_energy)
        else:
            mask = (density > density.max() * 0.001).astype(np.float32)
            count = mask.sum(dtype=np.float64)
            # Use the same support for the source moments and the target window.
            density *= mask
            total = density.sum(dtype=np.float64)
            norm = np.einsum("ijk,ijk->", density, density, dtype=np.float64) - total * total / count
            mask_fft = rfftn(mask[::-1, ::-1, ::-1], fft_shape)
            cross = irfftn(target_fft * rfftn(density[::-1, ::-1, ::-1], fft_shape), fft_shape)
            local_sum = irfftn(target_fft * mask_fft, fft_shape)
            local_squared = irfftn(target_squared_fft * mask_fft, fft_shape)
            cross -= total * local_sum / count
            variance = local_squared - local_sum * local_sum / count
            denominator = np.sqrt(np.maximum(variance * norm, 1e-30))
            cross /= denominator
            cross[variance < max(float(local_squared.max()) * 1e-5, 1e-15)] = -1
            np.clip(cross, -1, 1, out=cross)
        # Linear convolution, not circular wraparound. This slice indexes
        # template centers at actual target-grid coordinates.
        view = cross[half:half + target.shape[0], half:half + target.shape[1],
                     half:half + target.shape[2]]
        for _ in range(2):
            k, j, i = np.unravel_index(np.argmax(view), view.shape)
            center = target_grid.ijk_to_xyz_transform * np.array((i, j, k), dtype=float)
            pose = Place(matrix=np.column_stack((rotation, center)))
            peaks.append((float(view[k, j, i]), pose))
            distance = max(2, int(radius / (3 * step)))
            view[max(0, k-distance):k+distance+1, max(0, j-distance):j+distance+1,
                 max(0, i-distance):i+distance+1] = -1

    for index, rotation in enumerate(rotations):
        _shape_checkpoint()
        if status and index % 16 == 0:
            status(f"{'Gradient' if gradient else 'Partial'} shape search: orientation {index + 1}/{len(rotations)} "
                   f"at {resolution:g} A (all translations)...")
        evaluate(rotation)
    peaks.sort(key=lambda item: item[0], reverse=True)
    probes = mx[np.linspace(0, len(mx) - 1, min(64, len(mx)), dtype=int)]
    starts = []
    for _, pose in peaks:
        placed = pose.transform_points(probes)
        if any(np.sqrt(np.mean(np.sum((placed - kept.transform_points(probes)) ** 2, axis=1)))
               < resolution / 2 for kept in starts):
            continue
        starts.append(pose)
        if len(starts) == max_starts:
            break
    return starts, {"partial_orientations": len(evaluated),
                    "partial_neighborhood_frames": len(anchors),
                    "partial_fft_voxels": int(np.prod(fft_shape))}


def _refine_shape_level(moving, reference, blur, candidates, max_points,
                        partial=False, coarse=False, objective="density",
                        moving_weights=None, reference_weights=None):
    """Refine and rank poses against the same complete input densities."""
    import numpy as np
    from chimerax.map_fit.fitmap import locate_maximum, volume_values, overlap_and_correlation

    # Cell averaging preserves atom counts and hence relative density mass.
    mx, mw = _shape_points(moving, blur / 6, moving_weights)
    mg = _shape_grid(mx, mw, blur)
    _, _, rg = _shape_reference_atom_grid(reference, blur, reference_weights)
    points, weights, sampling, mnorm = _shape_samples(mg, max_points)
    rmatrix, reference_energy = _shape_reference_field(rg, objective)
    if objective == "gradient":
        _, moving_energy = _gradient_density_field(mg)
        denominator = sqrt(moving_energy * reference_energy)
    else:
        denominator = sqrt(mnorm * reference_energy)
    results = []
    for seed_index, initial in enumerate(candidates):
        _shape_checkpoint()
        placed = initial.transform_points(points)
        delta, _ = locate_maximum(
            placed, weights, rmatrix, rg.xyz_to_ijk_transform,
            max_steps=160 if coarse and not partial else 240,
            ijk_step_size_min=0.01 if partial else 0.02, ijk_step_size_max=0.5,
            metric="correlation about mean" if partial else "sum product",
        )
        transform = delta * initial
        values = volume_values(points, rg.xyz_to_ijk_transform * transform, rmatrix)
        score = (overlap_and_correlation(weights, values)[2] if partial else
                 float(np.dot(weights.astype(np.float64), values)) * sampling / denominator)
        results.append((score, transform, seed_index))
    results.sort(key=lambda item: item[0], reverse=True)
    return results, points


def _fit_gradient_points(moving, reference, resolution, max_points, extra_starts,
                          status, initial_transform, local, fallback_correlation,
                          moving_weights=None, reference_weights=None):
    """Experimental whole-density gradient fit; the objective never switches.

    Try the cheap principal-axis search first. A weak score can trigger a wider
    FFT search at the next finer scale, where secondary-structure features are
    less blurred. Its candidates bypass coarser stages and retain a wider beam.
    """
    from time import perf_counter
    import numpy as np
    from chimerax.geometry import Place, translation

    started = perf_counter()
    moving = np.asarray(moving, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    cm, cr = _shape_center(moving, moving_weights), _shape_center(reference, reference_weights)
    moving, reference = moving - cm, reference - cr
    levels, schedule_stats = _shape_resolution_schedule(
        moving, reference, resolution, local, moving_weights, reference_weights)
    if local:
        initial = Place() if initial_transform is None else initial_transform
        starts = [translation(-cr) * initial * translation(cm)]
    else:
        starts = _shape_starts(moving, reference, extra_starts, moving_weights, reference_weights)
        if initial_transform is not None:
            starts.append(translation(-cr) * initial_transform * translation(cm))

    def refine(candidates, first_level=0):
        counts = []
        for level, blur in enumerate(levels[first_level:], first_level):
            if status:
                status(f"Whole-structure gradient fit: {blur:.3g} A ({len(candidates)} starts)...")
            results, points = _refine_shape_level(
                moving, reference, blur, candidates, max_points, coarse=level == 0,
                objective="gradient",
                moving_weights=moving_weights, reference_weights=reference_weights,
            )
            counts.append(len(points))
            candidates = []
            probes = points[np.linspace(0, len(points) - 1, min(64, len(points)), dtype=int)]
            for _, pose, _ in results:
                placed = pose.transform_points(probes)
                if any(np.sqrt(np.mean(np.sum((placed - kept.transform_points(probes)) ** 2, axis=1)))
                       < blur / 3 for kept in candidates):
                    continue
                candidates.append(pose)
                if len(candidates) == 48:
                    break
        return results[0][:2], counts

    (score, transform), counts = refine(starts)
    stats = {"objective": "gradient", "partial": False, "local": local,
             "gradient_search": False, "resolutions": levels, "fit_points": counts,
             "starts": len(starts), "moving_atoms": len(moving), "reference_atoms": len(reference),
             **schedule_stats}
    if not local and fallback_correlation > 0 and score < fallback_correlation:
        # Choosing the smaller numerical template reduces FFT work. Its norm
        # still includes all its density, as does the target's norm: not a
        # partial-template or chain fit, despite sharing the search machinery.
        reverse = len(moving) > len(reference)
        source, target = (reference, moving) if reverse else (moving, reference)
        sw, tw = (reference_weights, moving_weights) if reverse else (moving_weights, reference_weights)
        first_level = 1 if len(levels) > 2 else 0
        search_resolution = levels[first_level]
        seeds, search_stats = _partial_shape_starts(
            source, target, search_resolution, extra_starts, status,
            objective="gradient", max_starts=48,
            moving_weights=sw, reference_weights=tw,
        )
        if reverse:
            seeds = [pose.inverse() for pose in seeds]
        (candidate_score, candidate), _ = refine([transform, *seeds], first_level)
        improved = candidate_score > score
        stats.update(gradient_search=True, gradient_search_improved=improved,
                     gradient_initial_score=score, gradient_search_resolution=search_resolution,
                     gradient_search_starts=len(seeds),
                     gradient_search_orientations=search_stats["partial_orientations"],
                     gradient_search_fft_voxels=search_stats["partial_fft_voxels"])
        if improved:
            score, transform = candidate_score, candidate
    stats.update(score=score, seconds=perf_counter() - started)
    return translation(cr) * transform * translation(-cm), stats


def fit_shape_points(moving, reference, resolution=None, max_points=12000,
                     extra_starts=0, status=None, initial_transform=None, local=False,
                     partial=None, fallback_correlation=None, objective="density",
                     moving_weights=None, reference_weights=None):
    """Return a rigid transform and diagnostics; input arrays are never changed.

    local=True refines only initial_transform (identity if omitted). Centering
    never translates either cloud into alignment as a separate starting step.
    In density mode, partial=None uses partial-template search for unequal atom
    counts or as a backup below fallback_correlation. True forces that search;
    False disables it. Local refinement never restarts globally.
    resolution=None derives a search pyramid from the smaller fitted ellipsoid.
    An explicit resolution overrides only the final endpoint, not the size-based
    coarse scale. Correlation-coloring scales are independent of this pyramid.
    objective="gradient" instead uses full-density gradient normalization and
    never enters partial fitting. Its same-objective FFT backup defaults to a
    0.5 cutoff, versus 0.7 for the existing density mode's partial backup.
    Optional per-atom weights affect geometry and density at every scale;
    zero-weight atoms are excluded before sizing grids. Arrays are untouched.
    """
    from time import perf_counter
    import numpy as np
    from chimerax.geometry import Place, translation
    from chimerax.core.errors import UserError

    if objective not in ("density", "gradient"):
        raise UserError("objective must be density or gradient.")
    if fallback_correlation is None:
        fallback_correlation = 0.5 if objective == "gradient" else 0.7
    def weighted_input(xyz, weights):
        xyz = np.asarray(xyz, dtype=np.float64)
        if weights is None:
            return xyz, None
        weights = np.asarray(weights, dtype=np.float64)
        if weights.shape != (len(xyz),) or not np.isfinite(weights).all() or np.any(weights < 0):
            raise UserError("Shape weights must be finite, nonnegative, and match the fitting atoms.")
        keep = weights > 0
        xyz, weights = xyz[keep], weights[keep]
        if len(xyz) < 4 or np.linalg.matrix_rank(xyz - _shape_center(xyz, weights)) < 2:
            raise UserError("Confidence weighting leaves fewer than four non-collinear fitting atoms.")
        return xyz, weights

    moving, moving_weights = weighted_input(moving, moving_weights)
    reference, reference_weights = weighted_input(reference, reference_weights)
    if objective == "gradient":
        if partial is True:
            raise UserError("objective gradient uses whole structures; do not combine it with partial true.")
        return _fit_gradient_points(moving, reference, resolution, max_points, extra_starts,
                                    status, initial_transform, local, fallback_correlation,
                                    moving_weights, reference_weights)

    started = perf_counter()
    moving = np.asarray(moving, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    moving_count = _shape_effective_count(moving, moving_weights)
    reference_count = _shape_effective_count(reference, reference_weights)
    use_partial = not local and (partial is True or (partial is None and
        min(moving_count, reference_count) / max(moving_count, reference_count) < _PARTIAL_SIZE_RATIO))
    if use_partial and moving_count > reference_count:
        # Always fit the smaller template into the larger target. Invert the
        # result, not the user-visible choice of which model is allowed to move.
        reverse, stats = fit_shape_points(
            reference, moving, resolution, max_points, extra_starts, status,
            initial_transform.inverse() if initial_transform is not None else None,
            partial=True,
            moving_weights=reference_weights, reference_weights=moving_weights,
        )
        stats.update(moving_atoms=len(moving), reference_atoms=len(reference),
                     partial_template_reversed=True, seconds=perf_counter() - started,
                     partial_trigger="size" if partial is None else "requested")
        stats["moving_size_diameter"], stats["reference_size_diameter"] = (
            stats["reference_size_diameter"], stats["moving_size_diameter"])
        return reverse.inverse(), stats
    cm, cr = _shape_center(moving, moving_weights), _shape_center(reference, reference_weights)
    moving, reference = moving - cm, reference - cr
    levels, schedule_stats = _shape_resolution_schedule(
        moving, reference, resolution, local, moving_weights, reference_weights)
    if status:
        label = "Automatic" if resolution is None else "Size-based, explicit final"
        status(f"{label} resolution range: " + " -> ".join(f"{value:.3g}" for value in levels) +
               f" A (fitted size {schedule_stats['object_size']:.3g} A)...")
    if local:
        # Centering is numerical bookkeeping, not a new initial alignment.
        # Preserve the existing pose even if the cropped clouds have different
        # centroids; no PCA, centroid-matching, or broad orientation search.
        initial = Place() if initial_transform is None else initial_transform
        starts = [translation(-cr) * initial * translation(cm)]
    elif use_partial:
        starts, partial_stats = _partial_shape_starts(
            moving, reference, levels[0], extra_starts, status,
            moving_weights=moving_weights, reference_weights=reference_weights)
        if initial_transform is not None:
            starts.append(translation(-cr) * initial_transform * translation(cm))
    else:
        starts = _shape_starts(moving, reference, extra_starts, moving_weights, reference_weights)
        if initial_transform is not None:
            starts.append(translation(-cr) * initial_transform * translation(cm))
    candidates = starts
    counts = []
    for level, blur in enumerate(levels):
        if status:
            mode = "Local shape alignment" if local else "Partial shape alignment" if use_partial else "Shape alignment"
            status(f"{mode}: {blur:g} A density ({len(candidates)} starts)...")
        results, points = _refine_shape_level(moving, reference, blur, candidates,
                                              max_points, partial=use_partial, coarse=level == 0,
                                              moving_weights=moving_weights, reference_weights=reference_weights)
        # Keep several basins through the intermediate scale; symmetries and
        # attached factors can make the coarsest-scale winner ambiguous.
        candidates = []
        probes = points[np.linspace(0, len(points) - 1, min(64, len(points)), dtype=int)]
        for _, transform, _ in results:
            placed = transform.transform_points(probes)
            if any(np.sqrt(np.mean(np.sum((placed - kept.transform_points(probes)) ** 2, axis=1)))
                   < blur / 3 for kept in candidates):
                continue
            candidates.append(transform)
            if len(candidates) == (12 if level == 0 else 6):
                break
        counts.append(len(points))
    score, transform, _ = results[0]
    if not local and partial is None and not use_partial and score < fallback_correlation:
        if status:
            status(f"Shape density score {score:.3f} is below {fallback_correlation:g}; "
                   "trying partial-template search...")
        refined, stats = fit_shape_points(
            moving, reference, resolution, max_points, extra_starts, status,
            initial_transform=transform, partial=True,
            moving_weights=moving_weights, reference_weights=reference_weights,
        )
        stats.update(seconds=perf_counter() - started, partial_trigger="low_correlation",
                     fallback_initial_score=score, fallback_correlation=fallback_correlation)
        return translation(cr) * refined * translation(-cm), stats
    transform = translation(cr) * transform * translation(-cm)
    stats = {
        "score": score, "seconds": perf_counter() - started,
        "resolutions": levels, "fit_points": counts, "starts": len(starts),
        "moving_atoms": len(moving), "reference_atoms": len(reference),
        "local": local,
        "partial": use_partial,
        **schedule_stats,
    }
    if use_partial:
        stats.update(partial_stats, partial_template_reversed=False,
                     partial_trigger="size" if partial is None else "requested")
    return transform, stats


def _shape_chain_regions(atoms):
    """Substantial selected polymer chains, not ligands or the whole input."""
    import numpy as np
    from chimerax.atomic import Residue

    residues = atoms.residues
    polymer = residues.polymer_types != Residue.PT_NONE
    chain_ids = residues.chain_ids
    coordinates = atoms.coords
    regions = []
    for chain_id in np.unique(chain_ids[polymer]):
        indices = np.flatnonzero(polymer & (chain_ids == chain_id))
        if max(4, _CHAIN_BACKUP_MIN_FRACTION * len(atoms)) <= len(indices) < 0.9 * len(atoms):
            xyz = coordinates[indices]
            if np.linalg.matrix_rank(xyz - xyz.mean(axis=0)) >= 2:
                regions.append((chain_id, indices))
    return regions


def _chain_shape_backup(moving, reference, source_regions, target_regions,
                        transform, stats, resolution, max_points, extra_starts, logger,
                        moving_weights=None, reference_weights=None):
    """Generate chain seeds, but accept them only on a common whole-input score.

    Non-nested complexes can contain unmatched parts on both sides. A chain
    seed must reach the final scale directly: whole-object coarse blurring can
    otherwise drag it back into a wrong basin before helix detail is available.
    """
    import numpy as np
    from chimerax.geometry import translation
    from chimerax.core.errors import UserError

    moving = np.asarray(moving, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    regions = [(len(indices), False, label, indices) for label, indices in source_regions]
    regions += [(len(indices), True, label, indices) for label, indices in target_regions]
    regions.sort(key=lambda item: (-item[0], item[1], item[2]))
    poses, labels = [transform], [None]
    for _, reverse, label, indices in regions[:_MAX_CHAIN_BACKUP_STARTS]:
        logger.status(f"Shape alignment backup: {'reference' if reverse else 'moving'} chain {label}...")
        source, target = (reference, moving) if reverse else (moving, reference)
        sw, tw = (reference_weights, moving_weights) if reverse else (moving_weights, reference_weights)
        try:
            pose, _ = fit_shape_points(source[indices], target, resolution, max_points, extra_starts,
                                       status=logger.status, partial=True, fallback_correlation=0,
                                       moving_weights=None if sw is None else sw[indices], reference_weights=tw)
        except UserError as error:
            # An optional seed must not invalidate an existing fit, for example
            # if its finer size-derived FFT exceeds the allocation guard.
            logger.warning(f"Skipped shape-alignment chain seed {label}: {error}")
            continue
        poses.append(pose.inverse() if reverse else pose)
        labels.append(f"{'reference' if reverse else 'moving'} chain {label}")
    stats.update(chain_backup_attempts=min(len(regions), _MAX_CHAIN_BACKUP_STARTS),
                 chain_backup_used=False)
    if len(poses) == 1:
        return transform

    # Match the original partial fit's template direction, grid, and sample set.
    # Comparing scores from separately fitted chains would favor small subsets.
    reverse = stats["partial_template_reversed"]
    source, target = (reference, moving) if reverse else (moving, reference)
    sw, tw = (reference_weights, moving_weights) if reverse else (moving_weights, reference_weights)
    cm, cr = _shape_center(source, sw), _shape_center(target, tw)
    starts = [translation(-cr) * (pose.inverse() if reverse else pose) * translation(cm)
              for pose in poses]
    results, _ = _refine_shape_level(source - cm, target - cr, stats["resolutions"][-1],
                                      starts, max_points, partial=True, moving_weights=sw, reference_weights=tw)
    score, best, seed_index = results[0]
    if score <= stats["score"] + 1e-6:
        return transform
    result = translation(cr) * best * translation(-cm)
    result = result.inverse() if reverse else result
    stats.update(chain_backup_used=True, chain_backup_initial_score=stats["score"],
                 chain_backup_seed=labels[seed_index] or "original pose", score=score)
    return result


def _prepare_shape_atomic_fit(session, atoms, to, resolution=None, max_points=12000,
                     extra_starts=0, move=True, radius=None, correlation=False,
                     color_fit=False, correlation_resolution=8.0, correlation_window=12.0,
                     sphere_center=None, partial=None, fallback_correlation=None, objective="density",
                     confidence="auto", confidence_moving=None, confidence_reference=None, copy=False):
    """Validate an atomic fit and snapshot its atoms, options, and scene frames.

    Exclude water/hydrogens and optionally crop both fitting subsets to the
    starting CoFR sphere. No density calculation or model movement happens here.
    Keep native model/atom handles on the main thread; workers receive only the
    numeric subset returned by _shape_atomic_snapshot().
    """
    from math import isfinite
    import numpy as np
    from chimerax.atomic import Residue
    from chimerax.core.errors import UserError

    if resolution is not None and (not isfinite(resolution) or resolution < 2):
        raise UserError("resolution must be finite and at least 2 Angstroms, or omitted for automatic scaling.")
    if max_points < 1000:
        raise UserError("max_points must be at least 1000.")
    if extra_starts < 0:
        raise UserError("extra_starts must be nonnegative.")
    if partial not in (None, True, False):
        raise UserError("partial must be true or false, or omitted for automatic selection.")
    if objective not in ("density", "gradient"):
        raise UserError("objective must be density or gradient.")
    if objective == "gradient" and partial is True:
        raise UserError("objective gradient uses whole structures; do not combine it with partial true.")
    if fallback_correlation is None:
        fallback_correlation = 0.5 if objective == "gradient" else 0.7
    if not isfinite(fallback_correlation) or not 0 <= fallback_correlation <= 1:
        raise UserError("fallback_correlation must be a finite value between 0 and 1.")
    score_fit = correlation or color_fit
    for name, value in (("correlation_resolution", correlation_resolution),
                        ("correlation_window", correlation_window)):
        if not isfinite(value) or value <= 0:
            raise UserError(f"{name} must be finite and greater than zero (Angstroms).")
    if score_fit and not move:
        raise UserError("correlation/color_fit require move true so scores describe the applied pose.")
    if radius is not None:
        if not isfinite(radius) or radius <= 0:
            raise UserError("radius must be finite and greater than zero (Angstroms).")
        if extra_starts:
            raise UserError("extra_starts cannot be used with radius_refine; local fitting refines only the current pose.")
        if partial:
            raise UserError("partial true cannot be used with radius_refine; local fitting must keep the current pose.")
    structures = atoms.unique_structures
    targets = to.unique_structures
    if len(structures) != 1 or len(targets) != 1 or structures[0] is targets[0]:
        raise UserError("Specify one moving model and one different reference model.")
    model, reference = structures[0], targets[0]
    if model in reference.drawing_lineage or reference in model.drawing_lineage:
        raise UserError("The moving and reference models must not contain one another.")
    for structure in (model, reference):
        parent = structure
        while parent is not None:
            if len(parent.positions) != 1:
                raise UserError("Instanced models are not supported; use a single model copy.")
            parent = parent.parent

    center = None
    if radius is not None:
        # Respect a manually chosen CoFR; do not recenter on the view or model.
        center = np.array(session.main_view.center_of_rotation if sphere_center is None else sphere_center,
                          dtype=np.float64, copy=True)
        if center.shape != (3,) or not np.isfinite(center).all():
            raise UserError("The current center of rotation must be a finite XYZ position.")

    def fitting_atoms(selected):
        keep = (selected.element_numbers > 1) & ~np.isin(
            selected.residues.names, list(Residue.water_res_names)
        )
        # Grid in each model's local frame: arbitrary scene rotations then
        # do not change voxelization, and stale atomic scene caches are avoided.
        selected = selected.filter(keep)
        xyz = selected.coords
        if not np.isfinite(xyz).all():
            raise UserError("Each fitting region needs finite atom coordinates.")
        if center is not None and len(selected):
            # Use the full model hierarchy, not cached atomic scene coordinates.
            offsets = selected[0].structure.scene_position.transform_points(xyz) - center
            inside = np.linalg.norm(offsets, axis=1) <= radius
            xyz, selected = xyz[inside], selected.filter(inside)
        if len(xyz) < 4:
            where = f" within {radius:g} A of the center of rotation" if radius is not None else ""
            raise UserError(f"Each fitting region needs at least four non-water heavy atoms{where}.")
        if np.linalg.matrix_rank(xyz - xyz.mean(axis=0)) < 2:
            raise UserError("A collinear fitting region cannot define a unique rigid orientation.")
        return selected, xyz

    source_scene, reference_scene = model.scene_position, reference.scene_position
    source_atoms, source_xyz = fitting_atoms(atoms)
    reference_atoms, reference_xyz = fitting_atoms(to)
    confidence_inputs = (source_atoms, reference_atoms)
    source_atoms, source_xyz, source_weights = _shape_weighted_atoms(
        source_atoms, source_xyz, confidence_moving or confidence)
    reference_atoms, reference_xyz, reference_weights = _shape_weighted_atoms(
        reference_atoms, reference_xyz, confidence_reference or confidence)
    _shape_log_confidence(session, model, source_weights, "moving")
    _shape_log_confidence(session, reference, reference_weights, "reference")
    options = dict(resolution=resolution, max_points=max_points, extra_starts=extra_starts,
                   move=move, radius=radius, correlation=correlation, color_fit=color_fit,
                   correlation_resolution=correlation_resolution, correlation_window=correlation_window,
                   partial=partial, fallback_correlation=fallback_correlation, objective=objective, copy=copy)
    backup = objective == "density" and radius is None and partial is not False and fallback_correlation > 0
    return dict(model=model, reference=reference, source_atoms=source_atoms, reference_atoms=reference_atoms,
                source_xyz=source_xyz.astype(np.float64), reference_xyz=reference_xyz.astype(np.float64),
                source_weights=source_weights, reference_weights=reference_weights,
                confidence_inputs=[(selected[0].structure, selected, selected.pointers.copy(), selected.bfactors.copy())
                                   for selected, weights in zip(confidence_inputs, (source_weights, reference_weights))
                                   if weights is not None],
                source_scene=source_scene, reference_scene=reference_scene, center=center, options=options,
                source_pointers=source_atoms.pointers.copy(), reference_pointers=reference_atoms.pointers.copy(),
                source_regions=_shape_chain_regions(source_atoms) if backup else [],
                reference_regions=_shape_chain_regions(reference_atoms) if backup else [])


def _compute_shape_atomic_fit(data, logger):
    """Numerical work only: no session, atom, model, or GUI objects are required."""
    from time import perf_counter
    from chimerax.geometry import Place

    fitting_started = perf_counter()
    source_xyz, reference_xyz = data["source_xyz"], data["reference_xyz"]
    source_weights, reference_weights = data.get("source_weights"), data.get("reference_weights")
    options = data["options"]
    resolution, max_points, extra_starts = (options[name] for name in ("resolution", "max_points", "extra_starts"))
    radius, partial, objective = (options[name] for name in ("radius", "partial", "objective"))
    fallback_correlation = options["fallback_correlation"]
    local_transform, stats = fit_shape_points(
        source_xyz, reference_xyz, resolution, max_points, extra_starts,
        status=logger.status,
        initial_transform=Place(matrix=data["initial_matrix"]),
        local=radius is not None,
        partial=partial, fallback_correlation=fallback_correlation, objective=objective,
        moving_weights=source_weights, reference_weights=reference_weights,
    )
    if (objective == "density" and radius is None and partial is not False and stats["partial"] and
            stats["score"] < min(_CHAIN_BACKUP_CORRELATION, fallback_correlation)):
        local_transform = _chain_shape_backup(
            source_xyz, reference_xyz, data["source_regions"],
            data["reference_regions"], local_transform, stats,
            resolution, max_points, extra_starts, logger,
            source_weights, reference_weights,
        )
    stats["seconds"] = perf_counter() - fitting_started
    stats.update(moving_confidence="off" if source_weights is None else "plddt",
                 reference_confidence="off" if reference_weights is None else "plddt")
    return local_transform, stats


def _shape_atomic_snapshot(prepared):
    return dict(source_xyz=prepared["source_xyz"], reference_xyz=prepared["reference_xyz"],
                source_weights=prepared["source_weights"], reference_weights=prepared["reference_weights"],
                initial_matrix=(prepared["reference_scene"].inverse() * prepared["source_scene"]).matrix,
                options=prepared["options"], source_regions=prepared["source_regions"],
                reference_regions=prepared["reference_regions"])


class _ShapeCopyPresence:
    """Remove/re-add the same native structure for Undo/Redo, without refitting."""

    def __init__(self, session, model):
        self.session, self.model = session, model
        self.model_id = model.id

    def restore(self, present):
        model, models = self.model, self.session.models
        if model.deleted:
            return
        if not present and models.have_model(model):
            models.remove([model])
        elif present and not models.have_model(model):
            # Preserve the ID if still available; otherwise use a new free ID.
            if not models.list(model_id=self.model_id):
                model.id = self.model_id
            models.add([model])


def _shape_place_model(session, model, atoms, transform, undo, copy=False):
    """Only create a full-model copy after numerical fitting has succeeded."""
    if not copy:
        old = model.position
        parent = model.parent.scene_position
        new = parent.inverse() * transform * parent * old
        undo.add(model, "position", old, new)
        model.position = new
        return model, atoms

    indices = model.atoms.indices(atoms)
    duplicate = model.copy(name=f"{model.name} fit")
    try:
        # Flatten the original hierarchy into the new independent owner's pose.
        duplicate.position = transform * model.scene_position
        for attribute in ("database", "filename"):
            if hasattr(model, attribute):
                setattr(duplicate, attribute, getattr(model, attribute))
        session.models.add([duplicate])
        presence = _ShapeCopyPresence(session, duplicate)
        undo.add(presence, "restore", False, True, option="M")
        return duplicate, duplicate.atoms[indices]
    except BaseException:
        session.models.close([duplicate])
        raise


def _finish_shape_atomic_fit(session, prepared, result):
    """Apply one computed result and its optional coloring on the main thread."""
    from time import perf_counter
    from chimerax.core.undo import UndoState

    model, reference = prepared["model"], prepared["reference"]
    source_atoms = prepared["source_atoms"]
    source_xyz, reference_xyz = prepared["source_xyz"], prepared["reference_xyz"]
    source_scene, reference_scene = prepared["source_scene"], prepared["reference_scene"]
    center, options = prepared["center"], prepared["options"]
    radius, move, objective = (options[name] for name in ("radius", "move", "objective"))
    score_fit, color_fit = options["correlation"] or options["color_fit"], options["color_fit"]
    correlation_resolution, correlation_window = options["correlation_resolution"], options["correlation_window"]
    fallback_correlation = options["fallback_correlation"]
    local_transform, stats = result
    if center is not None:
        stats.update(radius=radius, radius_refine=radius, center_of_rotation=center.copy())
    transform = reference_scene * local_transform * source_scene.inverse()
    scores = None
    if score_fit:
        from time import perf_counter

        scoring_started = perf_counter()
        session.logger.status("Shape alignment: calculating local density correlations...")
        scores, correlation_stats = shape_correlations(
            local_transform.transform_points(source_xyz), reference_xyz,
            correlation_resolution, correlation_window,
            prepared["source_weights"], prepared["reference_weights"],
        )
        stats.update(correlation_stats)
    undo = UndoState("shape align")
    original_model = model
    try:
        if move:
            model, source_atoms = _shape_place_model(
                session, model, source_atoms, transform, undo, options.get("copy", False))
        if scores is not None:
            _assign_shape_correlations(session, model, source_atoms, scores, undo, color_fit)
    except BaseException:
        undo.undo()
        if model is not original_model:
            session.models.close([model])
        raise
    if undo.state:
        session.undo.register(undo)
    if score_fit:
        stats["correlation_grid_seconds"] = stats["correlation_seconds"]
        stats["correlation_seconds"] = perf_counter() - scoring_started
    stats.update(moving_model_id=model.id_string, reference_model_id=reference.id_string)
    if model is not original_model:
        stats["source_model_id"] = original_model.id_string
        session.logger.info(f"Shape alignment: fitted copy #{model.id_string}; original #{original_model.id_string} unchanged.")
    mode = ("Local shape alignment" if radius is not None else
            "Partial shape alignment" if stats["partial"] else "Shape alignment")
    metric = ("whole-density gradient correlation" if objective == "gradient" else
              "partial density correlation" if stats["partial"] else "sampled density overlap")
    reason = ""
    if stats.get("partial_trigger") == "low_correlation":
        reason = (f"fallback from initial score {stats['fallback_initial_score']:.3f} "
                  f"< {fallback_correlation:g}; ")
    elif stats.get("partial_trigger") == "size":
        reason = "automatic size-mismatch search; "
    region = (f"{stats['moving_atoms']}/{stats['reference_atoms']} atoms in a {radius:g} A CoFR sphere; "
              if radius is not None else "")
    message = (
        f"{mode} #{model.id_string} to #{reference.id_string}: {region}"
        f"{metric} {stats['score']:.3f}; {reason}{stats['seconds']:.2f} s; "
        f"{stats['resolutions'][-1]:.3g} A final resolution" + (" (not moved)." if not move else ".")
    )
    session.logger.info(message)
    session.logger.status(message)
    session.logger.info(
        f"Resolution range ({stats['resolution_mode']}, fitted size {stats['object_size']:.3g} A): " +
        " -> ".join(f"{value:.3g}" for value in stats["resolutions"]) + " A."
    )
    if stats.get("chain_backup_used"):
        session.logger.info(
            f"Chain-seeded backup ({stats['chain_backup_seed']}): whole-input partial density "
            f"correlation improved from {stats['chain_backup_initial_score']:.3f} to {stats['score']:.3f}. "
            "The whole model moved rigidly."
        )
    if stats.get("gradient_search"):
        session.logger.info(
            f"Whole-structure gradient FFT search at {stats['gradient_search_resolution']:.3g} A: "
            f"{stats['gradient_search_starts']} seeds; score "
            f"{stats['gradient_initial_score']:.3f} -> {stats['score']:.3f}. No partial or chain fitting."
        )
    if score_fit:
        session.logger.info(
            f"Local density correlation: {stats['correlation_defined_atoms']}/{len(source_atoms)} "
            f"atoms defined; density resolution {correlation_resolution:g} A, "
            f"Gaussian window FWHM {correlation_window:g} A; "
            f"{stats['correlation_seconds']:.2f} s extra including attributes/coloring. "
            "Saved shape_correlation on atoms and residue means" +
            ("; colored fitted atoms and cartoons." if color_fit else "; colors unchanged.")
        )
    return transform, stats


def _shape_align_one(session, atoms, to, resolution=None, max_points=12000,
                     extra_starts=0, move=True, radius=None, correlation=False,
                     color_fit=False, correlation_resolution=8.0, correlation_window=12.0,
                     sphere_center=None, partial=None, fallback_correlation=None, objective="density",
                     confidence="auto", confidence_moving=None, confidence_reference=None, copy=False):
    prepared = _prepare_shape_atomic_fit(
        session, atoms, to, resolution, max_points, extra_starts, move, radius, correlation,
        color_fit, correlation_resolution, correlation_window, sphere_center, partial,
        fallback_correlation, objective, confidence, confidence_moving, confidence_reference, copy)
    with _ShapeMapInputGuard((*atoms.by_structure, *to.by_structure)) as guard:
        result = _compute_shape_atomic_fit(_shape_atomic_snapshot(prepared), session.logger)
        guard.validate()
    return _finish_shape_atomic_fit(session, prepared, result)


def _shape_available_cpus():
    """Respect process CPU limits where supported, not just host hardware size."""
    import os

    count = getattr(os, "process_cpu_count", lambda: None)()
    if count is not None:
        return max(1, count)
    if hasattr(os, "sched_getaffinity"):
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except (OSError, NotImplementedError):
            pass
    return os.cpu_count() or 1


def _shape_worker_count(workers, models, local=False):
    from chimerax.core.errors import UserError

    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 0:
        raise UserError("workers must be a nonnegative integer: 0 automatic, 1 serial, or a process count.")
    cpus = _shape_available_cpus()
    if workers == 0:
        return min(max(1, cpus - 1), models) if models >= _SHAPE_AUTO_PARALLEL_MIN and not local else 1
    return min(workers, cpus, models)


def _shape_atomic_worker_main(directory):
    """Private subprocess entry point; exchange only arrays and JSON, never models."""
    import json
    from pathlib import Path
    import sys
    import traceback
    import numpy as np

    directory = Path(directory)
    config = json.loads((directory / "config.json").read_text())
    reference = np.load(directory / "reference.npy", mmap_mode="r", allow_pickle=False)
    reference_weights = (np.load(directory / "reference_weights.npy", mmap_mode="r", allow_pickle=False)
                         if config.get("reference_weighted") else None)

    class Logger:
        def __init__(self):
            self.warnings = []

        def status(self, message):
            pass

        def warning(self, message):
            self.warnings.append(str(message))

    def json_value(value):
        return value.tolist() if isinstance(value, np.ndarray) else value.item()

    with _ShapeBatchCache():
        for line in sys.stdin:
            job = json.loads(line)
            logger = Logger()
            try:
                data = dict(source_xyz=np.load(directory / f"source_{job['index']}.npy", allow_pickle=False),
                            reference_xyz=reference, options=config["options"],
                            source_weights=(np.load(directory / f"weights_{job['index']}.npy", allow_pickle=False)
                                            if job.get("source_weighted") else None),
                            reference_weights=reference_weights,
                            initial_matrix=job["initial_matrix"], source_regions=job["source_regions"],
                            reference_regions=config["reference_regions"])
                pose, stats = _compute_shape_atomic_fit(data, logger)
                result = dict(index=job["index"], pose=pose.matrix, stats=stats, warnings=logger.warnings)
            except Exception:
                result = dict(index=job["index"], error=traceback.format_exc())
            # A prefix keeps occasional native-library stdout messages out of
            # the protocol. Worker diagnostics go to a separate temporary file.
            print("SHAPEALIGN:" + json.dumps(result, default=json_value), flush=True)


def _shape_atomic_workers(session, prepared, workers):
    """Run a bounded pool; only pipe I/O is threaded in the parent process."""
    import json
    import os
    from pathlib import Path
    import queue
    from shutil import copyfile
    import subprocess
    import sys
    import tempfile
    import threading
    from chimerax.core.errors import UserError
    from chimerax.geometry import Place
    import numpy as np

    try:
        from chimerax.core.python_utils import chimerax_python_executable

        executable = chimerax_python_executable()
    except ImportError:
        executable = ""
    if not Path(executable).is_file():
        if Path(sys.executable).name.lower().startswith("python"):
            executable = sys.executable
        else:
            raise UserError("Cannot locate ChimeraX's Python executable; use workers 1.")
    script = Path(__file__).resolve()
    if not script.is_file():
        raise UserError("Cannot locate chimerax_trimmings.py for worker startup; use workers 1.")
    environment = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        environment[variable] = "1"
    environment.pop("PYTHONHOME", None)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    stop = threading.Event()
    waiting, completed = queue.Queue(), queue.Queue()
    processes, threads, diagnostics = [], [], []
    results = [None] * len(prepared)

    def exchange(process, diagnostic_path):
        try:
            while not stop.is_set():
                try:
                    job = waiting.get_nowait()
                except queue.Empty:
                    break
                process.stdin.write(json.dumps(job) + "\n")
                process.stdin.flush()
                while True:
                    line = process.stdout.readline()
                    if not line:
                        raise RuntimeError("A shapealign worker exited unexpectedly; try workers 1.")
                    position = line.find("SHAPEALIGN:")
                    if position >= 0:
                        result = json.loads(line[position + len("SHAPEALIGN:"):])
                        if result["index"] != job["index"]:
                            raise RuntimeError("A shapealign worker returned an unexpected job index.")
                        completed.put(result)
                        break
        except Exception as error:
            details = Path(diagnostic_path).read_text(errors="replace")[-4000:]
            completed.put(dict(error=str(error) + ("\n" + details if details else "")))

    def regions(values):
        return [(label, indices.tolist()) for label, indices in values]

    with tempfile.TemporaryDirectory(prefix="shapealign_") as temporary:
        directory = Path(temporary)
        try:
            worker_source = directory / "trimmings.py"
            copyfile(script, worker_source)
            np.save(directory / "reference.npy", prepared[0]["reference_xyz"], allow_pickle=False)
            reference_weights = prepared[0]["reference_weights"]
            if reference_weights is not None:
                np.save(directory / "reference_weights.npy", reference_weights, allow_pickle=False)
            config = dict(script=str(worker_source), sys_path=[os.path.abspath(p) for p in sys.path],
                          options=prepared[0]["options"], reference_regions=regions(prepared[0]["reference_regions"]),
                          reference_weighted=reference_weights is not None)
            (directory / "config.json").write_text(json.dumps(config))
            for index, item in enumerate(prepared):
                data = _shape_atomic_snapshot(item)
                np.save(directory / f"source_{index}.npy", data["source_xyz"], allow_pickle=False)
                if data["source_weights"] is not None:
                    np.save(directory / f"weights_{index}.npy", data["source_weights"], allow_pickle=False)
                waiting.put(dict(index=index, initial_matrix=data["initial_matrix"].tolist(),
                                 source_regions=regions(data["source_regions"]),
                                 source_weighted=data["source_weights"] is not None))
            bootstrap = directory / "worker.py"
            bootstrap.write_text(
                "import json, pathlib, runpy, sys\n"
                "root = pathlib.Path(sys.argv[1])\n"
                "config = json.loads((root / 'config.json').read_text())\n"
                "sys.path[:] = config['sys_path']\n"
                "namespace = runpy.run_path(config['script'])\n"
                "namespace['_shape_atomic_worker_main'](root)\n")
            for number in range(workers):
                diagnostic = (directory / f"worker_{number}.log").open("w+")
                diagnostics.append(diagnostic)
                process = subprocess.Popen(
                    [executable, "-u", str(bootstrap), str(directory)], stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE, stderr=diagnostic, text=True, encoding="utf-8",
                    env=environment, **({"creationflags": subprocess.CREATE_NO_WINDOW} if sys.platform == "win32" else {}))
                processes.append(process)
                thread = threading.Thread(target=exchange, args=(process, diagnostic.name), daemon=True)
                threads.append(thread)
                thread.start()
            count = 0
            while count < len(prepared):
                session.logger.status(f"Shape alignment: {count}/{len(prepared)} fits complete ({workers} workers)...")
                _shape_checkpoint()
                try:
                    result = completed.get(timeout=.1)
                except queue.Empty:
                    continue
                if "error" in result:
                    raise UserError("Shape alignment worker failed; no results were applied.\n" + result["error"])
                index = result["index"]
                pose = Place(matrix=np.array(result["pose"], dtype=float))
                if not np.isfinite(pose.matrix).all():
                    raise UserError("A shapealign worker returned a nonfinite transform; no results were applied.")
                results[index] = pose, result["stats"]
                for warning in result["warnings"]:
                    session.logger.warning(f"Shape alignment #{prepared[index]['model'].id_string}: {warning}")
                count += 1
        except OSError as error:
            raise UserError(f"Could not start shapealign workers; use workers 1. {error}") from error
        finally:
            stop.set()
            for process in processes:
                if process.poll() is None:
                    process.terminate()
            for process in processes:
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            for thread in threads:
                thread.join(timeout=5)
            for process in processes:
                process.stdin.close()
                process.stdout.close()
            for diagnostic in diagnostics:
                diagnostic.close()
    return results


def _validate_shape_atomic_snapshot(prepared):
    import numpy as np
    from chimerax.core.errors import UserError

    for prefix, model_key, atoms_key in (("source", "model", "source_atoms"),
                                         ("reference", "reference", "reference_atoms")):
        model, atoms = prepared[model_key], prepared[atoms_key]
        if (model.deleted or not np.array_equal(model.scene_position.matrix, prepared[prefix + "_scene"].matrix)
                or not np.isin(prepared[prefix + "_pointers"], model.atoms.pointers).all()
                or not np.array_equal(atoms.pointers, prepared[prefix + "_pointers"])
                or not np.array_equal(atoms.coords, prepared[prefix + "_xyz"])):
            raise UserError("A fitting model changed while workers were running; no results were applied.")
    for model, atoms, pointers, values in prepared.get("confidence_inputs", ()):
        # Check liveness before accessing a collection that may contain deleted atoms.
        if (model.deleted or not np.isin(pointers, model.atoms.pointers).all()
                or not np.array_equal(atoms.bfactors, values)):
            raise UserError("Fitting confidence scores changed while workers were running; no results were applied.")


def _shape_align_atomic_models(session, atoms, to, resolution=None, max_points=12000,
                extra_starts=0, move=True, radius=None, correlation=False,
                color_fit=False, correlation_resolution=8.0, correlation_window=12.0,
                partial=None, fallback_correlation=None, objective="density", workers=0,
                confidence="auto", confidence_moving=None, confidence_reference=None, copy=False):
    """Independently fit each moving model to one fixed reference by shape.

    Atom specifications restrict fitting regions; each whole moving model gets
    its own rigid transform. Search resolutions default to a size-derived range;
    resolution overrides its fine endpoint. radius enables local refinement
    in a sphere around the starting CoFR. correlation optionally saves local
    density CC attributes; color_fit also colors fitted atoms and cartoons.
    Scoring is off by default. A batch is one undo action and rolls back if any
    fit fails. Single fits return (transform, stats); batches return those pairs.
    """
    from contextlib import nullcontext
    from time import perf_counter
    from chimerax.core.errors import UserError

    targets = to.unique_structures
    moving = atoms.unique_structures
    if len(targets) != 1 or not len(moving):
        raise UserError("Specify moving atomic models and exactly one reference model.")
    reference = targets[0]
    if reference in moving:
        from chimerax.atomic import concatenate

        items = _shape_exclude_reference(session, list(atoms.by_structure), reference)
        if not items:
            return []
        atoms = concatenate([selected for _, selected in items])
        moving = atoms.unique_structures
    worker_count = _shape_worker_count(workers, len(moving), radius is not None)
    for model in moving:
        if reference in model.drawing_lineage or model in reference.drawing_lineage:
            raise UserError("The moving and reference models must not contain one another.")
        if any(other in model.drawing_lineage[:-1] for other in moving):
            raise UserError("Moving models must not contain one another; specify independent models.")
    options = dict(resolution=resolution, max_points=max_points, extra_starts=extra_starts,
                   move=move, radius=radius, correlation=correlation, color_fit=color_fit,
                   correlation_resolution=correlation_resolution, correlation_window=correlation_window,
                   partial=partial, fallback_correlation=fallback_correlation, objective=objective,
                   confidence=confidence, confidence_moving=confidence_moving, confidence_reference=confidence_reference,
                   copy=copy)
    if len(moving) == 1:
        return _shape_align_one(session, atoms, to, **options)
    if radius is not None:
        # Read once: automatic CoFR modes must not change the sphere mid-batch.
        options["sphere_center"] = session.main_view.center_of_rotation.copy()
    ordered = sorted(atoms.by_structure, key=lambda item: item[0].id)
    prepared = computed = None
    started = perf_counter()
    if worker_count > 1:
        prepared = [_prepare_shape_atomic_fit(session, selected, to, **options) for _, selected in ordered]
        session.logger.info(f"Shape alignment: fitting {len(prepared)} atomic models with {worker_count} workers.")
        with _ShapeMapInputGuard((*atoms.by_structure, *to.by_structure)) as guard:
            computed = _shape_atomic_workers(session, prepared, worker_count)
            guard.validate()
        for item in prepared:
            _validate_shape_atomic_snapshot(item)
    results = []
    context = session.undo.aggregate("shape align models") if move else nullcontext()
    with _ShapeBatchCache(), context as batch:
        try:
            for index, (_, selected) in enumerate(ordered):
                if prepared is None:
                    results.append(_shape_align_one(session, selected, to, **options))
                else:
                    computed[index][1]["workers"] = worker_count
                    results.append(_finish_shape_atomic_fit(session, prepared[index], computed[index]))
        except BaseException:
            if batch is not None:
                for action in reversed(batch.actions):
                    action.undo()
                if batch.actions:
                    session.logger.warning("Shape alignment batch failed; earlier changes were undone.")
            raise
    if worker_count > 1:
        message = (f"Shape alignment: {len(results)} atomic models fitted with {worker_count} workers "
                   f"in {perf_counter() - started:.2f} s; " + ("one Undo action." if move else "not moved."))
        session.logger.info(message)
        session.logger.status(message)
    return results


def _shape_exclude_reference(session, items, reference):
    remaining = [(model, selected) for model, selected in items if model is not reference]
    if len(remaining) != len(items):
        message = f"Shape alignment: skipping reference #{reference.id_string} in the moving selection."
        if not remaining:
            message += " No other moving models remain."
            session.logger.status(message)
        session.logger.info(message)
    return remaining


def _shape_input_items(selection):
    """Retain atom subsets and map owners from a mixed command specification."""
    from chimerax.atomic import Atoms, Structure
    from chimerax.core.errors import UserError
    from chimerax.map import Volume

    if isinstance(selection, Volume):
        return [(selection, None)]
    if isinstance(selection, Structure):
        return [(selection, selection.atoms)]
    if isinstance(selection, Atoms):
        return sorted(selection.by_structure, key=lambda item: item[0].id)
    atoms = selection.atoms
    items = dict(atoms.by_structure)
    models = selection.models
    for model in models:
        if isinstance(model, Volume):
            items[model] = None
    for model in models:
        if not isinstance(model, (Structure, Volume)):
            # A selected group or a map's surface is not an independent input.
            if not any(model in owner.drawing_lineage or owner in model.drawing_lineage
                       for owner in items):
                raise UserError("shapealign accepts atomic models and density maps, not other drawings.")
    return sorted(items.items(), key=lambda item: item[0].id)


def _shape_grid_corners(grid):
    from itertools import product
    import numpy as np

    return grid.ijk_to_xyz_transform.transform_points(np.array(list(product(
        *[(0, size - 1) for size in grid.size])), dtype=np.float64))


def _resample_shape_density(grid, origin, size, step, pose=None):
    """Interpolate in bounded chunks; pose maps source XYZ to destination XYZ."""
    import numpy as np
    from chimerax.map_fit.fitmap import volume_values

    matrix = np.empty(tuple(reversed(size)), dtype=np.float32)
    transform = grid.xyz_to_ijk_transform
    if pose is not None:
        transform = transform * pose.inverse()
    nx, ny, _ = size
    for start in range(0, matrix.size, 65536):
        _shape_checkpoint()
        indices = np.arange(start, min(start + 65536, matrix.size))
        k, rem = np.divmod(indices, nx * ny)
        j, i = np.divmod(rem, nx)
        points = np.column_stack((i, j, k)).astype(np.float32)
        points *= step
        points += origin
        matrix.ravel()[start:start + len(indices)] = volume_values(points, transform, grid.matrix())
    return matrix


def _shape_cofr_crop_bounds(data, center, radius):
    """Native-grid box enclosing a physical CoFR sphere, clipped to the map."""
    from itertools import product
    import numpy as np

    corners = np.array(list(product(*[(v - radius, v + radius) for v in center])))
    ijk = data.xyz_to_ijk_transform.transform_points(corners)
    lower = np.maximum(np.floor(ijk.min(axis=0)).astype(int), 0)
    upper = np.minimum(np.ceil(ijk.max(axis=0)).astype(int) + 1, data.size)
    return lower, upper


def _shape_cofr_radius_limit(data, center, rejected_radius):
    """Translate the crop allocation limit into an actionable physical radius."""
    from math import floor, log10
    import numpy as np

    low, high = 0., rejected_radius
    for _ in range(48):
        middle = (low + high) / 2
        lower, upper = _shape_cofr_crop_bounds(data, center, middle)
        if np.prod(np.maximum(0, upper - lower)) <= _MAX_GRID_VOXELS:
            low = middle
        else:
            high = middle
    if low <= 0:
        return 0.
    # Round down, not to nearest: the suggested number must pass the same guard.
    quantum = 10 ** (floor(log10(low)) - 1)
    return floor(low / quantum) * quantum


class _ShapeDensityInput:
    """A selected atomic density or native full map, in the owner's local frame."""

    def __init__(self, model, atoms, center=None, radius=None, confidence="auto", allow_empty=False):
        import numpy as np
        from chimerax.atomic import Residue
        from chimerax.core.errors import UserError
        from chimerax.map_data import ArrayGridData

        self.model, self.atoms = model, atoms
        self.is_map = atoms is None
        self.weights = None
        self._cache = {}
        self._grid_cache_key = object()
        self._shared_reference = False
        if atoms is not None:
            keep = (atoms.element_numbers > 1) & ~np.isin(atoms.residues.names, list(Residue.water_res_names))
            self.atoms = atoms.filter(keep)
            xyz = self.atoms.coords
            if not np.isfinite(xyz).all():
                raise UserError("Each fitting region needs finite atom coordinates.")
            if center is not None:
                inside = np.linalg.norm(model.scene_position.transform_points(xyz) - center, axis=1) <= radius
                xyz, self.atoms = xyz[inside], self.atoms.filter(inside)
            if len(xyz) < 4 or np.linalg.matrix_rank(xyz - xyz.mean(axis=0)) < 2:
                raise UserError("Each atomic fitting region needs four non-collinear non-water heavy atoms.")
            self.atoms, xyz, self.weights = _shape_weighted_atoms(self.atoms, xyz, confidence)
            self.xyz = xyz
            self.center = _shape_center(xyz, self.weights)
            self.landmarks = xyz - self.center
            self.frame_landmarks = _shape_geometry(self.landmarks, self.weights)
            return

        data = model.data
        if len(data.size) != 3 or min(data.size) < 2:
            raise UserError("Each map must contain a three-dimensional density grid.")
        if center is not None:
            local_center = model.scene_position.inverse() * center
            lower, upper = _shape_cofr_crop_bounds(data, local_center, radius)
            if np.any(upper <= lower):
                raise UserError(f"The {radius:.3g} A CoFR sphere does not intersect the map grid; "
                                "move CoFR into the desired map region or increase the radius.")
            voxels = int(np.prod(upper - lower))
            if voxels > _MAX_GRID_VOXELS:
                limit = _shape_cofr_radius_limit(data, local_center, radius)
                example = min(50., limit)
                raise UserError(f"The CoFR radius is {radius:.3g} A; the input-crop limit here is about {limit:g} A. "
                                f"Set a smaller explicit radius, e.g. radius_search {example:g} "
                                f"(or radius_refine {example:g} for local refinement). "
                                "A bare radius_search chooses the radius automatically.")
            matrix = data.matrix(tuple(lower), tuple(upper - lower)).copy()
            if matrix.ndim != 3:
                raise UserError("shapealign requires scalar three-dimensional density, not color/vector maps.")
            origin = data.ijk_to_xyz_transform * lower
            data = ArrayGridData(matrix, origin, data.step, data.cell_angles, data.rotation)
            for k, plane in enumerate(matrix):
                j, i = np.indices(plane.shape)
                points = np.column_stack((i.ravel(), j.ravel(), np.full(plane.size, k))).astype(float)
                points = data.ijk_to_xyz_transform.transform_points(points)
                plane[np.linalg.norm(points - local_center, axis=1).reshape(plane.shape) > radius] = 0
        else:
            matrix = data.full_matrix()
        if matrix.ndim != 3:
            raise UserError("shapealign requires scalar three-dimensional density, not color/vector maps.")
        self.data = data
        self.matrix = matrix

        # All voxel values, including negative density, contribute through rho^2
        # to geometry. No displayed contour, surface mask, or display region is read.
        peak = 0.0
        for plane in matrix:
            values = plane.astype(np.float64)
            if not np.isfinite(values).all():
                raise UserError("Fitting maps must contain only finite density values.")
            peak = max(peak, float(np.max(np.abs(values))))
        if peak == 0:
            if allow_empty:
                self.empty = True
                return
            raise UserError("The map fitting region has no nonzero density.")
        self.peak = peak
        mass, first, second = 0.0, np.zeros(3), np.zeros((3, 3))
        lower, upper = np.array(matrix.shape[::-1]), np.zeros(3, dtype=int)
        i, j = np.arange(matrix.shape[2], dtype=float), np.arange(matrix.shape[1], dtype=float)
        for k, plane in enumerate(matrix):
            weights = (plane.astype(np.float64) / peak) ** 2
            total = weights.sum()
            if total == 0:
                continue
            wi, wj = weights.sum(axis=0), weights.sum(axis=1)
            si, sj = np.dot(wi, i), np.dot(wj, j)
            mass += total
            first += (si, sj, k * total)
            second[0, 0] += np.dot(wi, i * i)
            second[1, 1] += np.dot(wj, j * j)
            second[2, 2] += k * k * total
            second[0, 1] += np.dot(j, weights @ i)
            second[0, 2] += k * si
            second[1, 2] += k * sj
            ii, jj = np.flatnonzero(wi), np.flatnonzero(wj)
            lower = np.minimum(lower, (ii[0], jj[0], k))
            upper = np.maximum(upper, (ii[-1], jj[-1], k))
        second += np.triu(second, 1).T
        mean = first / mass
        affine = data.ijk_to_xyz_transform.matrix[:, :3]
        covariance = affine @ (second / mass - np.outer(mean, mean)) @ affine.T
        values, axes = np.linalg.eigh(covariance)
        if values[-2] <= max(values[-1], 1) * 1e-12:
            if allow_empty:
                self.empty = True
                return
            raise UserError("The map density is collinear and cannot define a rigid orientation.")
        self.empty = False
        self.center = data.ijk_to_xyz_transform * mean
        offsets = axes * np.sqrt(3 * np.maximum(values, 0))
        self.landmarks = np.concatenate((offsets.T, -offsets.T))
        self.frame_landmarks = self.landmarks
        self.lower, self.upper = lower, upper

    def grid(self, resolution):
        import numpy as np
        from chimerax.core.errors import UserError
        from chimerax.map_data import ArrayGridData

        cache = _SHAPE_BATCH_CACHE.get() if self._shared_reference else None
        cache_key = ("density_grid", self._grid_cache_key, resolution)
        if cache is not None:
            grid = cache.get(cache_key)
            if grid is not None:
                return grid
        elif resolution in self._cache:
            return self._cache[resolution]
        if not self.is_map:
            grid = _shape_grid(*_shape_points(self.landmarks, resolution / 6, self.weights), resolution)
        else:
            # Cell-average before decimation, retaining boundary bins through
            # zero extension. This avoids point-sampling a fine map at coarse steps.
            affine = self.data.ijk_to_xyz_transform.matrix[:, :3]
            factors = np.maximum(1, np.floor((resolution / 3) / np.linalg.norm(affine, axis=0))).astype(int)
            shape = np.ceil((self.upper - self.lower + 1) / factors).astype(int)
            sigma = resolution / (pi * sqrt(2))
            base_affine = affine @ np.diag(factors)
            inverse = np.linalg.inv(base_affine)
            halo = np.ceil(4 * sigma * np.linalg.norm(inverse, axis=1)).astype(int) + 1
            padded_size = shape + 2 * halo
            if np.prod(padded_size) > _MAX_GRID_VOXELS:
                raise UserError("The map fitting grid is too large; use a coarser resolution or smaller CoFR region.")
            bounds = self.data.ijk_to_xyz_transform.transform_points(np.array([
                (i, j, k) for i in (self.lower[0], self.upper[0])
                for j in (self.lower[1], self.upper[1]) for k in (self.lower[2], self.upper[2])], dtype=float))
            step = resolution / 3
            origin = bounds.min(axis=0) - self.center - 4 * sigma
            size = np.ceil((bounds.max(axis=0) - self.center + 4 * sigma - origin) / step).astype(int) + 1
            if np.prod(size) > _MAX_GRID_VOXELS:
                raise UserError("The Cartesian map fitting grid is too large; use a coarser resolution.")
            work = self.matrix[tuple(slice(a, b + 1) for a, b in zip(self.lower[::-1], self.upper[::-1]))]
            for axis in np.argsort(factors[::-1])[::-1]:
                factor = int(factors[2 - axis])
                work = np.add.reduceat(work, np.arange(0, work.shape[axis], factor), axis=axis, dtype=np.float64)
                work /= factor
            work /= self.peak
            work = np.pad(work.astype(np.float32), tuple((int(v), int(v)) for v in halo[::-1]))
            # Physical isotropic Gaussian filtering also handles non-orthogonal
            # cells. Padding prevents periodic wraparound in this temporary FFT.
            if getattr(self, "_finite_gaussian", not getattr(self, "_normalize_grid", True)):
                work = _shape_region_gaussian(work, base_affine, sigma)
            else:
                from scipy.fft import rfftn, irfftn, fftfreq, rfftfreq
                fz = fftfreq(work.shape[0])[:, None, None]
                fy = fftfreq(work.shape[1])[None, :, None]
                fx = rfftfreq(work.shape[2])[None, None, :]
                covariance = sigma * sigma * inverse @ inverse.T
                q = (covariance[0, 0] * fx * fx + covariance[1, 1] * fy * fy + covariance[2, 2] * fz * fz +
                     2 * covariance[0, 1] * fx * fy + 2 * covariance[0, 2] * fx * fz + 2 * covariance[1, 2] * fy * fz)
                fft_threads = _shape_map_fft_threads(work.size)
                spectrum = rfftn(work, workers=fft_threads)
                spectrum *= np.exp(-2 * pi * pi * q)
                work = irfftn(spectrum, s=work.shape, workers=fft_threads).astype(np.float32)
            base_origin = self.data.ijk_to_xyz_transform * (self.lower + (factors - 1) / 2 - factors * halo)
            filtered = ArrayGridData(work, base_origin - self.center,
                                     np.asarray(self.data.step) * factors,
                                     self.data.cell_angles, self.data.rotation)
            matrix = _resample_shape_density(filtered, origin, size, step)
            grid = ArrayGridData(matrix, origin, (step,) * 3)
        # Positive scaling does not change either objective and stabilizes fits
        # between experimental map units and simulated atomic densities.
        peak = np.max(np.abs(grid.matrix()))
        if not np.isfinite(peak) or peak == 0:
            raise UserError("No finite nonzero density remains at this fitting scale.")
        if getattr(self, "_normalize_grid", True):
            grid.matrix()[:] /= peak
        if cache is not None:
            cache.put(cache_key, grid)
        else:
            self._cache.clear()  # Moving inputs and single fits retain only one level.
            self._cache[resolution] = grid
        return grid


def _shape_reference_input(model, atoms, center, radius, confidence="auto"):
    def prepare():
        value = _ShapeDensityInput(model, atoms, center, radius, confidence)
        value._shared_reference = True
        return value

    if _SHAPE_BATCH_CACHE.get() is None:
        return prepare()
    selection = None if atoms is None else (_shape_array_key(atoms.pointers), _shape_array_key(atoms.coords))
    confidence_key = None if atoms is None or confidence == "off" else (
        confidence, _shape_has_plddt(model) if confidence == "auto" else True, _shape_array_key(atoms.bfactors))
    key = ("reference_input", id(model), selection,
           None if center is None else tuple(center), radius, model.scene_position.matrix.tobytes(),
           confidence_key)
    return _shape_cached(key, prepare)


def _shape_density_samples(grid, maximum):
    """Signed, importance-weighted samples of the complete density integral."""
    import numpy as np

    matrix = grid.matrix()
    if matrix.size <= maximum:
        indices = np.flatnonzero(matrix.ravel())
        weights = matrix.ravel()[indices].copy()
    else:
        masses = np.array([np.abs(plane).sum(dtype=np.float64) for plane in matrix])
        cumulative = np.cumsum(masses)
        total = cumulative[-1]
        levels = (np.arange(maximum) + np.random.default_rng(0).random()) * total / maximum
        planes = np.searchsorted(cumulative, levels, side="right")
        indices = np.empty(maximum, dtype=np.int64)
        for k in np.unique(planes):
            chosen = np.flatnonzero(planes == k)
            local = levels[chosen] - (cumulative[k - 1] if k else 0)
            cdf = np.cumsum(np.abs(matrix[k]).ravel(), dtype=np.float64)
            indices[chosen] = int(k) * matrix.shape[1] * matrix.shape[2] + np.searchsorted(cdf, local, side="right")
        weights = np.sign(matrix.ravel()[indices]) * (total / maximum)
    k, j, i = np.unravel_index(indices, matrix.shape)
    points = np.column_stack((i, j, k)).astype(np.float32)
    grid.ijk_to_xyz_transform.transform_points(points, in_place=True)
    return points, weights.astype(np.float32)


def _refine_density_grids(source, target, candidates, maximum, objective):
    import numpy as np

    points, weights = _shape_density_samples(source, maximum)
    field, second = _shape_reference_field(target, objective)
    if objective == "gradient":
        _, first = _gradient_density_field(source)
    else:
        first = np.einsum("ijk,ijk->", source.matrix(), source.matrix(), dtype=np.float64)
    denominator = sqrt(first * second)
    workers = _SHAPE_MAP_WORKERS.get()
    if workers is not None:
        results = workers.refine(points, weights, field, target.xyz_to_ijk_transform, candidates, denominator)
        if results is not None:
            return results, points
    results = []
    for index, initial in enumerate(candidates):
        _shape_checkpoint()
        score, pose = _shape_refine_grid_pose(points, weights, field, target.xyz_to_ijk_transform,
                                             initial, denominator)
        results.append((score, pose, index))
    return sorted(results, key=lambda item: item[0], reverse=True), points


def _shape_map_fft_starts(source, target, source_landmarks, target_landmarks, objective, extra_starts,
                          valid_pose=None, site_distance=None, peaks_per_orientation=2, _site_center=None):
    import numpy as np
    from scipy.fft import rfftn, irfftn, next_fast_len
    from chimerax.core.errors import UserError
    from chimerax.geometry import Place

    step = source.step[0]
    radius = np.linalg.norm(_shape_grid_corners(source), axis=1).max()
    target_radius = np.linalg.norm(_shape_grid_corners(target), axis=1).max()
    if target_radius < radius:
        # This only chooses the cheaper numerical template. Full-density norms
        # stay global, and inversion preserves which owner the user asked to move.
        center = _site_center or (lambda pose: pose.origin())
        return [pose.inverse() for pose in _shape_map_fft_starts(
            target, source, target_landmarks, source_landmarks, objective, extra_starts,
            None if valid_pose is None else lambda pose: valid_pose(pose.inverse()),
            site_distance, peaks_per_orientation, lambda pose: center(pose.inverse()))]
    side = int(np.ceil(2 * radius / step)) + 1
    side += 1 - side % 2
    fft_shape = tuple(next_fast_len(n + side - 1) for n in target.matrix().shape)
    if np.prod(fft_shape) > _MAX_PARTIAL_FFT_VOXELS:
        raise UserError("The map FFT search exceeds its memory limit; use a coarser resolution or smaller region.")
    half = side // 2
    origin = np.full(3, -half * step)
    field, second = _shape_reference_field(target, objective)
    if objective == "gradient":
        _, first = _gradient_density_field(source)
    else:
        first = np.sum(source.matrix().astype(float) ** 2)
    fft_threads = _shape_map_fft_threads(int(np.prod(fft_shape)))
    target_fft = _shape_reference_fft(target, field, fft_shape, objective, workers=fft_threads)
    peaks = []
    for initial in _shape_starts(source_landmarks, target_landmarks, 96 + extra_starts):
        _shape_checkpoint()
        rotation = Place(matrix=np.column_stack((initial.matrix[:, :3], np.zeros(3))))
        density = _resample_shape_density(source, origin, (side,) * 3, step, rotation)
        cross = irfftn(target_fft * rfftn(density[::-1, ::-1, ::-1], fft_shape, workers=fft_threads),
                        fft_shape, workers=fft_threads) / sqrt(first * second)
        # argmax otherwise copies this strided FFT slice for every peak.
        view = cross[half:half + field.shape[0], half:half + field.shape[1], half:half + field.shape[2]].copy()
        accepted = 0
        for _ in range(peaks_per_orientation if valid_pose is None else view.size):
            k, j, i = np.unravel_index(np.argmax(view), view.shape)
            value = float(view[k, j, i])
            if not np.isfinite(value):
                break
            center = target.ijk_to_xyz_transform * np.array((i, j, k), dtype=float)
            pose = Place(matrix=np.column_stack((rotation.matrix[:, :3], center)))
            if valid_pose is None or valid_pose(pose):
                peaks.append((value, pose))
                accepted += 1
            distance = max(2, int(max(radius / 3, site_distance or 0) / step))
            view[max(0, k-distance):k+distance+1, max(0, j-distance):j+distance+1,
                 max(0, i-distance):i+distance+1] = -np.inf
            if accepted == peaks_per_orientation:
                break
    peaks.sort(key=lambda item: item[0], reverse=True)
    kept, sites = [], []
    probes = _shape_grid_corners(source)
    for _, pose in peaks:
        if any(np.sqrt(np.mean(np.sum((pose.transform_points(probes) - old.transform_points(probes)) ** 2, axis=1)))
               < step for old in kept):
            continue
        center = (_site_center(pose) if _site_center is not None else pose.origin())
        # Reserve room for weaker spatial basins before fine-detail scoring.
        # A bright site must not occupy the entire tiled-search orientation beam.
        if site_distance is not None and sum(np.linalg.norm(center - old) < site_distance for old in sites) >= 4:
            continue
        kept.append(pose)
        sites.append(center)
        if len(kept) == 48:
            break
    return kept


def _shape_region_gaussian(matrix, affine, sigma):
    """Finite physical Gaussian: constant interiors stay constant even below a voxel.

    A truncated Fourier Gaussian has long ringing tails when its transfer at
    Nyquist is nonzero. Those tails cannot be removed with a small spatial halo.
    Use a normalized spatial kernel for local bands instead.
    """
    import numpy as np
    from scipy.ndimage import gaussian_filter
    from chimerax.core.errors import UserError

    metric = affine.T @ affine
    if np.allclose(metric, np.diag(np.diag(metric)), atol=float(np.max(np.diag(metric))) * 1e-8):
        gaussian_filter(matrix, (sigma / np.sqrt(np.diag(metric)))[::-1],
                         output=matrix, mode="constant", truncate=4)
        return matrix
    from scipy.signal import fftconvolve
    from scipy.fft import next_fast_len

    half = np.ceil(4 * sigma * np.linalg.norm(np.linalg.inv(affine), axis=1)).astype(int)
    kernel_shape = tuple(2 * half[::-1] + 1)
    fft_shape = tuple(next_fast_len(n + k - 1) for n, k in zip(matrix.shape, kernel_shape))
    if np.prod(fft_shape) > _MAX_GRID_VOXELS:
        raise UserError("The skew-map filtering grid is too large; restrict fitting atoms or use coarser settings.")
    k, j, i = np.indices(kernel_shape, dtype=np.float64)
    xyz = np.stack((i - half[0], j - half[1], k - half[2]), axis=-1) @ affine.T
    squared = np.sum(xyz * xyz, axis=-1)
    kernel = np.exp(-squared / (2 * sigma * sigma))
    kernel[squared > (4 * sigma) ** 2] = 0
    kernel /= kernel.sum()
    return fftconvolve(matrix, kernel.astype(np.float32), mode="same").astype(np.float32)


def _shape_map_box(target, points, pad, resolution=None):
    """Numerical rectangular view with the original map frame and intensity scale.

    Padding is physical (Angstroms), so rotated/skewed map grids work too.
    This never creates a Volume or modifies the input matrix or surface.
    """
    from itertools import product
    import numpy as np
    from chimerax.map_data import ArrayGridData

    lower = np.min(points, axis=0) + target.center - pad
    upper = np.max(points, axis=0) + target.center + pad
    corners = np.array(list(product(*zip(lower, upper))), dtype=float)
    ijk = target.data.xyz_to_ijk_transform.transform_points(corners)
    lo = np.maximum(0, np.floor(ijk.min(axis=0)).astype(int))
    hi = np.minimum(target.data.size, np.ceil(ijk.max(axis=0)).astype(int) + 1)
    if resolution is not None:
        affine = target.data.ijk_to_xyz_transform.matrix[:, :3]
        factors = np.maximum(1, np.floor((resolution / 3) / np.linalg.norm(affine, axis=0))).astype(int)
        # Anchor decimation bins to the original grid, not the tile boundary.
        lo = (lo // factors) * factors
        hi = np.minimum(target.data.size, ((hi + factors - 1) // factors) * factors)
    if np.any(hi <= lo):
        return None
    matrix = target.matrix[tuple(slice(a, b) for a, b in zip(lo[::-1], hi[::-1]))]
    if not np.any(matrix):
        return None
    data = ArrayGridData(matrix, target.data.ijk_to_xyz_transform * lo,
                         target.data.step, target.data.cell_angles, target.data.rotation)
    # Bypass moment estimation: a working box must not recenter each candidate
    # or renormalize its density independently of neighbouring boxes.
    region = _ShapeDensityInput.__new__(_ShapeDensityInput)
    region.model, region.atoms, region.is_map, region.weights = target.model, None, True, None
    region.data, region.matrix, region.center, region.peak = data, matrix, target.center, target.peak
    region.lower, region.upper = np.zeros(3, dtype=int), np.array(data.size) - 1
    region._cache, region._grid_cache_key = {}, object()
    region._shared_reference, region._normalize_grid = False, False
    return region


def _shape_map_contains(target, points):
    import numpy as np

    # An explicit CoFR sphere is a scoring selection, not the physical map edge.
    data = target.model.data
    ijk = data.xyz_to_ijk_transform.transform_points(points + target.center)
    return np.all((ijk >= 0) & (ijk <= np.array(data.size) - 1), axis=1)


def _shape_map_grid_voxels(region, resolution):
    """Preflight the original map filter and Cartesian grid without allocating."""
    from itertools import product
    import numpy as np

    affine = region.data.ijk_to_xyz_transform.matrix[:, :3]
    factors = np.maximum(1, np.floor((resolution / 3) / np.linalg.norm(affine, axis=0))).astype(int)
    shape = np.ceil((region.upper - region.lower + 1) / factors).astype(int)
    sigma = resolution / (pi * sqrt(2))
    inverse = np.linalg.inv(affine @ np.diag(factors))
    halo = np.ceil(4 * sigma * np.linalg.norm(inverse, axis=1)).astype(int) + 1
    corners = np.array(list(product(*zip(region.lower, region.upper))), dtype=float)
    corners = region.data.ijk_to_xyz_transform.transform_points(corners)
    cartesian = np.ceil((np.ptp(corners, axis=0) + 8 * sigma) / (resolution / 3)).astype(int) + 1
    return max(int(np.prod(shape + 2 * halo)), int(np.prod(cartesian)))


def _shape_map_fft_voxels(source, target):
    """The allocation used by the original global FFT search, including reversal."""
    import numpy as np
    from scipy.fft import next_fast_len

    radius = float(np.linalg.norm(_shape_grid_corners(source), axis=1).max())
    target_radius = float(np.linalg.norm(_shape_grid_corners(target), axis=1).max())
    if target_radius < radius:
        source, target, radius = target, source, target_radius
    side = int(np.ceil(2 * radius / source.step[0])) + 1
    side += 1 - side % 2
    return int(np.prod([next_fast_len(n + side - 1) for n in target.matrix().shape]))


def _shape_submap_input(target, lower, upper, padding):
    """A numerical submap with its own density geometry, using the original engine."""
    from types import SimpleNamespace

    box = _shape_map_box(target, (lower, upper), padding)
    if box is None or min(box.data.size) < 2:
        return None
    region = _ShapeDensityInput(SimpleNamespace(data=box.data), None, allow_empty=True)
    if region.empty:
        return None
    # All coordinates remain in the original map's local coordinate system.
    # No Volume is created and the original owner/scene transform is retained.
    region.model = target.model
    return region


def _shape_map_tiles(source, target, levels, tile_size=None):
    """Partition the map into overlapping submaps sized from the fitted body.

    The confidence-weighted ellipsoid plans the regions; it does not remove
    fitting atoms. The original global search still fits the complete selection.
    Submaps share a fitting pyramid and are shrunk if either original allocation
    guard would be exceeded. Only numerical views/grids are retained.
    """
    from itertools import product
    import numpy as np
    from chimerax.core.errors import UserError

    bounds = _shape_grid_corners(target.data) - target.center
    lower, upper = bounds.min(axis=0), bounds.max(axis=0)
    geometry = source.frame_landmarks
    covariance = geometry.T @ geometry / len(geometry)
    padding = max(float(np.sqrt(5 * np.linalg.eigvalsh(covariance)[-1])), levels[0])
    core = float(tile_size) if tile_size is not None else 2.5 * _shape_spatial_size(geometry)
    source_grid = source.grid(levels[0])
    for _ in range(32):
        counts = np.maximum(1, np.ceil((upper - lower) / core).astype(int))
        regions, largest_fft, fits = [], 0, True
        for index in product(*(range(int(n)) for n in counts)):
            lo = lower + np.array(index) * core
            hi = np.minimum(lo + core, upper)
            region = _shape_submap_input(target, lo, hi, padding)
            if region is None:
                continue
            if any(_shape_map_grid_voxels(region, level) > _MAX_GRID_VOXELS for level in levels):
                fits = False
                break
            voxels = _shape_map_fft_voxels(source_grid, region.grid(levels[0]))
            largest_fft = max(largest_fft, voxels)
            if voxels > _MAX_PARTIAL_FFT_VOXELS:
                fits = False
                break
            regions.append((region, lo, hi))
        if fits:
            return regions, dict(tiles=int(np.prod(counts)), nonempty_tiles=len(regions),
                                 tile_core_angstroms=core, tile_padding_angstroms=padding,
                                 tile_fft_voxels=largest_fft)
        for region, _, _ in regions:
            region._cache.clear()
        core *= .8
        if core < levels[0] / 3:
            break
    raise UserError("The moving model/submap search exceeds the working-grid limits; "
                    "restrict fitting atoms or use a coarser resolution.")


def _shape_keep_map_poses(scored, probes, distance, maximum=_SHAPE_MAP_BEAM):
    """Keep distinct rigid placements by deterministic RMS-displacement checks."""
    import numpy as np

    kept, placed = [], []
    for entry in sorted(scored, key=lambda item: item[0], reverse=True):
        if not np.isfinite(entry[0]):
            continue
        pose = entry[1]
        trial = pose.transform_points(probes)
        if any(np.sqrt(np.mean(np.sum((trial - old) ** 2, axis=1))) < distance for old in placed):
            continue
        kept.append(entry)
        placed.append(trial)
        if len(kept) == maximum:
            break
    return kept


class _ShapeBandScorer:
    """Smooth Gaussian-difference band, measured inside a fixed model-local mask.

    Both signs of the filtered density are retained. A deterministic importance
    sample integrates the same soft mask for every candidate, independent of
    tile size, background, map threshold, or candidate-specific peak values.
    """

    def __init__(self, source, coarse, fine, maximum):
        import numpy as np
        from chimerax.map.molmap import add_gaussians, bounding_grid
        from chimerax.core.errors import UserError

        step, pad = fine / 3, 4 * coarse / (pi * sqrt(2))
        dimensions = np.ceil((np.ptp(source.landmarks, axis=0) + 2 * pad) / step)
        if np.prod(dimensions) > _MAX_GRID_VOXELS:
            raise UserError("The model-local scoring grid is too large; use a coarser map_resolution "
                            "or restrict fitting atoms.")
        grid = bounding_grid(source.landmarks, step, pad)
        xyz, weights = _shape_points(source.landmarks, fine / 6, source.weights)
        values = []
        for blur in (fine, coarse, max(4, fine)):
            grid.matrix()[:] = 0
            add_gaussians(grid, xyz, weights, blur / (pi * sqrt(2)), 4)
            values.append(grid.matrix().copy())
        mask = values.pop()
        positive = mask[mask > float(mask.max()) * .001]
        mask /= max(float(np.median(positive)), 1e-30)
        np.minimum(mask, 30, out=mask)
        np.negative(mask, out=mask)
        np.expm1(mask, out=mask)
        mask *= -1
        mask[mask < .001] = 0
        cumulative = np.cumsum(mask.ravel(), dtype=np.float64)
        indices = np.searchsorted(cumulative, (np.arange(maximum) + .5) * cumulative[-1] / maximum)
        k, j, i = np.unravel_index(indices, mask.shape)
        self.points = grid.ijk_to_xyz_transform.transform_points(np.column_stack((i, j, k)).astype(np.float32))
        self.values = (values[0].ravel()[indices] - values[1].ravel()[indices]).astype(np.float64)
        self.values -= self.values.mean()
        self.norm = float(np.dot(self.values, self.values))
        self.coarse, self.fine = coarse, fine

    def score(self, target, pose):
        import numpy as np
        from chimerax.map_fit.fitmap import volume_values

        placed = pose.transform_points(self.points)
        # Missing map voxels must not improve a fit by shrinking its mask.
        if np.mean(_shape_map_contains(target, placed)) < .99:
            return -np.inf
        values = []
        for blur in (self.fine, self.coarse):
            # Six sigmas keep artificial crop-edge tails below float32 noise.
            region = _shape_map_box(target, placed, 6 * self.coarse / (pi * sqrt(2)), blur)
            if region is None:
                return -np.inf
            grid = region.grid(blur)
            values.append(volume_values(placed, grid.xyz_to_ijk_transform, grid.matrix()).astype(np.float64))
            del grid, region
        band = values[0] - values[1]
        band -= band.mean()
        variance = float(np.dot(band, band))
        if variance <= max(1e-20, float(np.dot(values[0], values[0])) * 1e-10):
            return -np.inf  # Reject constant/undefined bands, not float32 FFT roundoff.
        denominator = sqrt(self.norm * variance)
        return float(np.clip(np.dot(self.values, band) / denominator, -1, 1)) if denominator > 1e-20 else -np.inf


class _ShapeEmptyMapSearch(Exception):
    """A submap has no physically valid placement; other submaps may still fit."""


def _shape_fit_levels(source, target, options):
    final = options["resolution"]
    if final is None and options.get("map_resolution") is not None:
        automatic, _ = _shape_resolution_schedule(source.frame_landmarks, target.frame_landmarks)
        final = max(automatic[-1], options["map_resolution"])
    levels, stats = _shape_resolution_schedule(source.frame_landmarks, target.frame_landmarks,
                                              final, options["radius"] is not None)
    stats["resolution_mode"] = "automatic" if options["resolution"] is None else "explicit"
    return levels, stats


def _shape_refine_regional(source, target, pose, blur, maximum, objective):
    """Use the original native sum-product refiner in a fresh padded map box."""
    grid = source.grid(blur)
    points, _ = _shape_density_samples(grid, maximum)
    region = _shape_map_box(target, pose.transform_points(points),
                            6 * blur / (pi * sqrt(2)) + 2 * blur, blur)
    if region is None:
        return float("-inf"), pose
    region._normalize_grid = True  # Original map filtering and normalization.
    results, _ = _refine_density_grids(grid, region.grid(blur), [pose], maximum, objective)
    return results[0][:2]


def _shape_common_map_score(source, target, pose, blur, maximum, objective):
    """The original overlap objective on one common map scale, without tile norms.

    The full reference energy is a constant for all placements. Omitting that
    common denominator preserves ranking without allocating a whole-map fine
    grid. This score is an overlap, not a normalized correlation coefficient.
    """
    import numpy as np
    from chimerax.map_fit.fitmap import volume_values

    grid = source.grid(blur)
    points, weights = _shape_density_samples(grid, maximum)
    region = _shape_map_box(target, pose.transform_points(points),
                            6 * blur / (pi * sqrt(2)), blur)
    if region is None:
        return -np.inf
    region._finite_gaussian = False  # Same filter as the original search.
    reference = region.grid(blur)
    field, _ = _shape_reference_field(reference, objective)
    if objective == "gradient":
        _, energy = _gradient_density_field(grid)
    else:
        energy = float(np.einsum("ijk,ijk->", grid.matrix(), grid.matrix(), dtype=np.float64))
    values = volume_values(points, reference.xyz_to_ijk_transform * pose, field)
    return float(np.dot(weights.astype(np.float64), values) / sqrt(energy))


def _shape_symmetry_operators(session, model, symmetry, center, axis, copies=1):
    """Native symmetry operators in reference-map coordinates, never scene axes."""
    import numpy as np
    from chimerax.atomic.args import parse_symmetry
    from chimerax.core.errors import UserError

    group = symmetry if isinstance(symmetry, str) else symmetry.group
    requested_group = group
    fields = group.lower().split(',')
    kind = fields[0]
    automatic_helix = kind == 'h' and len(fields) == 3 and '*' not in group
    if automatic_helix:
        fields += [str(copies), str(-(copies // 2))]
        group = ','.join(fields)
    if '*' in group:
        count = 1
        for component in group.split('*'):
            if component.lower().startswith('h,') and len(component.split(',')) == 3:
                raise UserError('Within symmetry products, use h,rise,twist,count[,offset].')
            ops, _ = _shape_symmetry_operators(session, model, component, center, axis)
            count *= len(ops)
            if count > 256:
                raise UserError('A symmetry product must contain at most 256 operators.')
    elif kind.startswith(('c', 'd')):
        try:
            count = int(kind[1:]) * (2 if kind[0] == 'd' else 1)
        except ValueError:
            raise UserError('Use a native point group, such as C6, D4, T, O, or I.') from None
        if not 1 <= count <= 256:
            raise UserError('Symmetry must contain between 1 and 256 operators.')
    elif kind == 'h':
        try:
            rise, twist, count = map(float, fields[1:4])
            offset = float(fields[4]) if len(fields) == 5 else 0.
            valid = (len(fields) in (4, 5) and np.isfinite([rise, twist, count, offset]).all()
                     and rise != 0 and count == int(count) and 1 <= count <= 512
                     and offset == int(offset) and offset <= 0 < offset + count)
        except (ValueError, OverflowError):
            valid = False
        if not valid:
            raise UserError('Use symmetry h,rise,twist[,count,offset]; rise must be nonzero, '
                            'count must be 1-512, and an explicit repeat range must include zero.')
    elif kind not in ('t', 'o', 'i'):
        raise UserError('Symmetry supports native point groups Cn, Dn, T, O, I, or h,rise,twist,count[,offset].')

    frame = model.scene_position
    c = (model.data.ijk_to_xyz_transform * (np.array(model.data.size, dtype=float) / 2) if center is None else
         frame.inverse() * center.scene_coordinates(frame) if hasattr(center, 'scene_coordinates')
         else np.array(center, dtype=float))
    a = ((0, 0, 1) if axis is None else
         frame.inverse().transform_vectors(axis.scene_coordinates(frame)) if hasattr(axis, 'scene_coordinates')
         else np.array(axis, dtype=float))
    a = np.array(a, dtype=float)
    if (np.shape(c) != (3,) or a.shape != (3,) or not np.isfinite(c).all()
            or not np.isfinite(a).all() or np.linalg.norm(a) == 0):
        raise UserError('symmetry_center and symmetry_axis must be finite XYZ values; the axis must be nonzero.')
    a /= np.linalg.norm(a)
    operators = []
    for operator in parse_symmetry(session, group, tuple(c), tuple(a)):
        if not any(np.allclose(operator.matrix, old.matrix, atol=1e-8, rtol=0) for old in operators):
            operators.append(operator)
    operators.sort(key=lambda op: not op.is_identity())
    metadata = dict(symmetry=requested_group, symmetry_center=np.asarray(c).tolist(),
                    symmetry_axis=a.tolist(), symmetry_operators=len(operators))
    if kind == 'h' and '*' not in group:
        metadata.update(helical_rise=rise, helical_twist=twist, helical_auto_range=automatic_helix)
    return operators, metadata


def _shape_pose_inside(source, target, raw):
    import numpy as np

    points = raw.transform_points(source.xyz) - target.center
    return np.average(_shape_map_contains(target, points), weights=source.weights) >= .99


def _shape_symmetry_clashes(source, context, raw, maximum, distance):
    """Test implicit neighbours using one confidence-weighted source KD-tree."""
    import numpy as np
    from chimerax.geometry import Place, helical_symmetry_matrix

    if '_body' not in context:
        context['_body'] = _ShapeOccupiedCopy(source, Place())
    body = context['_body']
    parameters = context.get('parameters', {})
    if 'helical_rise' in parameters:
        axis = np.array(parameters['symmetry_axis'])
        span = np.ptp(raw.transform_points(source.xyz) @ axis)
        count = int(np.ceil((span + distance) / abs(parameters['helical_rise'])))
        # Inverse repeats have identical symmetric overlap. Positive repeats
        # suffice, and axial separation bounds all possible contacting neighbours.
        operators = (helical_symmetry_matrix(parameters['helical_rise'], parameters['helical_twist'],
                     axis=axis, center=parameters['symmetry_center'], n=n) for n in range(1, count + 1))
    else:
        operators = (op for op in context['operators'] if not op.is_identity())
    inverse = raw.inverse()
    for operator in operators:
        relative = inverse * operator * raw
        for pose in (relative, relative.inverse()):
            points = pose.transform_points(source.xyz)
            if np.any(points.max(axis=0) + distance < body.lower) or np.any(body.upper + distance < points.min(axis=0)):
                continue
            neighbours = body.tree.query_ball_point(points, distance)
            weights = np.array([body.weights[indices].max() if len(indices) else 0 for indices in neighbours])
            if np.dot(body.weights, weights) / body.weights.sum() > maximum:
                return True
    return False


def _shape_expansion_operators(source, context, base):
    """Finite point/range operators, or an automatic helix within the map box."""
    import numpy as np
    from chimerax.geometry import Place, helical_symmetry_matrix

    parameters = context.get('parameters', {})
    if not parameters.get('helical_auto_range'):
        return iter(context['operators'])
    target = context['target']
    axis = np.array(parameters['symmetry_axis'])
    bounds = _shape_grid_corners(target.model.data) @ axis
    projected = base.transform_points(source.xyz) @ axis
    rise = parameters['helical_rise']
    interval = sorted(((bounds.min() - projected.max()) / rise, (bounds.max() - projected.min()) / rise))
    lower, upper = int(np.ceil(interval[0])), int(np.floor(interval[1]))
    def generate():
        yield Place()
        for n in range(1, max(abs(lower), abs(upper)) + 1):
            for repeat in (-n, n):
                if lower <= repeat <= upper:
                    yield helical_symmetry_matrix(rise, parameters['helical_twist'], axis=axis,
                                                  center=parameters['symmetry_center'], n=repeat)
    return generate()


def _shape_expand_seed(source, target, stats, options, best):
    """Choose one seed, then expand it; no joint or repeated global fitting."""
    import numpy as np
    from chimerax.geometry import Place, translation
    from chimerax.core.errors import UserError

    context = options['_symmetry']
    root = context['target']
    rows = stats.get('candidates', [dict(score=stats['score'], matrix=best.matrix.tolist())])
    base, score = None, None
    for entry in sorted(rows, key=lambda row: row['score'], reverse=True):
        _shape_checkpoint()
        candidate = Place(matrix=np.array(entry['matrix']))
        if (_shape_pose_inside(source, root, candidate) and not
                _shape_symmetry_clashes(source, context, candidate, options['max_overlap'], options['overlap_distance'])):
            base, score = candidate, entry['score']
            break
    if base is None:
        raise UserError('No seed fit is compatible with its implicit symmetry neighbours.')
    accepted, occupied = [], []
    _shape_map_status(source.model.session.logger, 'Expanding the fitted seed by symmetry...')
    for operator in _shape_expansion_operators(source, context, base):
        _shape_checkpoint()
        ideal = operator * base
        if not _shape_pose_inside(source, root, ideal):
            continue
        pose = ideal
        if options.get('refine_copies'):
            _shape_map_status(source.model.session.logger,
                              f"Local refinement of symmetry copy {len(accepted) + 1}/{options.get('copies', 1)}...")
            centered = translation(-root.center) * ideal * translation(source.center)
            _, refined = _shape_refine_regional(source, root, centered, stats['resolutions'][-1],
                                                options['max_points'], options['objective'])
            proposed = translation(root.center) * refined * translation(-source.center)
            if _shape_pose_inside(source, root, proposed):
                pose = proposed
        body = _ShapeOccupiedCopy(source, pose)
        if any(body.overlap(old, options['overlap_distance']) > options['max_overlap'] for old in occupied):
            pose, body = ideal, _ShapeOccupiedCopy(source, ideal)
            if any(body.overlap(old, options['overlap_distance']) > options['max_overlap'] for old in occupied):
                continue
        accepted.append((pose, score))
        occupied.append(body)
        if len(accepted) == options.get('copies', 1):
            break
    if not accepted:
        raise UserError('No symmetry-expanded copy fits inside the map.')
    stats.update(symmetry_mode='seed_then_expand', symmetry_search_copies=1,
                 refine_copies=bool(options.get('refine_copies')), seed_matrix=base.matrix.tolist())
    return accepted


class _ShapeOccupiedCopy:
    """Cache an accepted placement's atoms, confidence, and spatial index."""

    def __init__(self, source, pose):
        import numpy as np
        from scipy.spatial import cKDTree

        self.points = pose.transform_points(source.xyz)
        self.weights = np.ones(len(self.points)) if source.weights is None else source.weights / source.weights.max()
        self.lower, self.upper = self.points.min(axis=0), self.points.max(axis=0)
        self.tree = cKDTree(self.points)

    def overlap(self, other, distance):
        import numpy as np

        if np.any(self.upper + distance < other.lower) or np.any(other.upper + distance < self.lower):
            return 0.

        def fraction(first, second):
            neighbours = second.tree.query_ball_point(first.points, distance)
            weights = np.array([second.weights[indices].max() if len(indices) else 0 for indices in neighbours])
            return float(np.dot(first.weights, weights) / first.weights.sum())

        # Check both directions: a small confident core must not be hidden by
        # the larger flexible region of the other placement.
        return max(fraction(self, other), fraction(other, self))


def _shape_copy_site_candidates(source, candidates, maximum, options):
    """Keep different finished sites, not several orientations of one tile hit."""
    import numpy as np
    from chimerax.geometry import Place

    kept, occupied = [], []
    for entry in candidates:
        _shape_checkpoint()
        pose = Place(matrix=np.array(entry['matrix']))
        body = _ShapeOccupiedCopy(source, pose)
        if any(body.overlap(old, options['overlap_distance']) > options['max_overlap'] for old in occupied):
            continue
        kept.append(entry)
        occupied.append(body)
        if len(kept) == maximum:
            break
    return kept


def _shape_select_copies(source, target, stats, options, best):
    """Greedy score-ranked, confidence-weighted non-overlap across all tiles."""
    import numpy as np
    from chimerax.geometry import Place
    from chimerax.core.errors import UserError

    if options.get('_symmetry') is not None:
        return _shape_expand_seed(source, target, stats, options, best)
    candidates = stats.get('candidates', [dict(score=stats['score'], matrix=best.matrix.tolist())])
    candidates = sorted(candidates, key=lambda row: row['score'], reverse=True)
    operators = [Place()]
    requested = options.get('copies', 1)
    maximum, distance = options['max_overlap'], options['overlap_distance']
    accepted, occupied = [], []
    cutoff = candidates[0]['score'] * options['score_fraction']
    if requested > 1 and candidates[0]['score'] <= 0:
        raise UserError('No positive-scoring placements were found for multiple-copy fitting.')
    for entry in candidates:
        if accepted and entry['score'] < cutoff:
            break
        base = Place(matrix=np.array(entry['matrix'], dtype=float))
        orbit, blocked = [], False
        for operator in operators:
            pose = operator * base
            candidate = _ShapeOccupiedCopy(source, pose)
            for old_pose, old in [*occupied, *orbit]:
                displacement = (pose.transform_points(source.frame_landmarks + source.center) -
                                old_pose.transform_points(source.frame_landmarks + source.center))
                if (np.sqrt(np.mean(np.sum(displacement ** 2, axis=1))) < stats['resolutions'][-1] / 3
                        or candidate.overlap(old, distance) > maximum):
                    blocked = True
                    break
            if blocked:
                break
            orbit.append((pose, candidate))
        if blocked:
            continue
        for pose, candidate in orbit:
            accepted.append((pose, entry['score']))
            occupied.append((pose, candidate))
            if len(accepted) == requested:
                return accepted
    if not accepted:
        raise UserError('No non-overlapping placement was found.')
    return accepted


def _fit_shape_tiled_maps(source, target, options, center):
    """Run the original global-search engine independently on padded submaps."""
    from contextlib import nullcontext
    import numpy as np
    from chimerax.geometry import Place, translation
    from chimerax.core.errors import UserError

    levels, stats = _shape_fit_levels(source, target, options)
    windows, tile_stats = _shape_map_tiles(source, target, levels, options.get("tile_size"))
    if not windows:
        raise UserError("No submap has enough nonzero three-dimensional density for fitting.")
    logger = source.model.session.logger
    progress = [0, len(windows)]
    progress_token = _SHAPE_TILE_PROGRESS.set(progress)
    scorers, pooled = {}, []
    context = _ShapeBatchCache() if _SHAPE_BATCH_CACHE.get() is None else nullcontext()
    previous_shared = source._shared_reference
    source._shared_reference = True
    try:
        with context:
            for index, (window, _, _) in enumerate(windows):
                _shape_map_status(logger, f"Tiled shape search: tile {index + 1}/{len(windows)} (original global search)...")
                shift = translation(window.center - target.center)
                def valid(pose):
                    inside = _shape_map_contains(target, (shift * pose).transform_points(source.landmarks))
                    weights = source.weights if options.get('_symmetry') is not None else None
                    return np.average(inside, weights=weights) >= .99
                finished = False
                try:
                    _, local = _fit_shape_original_map_search(
                        source, window, options, None, force_global=True, collect=True,
                        ranking_target=target, shared_levels=levels, scorers=scorers, valid_pose=valid)
                    finished = True
                    # Compare finished fits, not prematurely normalized tile
                    # peaks. Each submap contributes several independent basins.
                    contribution = max(_SHAPE_SUBMAP_FITS, options.get('copies', 1))
                    candidates = local['candidates']
                    if options.get('copies', 1) > 1:
                        candidates = _shape_copy_site_candidates(source, candidates, contribution, options)
                    for candidate in candidates[:contribution]:
                        raw = Place(matrix=np.array(candidate["matrix"], dtype=float))
                        pooled.append((candidate["score"], translation(-target.center) * raw *
                                       translation(source.center), candidate["fit_score"]))
                except _ShapeEmptyMapSearch:
                    finished = True
                    continue
                finally:
                    window._cache.clear()
                    if finished:
                        progress[0] += 1
                        _shape_map_status(logger, f"Finished tile {index + 1}/{len(windows)}.")
            if not pooled:
                raise UserError("No complete moving-model placement was found inside the reference map.")
            merged = []
            _shape_map_status(logger, 'Comparing completed submap fits...')
            for finished_score, pose, local_score in pooled:
                if options.get("bandpass"):
                    # Each submap already used this same final-band scorer
                    # against the original map; do not filter it a second time.
                    score = finished_score
                else:
                    score = _shape_common_map_score(source, target, pose, levels[-1],
                                                    options["max_points"], options["objective"])
                merged.append((score, pose, local_score))
            pool_size = max(_SHAPE_MAP_BEAM, options.get('copies', 1) * 8)
            kept = _shape_keep_map_poses(merged, source.frame_landmarks, levels[-1] / 3, maximum=pool_size)
    finally:
        source._shared_reference = previous_shared
        _SHAPE_TILE_PROGRESS.reset(progress_token)
    if not kept:
        raise UserError("No submap fit has a defined common-reference score.")
    score, pose, fit_score = kept[0]
    bands = [(max(blur, 1.5 * options["map_resolution"]), max(options["map_resolution"], blur / 2))
             for blur in levels] if options.get("bandpass") else []
    stats.update(tile_stats, score=score, fit_score=fit_score, objective=options["objective"], partial=False,
                 local=False, map_search=True, tiled=True, bandpass=bool(options.get("bandpass")),
                 bandpasses=bands, map_resolution=options.get("map_resolution"), resolutions=levels,
                 moving_input_type="atoms", reference_input_type="map",
                 search_engine="original_global", submap_candidates=len(pooled),
                 search_orientations=122 + options["extra_starts"],
                 score_metric=("model-local bandpass correlation" if bands else
                               f"common-map {options['objective']} overlap"),
                 candidates=[dict(score=s, fit_score=f, matrix=(translation(target.center) * p *
                                  translation(-source.center)).matrix.tolist()) for s, p, f in
                             (kept if options.get('copies', 1) > 1 or options.get('_symmetry') else kept[:10])])
    if options.get("radius_search") is not None:
        stats["radius_search"] = options["radius_search"]
    if center is not None:
        stats["center_of_rotation"] = center.copy()
    return translation(target.center) * pose * translation(-source.center), stats


def _fit_shape_density_inputs(source, target, options, center):
    from chimerax.core.errors import UserError

    if options.get("tiled"):
        return _fit_shape_tiled_maps(source, target, options, center)
    try:
        return _fit_shape_original_map_search(source, target, options, center,
                                              force_global=options.get('copies', 1) > 1 and options['radius'] is None,
                                              collect=bool(options.get("bandpass") or options.get('_symmetry')
                                                           or options.get('copies', 1) > 1))
    except _ShapeEmptyMapSearch:
        raise UserError("No valid placement with a defined fitting score was found inside the map.") from None


def _fit_shape_original_map_search(source, target, options, center, force_global=False,
                                    collect=False, ranking_target=None, shared_levels=None,
                                    scorers=None, valid_pose=None):
    """The original PCA/FFT search and sum-product refiner, shared by all submaps.

    With bandpass off this is the existing algorithm, not a local Pearson fit.
    Optional band scores change candidate ranking only, never optimization.
    """
    import numpy as np
    from chimerax.geometry import translation

    context = options.get('_symmetry')
    if context is not None:
        original_valid = valid_pose
        raw_frame = translation(target.center)
        valid_pose = lambda pose: ((original_valid is None or original_valid(pose)) and
                                  _shape_pose_inside(source, context['target'],
                                                     raw_frame * pose * translation(-source.center)))
    objective, radius = options["objective"], options["radius"]
    levels, stats = _shape_fit_levels(source, target, options)
    if shared_levels is not None:
        levels = shared_levels
    ranking_target = target if ranking_target is None else ranking_target
    scorers = {} if scorers is None else scorers
    fft_options = {} if valid_pose is None else dict(valid_pose=valid_pose)
    if force_global:
        # Keep several translation basins at model-size spacing. Fine-band
        # ranking cannot recover sites removed by a bright coarse-density peak.
        fft_options["site_distance"] = _shape_spatial_size(source.frame_landmarks)
        fft_options["peaks_per_orientation"] = 8
    last_kept = []
    initial = (translation(-target.center) * target.model.scene_position.inverse() *
               source.model.scene_position * translation(source.center))
    starts = [initial] if radius is not None else [*_shape_starts(
        source.frame_landmarks, target.frame_landmarks, options["extra_starts"]), initial]
    search_radius = options.get("radius_search")
    if search_radius is not None or force_global:
        # A regional global search must explore translations even when PCA
        # already gives a good score. Only the reference density was cropped.
        _shape_map_status(source.model.session.logger, "Shape alignment: global search of " +
                          ("a numerical submap..." if force_global else "the reference CoFR sphere..."))
        starts.extend(_shape_map_fft_starts(
            source.grid(levels[0]), target.grid(levels[0]), source.frame_landmarks,
            target.frame_landmarks, objective, options["extra_starts"], **fft_options))

    def refine(candidates, first=0):
        nonlocal last_kept
        for blur in levels[first:]:
            _shape_map_status(source.model.session.logger,
                              f"Shape alignment: map density at {blur:.3g} A ({len(candidates)} starts)...")
            results, points = _refine_density_grids(source.grid(blur), target.grid(blur), candidates,
                                                   options["max_points"], objective)
            ranked = []
            scorer = None
            if options.get("bandpass"):
                floor = options["map_resolution"]
                key = (max(blur, 1.5 * floor), max(floor, blur / 2))
                if key not in scorers:
                    scorers[key] = _ShapeBandScorer(source, *key, options["max_points"])
                scorer = scorers[key]
            shift = translation(target.center - ranking_target.center)
            for fit_score, pose, _ in results:
                _shape_checkpoint()
                if valid_pose is not None and not valid_pose(pose):
                    continue
                if context is not None:
                    raw = translation(target.center) * pose * translation(-source.center)
                    if _shape_symmetry_clashes(source, context, raw, options['max_overlap'], options['overlap_distance']):
                        continue
                score = scorer.score(ranking_target, shift * pose) if scorer is not None else fit_score
                if np.isfinite(score):
                    ranked.append((score, pose, fit_score))
            if scorer is not None:
                ranked.sort(key=lambda row: row[0], reverse=True)
            candidates, last_kept = [], []
            probes = points[np.linspace(0, len(points) - 1, min(64, len(points)), dtype=int)]
            for entry in ranked:
                _, pose, _ = entry
                if any(np.sqrt(np.mean(np.sum((pose.transform_points(probes) - old.transform_points(probes)) ** 2, axis=1)))
                       < blur / 3 for old in candidates):
                    continue
                candidates.append(pose)
                last_kept.append(entry)
                if len(candidates) == 48:
                    break
            if not last_kept:
                raise _ShapeEmptyMapSearch()
        return last_kept[0][:2]

    score, pose = refine(starts)
    stats.update(objective=objective, partial=False, local=radius is not None, map_search=False,
                 moving_input_type="map" if source.is_map else "atoms",
                 reference_input_type="map" if target.is_map else "atoms", resolutions=levels)
    if search_radius is not None:
        stats.update(map_search=True, map_search_resolution=levels[0], radius_search=search_radius)
    elif force_global:
        stats.update(map_search=True, map_search_resolution=levels[0])
    elif radius is None and options["fallback_correlation"] > 0 and last_kept[0][2] < options["fallback_correlation"]:
        first = 1 if len(levels) > 2 else 0
        blur = levels[first]
        seeds = _shape_map_fft_starts(source.grid(blur), target.grid(blur), source.frame_landmarks,
                                      target.frame_landmarks, objective, options["extra_starts"], **fft_options)
        previous_kept = last_kept
        proposed_score, proposed = refine([pose, *seeds], first)
        stats.update(map_search=True, map_search_resolution=blur, map_initial_score=score)
        if proposed_score > score:
            score, pose = proposed_score, proposed
        else:
            last_kept = previous_kept
    stats["score"] = score
    if options.get("bandpass"):
        stats.update(bandpass=True, tiled=False, fit_score=last_kept[0][2],
                     score_metric="model-local bandpass correlation",
                     bandpasses=[(max(blur, 1.5 * options["map_resolution"]),
                                  max(options["map_resolution"], blur / 2)) for blur in levels])
    if collect:
        stats["candidates"] = [dict(score=s, fit_score=f, matrix=(translation(target.center) * p *
                                translation(-source.center)).matrix.tolist()) for s, p, f in last_kept]
    if center is not None:
        stats["center_of_rotation"] = center.copy()
        if radius is not None:
            stats.update(radius=radius, radius_refine=radius)
    return translation(target.center) * pose * translation(-source.center), stats


def _shape_map_correlations(source, target, pose, resolution, window, regional=False):
    import numpy as np
    from chimerax.core.errors import UserError
    from chimerax.geometry import translation
    from chimerax.map_data import ArrayGridData
    from chimerax.map_fit.fitmap import volume_values

    placed = pose.transform_points(source.xyz) - target.center
    step = resolution / 3
    pad = 4 * resolution / (pi * sqrt(2)) + 3 * window / sqrt(8 * np.log(2)) + step
    origin = placed.min(axis=0) - pad
    size = np.ceil((placed.max(axis=0) + pad - origin) / step).astype(int) + 1
    if np.prod(size) > _MAX_CORRELATION_VOXELS:
        raise UserError("The local correlation grid is too large; use a coarser correlation resolution.")
    if regional:
        target = _shape_map_box(target, placed, pad)
        if target is None:
            raise UserError("The fitted model does not overlap the map correlation region.")
    centered_pose = translation(-target.center) * pose * translation(source.center)
    first = _resample_shape_density(source.grid(resolution), origin, size, step, centered_pose)
    second = _resample_shape_density(target.grid(resolution), origin, size, step)
    field, valid = _local_correlation_field(first, second, window / sqrt(8 * np.log(2)) / step)
    grid = ArrayGridData(field, origin, (step,) * 3)
    scores = volume_values(placed, grid.xyz_to_ijk_transform, field).copy()
    scores[volume_values(placed, grid.xyz_to_ijk_transform, valid) < .999] = np.nan
    return scores, {"correlation_resolution": resolution, "correlation_window": window,
                    "correlation_grid_voxels": int(field.size),
                    "correlation_defined_atoms": int(np.isfinite(scores).sum())}


def _shape_align_with_map(session, source_item, target_item, options, center, target_input=None, started=None):
    from time import perf_counter
    from chimerax.core.undo import UndoState

    started = perf_counter() if started is None else started
    radius = options["radius"]
    # Local refinement crops both inputs; global search keeps the moving input
    # intact, including when it starts completely outside the reference sphere.
    source = _ShapeDensityInput(*source_item, center if radius is not None else None, radius,
                                options.get("confidence_moving") or options.get("confidence", "auto"))
    auto_search = options.get("radius_search") is True
    if auto_search:
        search_radius = _shape_search_radius(source)
        options = dict(options, radius_search=search_radius)
        session.logger.info(f"Automatic radius_search for #{source.model.id_string}: {search_radius:.3g} A "
                            "(1.5 times the moving input's enclosing radius).")
    target = target_input if target_input is not None else _shape_reference_input(
        *target_item, center, radius if radius is not None else options.get("radius_search"),
        options.get("confidence_reference") or options.get("confidence", "auto"))
    symmetry_spec = options.get('_symmetry_spec')
    if symmetry_spec is not None:
        # The seed search sphere does not crop the density of symmetry mates.
        full_target = _shape_reference_input(target.model, None, None, None)
        options = dict(options, _symmetry=dict(target=full_target, operators=symmetry_spec[0],
                                              parameters=symmetry_spec[1]))
        description = symmetry_spec[1]
        axis_text = ','.join(f'{v:.4g}' for v in description['symmetry_axis'])
        center_text = ','.join(f'{v:.4g}' for v in description['symmetry_center'])
        session.logger.info(f"Symmetry {description['symmetry']}: one seed, then expansion; "
                            f"reference-map axis {axis_text}, center {center_text}.")
    _shape_log_confidence(session, source.model, source.weights, "moving")
    _shape_log_confidence(session, target.model, target.weights, "reference")
    workers = _SHAPE_MAP_WORKERS.get()
    multiple = options.get('copies', 1) > 1
    with _ShapeMapInputGuard((source_item, target_item)) as guard:
        search_options = dict(options, copies=1) if symmetry_spec is not None else options
        pose, stats = _fit_shape_density_inputs(source, target, search_options, center)
        placements = (_shape_select_copies(source, target, stats, options, pose)
                      if multiple or symmetry_spec is not None else [(pose, stats['score'])])
        guard.validate()
        if workers is not None and workers.used:
            stats['map_workers'] = len(workers.processes)
    stats.update(moving_confidence="off" if source.weights is None else "plddt",
                 reference_confidence="off" if target.weights is None else "plddt")
    if symmetry_spec is not None:
        stats.update(symmetry_spec[1])
    if options.get("radius_search") is not None:
        stats["radius_search_automatic"] = auto_search
    stats['score'] = placements[0][1]
    stats["seconds"] = perf_counter() - started
    model, reference = source.model, target.model
    transform = reference.scene_position * placements[0][0] * model.scene_position.inverse()
    undo = UndoState("shape align")
    original_model = model
    created, records = [], []
    scores = None
    scoring_started = perf_counter()
    try:
        for index, (local_pose, score) in enumerate(placements):
            extra = {}
            if options['correlation'] or options['color_fit']:
                scoring_target = options['_symmetry']['target'] if symmetry_spec is not None else target
                t0 = perf_counter()
                scores, extra = _shape_map_correlations(
                    source, scoring_target, local_pose, options['correlation_resolution'],
                    options['correlation_window'],
                    regional=bool(stats.get('tiled') or stats.get('bandpass') or symmetry_spec is not None))
                extra['correlation_seconds'] = perf_counter() - t0
            current, scoring_atoms = original_model, source.atoms
            if options['move']:
                motion = reference.scene_position * local_pose * original_model.scene_position.inverse()
                current, scoring_atoms = _shape_place_model(
                    session, original_model, source.atoms, motion, undo, options.get('copy', False) or multiple)
                if current is not original_model:
                    created.append(current)
                    if multiple:
                        current.name = f'{original_model.name} fit {index + 1}'
            if scores is not None:
                _assign_shape_correlations(session, current, scoring_atoms, scores, undo, options['color_fit'])
            if index == 0:
                model = current
                stats.update(extra)
            records.append(dict(score=score, matrix=local_pose.matrix.tolist(),
                                model_id=current.id_string if options['move'] else None))
    except BaseException:
        undo.undo()
        if created:
            session.models.close(created)
        raise
    if undo.state:
        session.undo.register(undo)
    if scores is not None:
        stats["correlation_grid_seconds"] = stats["correlation_seconds"]
        stats["correlation_seconds"] = perf_counter() - scoring_started
    stats.update(moving_model_id=model.id_string, reference_model_id=reference.id_string)
    if multiple or symmetry_spec is not None:
        stats.update(copies_requested=options.get('copies', 1), copies_found=len(records), placements=records,
                     max_overlap=options['max_overlap'], overlap_distance=options['overlap_distance'],
                     score_fraction=options['score_fraction'])
    if multiple:
        stats['moving_model_ids'] = [record['model_id'] for record in records if record['model_id'] is not None]
        if len(records) < options['copies']:
            if symmetry_spec is not None:
                session.logger.warning(f"Expanded {len(records)}/{options['copies']} copies; the symmetry range, "
                                       "map bounds, or overlap checks limit expansion.")
            else:
                session.logger.warning(f"Found {len(records)}/{options['copies']} non-overlapping copies among the "
                                       "retained candidates above the score-fraction cutoff.")
    if model is not original_model:
        stats["source_model_id"] = original_model.id_string
        fitted = ', '.join(f'#{item.id_string}' for item in created)
        session.logger.info(f"Shape alignment: fitted {'copies' if len(created) > 1 else 'copy'} {fitted}; "
                            f"original #{original_model.id_string} unchanged.")
    metric = stats.get("score_metric") or ("whole gradient correlation" if options["objective"] == "gradient" else
                                            "whole density overlap")
    region = (f"reference density in a {options['radius_search']:g} A CoFR search sphere."
              if options.get("radius_search") is not None else
              f"density in a {radius:g} A CoFR refinement sphere." if radius is not None else
              "full map density, independent of contour.")
    message = (f"Shape alignment #{model.id_string} ({stats['moving_input_type']}) to "
               f"#{reference.id_string} ({stats['reference_input_type']}): {metric} "
               f"score {stats['score']:.3f}; " + (f"{len(records)}/{options['copies']} copies; " if multiple else '') +
               f"{stats['seconds']:.2f} s; {region}" +
               (" (not moved)" if not options["move"] else ""))
    session.logger.info(message)
    session.logger.status(message)
    session.logger.info("Fitting scales: " + " -> ".join(f"{v:.3g}" for v in stats["resolutions"]) +
                        " A; Gaussian blur is additional to the map's existing resolution.")
    if stats.get("tiled"):
        session.logger.info(f"Tiled search: {stats['tiles']} translation cores, {stats['nonempty_tiles']} nonempty; "
                            f"up to {stats['search_orientations']} orientations per tile; "
                            f"largest FFT {stats['tile_fft_voxels']:,} voxels.")
    if stats.get("bandpass"):
        session.logger.info("Model-local validation bands (smooth Gaussian differences): " +
                            " -> ".join(f"{coarse:.3g}-{fine:.3g}" for coarse, fine in stats["bandpasses"]) +
                            " A. These are ranking scores, not probabilities or independent validation data.")
    if scores is not None:
        session.logger.info(f"Local density CC: {stats['correlation_defined_atoms']}/{len(source.atoms)} atoms defined; "
                            f"{stats['correlation_seconds']:.2f} s extra; saved shape_correlation" +
                            (" and colored fitted atoms/cartoons." if options["color_fit"] else "; colors unchanged."))
    return transform, stats


def _shape_search_radius(source):
    """Leave 50% padding around a rotation-independent enclosing sphere."""
    import numpy as np

    if source.is_map:
        # Use full nonzero density bounds, never a contour or displayed surface.
        from itertools import product

        corners = np.array(list(product(*zip(source.lower, source.upper))), dtype=float)
        points = source.data.ijk_to_xyz_transform.transform_points(corners) - source.center
    else:
        points = source.landmarks  # Fitting atoms relative to their centroid.
    return 1.5 * float(np.linalg.norm(points, axis=1).max())


def _shape_align_help(session):
    """Show a compact quick reference; native usage supplies the full option list."""
    rows = (
        ("radius_search [A]", "Global search in a reference-map CoFR sphere; omit A for 1.5x model radius."),
        ("radius_refine A", "Local refinement at CoFR. Use only one radius; both are in Angstroms."),
        ("copies N", "Model/map: up to N non-overlapping copies; N > 1 leaves the source unchanged."),
        ("symmetry C6 / h,rise,twist", "Fit one seed, check implied clashes, then expand; default axis: map Z."),
        ("resolution A / objective density|gradient", "Final scale: automatic; objective: density."),
        ("confidence auto|off|plddt", "Weight AlphaFold confidence automatically, disable, or force B-factor pLDDT."),
        ("workers 0|N / color_fit true", "CPUs: 0 automatic, 1 serial; or color atoms by local density agreement."),
    )
    description = ("Rigidly fit models/maps; copy true fits a model copy. "
                   "Model/map: tiled true; bandpass true map_resolution A. Supports Undo.")
    example = "shapealign #2 to #1 radius_search"
    if session.ui.is_gui:
        from html import escape
        from chimerax.log.cmd import log_show

        log_show(session)
        table = "".join(f"<tr><td><b>{escape(option)}</b></td><td>{escape(text)}</td></tr>"
                        for option, text in rows)
        session.logger.info(
            "<b>shapealign</b> <i>moving-spec</i> <b>to</b> <i>reference-spec</i> [options]"
            f"<br>{description}<table cellpadding=\"2\">{table}</table>"
            f"Example: <code>{example}</code><br>"
            '<a href="cxcmd:usage%20shapealign">All options</a> '
            "(<b>usage shapealign</b>); <b>radius</b> is the legacy local-refinement alias.",
            is_html=True,
        )
        session.logger.status("shapealign help is shown in the Log.")
    else:
        session.logger.info(
            "Usage: shapealign <moving-spec> to <reference-spec> [options]\n" + description + "\n" +
            "\n".join(f"  {option}: {text}" for option, text in rows) +
            f"\nExample: {example}. All options: usage shapealign (radius is the legacy local-refinement alias)."
        )


def shape_align(session, moving=None, to=None, resolution=None, max_points=12000,
                extra_starts=0, move=True, radius=None, correlation=False,
                color_fit=False, correlation_resolution=8.0, correlation_window=12.0,
                partial=None, fallback_correlation=None, objective="density",
                radius_search=None, radius_refine=None, workers=0,
                confidence="auto", confidence_moving=None, confidence_reference=None, copy=False,
                tiled=False, bandpass=False, map_resolution=None, tile_size=None,
                copies=1, max_overlap=.1, overlap_distance=4., score_fraction=.8,
                symmetry=None, symmetry_center=None, symmetry_axis=None, refine_copies=False):
    """Align independent inputs globally, regionally, or by local refinement.

    radius_search globally searches a reference-map CoFR sphere, keeping the
    moving selection intact. True requests 1.5 times its enclosing radius.
    radius_refine locally refines the current pose
    using both inputs' CoFR spheres. radius is the legacy radius_refine alias.
    """
    from contextlib import nullcontext
    from math import isfinite
    from time import perf_counter
    import numpy as np
    from chimerax.core.errors import UserError

    if moving is None and to is None:
        return _shape_align_help(session)
    if getattr(session, '_shapealign_running', False):
        raise UserError('shapealign is already running; wait for it to finish.')
    if moving is None or to is None:
        raise UserError("Specify both moving and reference inputs: shapealign <moving-spec> to <reference-spec>. "
                        "Type shapealign help for options.")
    sources, targets = _shape_input_items(moving), _shape_input_items(to)
    if not sources or len(targets) != 1:
        raise UserError("Specify moving atomic models/maps and exactly one reference model or map.")
    _shape_worker_count(workers, len(sources))
    for value in (confidence, confidence_moving, confidence_reference):
        if value is not None and value not in ("auto", "off", "plddt"):
            raise UserError("confidence, confidence_moving, and confidence_reference must be auto, off, or plddt.")
    if sum(value is not None for value in (radius, radius_refine, radius_search)) > 1:
        raise UserError("Use only one of radius_search (global), radius_refine (local), or its legacy alias radius.")
    if radius_refine is not None:
        radius = radius_refine
    for name, value in (("radius_refine", radius), ("radius_search", radius_search)):
        if value is not None and (not isfinite(value) or value <= 0):
            raise UserError(f"{name} must be finite and greater than zero (Angstroms).")
    if radius_search is not None and targets[0][1] is not None:
        raise UserError("radius_search requires a density map as the reference; use radius_refine for local model/model fitting.")
    target = targets[0][0]
    sources = _shape_exclude_reference(session, sources, target)
    if not sources:
        return []
    if isinstance(copies, bool) or not isinstance(copies, (int, np.integer)) or not 1 <= copies <= 512:
        raise UserError('copies must be an integer from 1 to 512 (the maximum number of fitted copies).')
    for name, value in (('max_overlap', max_overlap), ('score_fraction', score_fraction)):
        if not isfinite(value) or not 0 <= value <= 1:
            raise UserError(f'{name} must be between zero and one.')
    if not isfinite(overlap_distance) or overlap_distance <= 0:
        raise UserError('overlap_distance must be finite and positive, in Angstroms.')
    if (copies > 1 or symmetry is not None) and (targets[0][1] is not None or any(a is None for _, a in sources)):
        raise UserError('copies and symmetry currently require moving atomic models and a reference map.')
    if copies > 1 and radius is not None and symmetry is None:
        raise UserError('Multiple-copy search is global; use radius_search or tiled, not radius_refine.')
    if symmetry is None and (symmetry_center is not None or symmetry_axis is not None):
        raise UserError('symmetry_center and symmetry_axis require symmetry.')
    if refine_copies and symmetry is None:
        raise UserError('refine_copies is used after symmetry expansion; supply symmetry.')
    symmetry_spec = None
    if symmetry is not None:
        symmetry_spec = _shape_symmetry_operators(session, target, symmetry, symmetry_center, symmetry_axis, copies)
    if copy and not move:
        raise UserError("copy true requires move true; move false reports a fit without creating a copy.")
    if copy and any(atoms is None for _, atoms in sources):
        raise UserError("copy true supports moving atomic models only; moving maps are never duplicated.")
    regional_search = tiled or bandpass
    if regional_search and (targets[0][1] is not None or any(atoms is None for _, atoms in sources)):
        raise UserError("tiled and bandpass currently require moving atomic models and a reference map.")
    if tiled and radius is not None:
        raise UserError("tiled is a global search; do not combine it with radius_refine.")
    if bandpass and map_resolution is None:
        raise UserError("bandpass true requires map_resolution in Angstroms, for example map_resolution 2.2.")
    if map_resolution is not None:
        if not regional_search:
            raise UserError("map_resolution is used with tiled true or bandpass true.")
        if not isfinite(map_resolution) or map_resolution <= 0:
            raise UserError("map_resolution must be finite and positive (Angstroms).")
        if map_resolution + 1e-6 < 2 * max(target.data.step):
            raise UserError("map_resolution cannot be finer than the reference map's Nyquist sampling limit.")
        if resolution is not None and resolution < map_resolution:
            raise UserError("resolution cannot be finer than map_resolution.")
    if tile_size is not None and (not tiled or not isfinite(tile_size) or tile_size <= 0):
        raise UserError("tile_size requires tiled true and a finite positive length in Angstroms.")
    for model, _ in sources:
        if target in model.drawing_lineage or model in target.drawing_lineage:
            raise UserError("Moving and reference models must not contain one another.")
        if any(other in model.drawing_lineage[:-1] for other, _ in sources):
            raise UserError("Moving inputs must be independent models, not ancestor/descendant pairs.")
    options = dict(resolution=resolution, max_points=max_points, extra_starts=extra_starts, move=move,
                   radius=radius, correlation=correlation, color_fit=color_fit,
                   correlation_resolution=correlation_resolution, correlation_window=correlation_window,
                   partial=partial, fallback_correlation=fallback_correlation, objective=objective,
                   confidence=confidence, confidence_moving=confidence_moving, confidence_reference=confidence_reference,
                   copy=copy)
    if all(atoms is not None for _, atoms in [*sources, *targets]):
        from chimerax.atomic import concatenate

        with _ShapeUI(session):
            return _shape_align_atomic_models(session, concatenate([selected for _, selected in sources]),
                                              targets[0][1], workers=workers, **options)
    options.update(tiled=tiled, bandpass=bandpass, map_resolution=map_resolution, tile_size=tile_size)
    options.update(copies=copies, max_overlap=max_overlap, overlap_distance=overlap_distance,
                   score_fraction=score_fraction, _symmetry_spec=symmetry_spec, refine_copies=refine_copies)
    if radius_search is not None:
        options["radius_search"] = radius_search
    if objective not in ("density", "gradient") or partial not in (None, False):
        raise UserError("Map fitting supports objective density or gradient, with whole densities only; partial true is not supported.")
    if resolution is not None and (not isfinite(resolution) or resolution < 2):
        raise UserError("resolution must be finite and at least 2 Angstroms.")
    if max_points < 1000 or extra_starts < 0:
        raise UserError("max_points must be at least 1000 and extra_starts must be nonnegative.")
    cutoff = (.5 if objective == "gradient" else .7) if fallback_correlation is None else fallback_correlation
    if not isfinite(cutoff) or not 0 <= cutoff <= 1:
        raise UserError("fallback_correlation must be between 0 and 1.")
    options["fallback_correlation"] = cutoff
    if correlation or color_fit:
        if not move or any(atoms is None for _, atoms in sources):
            raise UserError("correlation/color_fit require moving atomic models and move true, not moving maps.")
    for value in (correlation_resolution, correlation_window):
        if not isfinite(value) or value <= 0:
            raise UserError("Correlation resolution and window must be finite and positive.")
    center = None
    if radius is not None and extra_starts:
        raise UserError("extra_starts is not supported with radius_refine (local refinement).")
    region_radius = radius if radius is not None else radius_search
    if region_radius is not None:
        center = np.array(session.main_view.center_of_rotation, dtype=float, copy=True)
        if center.shape != (3,) or not np.isfinite(center).all():
            raise UserError("The current center of rotation must be a finite XYZ position.")
    for model, _ in [*sources, *targets]:
        for owner in model.drawing_lineage:
            if len(owner.positions) != 1:
                raise UserError("Instanced models are not supported; use a single model copy.")
    results = []
    target_input = None
    cache_context = _ShapeBatchCache() if len(sources) > 1 else nullcontext()
    undo_context = session.undo.aggregate("shape align models") if move and len(sources) > 1 else nullcontext()
    with _ShapeUI(session), cache_context, undo_context as batch, _ShapeMapWorkers(session, workers):
        try:
            for item in sources:
                if item[1] is not None and targets[0][1] is not None:
                    result = _shape_align_one(session, item[1], targets[0][1], sphere_center=center, **options)
                else:
                    started = perf_counter()
                    # Automatic search sizes are per moving input, not inherited
                    # from the first member of a batch.
                    if target_input is None and radius_search is not True:
                        target_input = _shape_reference_input(
                            *targets[0], center, region_radius, confidence_reference or confidence)
                    result = _shape_align_with_map(session, item, targets[0], options, center, target_input, started)
                results.append(result)
        except BaseException:
            if batch is not None:
                for action in reversed(batch.actions):
                    action.undo()
                session.logger.warning("Shape alignment batch failed; earlier changes were undone.")
            raise
    return results[0] if len(results) == 1 else results


def _register_shape_align_command(logger):
    from chimerax.core.commands import (AxisArg, BoolArg, CenterArg, CmdDesc, EnumOf, FloatArg,
                                       IntArg, NoArg, ObjectsArg, Or, register)
    from chimerax.atomic import SymmetryArg

    class SearchRadiusArg(Or):
        # ChimeraX rejects empty keyword values before calling parse(), except
        # for NoArg. Opt into that flag check while retaining numeric parsing.
        def __eq__(self, other):
            return other is NoArg or other is self

        def unparse(self, value, session=None):
            return "" if value is True else FloatArg.unparse(value, session)

    register("shapealign", CmdDesc(
        optional=[("moving", ObjectsArg)],
        keyword=[("to", ObjectsArg), ("resolution", FloatArg),
                 ("max_points", IntArg), ("extra_starts", IntArg), ("move", BoolArg), ("copy", BoolArg),
                 ("radius_search", SearchRadiusArg(FloatArg, NoArg, name="radius in Angstroms or no value (automatic)")),
                 ("radius_refine", FloatArg), ("radius", FloatArg),
                 ("correlation", BoolArg), ("color_fit", BoolArg),
                 ("correlation_resolution", FloatArg), ("correlation_window", FloatArg),
                 ("partial", BoolArg), ("fallback_correlation", FloatArg),
                 ("workers", IntArg),
                 ("tiled", BoolArg), ("bandpass", BoolArg), ("map_resolution", FloatArg), ("tile_size", FloatArg),
                 ("copies", IntArg), ("max_overlap", FloatArg), ("overlap_distance", FloatArg),
                 ("score_fraction", FloatArg), ("symmetry", SymmetryArg),
                 ("symmetry_center", CenterArg), ("symmetry_axis", AxisArg), ("refine_copies", BoolArg),
                 ("confidence", EnumOf(("auto", "off", "plddt"))),
                 ("confidence_moving", EnumOf(("auto", "off", "plddt"))),
                 ("confidence_reference", EnumOf(("auto", "off", "plddt"))),
                 ("objective", EnumOf(("density", "gradient")))],
        synopsis="Align models/maps by shape; radius_search searches a map region, radius_refine refines locally",
    ), shape_align, logger=logger)
    register("shapealign help", CmdDesc(synopsis="Show shapealign usage, options, and examples"),
             _shape_align_help, logger=logger)


# BEGIN INTEGRATED DOMAIN RUNTIME

class _RuntimeNamespace:
    """Expose integrated routines without mixing the backends' public names."""
    def __init__(self, prefix):
        object.__setattr__(self, '_prefix', prefix)
    def __getattr__(self, name):
        try:
            return globals()[self._prefix + name]
        except KeyError:
            raise AttributeError(name) from None
    def __setattr__(self, name, value):
        globals()[self._prefix + name] = value
    def __delattr__(self, name):
        try:
            del globals()[self._prefix + name]
        except KeyError:
            raise AttributeError(name) from None
    def __dir__(self):
        return sorted(k[len(self._prefix):] for k in globals() if k.startswith(self._prefix))


# sequence_features

from collections import Counter as _sf_Counter
from functools import lru_cache as _sf_lru_cache
from hashlib import sha256 as _sf_sha256
import json as _sf_json
import math as _sf_math
from pathlib import Path as _sf_Path
import re as _sf_re
from time import monotonic as _sf_monotonic, time as _sf_time
from urllib.parse import urlencode as _sf_urlencode
from urllib.request import Request as _sf_Request, urlopen as _sf_urlopen
from uuid import uuid4 as _sf_uuid4
_sf_VERSION = 1
_sf_PROSITE_URL = 'https://ftp.expasy.org/databases/prosite/prosite.dat'
_sf_SCAN_URL = 'https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi'
_sf_AMINO_ACIDS = frozenset('ACDEFGHIKLMNPQRSTVWY')
_sf_CACHE_SECONDS = 14 * 86400

def _sf__check(cancel, deadline=None):
    if cancel.is_set():
        raise InterruptedError('Sequence-feature lookup cancelled.')
    if deadline is not None and _sf_monotonic() >= deadline:
        raise TimeoutError('Sequence-feature lookup timed out.')

def _sf__atomic_write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.' + _sf_uuid4().hex + '.tmp')
    try:
        temp.write_bytes(data)
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)

def _sf__request(url, timeout, data=None):
    request = _sf_Request(url, data=data, headers={'User-Agent': 'ChimeraX-Trimmings/1.0', 'Content-Type': 'application/x-www-form-urlencoded'})
    with _sf_urlopen(request, timeout=timeout) as response:
        content = response.read(40 * 1024 * 1024 + 1)
    if len(content) > 40 * 1024 * 1024:
        raise ValueError('Sequence-feature response exceeded the download limit.')
    return content

@_sf_lru_cache(maxsize=16384)
def _sf__lnfact(n):
    return round(_sf_math.lgamma(n + 1), 6) if n <= 10000 else (n + 0.5) * _sf_math.log(n) - n + 0.9189385332

def _sf__seg_probability(counts, length):
    state = sorted((n for n in counts.values() if n), reverse=True)
    classes = _sf_Counter(state + [0] * (20 - len(state)))
    assignments = _sf__lnfact(20) - sum((_sf__lnfact(n) for n in classes.values()))
    permutations = _sf__lnfact(length) - sum((_sf__lnfact(n) for n in state))
    return assignments + permutations - length * 2.995732273553991

def _sf__seg_trim(sequence, cancel):
    size, best, left, right = (len(sequence), 1.0, 0, len(sequence) - 1)
    for length in range(size, max(1, size - 50), -1):
        _sf__check(cancel)
        counts = _sf_Counter((aa for aa in sequence[:length] if aa in _sf_AMINO_ACIDS))
        for start in range(size - length + 1):
            probability = _sf__seg_probability(counts, length)
            if probability < best:
                best, left, right = (probability, start, start + length - 1)
            if start + length < size:
                old, new = (sequence[start], sequence[start + length])
                if old in _sf_AMINO_ACIDS:
                    counts[old] -= 1
                if new in _sf_AMINO_ACIDS:
                    counts[new] += 1
    return (left, right)

def _sf_seg(sequence, cancel):
    """Return inclusive one-based SEG intervals (12, 2.2, 2.5; maxtrim=50).

    Ambiguous residues are excluded from the 20-letter composition; windows
    containing >2 ambiguous residues cannot trigger or extend a segment.
    An explicit stack replaces the reference implementation's left-trim recursion.
    """
    pending, regions = ([(sequence, 0)], [])
    while pending:
        seq, offset = pending.pop()
        if len(seq) < 12:
            continue
        _sf__check(cancel)
        entropy = [-1.0] * len(seq)
        for start in range(len(seq) - 11):
            counts = _sf_Counter((aa for aa in seq[start:start + 12] if aa in _sf_AMINO_ACIDS))
            total = sum(counts.values())
            if total >= 10:
                entropy[start + 5] = abs(sum((n * _sf_math.log(n / total) / _sf_math.log(2) for n in sorted(counts.values(), reverse=True))) / total)
        i, last, lowlim = (5, len(seq) - 7, 5)
        while i <= last:
            if 0 <= entropy[i] <= 2.2:
                lo, hi = (i, i)
                while lo >= lowlim and 0 <= entropy[lo] <= 2.5:
                    lo -= 1
                lo += 1
                while hi <= last and 0 <= entropy[hi] <= 2.5:
                    hi += 1
                hi -= 1
                left, right = (lo - 5, hi + 6)
                trim_left, trim_right = _sf__seg_trim(seq[left:right + 1], cancel)
                left, right = (left + trim_left, left + trim_right)
                if i + 6 < left:
                    pending.append((seq[lo - 5:left], offset + lo - 5))
                regions.append((left + offset + 1, right + offset + 1))
                i = min(hi, right + 5)
                lowlim = i + 1
            i += 1
    return sorted(set(regions))

def _sf_low_complexity(sequence, cancel):
    hits = []
    for start, end in _sf_seg(sequence, cancel):
        composition = _sf_Counter(sequence[start - 1:end])
        summary = ', '.join((f'{aa} {n / (end - start + 1):.0%}' for aa, n in composition.most_common(3)))
        hits.append({'origin': 'sequence', 'feature_kind': 'low_complexity', 'type': 'Low complexity', 'name': 'Low-complexity region', 'accession': 'SEG', 'source': 'SEG', 'start': start, 'end': end, 'evalue': None, 'bitscore': None, 'outcompeted': False, 'entry_url': 'https://www.ncbi.nlm.nih.gov/books/NBK279684/', 'structures_url': '', 'description': 'Composition: ' + summary + '.', 'evidence': 'SEG prediction; window 12, trigger 2.2, extension 2.5, maxtrim 50. Low complexity does not establish intrinsic disorder.'})
    return hits

def _sf_compile_pattern(pattern):
    """Translate PROSITE's constrained pattern language; reject arbitrary regex.

    Character sets containing < or > allow that terminus instead of a residue.
    Positive lookahead preserves overlapping occurrences. Conserved/excluded
    positions never accept an ambiguous residue unless explicitly listed.
    """
    pattern = pattern.strip().removesuffix('.')
    tokens = _sf_re.findall('\\[[A-Z<>]+\\]|\\{[A-Z]+\\}|[A-Zx<>]|\\(\\d+(?:,\\d+)?\\)|-', pattern)
    if ''.join(tokens) != pattern:
        raise ValueError('Unsupported PROSITE pattern: ' + pattern)
    pieces = []
    repeatable = False
    for token in tokens:
        if token == '-':
            repeatable = False
        elif token.startswith('('):
            if not repeatable:
                raise ValueError('Invalid PROSITE repetition.')
            bounds = [int(n) for n in token[1:-1].split(',')]
            if not 0 <= bounds[0] <= bounds[-1] <= 1000:
                raise ValueError('Invalid PROSITE repetition bounds.')
            pieces.append('{' + token[1:-1] + '}')
            repeatable = False
        else:
            repeatable = token not in ('<', '>')
            if token.startswith('['):
                letters = token[1:-1].replace('<', '').replace('>', '')
                alternatives = ['[' + letters + ']'] if letters else []
                alternatives += (['\\A'] if '<' in token else []) + (['\\Z'] if '>' in token else [])
                pieces.append('(?:' + '|'.join(alternatives) + ')')
            elif token.startswith('{'):
                allowed = ''.join(sorted(_sf_AMINO_ACIDS - set(token[1:-1])))
                pieces.append('[' + allowed + ']')
            else:
                pieces.append({'x': '[A-Z]', '<': '\\A', '>': '\\Z'}.get(token, token))
    regex = _sf_re.compile('(?=(' + ''.join(pieces) + '))')
    if regex.match(''):
        raise ValueError('A PROSITE pattern must consume at least one residue.')
    return regex

def _sf_parse_prosite(text):
    release = _sf_re.search('^CC\\s+Release (.+?)\\.$', text, _sf_re.M)
    entries = {}
    for block in text.split('//'):
        identifier = _sf_re.search('^ID\\s+(\\w+); (PATTERN|MATRIX)\\.', block, _sf_re.M)
        if not identifier:
            continue
        accession = _sf_re.search('^AC\\s+(PS\\d{5});', block, _sf_re.M)
        if not accession:
            raise ValueError('Invalid PROSITE accession.')
        acc = accession[1]
        if acc in entries:
            raise ValueError('Duplicate PROSITE accession.')
        pattern = ''.join(_sf_re.findall('^PA\\s+(.+)$', block, _sf_re.M))
        if identifier[2] == 'PATTERN' and (not pattern):
            raise ValueError('Missing PROSITE pattern.')
        entries[acc] = {'accession': acc, 'name': identifier[1], 'kind': identifier[2], 'description': ' '.join(_sf_re.findall('^DE\\s+(.+)$', block, _sf_re.M)), 'common': bool(_sf_re.search('/SKIP-FLAG=TRUE;', block)), 'pattern': pattern, 'prosite_release': release[1] if release else 'local file', 'regex': _sf_compile_pattern(pattern) if pattern else None}
    if not entries:
        raise ValueError('No PROSITE entries found; the response may be an error page.')
    return entries

@_sf_lru_cache(maxsize=2)
def _sf__local_prosite(path, modified, size):
    """Keep compiled definitions between commands; file edits invalidate them."""
    data = _sf_Path(path).read_bytes()
    return (_sf_parse_prosite(data.decode('utf-8')), _sf_sha256(data).hexdigest())

def _sf_load_prosite(cache_dir, cancel, deadline, database=None, refresh=False, report=lambda text: None):
    path = _sf_Path(database).expanduser() if database else _sf_Path(cache_dir) / 'prosite.dat'
    _sf__check(cancel, deadline)
    if database or (path.is_file() and (not refresh)):
        stamp = path.stat()
        entries, digest = _sf__local_prosite(str(path.resolve()), stamp.st_mtime_ns, stamp.st_size)
    else:
        report('downloading PROSITE definitions (once; sequences stay local)')
        _sf__check(cancel, deadline)
        data = _sf__request(_sf_PROSITE_URL, min(60, max(0.1, deadline - _sf_monotonic())))
        entries, digest = (_sf_parse_prosite(data.decode('utf-8')), _sf_sha256(data).hexdigest())
        _sf__check(cancel, deadline)
        _sf__atomic_write(path, data)
    _sf__check(cancel, deadline)
    return (entries, digest)

def _sf_pattern_hits(sequence, entries, common, cancel):
    results = []
    for entry in entries.values():
        _sf__check(cancel)
        if entry['kind'] != 'PATTERN' or (entry['common'] and (not common)):
            continue
        for match in entry['regex'].finditer(sequence):
            start, end = (match.start(1) + 1, match.end(1))
            results.append(_sf__prosite_hit(entry, sequence, start, end, 'pattern'))
    return results

def _sf__prosite_hit(entry, sequence, start, end, kind):
    return {**{key: value for key, value in entry.items() if key not in ('regex', 'kind')}, 'origin': 'prosite', 'feature_kind': 'motif' if kind == 'pattern' else 'profile', 'database': 'prosite' if kind == 'pattern' else 'profile', 'source': 'PROSITE', 'type': 'Pattern match' if kind == 'pattern' else 'Profile-supported region', 'start': start, 'end': end, 'matched_sequence': sequence[start - 1:end], 'evalue': None, 'bitscore': None, 'outcompeted': False, 'entry_url': 'https://prosite.expasy.org/' + entry['accession'], 'structures_url': '', 'evidence': 'Sequence pattern match; function is not established by this match alone.' if kind == 'pattern' else "ScanProsite profile match at the standard level-0 cutoff. A profile region alone does not establish a particular residue's function."}

def _sf_parse_scanprosite(payload, sequence, identifier, entries, common=False):
    try:
        if not isinstance(payload, dict) or type(payload['n_match']) is not int or type(payload['n_seq']) is not int or (payload['n_seq'] not in (0, 1)) or (not isinstance(payload['matchset'], list)) or (payload['n_match'] != len(payload['matchset'])):
            raise ValueError()
        results = []
        for match in payload['matchset']:
            if match.get('sequence_id') != identifier and match.get('sequence_ac') != identifier:
                raise ValueError()
            acc, start, end = (match['signature_ac'], match['start'], match['stop'])
            if acc not in entries or type(start) is not int or type(end) is not int or (not 1 <= start <= end <= len(sequence)):
                raise ValueError()
            entry = entries[acc]
            if entry['common'] and (not common):
                continue
            if entry['kind'] == 'MATRIX':
                score, level = (float(match['score']), int(match['level']))
                if not _sf_math.isfinite(score) or level < 0:
                    raise ValueError()
                hit = _sf__prosite_hit(entry, sequence, start, end, 'profile')
                hit.update(score=score, profile_level=level, source='ScanProsite')
                results.append(hit)
        return results
    except (KeyError, TypeError, ValueError, OverflowError):
        raise ValueError('Invalid ScanProsite results or sequence coordinates; local features were retained.') from None

def _sf_search(sequences, cache_dir, timeout, cancel, report, *, database=None, common=False, profiles=False, refresh=False):
    """Local SEG + PROSITE, optionally enriched online; successful zeroes cache too.

    Local features survive a failed definition download or optional online search.
    Failures are warnings, never cached as valid empty remote responses.
    """
    deadline = _sf_monotonic() + timeout
    cache_dir = _sf_Path(cache_dir)
    results, warnings = ({}, [])
    try:
        entries, digest = _sf_load_prosite(cache_dir, cancel, deadline, database, refresh, report)
    except Exception as error:
        _sf__check(cancel, deadline)
        if database:
            raise
        entries, digest = ({}, 'unavailable')
        warnings.append('PROSITE patterns unavailable; low-complexity features retained: ' + str(error))
    for index, sequence in enumerate(dict.fromkeys(sequences), 1):
        _sf__check(cancel, deadline)
        report(f'local sequence {index}/{len(sequences)}')
        fingerprint = _sf_json.dumps([_sf_VERSION, sequence, digest, common])
        path = cache_dir / 'results' / (_sf_sha256(fingerprint.encode()).hexdigest() + '.json')
        local = _sf__read_cache(path, len(sequence)) if entries and (not refresh) else None
        if local is None:
            local = _sf_low_complexity(sequence, cancel) + _sf_pattern_hits(sequence, entries, common, cancel)
            if entries:
                _sf__atomic_write(path, _sf_json.dumps({'fetched': _sf_time(), 'hits': local}).encode())
        results[sequence] = list(local)
    if profiles and entries:
        for index, sequence in enumerate(results, 1):
            _sf__check(cancel)
            if _sf_monotonic() >= deadline:
                warnings.append('ScanProsite profile time budget exhausted; local features retained.')
                break
            fingerprint = _sf_json.dumps([_sf_VERSION, sequence, digest, common])
            path = cache_dir / 'results' / (_sf_sha256(fingerprint.encode()).hexdigest() + '.json')
            online_path = path.with_name(path.stem + '-profiles.json')
            online = _sf__read_cache(online_path, len(sequence)) if not refresh else None
            if online is None:
                report(f'ScanProsite profiles {index}/{len(sequences)} (sending sequence)')
                identifier = 'trimmings_' + _sf_sha256(sequence.encode()).hexdigest()
                data = _sf_urlencode({'seq': '>' + identifier + '\n' + sequence, 'output': 'json', 'skip': '0' if common else '1', 'lowscore': '0'}).encode()
                try:
                    payload = _sf_json.loads(_sf__request(_sf_SCAN_URL, min(60, max(0.1, deadline - _sf_monotonic())), data))
                    online = _sf_parse_scanprosite(payload, sequence, identifier, entries, common)
                    _sf__check(cancel, deadline)
                    _sf__atomic_write(online_path, _sf_json.dumps({'fetched': _sf_time(), 'hits': online}).encode())
                except Exception as error:
                    _sf__check(cancel)
                    warnings.append('ScanProsite profiles unavailable; local features retained: ' + str(error))
                    online = []
            results[sequence].extend(online)
    return (results, warnings)

def _sf__read_cache(path, sequence_length):
    try:
        payload = _sf_json.loads(path.read_text())
        if 0 <= _sf_time() - payload['fetched'] < _sf_CACHE_SECONDS and isinstance(payload['hits'], list) and all((isinstance(hit, dict) and {'start', 'end', 'accession', 'origin', 'feature_kind', 'name', 'type', 'source', 'outcompeted', 'evalue', 'bitscore'} <= hit.keys() and (type(hit['start']) is int) and (type(hit['end']) is int) and (1 <= hit['start'] <= hit['end'] <= sequence_length) for hit in payload['hits'])):
            return payload['hits']
    except (OSError, ValueError, TypeError, KeyError):
        pass
    return None

_SEQUENCE_FEATURE_RUNTIME = _RuntimeNamespace('_sf_')

# domain_pfam

from functools import lru_cache as _pl_lru_cache
import gzip as _pl_gzip
import hashlib as _pl_hashlib
import json as _pl_json
import math as _pl_math
import os as _pl_os
from pathlib import Path as _pl_Path
import re as _pl_re
import shutil as _pl_shutil
import signal as _pl_signal
import subprocess as _pl_subprocess
import tempfile as _pl_tempfile
import time as _pl_time
_pl_VERSION = 1

def _pl_executable(name):
    return _pl_shutil.which(name) or next((str(path) for path in (_pl_Path('/opt/homebrew/bin') / name, _pl_Path('/usr/local/bin') / name) if path.is_file() and _pl_os.access(path, _pl_os.X_OK)), None)

def _pl_discover(root, database=None, library=None):
    """Find an installed, pressed Pfam bundle; never download or index it here."""
    configured = database or _pl_os.environ.get('CHIMERAX_PFAM_DB')
    paths = []
    if configured:
        path = _pl_Path(configured).expanduser()
        paths = [path / 'Pfam-A.hmm' if path.is_dir() else path]
    else:
        base = _pl_Path(root) / '.domain_data'
        libraries = [_pl_Path(library).expanduser()] if library else []
        settings = base / 'default-library.json'
        if not library and settings.is_file():
            value = _pl_json.loads(settings.read_text()).get('library')
            if isinstance(value, str) and value:
                libraries.append(base / value)
        for folder in libraries:
            manifest_path = folder / 'manifest.json'
            if manifest_path.is_file():
                manifest = _pl_json.loads(manifest_path.read_text())
                if manifest.get('profile_database'):
                    paths.append(folder / manifest['profile_database'])
        candidates = list(base.glob('*/profiles/Pfam*/Pfam-A.hmm'))

        def release_key(path):
            return tuple((int(p) for p in _pl_re.findall('\\d+', path.parent.name)))
        paths.extend(sorted(candidates, key=release_key, reverse=True))
    hmmscan = _pl_executable('hmmscan')
    if hmmscan is None:
        if configured:
            raise ValueError('HMMER hmmscan is not installed or available on PATH.')
        return None
    for database_path in paths:
        database_path = database_path.resolve()
        metadata_path = _pl_Path(str(database_path) + '.dat.gz')
        files = [database_path, metadata_path, *(_pl_Path(str(database_path) + suffix) for suffix in ('.h3f', '.h3i', '.h3m', '.h3p'))]
        if not all((path.is_file() for path in files)):
            continue
        if any((path.stat().st_mtime_ns < database_path.stat().st_mtime_ns for path in files[2:])):
            continue
        version_file = database_path.parent / 'Pfam.version.gz'
        release = None
        if version_file.is_file():
            with _pl_gzip.open(version_file, 'rt') as stream:
                match = _pl_re.search('Pfam release\\s*:\\s*([\\d.]+)', stream.read())
                release = match[1] if match else None
            files.append(version_file)
        if release is None:
            match = _pl_re.fullmatch('Pfam([\\d.]+)', database_path.parent.name)
            release = match[1] if match else 'unknown release'
        programs = {'hmmscan': hmmscan, 'hmmfetch': _pl_executable('hmmfetch'), 'hmmpress': _pl_executable('hmmpress')}
        signature = [_pl_VERSION, programs, [(str(path), path.stat().st_size, path.stat().st_mtime_ns) for path in [*files, *(_pl_Path(p) for p in programs.values() if p)]]]
        identity = 'local-' + _pl_hashlib.sha256(_pl_json.dumps(signature).encode()).hexdigest()
        return {'kind': 'local', 'identity': identity, 'database': str(database_path), 'metadata': str(metadata_path), 'hmmscan': hmmscan, 'release': release, 'hmmfetch': programs['hmmfetch'], 'hmmpress': programs['hmmpress']}
    if configured:
        raise ValueError('The configured Pfam database needs matching .dat.gz metadata and current hmmpress indexes.')
    return None

@_pl_lru_cache(maxsize=4)
def _pl__metadata(path, size, modified):
    families, current = ({}, {})
    with _pl_gzip.open(path, 'rt', encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#=GF '):
                key, value = line[5:].split(None, 1)
                if key == 'NE':
                    current.setdefault('NE', []).append(value.strip())
                else:
                    current[key] = value.strip()
            elif line.strip() == '//':
                accession = current['AC'].split('.')[0]
                ga = [float(v.strip(';')) for v in current['GA'].split()]
                if not _pl_re.fullmatch('PF\\d{5,}', accession) or accession in families or len(ga) != 2 or (not all((_pl_math.isfinite(v) for v in ga))):
                    raise ValueError('Invalid Pfam family metadata.')
                families[accession] = {'accession': accession, 'version': current['AC'], 'ga': ga, 'name': current['ID'], 'description': current.get('DE', ''), 'type': current['TP'], 'length': int(current['ML']), 'clan': current.get('CL'), 'nested': current.get('NE', [])}
                current = {}
    if current or not families:
        raise ValueError('Incomplete Pfam family metadata.')
    return families

def _pl_metadata(backend):
    path = _pl_Path(backend['metadata'])
    stat = path.stat()
    return _pl__metadata(str(path), stat.st_size, stat.st_mtime_ns)

def _pl__run(arguments, output, deadline, cancel, progress):
    if cancel.is_set():
        raise InterruptedError('Pfam search cancelled.')
    if _pl_time.monotonic() >= deadline:
        raise TimeoutError('Local Pfam search timed out.')
    with output.open('w') as log:
        process = _pl_subprocess.Popen([str(a) for a in arguments], stdout=log, stderr=_pl_subprocess.STDOUT, start_new_session=_pl_os.name != 'nt')
        last_report = _pl_time.monotonic()
        try:
            while process.poll() is None:
                if cancel.is_set():
                    raise InterruptedError('Pfam search cancelled.')
                if _pl_time.monotonic() >= deadline:
                    raise TimeoutError('Local Pfam search timed out.')
                if _pl_time.monotonic() - last_report >= 1:
                    progress()
                    last_report = _pl_time.monotonic()
                cancel.wait(0.05)
        finally:
            if process.poll() is None:
                try:
                    if _pl_os.name == 'nt':
                        process.terminate()
                    else:
                        _pl_os.killpg(process.pid, _pl_signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=2)
                except _pl_subprocess.TimeoutExpired:
                    try:
                        if _pl_os.name == 'nt':
                            process.kill()
                        else:
                            _pl_os.killpg(process.pid, _pl_signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
    if process.returncode:
        raise ValueError('Local HMMER failed: ' + output.read_text(errors='replace')[-1500:].strip())

def _pl_parse_table(path, families, queries, release):
    hits = []
    for line in _pl_Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        row = line.split(maxsplit=22)
        if len(row) < 22:
            raise ValueError('Truncated local HMMER domain table.')
        accession = row[1].split('.')[0]
        family, query = (families[accession], queries[row[3]])
        if row[1] != family['version'] or int(row[2]) != family['length'] or int(row[5]) != len(query['sequence']):
            raise ValueError('Local Pfam profiles and metadata do not match.')
        first, last = (int(row[17]), int(row[18]))
        hf, hl = (int(row[15]), int(row[16]))
        seq_score, score, evalue, bias = (float(row[i]) for i in (7, 13, 12, 14))
        if not 1 <= first <= last <= len(query['sequence']) or not 1 <= hf <= hl <= family['length'] or (not all((_pl_math.isfinite(v) for v in (seq_score, score, evalue, bias)))) or (evalue < 0):
            raise ValueError('Invalid local Pfam alignment coordinates or scores.')
        offset = query['start'] - 1
        hits.append({'origin': 'pfam', 'source': 'Pfam local (' + release + ')', 'database': 'pfam', 'accession': accession, 'name': family['name'], 'type': family['type'], 'description': family['description'], 'start': first + offset, 'end': last + offset, 'evalue': evalue, 'bitscore': score, 'outcompeted': False, 'tentative': seq_score < family['ga'][0] or score < family['ga'][1], 'gathering_cutoff': family['ga'][1], 'sequence_gathering_cutoff': family['ga'][0], 'sequence_score': seq_score, 'profile_coverage': (hl - hf + 1) / family['length'], 'composition_bias': bias, 'pfam_backend': 'local', 'pfam_release': release, 'query_region': [query['start'], query['start'] + len(query['sequence']) - 1], '_owner': query['owner']})
    return hits

def _pl_resolve(hits, families):
    """Deduplicate retries, then apply same-clan competition with nested exceptions."""
    selected = []
    for original in sorted(hits, key=lambda h: (h['tentative'], h['evalue'], -h['bitscore'])):
        hit = {**original, 'outcompeted': False}

        def overlap(other):
            return max(0, min(hit['end'], other['end']) - max(hit['start'], other['start']) + 1)
        if any((old['accession'] == hit['accession'] and overlap(old) >= 0.5 * min(old['end'] - old['start'] + 1, hit['end'] - hit['start'] + 1) for old in selected)):
            continue
        family = families[hit['accession']]
        for old in selected:
            previous = families[old['accession']]
            if not old['outcompeted'] and old['accession'] != hit['accession'] and family['clan'] and (family['clan'] == previous['clan']) and overlap(old) and (old['name'] not in family['nested']) and (hit['name'] not in previous['nested']):
                hit['outcompeted'] = True
                break
        selected.append(hit)
    return selected

def _pl_retry_queries(sequence, hits, owner):
    curated = {p for h in hits if not h['tentative'] and (not h['outcompeted']) for p in range(h['start'], h['end'] + 1)}
    regions = set()
    for hit in hits:
        if hit['tentative'] and (not hit['outcompeted']):
            regions.add((max(1, hit['start'] - 12), min(len(sequence), hit['end'] + 12)))
    first = None
    for p in range(1, len(sequence) + 2):
        if p <= len(sequence) and p not in curated:
            if first is None:
                first = p
        elif first is not None:
            if p - first >= 20:
                regions.add((max(1, first - 12), min(len(sequence), p + 11)))
            first = None
    return [{'sequence': sequence[a - 1:b], 'start': a, 'owner': owner} for a, b in sorted(regions) if (a, b) != (1, len(sequence))]

def _pl_search(sequences, backend, timeout, cancel, report):
    """Curated full scan, reported tentative hits, and ordinary targeted retries."""
    if not _pl_math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Local Pfam timeout must be finite and positive.')
    if cancel.is_set():
        raise InterruptedError('Pfam search cancelled.')
    deadline = _pl_time.monotonic() + timeout
    sequences = list(dict.fromkeys(sequences))
    if not sequences:
        return {}
    families = _pl_metadata(backend)
    full = {f'chain{i}': {'sequence': seq, 'start': 1, 'owner': i} for i, seq in enumerate(sequences)}
    release = backend['release']
    with _pl_tempfile.TemporaryDirectory(prefix='trimmings-pfam-') as temporary:
        folder = _pl_Path(temporary)

        def scan(database, queries, name, curated, expected_models):
            report('Local Pfam ' + release + ': ' + name)
            fasta, table, log = (folder / (name + suffix) for suffix in ('.fasta', '.tsv', '.log'))
            fasta.write_text(''.join((f">{key}\n{q['sequence']}\n" for key, q in queries.items())))
            options = ['--cut_ga'] if curated else ['-E', '10', '--domE', '10']
            _pl__run([backend['hmmscan'], '--noali', '--cpu', '4', '-Z', len(families), '--domtblout', table, *options, database, fasta], log, deadline, cancel, lambda: report('Local Pfam ' + release + ': ' + name))
            counts = _pl_re.findall('Target model\\(s\\):\\s+(\\d+)', log.read_text())
            if len(counts) != len(queries) or any((int(n) != expected_models for n in counts)):
                raise ValueError('Local Pfam database inventory does not match its metadata.')
            hits = _pl_parse_table(table, families, queries, release)
            if curated and any((h['tentative'] for h in hits)):
                raise ValueError('A local curated match failed its recorded gathering thresholds.')
            return hits
        all_hits = scan(backend['database'], full, 'curated scan', True, len(families))
        reported = scan(backend['database'], full, 'tentative scan', False, len(families))
        all_hits.extend((h for h in reported if h['tentative'] and h['evalue'] <= 0.01))
        if backend.get('hmmfetch') and backend.get('hmmpress'):
            candidates_by_owner = {owner: {h['accession'] for h in all_hits if h['_owner'] == owner} for owner in range(len(sequences))}
            candidates = set().union(*candidates_by_owner.values())
            regions = {}
            for owner, seq in enumerate(sequences):
                hits = _pl_resolve([h for h in all_hits if h['_owner'] == owner], families)
                for query in _pl_retry_queries(seq, hits, owner):
                    regions[f'region{len(regions)}'] = query
            if candidates and regions:
                keys, subset = (folder / 'keys', folder / 'candidates.hmm')
                keys.write_text(''.join((families[a]['version'] + '\n' for a in sorted(candidates))))
                report('Local Pfam ' + release + ': preparing region retries')
                _pl__run([backend['hmmfetch'], '-o', subset, '-f', backend['database'], keys], folder / 'fetch.log', deadline, cancel, lambda: report('Preparing local Pfam profiles'))
                _pl__run([backend['hmmpress'], '-f', subset], folder / 'press.log', deadline, cancel, lambda: report('Indexing local Pfam profiles'))
                retries = scan(subset, regions, 'padded region retries', True, len(candidates))
                all_hits.extend((h for h in retries if h['accession'] in candidates_by_owner[h['_owner']]))
        results = {}
        for owner, seq in enumerate(sequences):
            results[seq] = [{key: value for key, value in h.items() if key != '_owner'} for h in _pl_resolve([h for h in all_hits if h['_owner'] == owner], families)]
        return results

_PFAM_LOCAL_RUNTIME = _RuntimeNamespace('_pl_')

# domain_exemplars

import json as _ex_json
import math as _ex_math
from pathlib import Path as _ex_Path
import re as _ex_re
import tempfile as _ex_tempfile
import time as _ex_time
from types import SimpleNamespace as _ex_SimpleNamespace
from urllib.request import Request as _ex_Request, urlopen as _ex_urlopen
from urllib.error import HTTPError as _ex_HTTPError, URLError as _ex_URLError

def _ex_source_spec(record, structure):
    """Return only explicit source identities and declared domain boundaries."""
    match = structure.get('domain_match', {})
    if match.get('locations'):
        return {'provider': 'pdb', 'id': structure['accession'].lower(), 'database': match['database'], 'accession': match['accession'], 'regions': match.get('regions', []), 'mapped': True}
    reference = record.get('structural_match', {}).get('reference', {})
    same_url = bool(structure.get('url')) and structure['url'] == reference.get('structure_url')
    if not same_url and structure.get('accession') not in (reference.get('pdb'), reference.get('ecod_domain'), reference.get('id'), reference.get('uniprot')):
        return None

    def source(provider, identifier, **extra):
        spec = {'provider': provider, 'id': identifier, 'regions': regions, 'mapped': False, **extra}
        uid = reference.get('ecod_uid')
        if reference.get('catalog') == 'ECOD' and isinstance(uid, int) and (uid > 0):
            spec['ecod_uid'] = uid
        return spec

    def domain_only():
        uid = reference.get('ecod_uid')
        return source('ecod', str(uid), domain_only=True) if reference.get('catalog') == 'ECOD' and isinstance(uid, int) and (uid > 0) else None
    predicted = None
    if reference.get('catalog') == 'ECOD' and reference.get('structure_source') == 'predicted_or_other':
        predicted = _ex_re.fullmatch('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})_F1_nD[1-9][0-9]*', reference.get('ecod_domain') or '')
    chain_id = reference.get('chain') or ('A' if predicted else None)
    regions = []
    text = reference.get('ecod_seqid_range', '')
    if text:
        pieces = {}
        for item in text.split(','):
            parsed = _ex_re.fullmatch('(?:(.+):)?(\\d+)-(\\d+)', item.strip())
            if not parsed or int(parsed[2]) > int(parsed[3]):
                return domain_only()
            chain = parsed[1] or chain_id
            if not chain:
                return domain_only()
            pieces.setdefault(chain, []).append({'start': int(parsed[2]), 'end': int(parsed[3])})
        regions = [{'chain': chain, 'segments': segments} for chain, segments in pieces.items()]
    elif reference.get('segments') and reference.get('chain') and (not reference.get('parent_target')):
        regions = [{'chain': reference['chain'], 'segments': [{'start': a, 'end': b} for a, b in reference['segments']]}]
    if not regions:
        return domain_only()
    if reference.get('pdb'):
        return source('pdb', reference['pdb'].lower())
    if reference.get('structure_source') == 'alphafold' and reference.get('uniprot'):
        return source('alphafold', reference['uniprot'])
    if predicted:
        return source('alphafold', predicted[1], fragment=1)
    return domain_only()

def _ex__coordinate_file_valid(path):
    with _ex_Path(path).open('rb') as stream:
        data = stream.read(65536).lstrip()
    if _ex_Path(path).suffix == '.pdb':
        return data.startswith(b'ATOM  ') or b'\nATOM  ' in data
    return data.startswith(b'data_')

def _ex__download(url, path, progress):
    request = _ex_Request(url, headers={'User-Agent': 'ChimeraX-Trimmings/1.0'})
    temporary = None
    try:
        with _ex_urlopen(request, timeout=45) as response, _ex_tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            temporary = _ex_Path(stream.name)
            total, size, updated = (int(response.headers.get('Content-Length', 0)), 0, 0)
            while (block := response.read(1024 * 1024)):
                size += len(block)
                if size > 1024 * 1024 * 1024:
                    raise ValueError('Structure download exceeds 1 GB.')
                stream.write(block)
                if progress and _ex_time.monotonic() - updated > 1:
                    progress(f'Downloading {size / 1048576:.1f}' + (f'/{total / 1048576:.1f}' if total else '') + ' MB')
                    updated = _ex_time.monotonic()
        if total and size != total:
            raise ValueError('Incomplete structure download.')
        with temporary.open('rb') as stream:
            data = stream.read(65536).lstrip()
        valid = data.startswith(b'ATOM  ') or b'\nATOM  ' in data if path.suffix == '.pdb' else data.startswith(b'data_')
        if not valid:
            raise ValueError('The source did not return a PDB structure.' if path.suffix == '.pdb' else 'The source did not return an mmCIF structure.')
        temporary.replace(path)
    finally:
        if temporary:
            temporary.unlink(missing_ok=True)

def _ex_fetch_source(spec, cache, progress=None):
    """Keep the exact ECOD domain available if its optional parent cannot be fetched."""
    try:
        return _ex__fetch_source(spec, cache, progress)
    except (_ex_HTTPError, _ex_URLError, TimeoutError, ValueError) as error:
        uid = spec.get('ecod_uid')
        if spec.get('provider') == 'ecod' or not isinstance(uid, int) or uid <= 0 or (isinstance(error, ValueError) and str(error) != 'No full AlphaFold model is available for this reference.'):
            raise
        if progress:
            progress('Full source unavailable; retrieving the exact ECOD domain')
        fallback = {**spec, 'provider': 'ecod', 'domain_only': True, 'full_source_error': str(error)}
        try:
            return _ex__fetch_source(fallback, cache, progress)
        except (_ex_HTTPError, _ex_URLError, TimeoutError) as fallback_error:
            raise ValueError(f'Could not retrieve the full source ({error}) or the ECOD domain ({fallback_error}). Try again later.') from fallback_error

def _ex__fetch_source(spec, cache, progress=None):
    """Network/cache work only; safe to run outside the GUI thread."""
    spec = dict(spec)
    cache = _ex_Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    provider, identifier = (spec['provider'], spec['id'])
    if provider == 'ecod':
        uid = spec.get('ecod_uid')
        if not isinstance(uid, int) or uid <= 0:
            raise ValueError('Invalid ECOD domain UID.')
        path = cache / f'ecod-{uid}.pdb'
        if not (path.is_file() and 0 <= _ex_time.time() - path.stat().st_mtime < 14 * 86400 and _ex__coordinate_file_valid(path)):
            _ex__download(f'http://prodata.swmed.edu/ecod/api/domain/{uid}/pdb', path, progress)
        return (path, spec)
    if provider == 'pdb':
        if not _ex_re.fullmatch('[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+', identifier):
            raise ValueError('Invalid PDB accession.')
        url = f'https://files.rcsb.org/download/{identifier.upper()}.cif'
    elif provider == 'alphafold':
        if not _ex_re.fullmatch('[A-Z0-9]{6,10}', identifier):
            raise ValueError('Invalid AlphaFold accession.')
        with _ex_urlopen(_ex_Request(f'https://alphafold.ebi.ac.uk/api/prediction/{identifier}', headers={'User-Agent': 'ChimeraX-Trimmings/1.0'}), timeout=30) as response:
            predictions = _ex_json.load(response)
        first = min((p['start'] for r in spec['regions'] for p in r['segments']))
        last = max((p['end'] for r in spec['regions'] for p in r['segments']))
        expected = f"AF-{identifier}-F{spec['fragment']}" if spec.get('fragment') else None
        model = next((m for m in predictions if m.get('uniprotAccession') == identifier and m.get('cifUrl') and (not m.get('isComplex')) and (expected is None or m.get('modelEntityId', m.get('entryId')) == expected) and (int(m.get('uniprotStart', m.get('sequenceStart', 1))) <= first) and (int(m.get('uniprotEnd', m.get('sequenceEnd', 0))) >= last)), None)
        if model is None:
            raise ValueError('No full AlphaFold model is available for this reference.')
        url = model['cifUrl']
        if not url.startswith('https://alphafold.ebi.ac.uk/'):
            raise ValueError('Unexpected AlphaFold coordinate source.')
        spec['sequence_offset'] = int(model.get('uniprotStart', model.get('sequenceStart', 1))) - 1
        if model.get('chainId'):
            spec['regions'] = [{**region, 'chain': model['chainId']} for region in spec['regions']]
    else:
        raise ValueError('Unsupported structure source.')
    suffix = f"-F{spec['fragment']}" if spec.get('fragment') else ''
    if provider == 'alphafold':
        suffix += '-v' + str(model.get('latestVersion', 'current'))
    path = cache / f'{provider}-{identifier}{suffix}.cif'
    usable = path.is_file() and 0 <= _ex_time.time() - path.stat().st_mtime < 14 * 86400
    if usable:
        usable = _ex__coordinate_file_valid(path)
    if not usable:
        _ex__download(url, path, progress)
    if spec.get('ecod_uid'):
        uid = spec['ecod_uid']
        domain_path = cache / f'ecod-{uid}.pdb'
        usable = domain_path.is_file() and 0 <= _ex_time.time() - domain_path.stat().st_mtime < 14 * 86400 and _ex__coordinate_file_valid(domain_path)
        if not usable:
            _ex__download(f'http://prodata.swmed.edu/ecod/api/domain/{uid}/pdb', domain_path, progress)
        spec['domain_path'] = str(domain_path)
    return (path, spec)

def _ex_region_residues(model, region, offset=0):
    from chimerax.atomic import Residues
    candidates = [c for c in model.chains if c.chain_id == region['chain']]
    if not candidates:
        candidates = [c for c in model.chains if any((r.mmcif_chain_id == region['chain'] for r in c.existing_residues))]
    selected = []
    for chain in candidates:
        for segment in region['segments']:
            a, b = (segment['start'] - offset, segment['end'] - offset)
            if not 1 <= a <= b <= len(chain.residues):
                continue
            members = [r for r in chain.residues[a - 1:b] if r is not None]
            for index, field in [(a - 1, 'author_start'), (b - 1, 'author_end')]:
                expected = segment.get(field)
                residue = chain.residues[index]
                if expected is not None and residue is not None and (str(residue.number) + residue.insertion_code.strip() != str(expected)):
                    raise ValueError('The exemplar residue numbering disagrees with its domain mapping.')
            selected.extend(members)
    return Residues(list(dict.fromkeys(selected)))

def _ex_copy_region(session, model, residues, name):
    """Copy only selected atoms/bonds, avoiding a whole-complex temporary copy."""
    from chimerax.atomic import AtomicStructure
    result = AtomicStructure(session, name=name, auto_style=False, log_info=False)
    atom_map = {}
    try:
        for residue in residues:
            new = result.new_residue(residue.name, residue.chain_id, residue.number, residue.insertion_code)
            new.ss_type, new.ss_id = (residue.ss_type, residue.ss_id)
            new.ribbon_color = residue.ribbon_color
            for atom in residue.atoms:
                copied = result.new_atom(atom.name, atom.element)
                new.add_atom(copied)
                copied.coord, copied.bfactor, copied.occupancy = (atom.coord, atom.bfactor, atom.occupancy)
                copied.color = atom.color
                atom_map[atom] = copied
        for bond in residues.atoms.intra_bonds:
            a, b = bond.atoms
            result.new_bond(atom_map[a], atom_map[b])
        result.scene_position = model.scene_position
        result.atoms.displays = False
        result.residues.ribbon_displays = True
        return result
    except BaseException:
        result.delete()
        raise

def _ex_align_full(session, record, full, spec):
    from chimerax.match_maker import match, defaults, CP_SPECIFIC_SPECIFIC
    baseline, best = (full.scene_position, None)
    failures = []
    for region in spec['regions']:
        full.scene_position = baseline
        try:
            selected = _ex_region_residues(full, region, spec.get('sequence_offset', 0))
            if len(selected) < 3:
                continue
            chains = selected.unique_chains
            if len(chains) != 1:
                raise ValueError('A domain occurrence must map to one polymer chain.')
            fits = match(session, CP_SPECIFIC_SPECIFIC, [(record['chain'], chains[0])], defaults['matrix'], defaults['alignment_algorithm'], defaults['gap_open'], defaults['gap_extend'], domain_residues=(record['residues'], selected), cutoff_distance=defaults['iter_cutoff'], always_raise_errors=True, log_parameters=False, verbose=None, compute_ss=False)
            if not fits:
                continue
            fit = fits[0]
            count, rmsd = (len(fit['final match atoms']), fit['final RMSD'])
            if count < 3 or not _ex_math.isfinite(rmsd):
                continue
            rank = (-count, rmsd)
            if best is None or rank < best[0]:
                best = (rank, full.scene_position, selected, count, rmsd)
        except Exception as error:
            failures.append(str(error))
    if best is None:
        full.scene_position = baseline
        raise ValueError('No mapped occurrence could be aligned: ' + (failures[0] if failures else 'fewer than three matching residues.'))
    _, full.scene_position, selected, count, rmsd = best
    return (selected, count, rmsd)

def _ex_open_aligned(session, record, spec, path):
    from chimerax.core.models import Model
    from chimerax.mmcif import open_mmcif
    group = Model(f"{spec['id'].upper()} aligned to {record['chain_spec']} {record['name']}", session)
    models = []
    try:
        full = None
        if not spec.get('domain_only'):
            models, _ = open_mmcif(session, str(path), auto_style=False, log_info=False, max_models=1)
            if not models or models[0].num_atoms == 0:
                raise ValueError('The coordinate file contains no atoms.')
            full = models[0]
            full.name = f"{spec['id'].upper()} full context"
            full.display = False
        session.models.add([group])
        if models:
            session.models.add(models, parent=group)
        if full is not None and (not full.ss_assigned):
            from chimerax.dssp import compute_ss
            compute_ss(full)
        if spec.get('domain_path') or spec.get('domain_only'):
            from chimerax.pdb import open_pdb
            domains, _ = open_pdb(session, str(path) if spec.get('domain_only') else spec['domain_path'], auto_style=False, log_info=False, max_models=1)
            if not domains or not domains[0].num_atoms:
                raise ValueError('ECOD returned an empty domain structure.')
            domain = domains[0]
            domain.display = False
            session.models.add(domains, parent=group)
            if not domain.ss_assigned:
                from chimerax.dssp import compute_ss
                compute_ss(domain)
            domain_spec = {'regions': [{'chain': c.chain_id, 'segments': [{'start': 1, 'end': len(c.residues)}]} for c in domain.chains if len(c.residues)]}
            members, count, rmsd = _ex_align_full(session, record, domain, domain_spec)
            context_rmsd = 0
            if full is not None:
                selected, _, context_rmsd = _ex_align_full(session, {'chain': members.unique_chains[0], 'residues': members}, full, spec)
                _ex_style_reference(full, selected)
            _ex_style_reference(domain, members)
            crop = _ex_copy_region(session, domain, members, f"ECOD {spec['ecod_uid']} aligned domain")
            session.models.add([crop], parent=group)
            crop.scene_position = domain.scene_position
            session.models.close(domains)
            relation = full.scene_position.inverse() * crop.scene_position if full is not None else None
        else:
            selected, count, rmsd = _ex_align_full(session, record, full, spec)
            _ex_style_reference(full, selected)
            crop = _ex_copy_region(session, full, selected, f"{spec['id'].upper()} aligned domain")
            session.models.add([crop], parent=group)
            crop.scene_position = full.scene_position
            relation, context_rmsd = (None, 0)
        return _ex_SimpleNamespace(group=group, domain=crop, full=full, count=count, rmsd=rmsd, source=spec, context_relation=relation, context_rmsd=context_rmsd)
    except BaseException:
        if group.id is not None:
            session.models.close([group])
        else:
            for model in models:
                model.delete()
            group.delete()
        raise

def _ex_style_reference(full, selected):
    import numpy as np
    full.atoms.displays = False
    full.residues.ribbon_displays = True
    full.atoms.colors = np.tile(np.array([180, 180, 180, 255], dtype=np.uint8), (len(full.atoms), 1))
    full.residues.ribbon_colors = np.tile(np.array([180, 180, 180, 255], dtype=np.uint8), (len(full.residues), 1))
    selected.atoms.colors = np.tile(np.array([255, 190, 70, 255], dtype=np.uint8), (len(selected.atoms), 1))
    selected.ribbon_colors = np.tile(np.array([255, 190, 70, 255], dtype=np.uint8), (len(selected), 1))

def _ex_realign_pair(session, record, pair):
    relation = getattr(pair, 'context_relation', None)
    if relation is not None or pair.full is None:
        original = pair.domain.scene_position
        try:
            regions = [{'chain': c.chain_id, 'segments': [{'start': 1, 'end': len(c.residues)}]} for c in pair.domain.chains if len(c.residues)]
            _, count, rmsd = _ex_align_full(session, record, pair.domain, {'regions': regions})
            if pair.full is not None:
                pair.full.scene_position = pair.domain.scene_position * relation.inverse()
            pair.count, pair.rmsd = (count, rmsd)
        except BaseException:
            pair.domain.scene_position = original
            raise
        return
    showing, original = (pair.full.display, pair.full.scene_position)
    crop = None
    pair.full.display = False
    try:
        selected, count, rmsd = _ex_align_full(session, record, pair.full, pair.source)
        _ex_style_reference(pair.full, selected)
        crop = _ex_copy_region(session, pair.full, selected, pair.domain.name)
        session.models.add([crop], parent=pair.group)
        crop.scene_position = pair.full.scene_position
        session.models.close([pair.domain])
        pair.domain, pair.count, pair.rmsd = (crop, count, rmsd)
        pair.full.display, pair.domain.display = (showing, not showing)
    except BaseException:
        if crop is not None and crop is not pair.domain and (not crop.deleted):
            if crop.id is not None:
                session.models.close([crop])
            else:
                crop.delete()
        pair.full.scene_position, pair.full.display = (original, showing)
        raise

def _ex_show_context(pair, show):
    if pair.full is None:
        raise ValueError('ECOD provides this domain, but its full source structure is not identified.')
    relation = getattr(pair, 'context_relation', None)
    if show:
        pair.full.scene_position = pair.domain.scene_position * relation.inverse() if relation is not None else pair.domain.scene_position
    else:
        pair.domain.scene_position = pair.full.scene_position * relation if relation is not None else pair.full.scene_position
    pair.full.display, pair.domain.display = (bool(show), not show)

def _ex_focus_pair(session, pair):
    """Frame the visible exemplar, updating clipping after a context switch."""
    model = pair.full if pair.full is not None and pair.full.display else pair.domain
    bounds = model.atoms.scene_bounds
    if bounds is None:
        return
    view, center = (session.main_view, bounds.center())
    view.view_all(bounds, pad=0.15)
    normal, depth = (view.camera.view_direction(), max(bounds.radius(), 1.0))
    view.clip_planes.set_clip_position('near', center - depth * normal, view)
    view.clip_planes.set_clip_position('far', center + depth * normal, view)
    view.center_of_rotation = center

_DOMAIN_EXEMPLAR_RUNTIME = _RuntimeNamespace('_ex_')

# domain_foldseek

import json as _fs_json
import math as _fs_math
import os as _fs_os
from pathlib import Path as _fs_Path
import re as _fs_re
import signal as _fs_signal
import subprocess as _fs_subprocess
import tempfile as _fs_tempfile
import time as _fs_time
_fs_FIELDS = ('query', 'target', 'evalue', 'bits', 'qstart', 'qend', 'qlen', 'tstart', 'tend', 'tlen', 'qcov', 'tcov', 'qtmscore', 'ttmscore', 'alntmscore', 'qaln', 'taln')
_fs_SCREEN_FIELDS = tuple((field for field in _fs_FIELDS if field not in ('qtmscore', 'ttmscore')))
_fs_ALIGNMENT_FIELDS = tuple((field for field in _fs_SCREEN_FIELDS if field != 'alntmscore'))
_fs_DEFAULT_THREADS = 4
_fs_SEARCH_VERSION = 12

def _fs_parse_segments(text):
    if not _fs_re.fullmatch('\\d+-\\d+(?:_\\d+-\\d+)*(?:,\\d+-\\d+(?:_\\d+-\\d+)*)*', text):
        raise ValueError('Use comma-separated domains, with underscore-separated pieces for discontinuous domains.')
    domains = [[tuple(map(int, pair.split('-'))) for pair in domain.split('_')] for domain in text.split(',')]
    if any((a < 1 or b < a for domain in domains for a, b in domain)):
        raise ValueError('Domain ranges must be positive and increasing.')
    if any((b >= c for domain in domains for (_, b), (c, _) in zip(domain, domain[1:]))):
        raise ValueError('Pieces within one domain must be ordered and non-overlapping.')
    return domains

def _fs_parse_hits(text, queries, manifest, fields=_fs_FIELDS):
    targets = {record['id']: record for record in manifest['records']}
    hits = []
    for line in text.splitlines():
        values = line.split('\t')
        if len(values) != len(fields):
            raise ValueError('Invalid Foldseek output columns.')
        row = dict(zip(fields, values))
        query, target = (row['query'].removesuffix('.pdb'), row['target'].removesuffix('.pdb'))
        if query not in queries or target not in targets:
            raise ValueError('Foldseek returned an unknown query or target identifier.')
        for field in ('evalue', 'bits', 'qcov', 'tcov', 'qtmscore', 'ttmscore', 'alntmscore'):
            if field not in fields and field in ('qtmscore', 'ttmscore'):
                row[field] = None
                continue
            row[field] = float(row[field])
            if not _fs_math.isfinite(row[field]) or row[field] < 0:
                raise ValueError('Invalid Foldseek score.')
        for field in ('qcov', 'tcov', 'qtmscore', 'ttmscore'):
            if row[field] is not None and row[field] > 1:
                raise ValueError(f"Invalid normalized Foldseek score: {field}={row[field]:g} for {queries[query].get('label', query)} against {target}.")
        for field in ('qstart', 'qend', 'qlen', 'tstart', 'tend', 'tlen'):
            row[field] = int(row[field])
        target_length = targets[target].get('length')
        if target_length is None:
            target_length = len(targets[target]['sequence'])
        if row['qlen'] != len(queries[query]['positions']) or not 1 <= row['qstart'] <= row['qend'] <= row['qlen'] or row['tlen'] != target_length or (not 1 <= row['tstart'] <= row['tend'] <= row['tlen']) or (len(row['qaln']) != len(row['taln'])) or (len(row['qaln'].replace('-', '')) != row['qend'] - row['qstart'] + 1) or (len(row['taln'].replace('-', '')) != row['tend'] - row['tstart'] + 1):
            raise ValueError('Invalid Foldseek alignment coordinates.')
        q, t, aligned = (row['qstart'] - 1, row['tstart'] - 1, [])
        for aa, bb in zip(row['qaln'], row['taln']):
            if aa != '-' and bb != '-':
                aligned.append(queries[query]['positions'][q])
            q += aa != '-'
            t += bb != '-'
        if q != row['qend'] or t != row['tend']:
            raise ValueError('Foldseek alignment strings do not match their endpoints.')
        if row['alntmscore'] > 1:
            span = min(row['qend'] - row['qstart'] + 1, row['tend'] - row['tstart'] + 1)
            if span <= 1 or len(aligned) != span or row['alntmscore'] > span / (span - 1) + 0.0005:
                raise ValueError(f"Invalid normalized Foldseek score: alntmscore={row['alntmscore']:g} for {queries[query].get('label', query)} against {target}.")
            row['alntmscore_raw'] = row['alntmscore']
            row['alntmscore'] = min(1.0, row['alntmscore'] * (span - 1) / span)
        hits.append({**row, 'query': query, 'target': target, 'reference': targets[target], 'aligned_positions': aligned})
    return sorted(hits, key=lambda hit: (hit['evalue'], -hit['bits']))

def _fs_qualifying_hit(hit):
    return hit['evalue'] <= 0.01 and hit['tcov'] >= 0.5 and (hit['alntmscore'] >= 0.5) and (len(hit['aligned_positions']) >= 20)

def _fs_select_matches(hits, per_family=True):
    """Keep distinct occurrences; collapse competing representatives of one occurrence."""
    selected, footprints = ([], {})
    for hit in sorted(hits, key=lambda h: (h['evalue'], -h['bits'])):
        if not _fs_qualifying_hit(hit):
            continue
        key = (hit['query'], hit['reference']['accession'] if per_family else None)
        positions = set(hit['aligned_positions'])
        if any((len(positions & previous) / min(len(positions), len(previous)) >= 0.5 for previous in footprints.get(key, []))):
            continue
        footprints.setdefault(key, []).append(positions)
        selected.append(hit)
    return selected

def _fs_library_reference_count(manifest):
    return manifest['reference_count'] if manifest.get('schema') == 2 else len(manifest['records'])

def _fs_library_file(library, name):
    path = (_fs_Path(library) / name).resolve()
    if not path.is_relative_to(_fs_Path(library).resolve()):
        raise ValueError('Reference-library files must be inside the library directory.')
    return path

def _fs_hit_metadata(library, manifest, text, fields=_fs_FIELDS):
    if manifest['schema'] == 1:
        return manifest
    import sqlite3
    from contextlib import closing
    targets = set()
    for line in text.splitlines():
        values = line.split('\t')
        if len(values) != len(fields):
            raise ValueError('Invalid Foldseek output columns.')
        targets.add(values[1].removesuffix('.pdb'))
    path = _fs_library_file(library, manifest['reference_metadata'])
    records = []
    with closing(sqlite3.connect(f'file:{path}?mode=ro', uri=True)) as connection:
        identifiers = sorted(targets)
        for first in range(0, len(identifiers), 500):
            batch = identifiers[first:first + 500]
            rows = connection.execute('SELECT record FROM refs WHERE id IN (' + ','.join(('?' for _ in batch)) + ')', batch)
            records.extend((_fs_json.loads(row[0]) for row in rows))
    return {**manifest, 'records': records}

def _fs__run_foldseek(command, log_path, timeout, cancel=None, progress=None):
    """Run one stage with cancellation, progress and a bounded lifetime."""
    if timeout <= 0:
        raise TimeoutError('Foldseek search timed out.')
    if cancel is not None and cancel.is_set():
        raise InterruptedError('Foldseek search cancelled.')
    started = _fs_time.monotonic()
    next_progress = started + 1
    with _fs_Path(log_path).open('a+') as log:
        process = _fs_subprocess.Popen(list(map(str, command)), stdout=log, stderr=_fs_subprocess.STDOUT, start_new_session=True)
        try:
            while process.poll() is None:
                if cancel is not None and cancel.is_set():
                    raise InterruptedError('Foldseek search cancelled.')
                if _fs_time.monotonic() - started > timeout:
                    raise TimeoutError('Foldseek search timed out.')
                if progress is not None and _fs_time.monotonic() >= next_progress:
                    progress()
                    next_progress = _fs_time.monotonic() + 1
                _fs_time.sleep(0.05)
            if process.returncode:
                log.seek(0)
                raise RuntimeError('Foldseek failed: ' + log.read()[-2500:])
        finally:
            if process.poll() is None:
                _fs_os.killpg(process.pid, _fs_signal.SIGTERM)
                try:
                    process.wait(timeout=3)
                except _fs_subprocess.TimeoutExpired:
                    _fs_os.killpg(process.pid, _fs_signal.SIGKILL)
                    process.wait()

def _fs__parallel_scoring(jobs, directory, deadline, cancel, progress):
    """Score independent alignment shards without exceeding the CPU budget."""
    from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
    from threading import Event
    from types import SimpleNamespace
    stopping = Event()
    token = SimpleNamespace(is_set=lambda: stopping.is_set() or (cancel is not None and cancel.is_set()))

    def worker(index, commands):
        for command in commands:
            _fs__run_foldseek(command, directory / f'score-{index}.log', deadline - _fs_time.monotonic(), token)
    pool = ThreadPoolExecutor(max_workers=len(jobs))
    pending = {pool.submit(worker, index, commands) for index, commands in enumerate(jobs)}
    next_progress = _fs_time.monotonic() + 1
    try:
        while pending:
            if token.is_set():
                raise InterruptedError('Foldseek search cancelled.')
            if _fs_time.monotonic() >= deadline:
                raise TimeoutError('Foldseek search timed out.')
            ready, pending = wait(pending, timeout=0.05, return_when=FIRST_COMPLETED)
            for future in ready:
                future.result()
            if progress is not None and _fs_time.monotonic() >= next_progress:
                progress()
                next_progress = _fs_time.monotonic() + 1
    finally:
        stopping.set()
        for future in pending:
            future.cancel()
        pool.shutdown(wait=True, cancel_futures=True)

def _fs__screen_alignments(query_db, database, alignment_db, result, directory, threads, deadline, cancel, progress, executable, convert):
    index = _fs_Path(str(alignment_db) + '.index')
    size = sum((int(line.split()[2]) for line in index.read_text().splitlines())) if index.is_file() else 0
    workers = min(threads, 8, max(1, size // 16384))
    if workers == 1:
        convert(alignment_db, result, _fs_SCREEN_FIELDS)
        return result.read_text()
    ordering = directory / 'ordering.tsv'
    convert(alignment_db, ordering, _fs_ALIGNMENT_FIELDS)
    rows = [line.split('\t') for line in ordering.read_text().splitlines()]
    if any((len(row) != len(_fs_ALIGNMENT_FIELDS) for row in rows)):
        raise ValueError('Invalid Foldseek alignment columns.')
    order = [tuple(row) for row in rows]
    if len(set(order)) != len(order):
        convert(alignment_db, result, _fs_SCREEN_FIELDS)
        return result.read_text()
    jobs, tables = ([], [])
    for shard in range(workers):
        source, table = (directory / f'shard-{shard}', directory / f'shard-{shard}.tsv')
        tables.append(table)
        jobs.append([[executable, 'filterdb', alignment_db, source, '--filter-expression', f'$1 % {workers} == {shard}', '--threads', '1'], [executable, 'convertalis', query_db, database, source, table, '--threads', '1', '--format-output', ','.join(_fs_SCREEN_FIELDS), '--exact-tmscore', '1']])
    _fs__parallel_scoring(jobs, directory, deadline, cancel, progress)
    scored = {}
    for table in tables:
        for line in table.read_text().splitlines():
            values = line.split('\t')
            if len(values) != len(_fs_SCREEN_FIELDS):
                raise ValueError('Invalid Foldseek score columns.')
            fields = dict(zip(_fs_SCREEN_FIELDS, values))
            key = tuple((fields[field] for field in _fs_ALIGNMENT_FIELDS))
            if key in scored:
                raise ValueError('Foldseek scored an alignment more than once.')
            scored[key] = line
    if set(scored) != set(order):
        raise ValueError('Foldseek changed or omitted alignments during parallel scoring.')
    return ''.join((scored[key] + '\n' for key in order))

def _fs__merge_alignment_shards(shards, destination, prefilter, deadline, cancel=None):
    """Restore structurealign's global score/length/target ordering per query."""
    expected = [int(line.split()[0]) for line in _fs_Path(str(prefilter) + '.index').read_text().splitlines()]
    rows = {key: [] for key in expected}
    seen, dbtype = (set(), None)
    for shard in shards:
        kind = _fs_Path(str(shard) + '.dbtype').read_bytes()
        if dbtype is not None and kind != dbtype:
            raise ValueError('Candidate alignment shards have different database types.')
        dbtype = kind
        with _fs_Path(shard).open('rb') as stream:
            for line in _fs_Path(str(shard) + '.index').read_text().splitlines():
                if cancel is not None and cancel.is_set():
                    raise InterruptedError('Foldseek search cancelled.')
                if _fs_time.monotonic() >= deadline:
                    raise TimeoutError('Foldseek search timed out.')
                key, offset, length = map(int, line.split())
                if key not in rows or offset < 0 or length < 1:
                    raise ValueError('Invalid candidate alignment shard index.')
                stream.seek(offset)
                data = stream.read(length)
                if len(data) != length or not data.endswith(b'\x00'):
                    raise ValueError('Truncated candidate alignment shard.')
                for record in data[:-1].splitlines(keepends=True):
                    fields = record.split()
                    if len(fields) < 10:
                        raise ValueError('Invalid candidate alignment record.')
                    target, score, target_length = (int(fields[0]), int(fields[1]), int(fields[9]))
                    if (key, target) in seen:
                        raise ValueError('Candidate alignment was calculated more than once.')
                    seen.add((key, target))
                    rows[key].append(((-score, target_length, target), record))
    with _fs_Path(destination).open('wb') as data, _fs_Path(str(destination) + '.index').open('w') as index:
        for key in expected:
            offset = data.tell()
            for _, record in sorted(rows[key], key=lambda row: row[0]):
                data.write(record)
            data.write(b'\x00')
            index.write(f'{key}\t{offset}\t{data.tell() - offset}\n')
    if dbtype is None:
        raise ValueError('No candidate alignment shards were produced.')
    _fs_Path(str(destination) + '.dbtype').write_bytes(dbtype)

def _fs__candidate_parallel_alignment(query_db, database, alignment_db, directory, threads, deadline, cancel, progress, executable, max_seqs, evalue):
    """Keep the full prefilter, then independently align target-candidate shards."""
    directory.mkdir()
    prefilter = directory / 'prefilter'
    _fs__run_foldseek([executable, 'ungappedprefilter', str(query_db) + '_ss', str(database) + '_ss', prefilter, '--threads', threads, '--max-seqs', max_seqs, '--comp-bias-corr', '1', '--comp-bias-corr-scale', '0.15', '--min-ungapped-score', '30', '--sub-mat', 'aa:3di.out,nucl:3di.out', '-c', '0', '--cov-mode', '0', '-e', '1.7976931348623157e+308', '--prefilter-mode', '1'], directory / 'prefilter.log', deadline - _fs_time.monotonic(), cancel, progress)
    workers = min(threads, 8)
    jobs, shards = ([], [])
    for worker in range(workers):
        selected, aligned = (directory / f'prefilter-{worker}', directory / f'alignments-{worker}')
        shards.append(aligned)
        jobs.append([[executable, 'filterdb', prefilter, selected, '--filter-expression', f'$1 % {workers} == {worker}', '--threads', '1'], [executable, 'structurealign', query_db, database, selected, aligned, '--threads', '1', '-a', '1', '--alignment-mode', '3', '--alignment-type', '2', '--sort-by-structure-bits', '1', '--comp-bias-corr', '1', '--comp-bias-corr-scale', '0.5', '--sub-mat', 'aa:3di.out,nucl:3di.out', '--min-aln-len', '20', '--exact-tmscore', '1', '-e', evalue, '--compressed', '0']])
    _fs__parallel_scoring(jobs, directory, deadline, cancel, progress)
    _fs__merge_alignment_shards(shards, alignment_db, prefilter, deadline, cancel)

def _fs__run_search_pass(query_dir, queries, library, timeout=120, cancel=None, executable='foldseek', threads=_fs_DEFAULT_THREADS, progress=None, manifest=None, max_seqs=1000, evalue=10, candidate_parallel=True):
    library = _fs_Path(library).resolve()
    if manifest is None:
        manifest = _fs_json.loads((library / 'manifest.json').read_text())
    if manifest.get('schema') not in (1, 2) or (manifest['schema'] == 1 and (not isinstance(manifest.get('records'), list))):
        raise ValueError('Unsupported domain library manifest.')
    database = _fs_library_file(library, manifest['database'])
    if timeout <= 0 or not _fs_math.isfinite(timeout) or threads < 1:
        raise ValueError('Search timeout and thread count must be positive.')
    if cancel is not None and cancel.is_set():
        raise InterruptedError('Foldseek search cancelled.')
    started = _fs_time.monotonic()
    command_threads = min(4, threads)
    with _fs_tempfile.TemporaryDirectory(prefix='domain-foldseek-') as directory:
        directory = _fs_Path(directory)
        query_db, alignment_db = (directory / 'query', directory / 'alignments')
        result = directory / 'hits.tsv'

        def run(*arguments):
            _fs__run_foldseek([executable, *arguments], directory / 'log', timeout - (_fs_time.monotonic() - started), cancel, progress)

        def convert(source, destination, fields):
            run('convertalis', query_db, database, source, destination, '--threads', command_threads, '--format-output', ','.join(fields), '--exact-tmscore', '1')
        run('createdb', query_dir, query_db, '--threads', command_threads, '--chain-name-mode', '0')
        query_lengths = [len(query['positions']) for query in queries.values()]
        parallel = candidate_parallel and command_threads > 1 and (max(query_lengths) >= 256) and (max(query_lengths) * command_threads >= 2 * sum(query_lengths)) and (_fs_library_reference_count(manifest) >= 512) and _fs_Path(str(query_db) + '_ss').is_file() and _fs_Path(str(database) + '_ss').is_file() and (not any((_fs_Path(str(database) + suffix).exists() for suffix in ('.idx', '_ss.idx'))))
        if parallel:
            _fs__candidate_parallel_alignment(query_db, database, alignment_db, directory / 'candidate-shards', command_threads, started + timeout, cancel, progress, executable, max_seqs, evalue)
        else:
            run('search', query_db, database, alignment_db, directory / 'work', '--threads', command_threads, '-a', '1', '--alignment-mode', '3', '--max-seqs', max_seqs, '-e', evalue, '-s', '9.5', '-k', '6', '--prefilter-mode', '1', '--min-aln-len', '20', '--exact-tmscore', '1')
        text = _fs__screen_alignments(query_db, database, alignment_db, result, directory, threads, started + timeout, cancel, progress, executable, convert)
        metadata = _fs_hit_metadata(library, manifest, text, _fs_SCREEN_FIELDS)
        hits = _fs_parse_hits(text, queries, metadata, _fs_SCREEN_FIELDS)
    return {'elapsed_seconds': _fs_time.monotonic() - started, 'pilot': manifest.get('pilot', False), 'references': _fs_library_reference_count(manifest), 'queries': queries, 'hits': hits, 'catalog': manifest.get('catalog', 'CDD')}

def _fs_review_windows(hits, queries, index_maps, tested):
    """Revisit provisional assignments without letting them block their own region."""
    stronger = {root: set() for root in queries}
    proposals = []
    for hit in _fs_select_matches(hits, per_family=False):
        root = hit['query']
        indices = {index_maps[root][p] for p in hit['aligned_positions']}
        first = max(0, min(indices) - 12)
        last = min(len(queries[root]['positions']), max(indices) + 13)
        runs, start = ([], None)
        for index in range(first, last + 1):
            available = index < last and index not in stronger[root]
            if available and start is None:
                start = index
            elif not available and start is not None:
                runs.append((start, index))
                start = None
        if runs:
            a, b = max(runs, key=lambda span: sum((span[0] <= p < span[1] for p in indices)))
            window = (root, a, b)
            if sum((a <= p < b for p in indices)) >= 20 and (*window, 5000) not in tested:
                proposals.append((hit['evalue'], hit['bits'], window))
        stronger[root].update(indices)
    counts, selected = ({}, [])
    for _, _, window in sorted(proposals, key=lambda row: (-row[0], row[1])):
        root = window[0]
        if counts.get(root, 0) < 4 and window not in selected:
            selected.append(window)
            counts[root] = counts.get(root, 0) + 1
    return selected

def _fs_search_windows(queries, index_maps, covered, pass_hits, tested, depth):
    """Propose unassigned stretches and bounded searches around plausible hits."""
    windows = {}
    for root, query in queries.items():
        start = None
        for index in range(len(query['positions']) + 1):
            unassigned = index < len(query['positions']) and index not in covered[root]
            if unassigned and start is None:
                start = index
            elif not unassigned and start is not None:
                window = (root, start, index)
                if index - start >= 20 and (*window, depth) not in tested:
                    windows[f'residual{len(windows)}'] = window
                start = None
    focused = {root: [] for root in queries}
    for hit in sorted(pass_hits, key=lambda h: (h['evalue'], -h['bits'])):
        if _fs_qualifying_hit(hit) or hit['evalue'] > 10 or hit['tcov'] < 0.5 or (hit['alntmscore'] < 0.5) or (len(hit['aligned_positions']) < 20):
            continue
        root = hit['query']
        indices = {index_maps[root][p] for p in hit['aligned_positions']}
        if len(indices & covered[root]) >= 0.5 * len(indices) or len(focused[root]) >= 4:
            continue
        first, last = (max(0, min(indices) - 12), min(len(queries[root]['positions']), max(indices) + 13))
        window = (root, first, last)
        if (*window, depth) in tested or window in windows.values():
            continue
        if any((max(0, min(last, b) - max(first, a)) / max(last - first, b - a) >= 0.7 for a, b in focused[root])):
            continue
        focused[root].append((first, last))
        windows[f'focus{len(windows)}'] = window
    return windows

def _fs_guided_windows(queries, guidance):
    """Pad confident sequence boundaries and cover every remaining modeled residue."""
    from bisect import bisect_left, bisect_right
    windows = {}
    useful = False
    for root, query in queries.items():
        positions = query['positions']
        spans = set()
        for a, b in (guidance or {}).get(root, []):
            first, last = (bisect_left(positions, a), bisect_right(positions, b))
            if last - first >= 20:
                spans.add((max(0, first - 12), min(len(positions), last + 12)))
        if not spans:
            spans.add((0, len(positions)))
        useful |= spans != {(0, len(positions))}
        covered = set((i for a, b in spans for i in range(a, b)))
        start = None
        for i in range(len(positions) + 1):
            if i < len(positions) and i not in covered:
                if start is None:
                    start = i
            elif start is not None:
                a, b = (max(0, start - 12), min(len(positions), i + 12))
                a = max(0, min(a, b - 20))
                b = min(len(positions), max(b, a + 20))
                spans.add((a, b))
                start = None
        for first, last in sorted(spans):
            windows[f'guided{len(windows)}'] = (root, first, last)
    return windows if useful else {}

def _fs_run_search(query_dir, queries, library, timeout=120, cancel=None, executable='foldseek', threads=_fs_DEFAULT_THREADS, progress=None, guidance=None, candidate_parallel=True):
    """Revisit unassigned stretches and focus plausible below-cutoff placements."""
    started = _fs_time.monotonic()
    manifest = _fs_json.loads((_fs_Path(library) / 'manifest.json').read_text())
    source_lines = {}
    for identifier, query in queries.items():
        lines = (_fs_Path(query_dir) / (identifier + '.pdb')).read_text().splitlines(keepends=True)
        source_lines[identifier] = [line for line in lines if line.startswith('ATOM  ') and line[12:16].strip() == 'CA']
        if len(source_lines[identifier]) != len(query['positions']):
            raise ValueError('Query coordinate file does not match its residue index map.')
    index_maps = {identifier: {p: i for i, p in enumerate(query['positions'])} for identifier, query in queries.items()}
    covered = {identifier: set() for identifier in queries}
    windows = {identifier: (identifier, 0, len(query['positions'])) for identifier, query in queries.items()}
    tested, hits, result, final_review = (set(), [], None, False)
    with _fs_tempfile.TemporaryDirectory(prefix='domain-repeat-') as directory:

        def prepare(next_windows, name):
            folder = _fs_Path(directory) / str(name)
            folder.mkdir()
            active = {}
            for identifier, (root, first, last) in next_windows.items():
                lines = [f'{line[:6]}{number:5d}{line[11:22]}{number:4d}{line[26:]}' for number, line in enumerate(source_lines[root][first:last], 1)]
                (folder / (identifier + '.pdb')).write_text(''.join(lines) + 'TER\nEND\n')
                active[identifier] = {'positions': queries[root]['positions'][first:last], 'label': queries[root].get('label', root)}
            return (folder, active)
        folder = _fs_Path(query_dir)
        active_queries = queries
        initial = _fs_guided_windows(queries, guidance)
        if initial:
            windows = initial
            folder, active_queries = prepare(windows, 'guided')
        review_after = 4 if initial else 3
        search_round = 0
        while True:
            search_round += 1
            guided = bool(initial) and search_round == 1
            max_seqs = 5000 if final_review else 1000
            remaining = timeout - (_fs_time.monotonic() - started)
            if remaining <= 0:
                raise TimeoutError('Foldseek search timed out.')
            phase = 'Pfam-guided search' if guided else 'final review' if final_review else 'search'
            message = f'{phase} pass {search_round}, {len(active_queries)} query region(s)'
            report = (lambda message=message: progress(message)) if progress is not None else None
            if report is not None:
                report()
            result = _fs__run_search_pass(folder, active_queries, library, remaining, cancel, executable, threads, report, manifest, max_seqs=max_seqs, evalue=0.01 if final_review or guided else 10, **{} if candidate_parallel else {'candidate_parallel': False})
            for identifier, window in windows.items():
                if not guided:
                    tested.add((*window, max_seqs))
            pass_hits = []
            for hit in result['hits']:
                root, first, last = windows[hit['query']]
                hit['search_query'] = hit['query']
                hit['query'] = root
                hit['search_round'] = search_round
                hit['prefilter_limit'] = max_seqs
                hit['search_region'] = [queries[root]['positions'][first], queries[root]['positions'][last - 1]]
                hit['original_query_length'] = len(queries[root]['positions'])
                if guided:
                    hit['search_strategy'] = 'pfam_guided'
                hits.append(hit)
                pass_hits.append(hit)
            covered = {identifier: set() for identifier in queries}
            for hit in _fs_select_matches(hits, per_family=False):
                root = hit['query']
                covered[root].update((index_maps[root][position] for position in hit['aligned_positions']))
            if final_review:
                next_windows = {key: window for key, window in _fs_search_windows(queries, index_maps, covered, pass_hits, tested, 5000).items() if key.startswith('residual')}
                if not next_windows:
                    break
                folder, active_queries = prepare(next_windows, search_round)
                windows = next_windows
                continue
            next_windows = _fs_search_windows(queries, index_maps, covered, pass_hits, tested, 1000)
            if not next_windows or search_round >= review_after:
                final_review = True
                next_windows = _fs_search_windows(queries, index_maps, covered, pass_hits, tested, 5000)
                for window in _fs_review_windows(hits, queries, index_maps, tested):
                    if window not in next_windows.values():
                        next_windows[f'review{len(next_windows)}'] = window
            if not next_windows:
                break
            folder, active_queries = prepare(next_windows, search_round)
            windows = next_windows
    return {**result, 'queries': queries, 'hits': sorted(hits, key=lambda h: (h['evalue'], -h['bits'])), 'elapsed_seconds': _fs_time.monotonic() - started, 'search_rounds': search_round}

_DOMAIN_SEARCH_RUNTIME = _RuntimeNamespace('_fs_')

# END INTEGRATED DOMAIN RUNTIME


_PFAM_API = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
_PFAM_POLL_SECONDS = 5
_PFAM_CACHE_VERSION = 3
_PFAM_TENTATIVE_EVALUE = 0.01
_PFAM_TENTATIVE_OVERLAP = 0.5
_DOMAIN_BROWSER_VERSION = 21
_DOMAIN_BROWSER_IMPLEMENTATION = object()
_DOMAIN_INFO_RENDER_VERSION = 3
_DOMAIN_INFO_ENRICH_VERSION = 2
_INTERPRO_CACHE_VERSION = 2


def _sequence_feature_kind(hit):
    kind = hit.get("feature_kind")
    if kind:
        return kind
    kind = hit.get("type", "").lower().replace("-", "_").replace(" ", "_")
    return kind if kind in {"motif", "disordered", "coiled_coil", "low_complexity",
                            "active_site", "binding_site", "conserved_site", "ptm", "site"} else None




def _sequence_feature_backend():
    return _SEQUENCE_FEATURE_RUNTIME


def _start_domain_features(entries, options, timeout, cancel, report):
    """Scan immutable sequences alongside domain work; never touch models here."""
    from concurrent.futures import Future
    from pathlib import Path
    from threading import Event, Thread
    from types import SimpleNamespace

    if options is None:
        return None
    sequences = list(dict.fromkeys(entry["sequence"] for entry in entries))
    future, stopped = Future(), Event()
    combined_cancel = SimpleNamespace(is_set=lambda: cancel.is_set() or stopped.is_set())

    def scan():
        try:
            module = _sequence_feature_backend()
            value = module.search(sequences, Path(__file__).parent / ".domain_data" / "sequence-features",
                                  timeout, combined_cancel, report, **options)
        except Exception as error:
            value = None, ["Sequence features unavailable; domain results retained: " + str(error)]
        future.set_result(value)

    Thread(target=scan, name="ChimeraX automatic sequence features", daemon=True).start()
    return SimpleNamespace(future=future, stopped=stopped, cancel=cancel)


def _await_domain_features(task, report):
    from concurrent.futures import TimeoutError as FutureTimeout
    from time import monotonic

    if task is None:
        return None, []
    next_report = monotonic() + 1
    while not task.future.done():
        if task.cancel.is_set():
            raise InterruptedError("Domain and feature lookup cancelled.")
        try:
            return task.future.result(timeout=.1)
        except FutureTimeout:
            pass
        if monotonic() >= next_report:
            report("waiting for sequence features")
            next_report = monotonic() + 1
    return task.future.result()


def _domain_completion_counts(records):
    features = sum(bool(_sequence_feature_kind(record)) for record in records)
    return f"{len(records) - features} matches" + (f"; {features} sequence features" if features else "")


def _pfam_local_backend(database=None, library=None):
    from pathlib import Path
    module = _PFAM_LOCAL_RUNTIME
    backend = module.discover(Path(__file__).parent, database, library)
    return {**backend, "module": module} if backend else None


def _pfam_backend(mode="auto", database=None, library=None):
    online = {"kind": "online", "identity": "online-v1", "requested": mode}
    if mode == "online":
        if database is not None:
            raise ValueError("pfam_database cannot be combined with pfam_backend online.")
        return online
    try:
        local = _pfam_local_backend(database, library)
    except (OSError, ValueError, KeyError) as error:
        if mode == "local":
            raise
        return {**online, "note": "Local Pfam unavailable; using the online service: " + str(error)}
    if local is None:
        if mode == "local":
            raise ValueError("No compatible local Pfam bundle and hmmscan installation were found.")
        return online
    return {**local, "requested": mode}


def _pfam_cached(session, backend):
    if (getattr(session, "_trimmings_pfam_cache_version", None) != _PFAM_CACHE_VERSION
            or getattr(session, "_trimmings_pfam_cache_backend", "online-v1") != backend["identity"]):
        return {}
    return getattr(session, "_trimmings_pfam_cache", {})


def _pfam_cache_results(session, results, backend):
    cache = dict(_pfam_cached(session, backend))
    cache.update(results)
    session._trimmings_pfam_cache = cache
    session._trimmings_pfam_cache_version = _PFAM_CACHE_VERSION
    session._trimmings_pfam_cache_backend = backend["identity"]


def _pfam_lookup(sequences, timeout, cancel, report, backend, all_sequences=None):
    """Use local profiles when available, with one deadline and coherent fallback."""
    from time import monotonic

    started = monotonic()
    if backend["kind"] == "online":
        return _pfam_search(sequences, timeout, cancel, report), backend, backend.get("note")
    try:
        budget = min(60, timeout / 2) if backend["requested"] == "auto" else timeout
        results = backend["module"].search(sequences, backend, budget, cancel, report)
        return results, backend, None
    except Exception as error:
        if cancel.is_set() or isinstance(error, InterruptedError) or backend["requested"] == "local":
            raise
        remaining = timeout - (monotonic() - started)
        if remaining <= 0:
            raise TimeoutError("Pfam lookup timed out.") from None
        note = "Local Pfam failed; using the online service: " + (str(error) or type(error).__name__)
        report(note)
        online = {"kind": "online", "identity": "online-v1", "requested": "auto"}
        # If some local results were cached, refetch those sequences as well;
        # one response must not mix profile releases or backends.
        results = _pfam_search(all_sequences or sequences, remaining, cancel, report)
        return results, online, note


def _pfam_request(path, timeout, payload=None):
    import json
    from urllib.error import HTTPError
    from urllib.request import Request, urlopen

    request = Request(_PFAM_API + path,
                      data=json.dumps(payload).encode("utf-8") if payload is not None else None, headers={
        "User-Agent": "ChimeraX-Trimmings/1.0", "Accept": "application/json",
        "Content-Type": "application/json",
    })
    try:
        with urlopen(request, timeout=timeout) as response:
            try:
                return json.load(response)
            except (ValueError, UnicodeError):
                raise ValueError("EMBL-EBI HMMER returned an invalid JSON response.") from None
    except HTTPError as error:
        raise ValueError(f"EMBL-EBI HMMER returned HTTP {error.code}: {error.reason}.") from None


def _pfam_parse_hits(payload, sequence):
    """Read hmmscan sequence coordinates, not the profile-HMM coordinates."""
    import math
    import re

    try:
        if payload["status"] != "SUCCESS" or not isinstance(payload["result"]["hits"], list):
            raise ValueError()
        results = []
        for hit in payload["result"]["hits"]:
            metadata = hit["metadata"]
            accession = metadata["accession"]
            if not isinstance(accession, str) or not re.fullmatch(r"PF\d{5,}", accession):
                raise ValueError()
            for domain in hit["domains"]:
                alignment = domain["alignment_display"]
                start, end = alignment["sqfrom"], alignment["sqto"]
                if (type(start) is not int or type(end) is not int
                        or not 1 <= start <= end <= len(sequence) or alignment["l"] != len(sequence)):
                    raise ValueError()
                evalue, score = float(domain["ievalue"]), float(domain["bitscore"])
                if not math.isfinite(evalue) or evalue < 0 or not math.isfinite(score):
                    raise ValueError()
                if not isinstance(metadata["identifier"], str) or not isinstance(metadata["type"], str):
                    raise ValueError()
                tentative = not domain["significant"]
                if tentative and evalue > _PFAM_TENTATIVE_EVALUE:
                    continue
                cutoff = metadata.get("dom_ga")
                if cutoff is not None:
                    cutoff = float(cutoff)
                    if not math.isfinite(cutoff):
                        raise ValueError()
                results.append({
                    "type": metadata["type"], "start": start, "end": end,
                    "evalue": evalue, "bitscore": score,
                    "accession": accession, "name": metadata["identifier"],
                    "description": metadata.get("description", ""),
                    "outcompeted": bool(domain["outcompeted"]),
                    "tentative": tentative, "gathering_cutoff": cutoff,
                })
        return results
    except (KeyError, TypeError, ValueError, OverflowError):
        raise ValueError("Malformed Pfam domain results or invalid sequence coordinates; no annotations were applied.") from None


def _pfam_search(sequences, timeout, cancel, report):
    from time import monotonic
    from uuid import UUID

    deadline = monotonic() + timeout
    results = {}

    def request(path, payload=None):
        remaining = deadline - monotonic()
        if remaining <= 0:
            raise ValueError(f"Pfam lookup timed out after {timeout:g} seconds.")
        return _pfam_request(path, min(60, remaining), payload)

    for index, sequence in enumerate(sequences, 1):
        if cancel.is_set():
            return None
        report(f"Searching Pfam sequence {index}/{len(sequences)} at EMBL-EBI HMMER.")
        submitted = request("/search/hmmscan", {
            "input": ">domain_query\n" + sequence, "database": "pfam", "threshold": "cut_ga",
        })
        try:
            rid = str(UUID(submitted["id"]))
        except (KeyError, TypeError, ValueError, AttributeError):
            raise ValueError("EMBL-EBI HMMER did not return a valid search ID.") from None
        page, page_count, hits = 1, 1, []
        while page <= page_count:
            if cancel.is_set():
                return None
            response = request(f"/result/{rid}?with_domains=true&page_size=100&page={page}")
            if not isinstance(response, dict):
                raise ValueError("EMBL-EBI HMMER returned an unrecognized search response.")
            status = response.get("status")
            if status in ("PENDING", "STARTED", "RUNNING", "RETRY", "RECEIVED"):
                report(f"Pfam search {index}/{len(sequences)} queued/running ({rid}).")
                cancel.wait(min(_PFAM_POLL_SECONDS, max(0, deadline - monotonic())))
                continue
            if status != "SUCCESS":
                raise ValueError(f"EMBL-EBI HMMER search {rid} failed (status {status!r}).")
            hits.extend(_pfam_parse_hits(response, sequence))
            page_count = response.get("page_count")
            if (type(page_count) is not int or page_count < 0
                    or (page_count == 0 and (page != 1 or response["result"]["hits"]))):
                raise ValueError("EMBL-EBI HMMER returned an invalid result page count.")
            page += 1
        results[sequence] = hits
    return results


_INTERPRO_DATABASES = {
    "pfam": "pfam", "cath-gene3d": "cathgene3d", "superfamily": "ssf", "cdd": "cdd",
    "smart": "smart", "prositepatterns": "prosite", "prositeprofiles": "profile",
    "panther": "panther", "pirsf": "pirsf", "hamap": "hamap", "prints": "prints",
    "sfld": "sfld", "ncbifam": "ncbifam", "ncbifams": "ncbifam", "interpro": "interpro",
}


def _interpro_matches_request(md5s, timeout):
    import json
    from urllib.error import HTTPError
    from urllib.request import Request, urlopen

    request = Request("https://www.ebi.ac.uk/interpro/matches/api/matches",
                      data=json.dumps({"md5": md5s}).encode("ascii"), headers={
                          "Accept": "application/json", "Content-Type": "application/json",
                          "User-Agent": "ChimeraX-Trimmings/1.0",
                      })
    try:
        with urlopen(request, timeout=timeout) as response:
            return json.load(response)
    except HTTPError as error:
        raise ValueError(f"InterPro Matches API HTTP {error.code}: {error.reason}.") from None


def _interpro_parse_matches(matches, sequence):
    """Retain independent member-database locations, including discontinuous domains."""
    from math import isfinite
    import re
    from urllib.parse import quote

    accepted_types = {"domain", "repeat", "homologous_superfamily", "family",
                      "active_site", "binding_site", "conserved_site", "ptm", "site",
                      "motif", "disordered", "coiled_coil", "low_complexity"}
    results = []
    try:
        if not isinstance(matches, list):
            raise ValueError()
        for match in matches:
            signature = match["signature"]
            provider = match.get("source") or signature["signatureLibraryRelease"]["library"]
            # Keep member-database evidence separate from InterPro-N predictions.
            if provider == "InterPro-N":
                continue
            entry = signature.get("entry") or {}
            signature_type = signature.get("type") or ""
            entry_type = (signature_type if _sequence_feature_kind({"type": signature_type})
                          else entry.get("type") or signature_type)
            if entry_type.lower().replace("-", "_").replace(" ", "_") not in accepted_types:
                continue
            accession = signature["accession"]
            if not isinstance(accession, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,100}", accession):
                raise ValueError()
            library = signature["signatureLibraryRelease"]["library"]
            database = _INTERPRO_DATABASES.get(library.lower())
            if database is None:
                continue
            integrated = entry.get("accession")
            if integrated is not None and not re.fullmatch(r"IPR\d+", integrated):
                raise ValueError()
            name = signature.get("name") or entry.get("name") or signature.get("description") or accession
            if not isinstance(name, str) or not isinstance(provider, str):
                raise ValueError()
            for location in match["locations"]:
                start, end = location["start"], location["end"]
                fragments = location.get("location-fragments") or [{"start": start, "end": end}]
                segments = [(f["start"], f["end"]) for f in fragments]
                if (type(start) is not int or type(end) is not int or not 1 <= start <= end <= len(sequence)
                        or any(type(a) is not int or type(b) is not int or not start <= a <= b <= end for a, b in segments)):
                    raise ValueError()
                evalue, score = location.get("evalue", match.get("evalue")), location.get("score", match.get("score"))
                for value in (evalue, score):
                    if value is not None and not isfinite(float(value)):
                        raise ValueError()
                if evalue is not None and float(evalue) < 0:
                    raise ValueError()
                hit = {
                    "origin": "interpro", "source": provider, "database": database,
                    "accession": accession, "name": name, "type": entry_type,
                    "description": signature.get("description") or entry.get("description") or "",
                    "start": start, "end": end, "segments": segments,
                    "evalue": float(evalue) if evalue is not None else None,
                    "score": float(score) if score is not None else None, "bitscore": None,
                    "integrated": integrated, "color_accession": integrated or accession,
                    "entry_url": f"https://www.ebi.ac.uk/interpro/entry/{database}/{quote(accession, safe='')}/",
                    "outcompeted": False,
                }
                results.append(hit)
                # Some providers attach conserved/catalytic sites to a domain
                # location instead of returning a separate site signature.
                for site in location.get("sites") or []:
                    pieces = site["siteLocations"]
                    site_segments = sorted(set((piece["start"], piece["end"]) for piece in pieces))
                    if (not site_segments or any(type(a) is not int or type(b) is not int
                                                or not start <= a <= b <= end for a, b in site_segments)):
                        raise ValueError()
                    description = site.get("description") or "Conserved site"
                    if not isinstance(description, str):
                        raise ValueError()
                    results.append({**hit, "feature_kind": "conserved_site", "type": "Conserved site",
                                    "name": description, "description": description,
                                    "evalue": None, "score": None,
                                    "parent_signature_name": name, "parent_evalue": hit["evalue"],
                                    "parent_score": hit["score"],
                                    "start": site_segments[0][0], "end": max(b for a, b in site_segments),
                                    "segments": site_segments,
                                    "evidence": "Member-database site prediction reported by InterPro; "
                                                "inspect the supporting family annotation."})
        return results
    except (KeyError, TypeError, ValueError, AttributeError):
        raise ValueError("Malformed InterPro matches or invalid sequence coordinates; no annotations were applied.") from None


def _interpro_search(sequences, timeout, cancel, report):
    from hashlib import md5
    from time import monotonic

    deadline, results = monotonic() + timeout, {}
    for offset in range(0, len(sequences), 100):
        if cancel.is_set():
            return None
        batch = sequences[offset:offset + 100]
        by_hash = {md5(seq.encode("ascii")).hexdigest().upper(): seq for seq in batch}
        report(f"Looking up {len(batch)} exact sequence(s) in InterPro Matches.")
        remaining = deadline - monotonic()
        if remaining <= 0:
            raise ValueError(f"InterPro lookup timed out after {timeout:g} seconds.")
        response = _interpro_matches_request(list(by_hash), min(60, remaining))
        try:
            items = response["results"]
            if not isinstance(items, list) or len(items) != len(by_hash):
                raise ValueError()
            seen = set()
            for item in items:
                checksum = item["md5"].upper()
                if checksum not in by_hash or checksum in seen or type(item["found"]) is not bool:
                    raise ValueError()
                seen.add(checksum)
                if not item["found"]:
                    raise LookupError("InterPro has no precomputed matches for one or more exact chain sequences. "
                                      "Use 'domains <atom-spec> source pfam' for a new sequence search.")
                sequence = by_hash[checksum]
                results[sequence] = _interpro_parse_matches(item["matches"], sequence)
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise ValueError(f"Invalid InterPro lookup response: {error}") from None
    return results


def _domain_hit_sort_key(hit):
    if hit.get("origin") == "foldseek":
        return (False, 0, -hit["bitscore"], hit["accession"], hit["start"], hit["end"])
    if hit.get("origin") == "interpro":
        # Scores from different member databases cannot be compared numerically.
        # Prefer structural superfamily coverage, then domain/repeat matches, then families.
        rank = {"ssf": 0, "cathgene3d": 1}.get(hit["database"], 2)
        if hit["type"].lower() == "family":
            rank = 3
        return (False, rank, -sum(b - a + 1 for a, b in hit["segments"]), hit["accession"], hit["start"], hit["end"])
    return (hit.get("tentative", False), hit["outcompeted"], -hit["bitscore"],
            hit["accession"], hit["start"], hit["end"])


def _domain_display_name(hit):
    return hit["name"] + (" (tentative)" if hit.get("tentative") else "")


def _domain_filter_hits(hits, include_overlaps, include_tentative):
    # Sequence features may legitimately lie inside a domain. Competition and
    # structural grouping only arbitrate domain assignments, never feature rows.
    features = [h for h in hits if _sequence_feature_kind(h)
                and (include_tentative or not h.get("tentative", False))]
    hits = [h for h in hits if not _sequence_feature_kind(h)]
    hits = [h for h in hits if (include_overlaps or not h["outcompeted"])
            and (include_tentative or not h.get("tentative", False))]
    if any(h.get("origin") == "foldseek" for h in hits):
        hits = _domain_collapse_structural_hits(hits)
    if not any(h.get("tentative") for h in hits):
        return hits + features

    def footprint(hit):
        return {position for a, b in hit.get("segments", [(hit["start"], hit["end"])])
                for position in range(a, b + 1)}

    # Compare full sequence footprints, independently of selected or modeled residues.
    curated = set()
    for hit in hits:
        if not hit.get("tentative"):
            curated.update(footprint(hit))
    visible = []
    for hit in hits:
        if hit.get("tentative"):
            positions = footprint(hit)
            if len(positions & curated) >= _PFAM_TENTATIVE_OVERLAP * len(positions):
                continue
        visible.append(hit)
    return visible + features


def _domain_collapse_structural_hits(hits):
    """One primary structural assignment per placement, retaining alternatives."""
    selected, footprints, families = [], [], []

    def family_keys(hit):
        match = hit.get("structural_match", {})
        reference = match.get("reference", {})
        accepted = (_foldseek_reference_annotation(match).get('pfam_accessions', [])
                    if reference.get('catalog') == 'ECOD' else [])
        # Cross-references may help group nearly identical placements; that is
        # weaker evidence than attaching a family name or functional description.
        keys = {('pfam', accession) for accession in reference.get('pfam_accessions', [])}
        keys.update(('pfam_clan', m['clan']) for m in reference.get('pfam_mappings', [])
                    if m.get('method') == 'hmmer' and m.get('accession') in accepted and m.get('clan'))
        ecod = reference.get("catalog_accession", hit["accession"])
        if ecod.startswith("ECOD:"):
            levels = ecod[5:].split(".")
            if len(levels) == 4 and all(level.isdigit() for level in levels):
                keys.add(("ecod_homology", ".".join(levels[:2])))
        if hit.get("database") == "pfam":
            keys.add(("pfam", hit["accession"]))
        return keys

    others = [hit for hit in hits if hit.get("origin") != "foldseek"]
    structural = sorted((hit for hit in hits if hit.get("origin") == "foldseek"),
                        key=lambda hit: (hit["evalue"], -hit["bitscore"], hit["accession"], hit["start"], hit["end"]))
    for hit in structural:
        keys = family_keys(hit)
        positions = {p for a, b in hit.get("aligned_segments", hit.get("segments", [(hit["start"], hit["end"])]))
                     for p in range(a, b + 1)}
        # Gaps in two alignments of the same domain can differ. Require actual
        # shared residues as well as comparable spans before grouping that case;
        # broad envelopes alone must not swallow inserted domains or repeats.
        overlaps = []
        for index, previous in enumerate(footprints):
            shared = len(positions & previous)
            overlap = shared / max(len(positions), len(previous))
            span_overlap = max(0, min(max(positions), max(previous)) - max(min(positions), min(previous)) + 1)
            mutual_span = span_overlap / max(max(positions) - min(positions) + 1, max(previous) - min(previous) + 1)
            # A weaker partial alignment may be completely contained yet fail
            # every mutual-overlap cutoff. Require real paired-residue support,
            # similar scale (or explicit shared family), and little extension.
            # Always compare against the primary, never an alternative's union.
            contained = (shared >= .9 * len(positions)
                         and span_overlap >= .9 * (max(positions) - min(positions) + 1)
                         and (len(positions) >= .5 * len(previous)
                              or bool(keys & families[index])))
            if overlap >= .7 or (overlap >= .5 and mutual_span >= .8) or contained:
                overlaps.append((overlap, index))
        if overlaps:
            _, index = max(overlaps)
            primary = selected[index]
            # Repeated queries/representatives of the same family are duplicates,
            # rather than additional family assignments in the details.
            for alternative in [hit, *hit.get("alternatives", [])]:
                if alternative["accession"] != primary["accession"] and not any(
                        previous["accession"] == alternative["accession"] for previous in primary["alternatives"]):
                    primary["alternatives"].append({key: value for key, value in alternative.items() if key != "alternatives"})
        else:
            selected.append({**hit, "alternatives": list(hit.get("alternatives", []))})
            footprints.append(positions)
            families.append(keys)
    return others + selected


def _domain_number(value):
    return "not reported" if value is None else f"{value:.3g}"


def _domain_entry_url(record):
    return record.get("entry_url", f'https://www.ebi.ac.uk/interpro/entry/pfam/{record["accession"]}/')


def _domain_chain_inputs(atoms):
    from chimerax.atomic import Residue
    from chimerax.core.errors import UserError

    targets = set(atoms.unique_residues)
    entries = []
    for structure in atoms.unique_structures:
        for chain in structure.chains:
            if chain.polymer_type != Residue.PT_PROTEIN:
                continue
            residues = tuple(chain.residues)
            selected = {r for r in residues if r is not None and r in targets}
            if not selected:
                continue
            sequence = chain.characters.upper()
            if len(sequence) != len(residues) or not sequence or not sequence.isascii() or not sequence.isalpha():
                raise UserError(f"Cannot map the sequence of {chain.atomspec} to its residues safely.")
            spec = chain.atomspec
            if not spec.startswith("#"):
                spec = f"#{structure.id_string}" + spec
            entries.append({"chain": chain, "sequence": sequence, "residues": residues,
                            "targets": selected, "spec": spec})
    if not entries:
        raise UserError("The atom specification must match at least one protein chain.")
    return entries


def _domain_color(accession):
    from colorsys import hsv_to_rgb
    from hashlib import sha256

    if not accession:
        return (0.65, 0.65, 0.65, 1)
    hue = int.from_bytes(sha256(accession.encode()).digest()[:4], "big") / 2**32
    return (*hsv_to_rgb(hue, 0.65, 0.95), 1)


def _domain_groups(records):
    """Group by actual signature, across chains/models; sites retain their names."""
    groups = {}
    for record in records:
        feature = _sequence_feature_kind(record)
        key = (feature or "domain", record.get("database", "pfam"), record["accession"],
               record["name"] if feature else "")
        groups.setdefault(key, []).append(record)
    return groups


def _domain_color_distance(a, b):
    """Perceptual separation in OKLab, plus circular hue separation."""
    from colorsys import rgb_to_hsv
    from math import sqrt
    def lab(color):
        r, g, b = (c / 12.92 if c <= .04045 else ((c + .055) / 1.055) ** 2.4 for c in color[:3])
        l = (.4122214708*r + .5363325363*g + .0514459929*b) ** (1/3)
        m = (.2119034982*r + .6806995451*g + .1073969566*b) ** (1/3)
        s = (.0883024619*r + .2817188376*g + .6299787005*b) ** (1/3)
        return (.2104542553*l + .793617785*m - .0040720468*s,
                1.9779984951*l - 2.428592205*m + .4505937099*s,
                .0259040371*l + .7827717662*m - .808675766*s)
    delta = abs(rgb_to_hsv(*a[:3])[0] - rgb_to_hsv(*b[:3])[0])
    return sqrt(sum((x-y)**2 for x,y in zip(lab(a),lab(b)))), min(delta, 1-delta)


def _domain_palette(records, previous):
    """Keep family colors consistent and separate sequence/spatial neighbors."""
    from colorsys import hsv_to_rgb
    from functools import lru_cache
    import numpy as np
    from chimerax.geometry import find_close_points

    chains, structures, neighbors = {}, {}, {}
    for record in records:
        if not _sequence_feature_kind(record):
            key = record["accession"]
            neighbors.setdefault(key, set())
            chains.setdefault(record["chain"], []).append(record)
            if 'residues' in record:
                atoms = record['residues'].atoms
                coords = atoms.filter(atoms.names == 'CA').coords
                if len(coords):
                    structures.setdefault(record['chain'].structure, []).append(
                        (key, coords, coords.min(axis=0), coords.max(axis=0)))
    for chain_records in chains.values():
        ordered = sorted(chain_records, key=lambda r: (r["start"], r["end"]))
        for first, second in zip(ordered, ordered[1:]):
            a, b = first["accession"], second["accession"]
            if a != b:
                neighbors[a].add(b)
                neighbors[b].add(a)
    # Domains can touch in 3D without being consecutive in the sequence.
    for domains in structures.values():
        for index, (a, xyz, lo, hi) in enumerate(domains):
            for b, other, other_lo, other_hi in domains[index+1:]:
                if a == b or b in neighbors[a]:
                    continue
                gap = np.maximum(0, np.maximum(lo-other_hi, other_lo-hi))
                if np.dot(gap, gap) > 64:
                    continue
                first, second = find_close_points(xyz, other, 8)
                if min(len(first), len(second)) >= 2:
                    neighbors[a].add(b)
                    neighbors[b].add(a)
    palette = dict(previous)
    candidates = [(*hsv_to_rgb(h / 36, saturation, value), 1)
                  for h in range(36) for saturation in (.55, .8) for value in (.82, .98)]
    @lru_cache(maxsize=None)
    def distance(a, b):
        return _domain_color_distance(a, b)
    def contrast(a, b):
        perceived, hue = distance(a, b)
        return min(perceived / .16, hue / .10)
    for key in sorted(neighbors, key=lambda k: (-len(neighbors[k]), k)):
        adjacent = [palette[n] for n in sorted(neighbors[key]) if n in palette]
        if key in palette and all(contrast(palette[key], other) >= 1 for other in adjacent):
            continue
        used = [palette[k] for k in neighbors if k in palette and k != key]
        preferred = _domain_color(key)
        available = [c for c in candidates if all(contrast(c, other) >= 1 for other in adjacent)] or candidates
        distinct = [c for c in available if all(distance(c, other)[0] >= .075 for other in used)]
        available = distinct or available
        palette[key] = max(available, key=lambda c: (
            min((contrast(c, other) for other in adjacent), default=0),
            min((distance(c, other)[0] for other in used), default=0),
            -distance(c, preferred)[0])) if used else preferred
    return palette


def _domain_record_color(session, record):
    if record.get("origin") == "tmsegments":
        return tuple(channel / 255 for channel in _TM_COLORS[record["tm_label"]])
    if not _sequence_feature_kind(record) and 'residues' in record:
        import numpy as np
        owned = np.fromiter((getattr(r, 'domain_accession', None) == record['accession']
                             for r in record['residues']), dtype=bool)
        if not owned.any():
            return None  # An overlapping alternative has no independently colored footprint.
        colors, counts = np.unique(record['residues'].ribbon_colors[owned], axis=0, return_counts=True)
        return tuple(float(c)/255 for c in colors[counts.argmax()])
    return getattr(session, "_trimmings_domain_palette", {}).get(record["accession"], _domain_color(record["accession"]))


def _set_domain_palette(session, value):
    session._trimmings_domain_palette = value
    tool = getattr(session, "_trimmings_domain_browser", None)
    if tool is not None and not tool._closed:
        tool.refresh()


def _domain_model_snapshots(session):
    snapshots = getattr(session, "_trimmings_domain_models", None)
    if snapshots is None:
        # Upgrade results produced before per-model browsing was available.
        snapshots = {}
        for record in getattr(session, "_trimmings_domains", []):
            snapshots.setdefault(record["structure"], {})[record["chain"]] = (
                record["sequence"], record["sequence_residues"],
            )
    return snapshots


def _set_domain_records(session, records, snapshots=None):
    session._trimmings_domains = records
    if snapshots is None:
        session._trimmings_domain_models = None
        snapshots = _domain_model_snapshots(session)
    session._trimmings_domain_models = snapshots
    tool = getattr(session, "_trimmings_domain_browser", None)
    if tool is not None and not tool._closed and not _domain_browser_is_current(tool):
        tool.delete()
        tool = None
    if tool is not None and not tool._closed:
        tool.refresh()


def _domain_chain_is_current(chain, sequence, residues):
    return (not chain.deleted and chain.structure is not None and not chain.structure.deleted
            and chain.structure.id is not None and chain.characters.upper() == sequence
            and tuple(chain.residues) == residues)


def _domain_analyzed_models(session):
    return sorted((model for model, chains in _domain_model_snapshots(session).items()
                   if not model.deleted and model.id is not None
                   and any(_domain_chain_is_current(chain, *snapshot) for chain, snapshot in chains.items())),
                  key=lambda model: model.id)


def _domain_record_is_current(record):
    return (record.get("origin", "pfam") in ("pfam", "interpro", "foldseek", "prosite", "sequence", "tmsegments")
            and bool(len(record["residues"])) and _domain_chain_is_current(
                record["chain"], record["sequence"], record["sequence_residues"]))


def _live_domain_records(session):
    # Records on one chain share a sequence snapshot. Validate it once per scan,
    # rather than rebuilding a potentially very long residue tuple for every hit.
    checked, live = {}, []
    for record in getattr(session, "_trimmings_domains", []):
        if record.get("origin", "pfam") not in ("pfam", "interpro", "foldseek", "prosite", "sequence", "tmsegments") or not len(record["residues"]):
            continue
        key = record["chain"], record["sequence"], id(record["sequence_residues"])
        if key not in checked:
            checked[key] = _domain_chain_is_current(record["chain"], record["sequence"], record["sequence_residues"])
        if checked[key]:
            live.append(record)
    return live


def _apply_domain_results(session, entries, results, color, include_overlaps, gui=True, source="pfam",
                          include_tentative=False, feature_results=None):
    from html import escape
    from types import SimpleNamespace
    from urllib.parse import quote
    from chimerax.atomic import Residue, Residues, concise_residue_spec
    from chimerax.core.colors import Color
    from chimerax.core.objects import Objects
    from chimerax.core.undo import UndoState
    from chimerax.std_commands.color import color as color_objects

    # One primary assignment controls display; structural alternatives stay in details.
    attributes = ("domain_name", "domain_accession", "domain_names", "domain_accessions", "domain_tentative")
    for name in attributes[:-1]:
        Residue.register_attr(session, name, "Trimmings domains", attr_type=str)
    Residue.register_attr(session, "domain_tentative", "Trimmings domains", attr_type=bool)
    feature_attributes = ("sequence_feature_names", "sequence_feature_accessions")
    for name in feature_attributes:
        Residue.register_attr(session, name, "Trimmings sequence features", attr_type=str)
    feature_only = source in ("features", "tmsegments")
    rows, records, count, updated_chains = [], [], 0, set()
    color_groups = {}
    previous_records = getattr(session, "_trimmings_domains", [])
    previous_models = _domain_model_snapshots(session)
    snapshots = {model: dict(chains) for model, chains in previous_models.items()}
    with session.undo.aggregate("annotate domains"):
        undo = UndoState("Domain attributes")
        for entry in entries:
            chain = entry["chain"]
            if (chain.deleted or chain.structure is None or chain.structure.deleted
                    or any(r.deleted for r in entry["targets"])
                    or chain.characters.upper() != entry["sequence"]
                    or tuple(chain.residues) != entry["residues"]):
                session.logger.warning(f"Skipped {entry['spec']}: chain changed during the domain lookup.")
                continue
            updated_chains.add(chain)
            snapshots.setdefault(chain.structure, {})[chain] = (entry["sequence"], entry["residues"])
            raw_hits = results[entry.get("lookup_key", entry["sequence"])]
            if feature_results is not None:
                # Keep raw domain caches independent of optional sequence features.
                raw_hits = [*raw_hits, *feature_results[entry["sequence"]]]
            display_hits = [_foldseek_display_hit(hit, entry) if hit.get("origin") == "foldseek" else hit
                            for hit in raw_hits]
            if entry.get("structural_representative") != entry["spec"] and entry.get("structural_representative"):
                display_hits = [{**hit, "structural_representative": entry["structural_representative"]}
                                for hit in display_hits]
            hits = _domain_filter_hits(display_hits, include_overlaps, include_tentative)
            hits = _domain_structural_display_boundaries(entry, hits)
            hits.sort(key=lambda hit: (bool(_sequence_feature_kind(hit)),
                                      (hit["start"], hit["end"], hit["accession"]) if _sequence_feature_kind(hit)
                                      else _domain_hit_sort_key(hit)))
            membership = {r: [] for r in entry["targets"]}
            for hit in hits:
                members = list(dict.fromkeys(
                    r for a, b in hit.get("segments", [(hit["start"], hit["end"])])
                    for r in entry["residues"][a - 1:b] if r in membership))
                if not _sequence_feature_kind(hit):
                    for residue in members:
                        membership[residue].append(hit)
                if members:
                    records.append({**hit, "structure": chain.structure, "chain": chain,
                                    "chain_spec": entry["spec"], "sequence": entry["sequence"],
                                    "sequence_residues": entry["residues"], "residues": Residues(members)})
                    spec = concise_residue_spec(session, Residues(members))
                    # Always retain a model ID in links, even in a one-model session.
                    if not spec.startswith("#"):
                        spec = f"#{chain.structure.id_string}" + spec
                    rows.append(
                        f'<tr><td>{escape(entry["spec"])}</td>'
                        f'<td><a href="cxcmd:{quote("select " + spec, safe="")}">{escape(_domain_display_name(hit))}</a></td>'
                        f'<td><a href="{_domain_entry_url(hit)}">'
                        f'{hit["accession"]}</a></td><td>{escape(hit["type"])}</td>'
                        f'<td>{hit["start"]}&ndash;{hit["end"]}</td><td>{escape(spec)}</td>'
                        f'<td>{_domain_number(hit["evalue"])}</td><td>{escape(hit.get("source", "Pfam"))}</td></tr>'
                    )
            groups = {}
            for residue, matches in membership.items():
                if feature_only:
                    continue
                primary = matches[0] if matches else None
                values = (
                    primary["name"] if primary else "", primary["accession"] if primary else "",
                    " | ".join(dict.fromkeys(h["name"] for h in matches)),
                    " | ".join(dict.fromkeys(h["accession"] for h in matches)),
                    bool(primary and primary.get("tentative")),
                )
                for name, value in zip(attributes, values):
                    # None is ChimeraX's unset attribute value and is undoable.
                    undo.add(residue, name, getattr(residue, name, None), value)
                    setattr(residue, name, value)
                groups.setdefault(values[1] if primary else "", []).append(residue)
                count += bool(primary)
            if color and not feature_only:
                for accession, members in groups.items():
                    color_groups.setdefault(accession, []).extend(members)
        records.sort(key=lambda r: (r["structure"].id, r["chain"].chain_id, r["start"], r["end"], r["accession"]))
        refreshed_features = {(record["chain"], record.get("origin", "pfam")) for record in records
                              if _sequence_feature_kind(record)}
        def preserve(record):
            if record["chain"] not in updated_chains:
                return True
            if source == "features":
                return record.get("origin", "pfam") not in ("sequence", "prosite")
            if source == "tmsegments":
                return record.get("origin", "pfam") != "tmsegments"
            if feature_results is not None and record.get("origin") in ("sequence", "prosite"):
                return False
            return (bool(_sequence_feature_kind(record)) and record.get("origin", "pfam") != source
                    and (record["chain"], record.get("origin", "pfam")) not in refreshed_features)
        combined = [record for record in _live_domain_records(session) if preserve(record)] + records
        combined.sort(key=lambda r: (r["structure"].id, r["chain"].chain_id, r["start"], r["end"], r["accession"]))
        if color and not feature_only:
            previous_palette = dict(getattr(session, "_trimmings_domain_palette", {}))
            palette = _domain_palette(combined, previous_palette)
            # Updating a family color must update its older copies too, so the
            # same annotation remains consistent across every analyzed model.
            for record in combined:
                key = record["accession"]
                if (_sequence_feature_kind(record) or record["chain"] in updated_chains
                        or palette.get(key) == previous_palette.get(key, _domain_color(key))):
                    continue
                color_groups.setdefault(key, []).extend(r for r in record["residues"]
                                                        if getattr(r, "domain_accession", None) == key)
            palette_state = SimpleNamespace(set_palette=lambda value: _set_domain_palette(session, value))
            undo.add(palette_state, "set_palette", previous_palette, palette, option="M")
            session._trimmings_domain_palette = palette
            for accession, members in color_groups.items():
                rgba = palette.get(accession, _domain_color(accession))
                color_objects(session, Objects(atoms=Residues(list(dict.fromkeys(members))).atoms), Color(rgba), target="ac")
        feature_membership = {residue: [] for entry in entries if entry["chain"] in updated_chains
                              for residue in entry["targets"]}
        for record in combined:
            if _sequence_feature_kind(record) and record["chain"] in updated_chains:
                for residue in record["residues"]:
                    if residue in feature_membership:
                        feature_membership[residue].append(record)
        for residue, matches in feature_membership.items():
            for name, key in zip(feature_attributes, ("name", "accession")):
                value = " | ".join(dict.fromkeys(hit[key] for hit in matches))
                undo.add(residue, name, getattr(residue, name, None), value)
                setattr(residue, name, value)
        # Keep the browser's results in step with annotation Undo/Redo.
        state = SimpleNamespace(set_records=lambda value: _set_domain_records(session, *value))
        undo.add(state, "set_records", (previous_records, previous_models), (combined, snapshots), option="M")
        state.set_records((combined, snapshots))
        session.undo.register(undo)
    source_name = {"pfam": "Pfam", "interpro": "InterPro", "foldseek": "Foldseek", "features": "Sequence features",
                   "tmsegments": "DeepTMHMM features"}[source]
    if rows:
        session.logger.info(
            f'<b>{source_name} annotations</b> (click a name to select its footprint)'
            '<table border="1" cellpadding="3"><tr><th>Chain</th><th>Annotation</th><th>Accession</th>'
            '<th>Hit type</th><th>Query positions</th><th>Structure residues</th><th>E-value</th><th>Source</th></tr>'
            + "".join(rows) + "</table>", is_html=True,
        )
    tentative_count = sum(bool(r.get("tentative")) for r in records)
    confidence = f"{tentative_count} tentative; " if tentative_count else ""
    session.logger.info(f"{source_name}: {len(records)} annotation matches; {confidence}"
                        + ("domain annotations retained." if source == "tmsegments" else
                           "domain colors and attributes retained." if feature_only else
                           f"annotated {count} target residues. Unassigned residues have empty domain_name/domain_accession."))
    if gui and session.ui.is_gui:
        domain_browser(session)
    return records


def _interpro_request(path):
    import json
    from urllib.error import HTTPError
    from urllib.request import Request, urlopen

    request = Request("https://www.ebi.ac.uk/interpro/api" + path,
                      headers={"Accept": "application/json", "User-Agent": "ChimeraX-Trimmings/1.0"})
    try:
        with urlopen(request, timeout=20) as response:
            if response.status == 204:
                return None
            return json.load(response)
    except HTTPError as error:
        if error.code in (204, 404):
            return None
        raise ValueError(f"InterPro HTTP {error.code}: {error.reason}") from None


def _domain_plain_text(markup):
    from html.parser import HTMLParser
    import re

    class Text(HTMLParser):
        def __init__(self):
            super().__init__()
            self.parts, self.hidden = [], 0

        def handle_starttag(self, tag, attrs):
            if tag in ("script", "style"):
                self.hidden += 1
            elif tag in ("p", "br", "li", "div"):
                self.parts.append(" ")

        def handle_endtag(self, tag):
            if tag in ("script", "style"):
                self.hidden = max(0, self.hidden - 1)
            elif tag in ("p", "li", "div"):
                self.parts.append(" ")

        def handle_data(self, text):
            if not self.hidden:
                self.parts.append(text)

    parser = Text()
    parser.feed(markup if isinstance(markup, str) else "")
    text = " ".join("".join(parser.parts).split())
    citation = r"\[+cite:[^\[\]]+\]+"
    text = re.sub(citation + r"(?:\s*[,;]\s*" + citation + r")*", "", text)
    # Also repair separators left in previously cached plain-text descriptions.
    text = re.sub(r"\s+[,;](?:\s*[,;])+(?=\s|[.!?]|$)", "", text)
    text = re.sub(r"[,;]\s*(?=[.!?])", "", text)
    text = re.sub(r"\s+([,.;:!?])", r"\1", text)
    return " ".join(text.split()).strip()


def _domain_structure_match(item, pdb, database, accession, representative=False):
    """Require a family match on the crystallized chain, not its parent protein."""
    import re

    metadata = item.get('metadata', {})
    if str(metadata.get('accession', '')).lower() != pdb.lower():
        return None
    locations, regions = [], []
    for entry in item.get('entries', []):
        identity = str(entry.get('accession', '')).casefold() == accession.casefold()
        if database == 'interpro':
            identity |= str(entry.get('entry_integrated', '')).casefold() == accession.casefold()
        if not identity or not entry.get('chain'):
            continue
        for location in entry.get('entry_structure_locations', []):
            segments = []
            for fragment in location.get('fragments', []):
                a, b = fragment.get('start'), fragment.get('end')
                if (not isinstance(a, int) or not isinstance(b, int) or not 1 <= a <= b):
                    continue
                locations.append({'chain': entry['chain'], 'start': a, 'end': b,
                                  'author_start': fragment.get('auth_start'), 'author_end': fragment.get('auth_end')})
                segments.append({'start': a, 'end': b, 'author_start': fragment.get('auth_start'),
                                 'author_end': fragment.get('auth_end')})
            if segments:
                regions.append({'chain': entry['chain'], 'segments': segments})
    if not locations or not re.fullmatch(r'[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+', pdb):
        return None
    return {'accession': pdb, 'name': metadata.get('name', ''), 'representative': representative,
            'resolution': metadata.get('resolution'), 'method': metadata.get('experiment_type'),
            'chains': sorted({p['chain'] for p in locations}),
            'domain_match': {'database': database, 'accession': accession, 'locations': locations, 'regions': regions}}


def _domain_family_structures(info, key):
    """Keep only exemplars explicitly mapped to the requested database entry."""
    return [s for s in info.get('structures', [])
            if s.get('domain_match', {}).get('database') == key[0]
            and str(s['domain_match'].get('accession', '')).casefold() == key[1].casefold()
            and s['domain_match'].get('locations')]


def _domain_structure_refresh_needed(info, key):
    return key[0] != 'pdb' and (len(_domain_family_structures(info, key)) != len(info.get('structures', []))
        or any(not s.get('domain_match', {}).get('regions') for s in info.get('structures', [])))


def _domain_checked_info(info, key):
    if not _domain_structure_refresh_needed(info, key):
        return info
    return {**{k: v for k, v in info.items() if not k.startswith('_')},
            'structures': _domain_family_structures(info, key)}


def _domain_descriptions(metadata):
    descriptions = metadata.get("description") or []
    if isinstance(descriptions, str):
        descriptions = [{"text": descriptions}]
    return [{"text": _domain_plain_text(d.get("text", "")), "generated": bool(d.get("llm"))}
            for d in descriptions if isinstance(d, dict) and d.get("text")]


def _fetch_domain_info(accession, database="pfam", integrated_hint=None):
    """Fetch entry-level descriptions and a database representative, without sequence submission."""
    import re
    from urllib.parse import quote

    if database == "pdb":
        return _fetch_structure_publication(accession)
    if database not in set(_INTERPRO_DATABASES.values()) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,100}", accession):
        raise ValueError("Invalid domain database or accession.")
    encoded_accession = quote(accession, safe="")
    info = {"descriptions": [], "structures": [], "go_terms": [], "errors": [], "integrated": None}

    def fetch(path, label):
        try:
            return _interpro_request(path)
        except Exception as error:
            info["errors"].append(f"{label}: {error}")
            return None

    source_entry = fetch(f"/entry/{database}/{encoded_accession}/", "Member database description")
    metadata = source_entry.get("metadata", {}) if isinstance(source_entry, dict) else {}
    if metadata.get('accession') and metadata['accession'].casefold() != accession.casefold():
        raise ValueError('InterPro returned a different member accession than requested.')
    info["descriptions"] = [{**d, "source": accession, "database": database} for d in _domain_descriptions(metadata)]
    name = metadata.get("name") or {}
    info["name"] = name.get("name", "") if isinstance(name, dict) else str(name)
    integrated = metadata.get('integrated') if 'integrated' in metadata else integrated_hint
    if isinstance(integrated, str) and re.fullmatch(r"IPR\d+", integrated):
        info["integrated"] = integrated
        entry = fetch(f"/entry/interpro/{integrated}/", "InterPro description")
        ipr = entry.get("metadata", {}) if isinstance(entry, dict) else {}
        if ipr.get('accession') and ipr['accession'] != integrated:
            raise ValueError('InterPro returned a different integrated accession than requested.')
        info["descriptions"].extend({**d, "source": integrated} for d in _domain_descriptions(ipr))
        info["go_terms"] = [term for term in ipr.get("go_terms") or []
                            if isinstance(term, dict) and re.fullmatch(r"GO:\d+", term.get("identifier", ""))]

    representative = metadata.get("representative_structure") or {}
    pdb = representative.get("accession", "") if isinstance(representative, dict) else ""
    if re.fullmatch(r"[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+", pdb):
        details = fetch(f"/structure/pdb/{pdb}/entry/{database}/{encoded_accession}/", "Representative structure details")
        if isinstance(details, dict):
            structure = _domain_structure_match(details, pdb, database, accession, representative=True)
            if structure:
                if not structure['name']:
                    structure['name'] = representative.get('name', '')
                info['structures'].append(structure)
    if not info['structures']:
        examples = fetch(f"/structure/pdb/entry/{database}/{encoded_accession}/?page_size=3", "Example structures")
        if isinstance(examples, dict):
            for item in examples.get("results", [])[:3]:
                data = item.get("metadata", {})
                pdb = data.get("accession", "")
                if re.fullmatch(r"[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+", pdb):
                    structure = _domain_structure_match(item, pdb, database, accession)
                    if structure:
                        info['structures'].append(structure)
    return info


def _fetch_structure_publication(pdb):
    """Fetch the deposited primary structure citation, without searching by name."""
    import json
    import re
    from urllib.error import HTTPError
    from urllib.request import Request, urlopen

    if not isinstance(pdb, str) or not re.fullmatch(r"[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+", pdb):
        raise ValueError("Invalid PDB accession for publication lookup.")
    request = Request(f"https://data.rcsb.org/rest/v1/core/entry/{pdb.lower()}",
                      headers={"Accept": "application/json", "User-Agent": "ChimeraX-Trimmings/1.0"})
    try:
        with urlopen(request, timeout=15) as response:
            data = json.load(response)
    except HTTPError as error:
        if error.code == 404:
            return {"publications": {}}
        raise
    citation = data.get("rcsb_primary_citation") or {}
    doi = citation.get("pdbx_database_id_DOI") or ""
    doi = doi.strip() if isinstance(doi, str) else ""
    if not re.fullmatch(r"10\.\d{4,9}/\S+", doi):
        doi = ""
    pubmed = str(citation.get("pdbx_database_id_PubMed") or "")
    if not pubmed.isdigit() or int(pubmed) < 1:
        pubmed = ""
    if not doi and not pubmed:
        return {"publications": {}}
    return {"publications": {pdb.lower(): {"title": citation.get("title") or "Structure publication",
            "doi": doi, "pubmed": pubmed, "year": citation.get("year"),
            "journal": citation.get("journal_abbrev") or ""}}}


def _domain_structure_pdb(structure, record=None):
    import re

    pdb = structure.get("pdb") or structure.get("accession")
    if isinstance(pdb, str) and re.fullmatch(r"[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+", pdb):
        return pdb.lower()
    if record is not None:
        reference = record.get("structural_match", {}).get("reference", {})
        if structure.get("url") and structure["url"] == reference.get("structure_url"):
            return _domain_structure_pdb({"accession": reference.get("pdb")})
    return None


def _domain_details_html(session, record, info=None, loading=False):
    from html import escape
    from chimerax.atomic import concise_residue_spec

    if _sequence_feature_kind(record):
        return _sequence_feature_details_html(session, record, info, loading)

    if record.get('origin') == 'foldseek':
        record = _foldseek_annotation_identity(record)
        # Caller-provided enrichment may have been assembled for an old label.
        info = _domain_info_for_record(session, record)

    accession = record["accession"]
    database = record.get("database", "pfam")
    score_text = (f'{escape(record["source"])} score {_domain_number(record.get("score"))}'
                  if record.get("origin") == "interpro" else f'bit score {_domain_number(record.get("bitscore"))}')
    structures_url = record.get("structures_url", f'https://www.ebi.ac.uk/interpro/structure/pdb/entry/{database}/{accession}/')
    structure_link = f' | <a href="{escape(structures_url, quote=True)}">All matching structures</a>' if structures_url else ""
    parts = [f'<h3>{escape(_domain_display_name(record))} ({escape(accession)})</h3>',
             f'<p>{escape(concise_residue_spec(session, record["residues"]))}<br>'
             f'Query positions {record["start"]}&ndash;{record["end"]}; {len(record["residues"])} modeled residues<br>'
             f'{escape(record["type"])}; E-value {_domain_number(record.get("evalue"))}; {score_text}</p>',
             f'<p><a href="{_domain_entry_url(record)}">{escape(record.get("source", "Pfam"))} entry</a>{structure_link}</p>']
    if record.get("tentative"):
        parts.append('<p><b>Tentative Pfam match:</b> below curated gathering thresholds. '
                     'This is a candidate assignment requiring supporting evidence.</p>')
    if record.get("structural_representative"):
        parts.append('<p>Structural results transferred from sequence-identical representative '
                     f'{escape(record["structural_representative"])} with the same modeled positions. '
                     'Use deduplicate false to search each copy separately.</p>')
    if record.get("structural_support"):
        parts.append('<p><b>Foldseek evidence for this region</b></p><ul>')
        for hit in record["structural_support"]:
            match = hit["structural_match"]
            parts.append(f'<li><a href="{escape(_domain_entry_url(hit), quote=True)}">'
                         f'{escape(hit["name"])} ({escape(hit["accession"])})</a>; '
                         f'positions {hit["start"]}–{hit["end"]}; '
                         f'E-value {_domain_number(hit["evalue"])}; '
                         f'alignment TM-score {match["alntmscore"]:.3f}</li>')
        parts.append('</ul>')
    if record.get("gathering_cutoff") is not None:
        parts.append(f'<p>Pfam domain gathering cutoff: {_domain_number(record["gathering_cutoff"])} bits.</p>')
    if record.get("pfam_backend") == "local":
        parts.append(f'<p>Local Pfam release: {escape(record["pfam_release"])}; '
                     f'profile coverage: {record["profile_coverage"]:.1%}; '
                     f'composition-bias score: {_domain_number(record["composition_bias"])} bits.</p>')
        if record["query_region"] != [1, len(record["sequence"])]:
            a, b = record["query_region"]
            parts.append(f'<p>Padded sequence retry: positions {a}&ndash;{b}; '
                         'curated gathering cutoffs are unchanged.</p>')
    if record.get("origin") == "foldseek":
        match = record["structural_match"]
        if match.get("search_round", 1) > 1 or match.get("search_strategy") == "pfam_guided":
            first, last = match["search_region"]
            label = "Pfam-guided search" if match.get("search_strategy") == "pfam_guided" else "Repeat search"
            parts.append(f'<p>{label}: query positions {first}&ndash;{last} '
                         f'({match["qlen"]} modeled residues). '
                         'Query coverage below refers to this searched segment.</p>')
        parts.append(f'<p><b>Structural match candidate</b><br>'
                     f'Query coverage: {match["qcov"]:.1%}; reference coverage: {match["tcov"]:.1%}<br>'
                     f'Alignment TM-score: {match["alntmscore"]:.3f}<br>'
                     f'Matched reference: {escape(match["target"])}<br>'
                     'The name describes the reference family. Structural similarity alone does not establish '
                     'the query\'s family membership or specific function.</p>')
        if "aligned_positions" in match:
            parts.append(f'<p>Structural alignment: {len(match["aligned_positions"])} paired residues.</p>')
        if record.get("display_boundary"):
            a, b = record["display_boundary"]["original_range"]
            parts.append(f'<p>Display boundary reconciled with neighboring structural domains. '
                         f'The scoring alignment is unchanged (original query span {a}&ndash;{b}).</p>')
        if match.get("alntmscore_raw") is not None:
            parts.append(f'<p>Alignment TM-score corrected for a Foldseek normalization error '
                         f'(reported value {match["alntmscore_raw"]:.4g}).</p>')
        if record.get("pilot"):
            parts.append('<p><b>Pilot library:</b> incomplete family coverage; scores depend on this reference collection.</p>')
        mapping = record.get("pfam_assignment")
        if mapping:
            parts.append(f'<p>Reference Pfam mapping: {escape(mapping["evidence"])}.<br>'
                         f'Original classification: <a href="{escape(mapping["entry_url"], quote=True)}">'
                         f'{escape(mapping["accession"])}</a> ({escape(mapping["name"])}).</p>')
            if record.get('pfam_mapping_status') == 'conflicting':
                parts.append('<p>The localized sequence match disagrees with an upstream family-wide cross-reference. '
                             'Descriptions use the sequence-supported reference family.</p>')
        elif record.get("pfam_mapping_status"):
            label = record["pfam_mapping_status"].replace("_", " ")
            parts.append(f'<p>Pfam mapping: {escape(label)}; the ECOD classification is retained.</p>')
            parts.append('<p>This is a fold-level classification. No Pfam functional description is assigned '
                         'without a supported signature in the matched reference region.</p>')
        if record.get("pfam_candidate_accessions"):
            parts.append('<p>ECOD family-level Pfam cross-references (not established for this matched region): ' + ', '.join(
                f'<a href="https://www.ebi.ac.uk/interpro/entry/pfam/{escape(a, quote=True)}/">{escape(a)}</a>'
                for a in record["pfam_candidate_accessions"]) + '</p>')
        reference = match.get("reference", {})
        if reference.get("structure_source") in ("alphafold", "prediction"):
            label = "AlphaFold" if reference["structure_source"] == "alphafold" else "Predicted"
            parts.append(f'<p>{label} reference; mean domain pLDDT {_domain_number(reference.get("mean_plddt"))}. '
                         f'{escape(reference.get("evidence", ""))}</p>')
        if info is None:
            info = _domain_info_for_record(session, record)
    if len(record.get("segments", [])) > 1:
        parts.append('<p>Segments: ' + ', '.join(f'{a}–{b}' for a, b in record["segments"]) + '</p>')
    if record.get("description"):
        parts.append(f'<p>{escape(_domain_plain_text(record["description"]))}</p>')
    if record.get("alternatives"):
        parts.append('<p><b>Alternative structural assignments for this region</b><br>'
                     'These overlap or are contained within the stronger primary match; '
                     'alignments can differ in gaps and domain boundaries.</p><ul>')
        for alternative in record["alternatives"]:
            parts.append(f'<li><a href="{escape(_domain_entry_url(alternative), quote=True)}">'
                         f'{escape(alternative["name"])} ({escape(alternative["accession"])})</a>; '
                         f'positions {alternative["start"]}–{alternative["end"]}; '
                         f'E-value {_domain_number(alternative.get("evalue"))}; '
                         f'bit score {_domain_number(alternative.get("bitscore"))}</li>')
        parts.append('</ul>')
    if info is None:
        parts.append("<p>Loading descriptions and structure links&hellip;</p>" if loading
                     else "<p>Details have not been fetched yet.</p>")
        return "".join(parts)
    parts.append(_domain_info_html(info))
    if record.get("origin") == "foldseek" and any(
            key[0] != "pdb" and _cached_domain_info(session, key[1], key[0], key[2]) is None
            for key in _domain_info_keys(record, session)):
        parts.append("<p>Loading mapped InterPro/Pfam descriptions and example structures&hellip;</p>")
    return "".join(parts)


def _sequence_feature_details_html(session, record, info=None, loading=False):
    from html import escape
    from chimerax.atomic import concise_residue_spec

    parts = [f'<h3>{escape(record["name"])} ({escape(record["accession"])})</h3>',
             f'<p>{escape(concise_residue_spec(session, record["residues"]))}<br>'
             f'Query positions {record["start"]}&ndash;{record["end"]}; {len(record["residues"])} modeled residues<br>'
             f'{escape(record["type"])}; source: {escape(record.get("source", "Pfam"))}</p>']
    if record.get("evalue") is not None:
        parts.append(f'<p>E-value {_domain_number(record["evalue"])}; '
                     f'score {_domain_number(record.get("bitscore") if record.get("bitscore") is not None else record.get("score"))}</p>')
    elif record.get("score") is not None:
        parts.append(f'<p>Profile score {_domain_number(record["score"])}</p>')
    if record.get("parent_signature_name"):
        parts.append(f'<p>Supporting family: {escape(record["parent_signature_name"])}; '
                     f'family E-value {_domain_number(record.get("parent_evalue"))}; '
                     f'family score {_domain_number(record.get("parent_score"))}. '
                     'These scores describe the family match, not a separate site assessment.</p>')
    for key, label in (("pattern", "Pattern"), ("matched_sequence", "Matched sequence"),
                       ("prosite_release", "Local PROSITE definitions")):
        if record.get(key):
            parts.append(f'<p>{label}: <code>{escape(record[key])}</code></p>')
    if record.get("origin") == "tmsegments":
        parts.append(f'<p>TM segment number within chain: {record["tm_segment_number"]}.</p>')
    if record.get("evidence"):
        parts.append(f'<p>{escape(record["evidence"])}</p>')
    if record.get("tentative"):
        parts.append('<p>Tentative Pfam match below the curated gathering threshold.</p>')
    if record.get("common"):
        parts.append('<p>PROSITE marks this signature as common or nonspecific (SKIP-FLAG). '
                     'It is shown because prosite_common is enabled.</p>')
    if record.get("description"):
        parts.append(f'<p>{escape(_domain_plain_text(record["description"]))}</p>')
    if len(record.get("segments", [])) > 1:
        parts.append('<p>Annotated site pieces: ' + ', '.join(
            str(a) if a == b else f'{a}–{b}' for a, b in record["segments"]) + '</p>')
    parts.append(f'<p><a href="{escape(_domain_entry_url(record), quote=True)}">Database documentation</a></p>')
    if record.get("feature_kind") == "motif" and record.get("origin") == "prosite":
        overlaps = [other for other in _live_domain_records(session)
                    if other["chain"] is record["chain"] and not _sequence_feature_kind(other)
                    and len(other["residues"].intersect(record["residues"]))]
        if overlaps:
            parts.append('<p>Overlapping domain annotations (context, not confirmation): '
                         + "; ".join(escape(other["name"]) for other in overlaps) + '</p>')
    if info is not None:
        parts.append(_domain_info_html(info))
    elif loading and record.get("origin") not in ("sequence", "tmsegments"):
        parts.append('<p>Loading descriptions and structure links&hellip;</p>')
    return "".join(parts)


def _domain_info_html(info):
    """Render plain metadata without touching molecular models or Qt widgets."""
    from html import escape

    if info.get('_rendered_html_version') == _DOMAIN_INFO_RENDER_VERSION and '_rendered_html' in info:
        return info["_rendered_html"]
    parts = []
    if info.get("mapped_entries"):
        parts.append('<p><b>Family details from InterPro/Pfam</b><br>Validated reference Pfam assignments: ' + ', '.join(
            f'<a href="https://www.ebi.ac.uk/interpro/entry/pfam/{accession}/">{accession}</a>'
            for accession in info["mapped_entries"]) + '</p>')
    for ipr in info.get("integrated_entries", [info["integrated"]] if info.get("integrated") else []):
        parts.append(f'<p>Integrated InterPro entry: '
                     f'<a href="https://www.ebi.ac.uk/interpro/entry/interpro/{ipr}/">{escape(ipr)}</a></p>')
    for description in info.get("descriptions", []):
        source = description["source"]
        description_database = description.get("database", "pfam" if source.startswith("PF") else "interpro")
        label = " (AI-generated database text)" if description.get("generated") else ""
        text = _domain_plain_text(description['text'])
        if len(text) > 1400:
            text = text[:1400].rsplit(" ", 1)[0] + "…"
        description_url = description.get("entry_url", f"https://www.ebi.ac.uk/interpro/entry/{description_database}/{source}/")
        parts.append(f'<p><b>{escape(source)} description{label}</b><br>{escape(text)} '
                     f'<a href="{escape(description_url, quote=True)}">Full entry</a></p>')
    terms = info.get("go_terms", [])
    if terms:
        parts.append("<p><b>Associated GO annotations</b><br>" + "<br>".join(
            f'<a href="https://www.ebi.ac.uk/QuickGO/term/{term["identifier"]}">{term["identifier"]}</a>: '
            f'{escape(term.get("name", ""))}' for term in terms) + "</p>")
    structures = info.get("structures", [])
    if structures:
        parts.append("<p><b>Structures</b></p><ul>")
        for structure_index, structure in enumerate(structures):
            pdb = structure["accession"]
            label = structure.get("label", "Database representative" if structure["representative"] else "Example match")
            details = [str(structure["method"])] if structure.get("method") else []
            if structure.get("resolution") is not None:
                details.append(f'{structure["resolution"]} Å')
            if structure.get("chains"):
                details.append("chains " + ", ".join(structure["chains"]))
            match = structure.get('domain_match')
            if match:
                details.append('Matched ' + match['accession'])
                for location in match['locations']:
                    author = location.get('author_start') is not None and location.get('author_end') is not None
                    a, b = (location['author_start'], location['author_end']) if author else (location['start'], location['end'])
                    details.append(f"{location['chain']}: {a}–{b}" + (' (author residues)' if author else ' (structure sequence positions)'))
            name = structure.get("name", "")
            if isinstance(name, dict):
                name = name.get("name", "")
            url = structure.get("url", f"https://www.rcsb.org/structure/{pdb}")
            alternate = (f' (<a href="https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb}">PDBe</a>)'
                         if "url" not in structure else "")
            parts.append(f'<li>{label}: <a href="{escape(url, quote=True)}">{pdb.upper()}</a>{alternate}<br>'
                         f'{escape(name)}<br>{escape("; ".join(details))}')
            parts.append(f'<!--trimmings-structure-actions:{structure_index}-->')
            publication = structure.get("publication")
            if publication:
                from urllib.parse import quote
                links = []
                if publication.get("doi"):
                    links.append(f'<a href="https://doi.org/{escape(quote(publication["doi"], safe="/"), quote=True)}">Paper / DOI</a>')
                if publication.get("pubmed"):
                    links.append(f'<a href="https://pubmed.ncbi.nlm.nih.gov/{escape(publication["pubmed"], quote=True)}/">PubMed</a>')
                parts.append('<br><b>Primary structure publication:</b> ' + escape(publication["title"]) +
                             '<br>' + ' | '.join(links))
            parts.append('</li>')
        parts.append("</ul>")
    elif not info.get("errors"):
        parts.append("<p>No verified experimental structure links are available for this entry.</p>")
    if not info.get("descriptions"):
        parts.append("<p>No additional functional description is available.</p>")
    if info.get("errors"):
        parts.append("<p>Some details could not be retrieved: " + escape("; ".join(info["errors"])) +
                     ". Use Reload details to retry.</p>")
    parts.append("<p><i>Annotations describe the domain family; exact function can differ between members.</i></p>")
    info["_rendered_html"] = "".join(parts)
    info['_rendered_html_version'] = _DOMAIN_INFO_RENDER_VERSION
    return info["_rendered_html"]


def _domain_info_key(accession, database="pfam", integrated=None):
    return database, accession, integrated or ""


def _domain_info_keys(record, session=None):
    """Use explicit database cross-references, never family-name guesses."""
    import re

    if record.get("origin") in ("sequence", "tmsegments"):
        return []

    if record.get("origin") != "foldseek":
        keys = [_domain_info_key(record["accession"], record.get("database", "pfam"), record.get("integrated"))]
    else:
        reference = record.get("structural_match", {}).get("reference", {})
        # A family-wide cross-reference or old cached accession is insufficient
        # to attach function text to this particular matched domain.
        if reference.get('catalog') == 'ECOD' and reference.get('name'):
            accepted = _foldseek_reference_annotation(record['structural_match']).get('pfam_accessions', [])
            keys = [_domain_info_key(accession) for accession in accepted]
        elif reference.get('catalog') == 'ECOD' or reference.get('pfam_accessions'):
            keys = []
        else:
            database, accession = reference.get('database', record.get('database')), reference.get('accession', record.get('accession'))
            keys = [_domain_info_key(accession, database)] if database in _INTERPRO_DATABASES.values() and accession else []
    structures = list(record.get("local_info", {}).get("structures", []))
    if session is not None:
        for key in keys:
            info = _cached_domain_info(session, key[1], key[0], key[2])
            if info is not None:
                structures.extend(_domain_family_structures(info, key))
    for structure in structures:
        pdb = _domain_structure_pdb(structure, record)
        if pdb:
            keys.append(_domain_info_key(pdb, "pdb"))
    return list(dict.fromkeys(keys))


def _domain_info_for_record(session, record):
    if record.get('origin') == 'foldseek':
        record = _foldseek_annotation_identity(record)
    keys = _domain_info_keys(record, session)
    if not keys and record.get("origin") != "foldseek":
        return None
    remote = tuple(_cached_domain_info(session, key[1], key[0], key[2]) for key in keys)
    structural = record.get("origin") == "foldseek"
    local = record["local_info"] if structural else (_domain_checked_info(remote[0], keys[0]) if remote[0] is not None else None)
    if local is None:
        return None
    previous = record.get("_enriched_info")
    if (previous is not None and len(previous) == 4
            and previous[3] == (_DOMAIN_INFO_ENRICH_VERSION, tuple(keys))
            and previous[0] is local and len(previous[1]) == len(remote)
            and all(old is new for old, new in zip(previous[1], remote))):
        return previous[2]
    available = [_domain_checked_info(info, key) for key, info in zip(keys, remote) if key[0] != "pdb" and info is not None]
    publications = {pdb: publication for key, info in zip(keys, remote) if key[0] == "pdb" and info is not None
                    for pdb, publication in info.get("publications", {}).items()}
    if not publications and (not structural or not available):
        return local
    descriptions = [d for info in available for d in info.get("descriptions", [])]
    structures, seen = [], set()
    for info in [*available, local]:
        for structure in info.get("structures", []):
            identifier = structure["accession"]
            key = structure.get("url", f"https://www.rcsb.org/structure/{identifier.lower()}")
            if key not in seen:
                seen.add(key)
                publication = publications.get(_domain_structure_pdb(structure, record))
                if publication:
                    structure = {**structure, "publication": publication}
                structures.append(structure)
    terms = {term["identifier"]: term for info in available for term in info.get("go_terms", [])}
    combined = {"descriptions": descriptions or local.get("descriptions", []), "structures": structures,
                "integrated_entries": list(dict.fromkeys(info["integrated"] for info in available if info.get("integrated"))),
                "go_terms": list(terms.values()), "errors": [error for info in available for error in info.get("errors", [])]}
    if descriptions or any(info.get("structures") for info in available):
        combined["mapped_entries"] = [key[1] for key in keys if key[0] == "pfam"
                                      and key[1] in record.get('pfam_accessions', [])]
    record["_enriched_info"] = local, remote, combined, (_DOMAIN_INFO_ENRICH_VERSION, tuple(keys))
    return combined


def _cached_domain_info(session, accession, database="pfam", integrated=None):
    cache = getattr(session, "_trimmings_domain_info", {})
    key = _domain_info_key(accession, database, integrated)
    if key in cache:
        return cache[key]
    # Accept entries created before database-specific cache keys were introduced.
    if not any(isinstance(k, tuple) and len(k) == 3 and k[1] == accession for k in cache):
        return cache.get(accession)
    return None


def _domain_info_cache_dir():
    from pathlib import Path

    return Path(__file__).parent / ".domain_data" / "domain-descriptions"


def _domain_info_disk_path(key):
    from hashlib import sha256
    import json

    return _domain_info_cache_dir() / (sha256(json.dumps(key).encode()).hexdigest() + ".json")


def _read_domain_info(key):
    import json
    from time import time

    try:
        with _domain_info_disk_path(key).open(encoding="utf-8") as stream:
            text = stream.read(1_000_001)
        if len(text) > 1_000_000:
            return None
        data = json.loads(text)
        if (not isinstance(data, dict) or data.get("schema") != 1 or data.get("key") != list(key)
                or not 0 <= time() - data["fetched"] < 14 * 86400 or not isinstance(data.get("info"), dict)):
            return None
        return {k: v for k, v in data["info"].items() if not k.startswith("_")}
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _write_domain_info(key, info, fetched):
    import json
    from pathlib import Path
    import tempfile

    if info.get("errors"):
        return
    temporary = None
    try:
        path = _domain_info_disk_path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        # An older request must not replace a newer forced refresh on disk.
        if path.is_file():
            try:
                previous = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                previous = {}
            if isinstance(previous, dict) and isinstance(previous.get("fetched"), (int, float)) and previous["fetched"] > fetched:
                return
        payload = {"schema": 1, "key": list(key), "fetched": fetched,
                   "info": {k: v for k, v in info.items() if not k.startswith("_")}}
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(payload, stream)
        temporary.replace(path)
    except (OSError, ValueError, TypeError):
        pass  # A read-only cache directory must not prevent description display.
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass


def _domain_info_state(session):
    from collections import OrderedDict
    from threading import Lock
    from types import SimpleNamespace

    state = getattr(session, "_trimmings_domain_info_state", None)
    if state is None:
        state = SimpleNamespace(queue=OrderedDict(), active={}, jobs={}, disk_lock=Lock())
        session._trimmings_domain_info_state = state
        session._trimmings_domain_info_pending = None  # Invalidate an older script's request.
    return state


def _pump_domain_info(session):
    from threading import Thread

    state = _domain_info_state(session)
    while state.queue and len(state.active) < 3:
        key = next(iter(state.queue))
        job = state.queue[key]
        # Keep a slot available for an entry selected during background prefetch.
        if not job["priority"] and len(state.active) >= 2:
            break
        state.queue.pop(key)
        state.active[job["token"]] = job

        def finish(job, info):
            state.active.pop(job["token"], None)
            key = job["key"]
            if state.jobs.get(key) is job:
                state.jobs.pop(key)
                cache = getattr(session, "_trimmings_domain_info", {})
                cache[key] = info
                session._trimmings_domain_info = cache
                tool = getattr(session, "_trimmings_domain_browser", None)
                if tool is not None and not tool._closed:
                    current = tool.current()
                    if current is not None and key in _domain_info_keys(current, session):
                        tool.show_current()
                for structure in _domain_family_structures(info, key):
                    pdb = _domain_structure_pdb(structure)
                    if pdb:
                        _ensure_domain_info(session, pdb, "pdb", priority=False)
            _pump_domain_info(session)

        def lookup(job=job, finish=finish):
            key = job["key"]
            info = None if job["refresh"] else _read_domain_info(key)
            if info is not None and _domain_structure_refresh_needed(info, key):
                info = None
            if info is not None:
                try:
                    _domain_info_html(info)
                except (KeyError, TypeError, ValueError, AttributeError):
                    info = None
            if info is None:
                try:
                    info = _fetch_domain_info(key[1], key[0], key[2] or None)
                    _domain_info_html(info)
                except Exception as error:
                    info = {"errors": [str(error)]}
                    _domain_info_html(info)
                with state.disk_lock:
                    if state.jobs.get(key) is job:
                        _write_domain_info(key, info, job["created"])
            session.ui.thread_safe(finish, job, info)

        Thread(target=lookup, name="ChimeraX domain details", daemon=True).start()


def _ensure_domain_info(session, accession, database="pfam", integrated=None, *, refresh=False, priority=True):
    from time import time

    cached = _cached_domain_info(session, accession, database, integrated)
    if not refresh and cached is not None and not _domain_structure_refresh_needed(cached, _domain_info_key(accession, database, integrated)):
        return
    state = _domain_info_state(session)
    key = _domain_info_key(accession, database, integrated)
    if key in state.jobs and not refresh:
        if priority and key in state.queue:
            state.queue[key]["priority"] = True
            state.queue.move_to_end(key, last=False)
    else:
        job = {"key": key, "token": object(), "refresh": refresh, "priority": priority, "created": time()}
        state.jobs[key] = state.queue[key] = job
        if priority:
            state.queue.move_to_end(key, last=False)
    _pump_domain_info(session)


def _prefetch_domain_info(session, records):
    from threading import Thread

    pending = getattr(session, "_trimmings_domain_local_info_pending", None)
    if pending is None:
        pending = session._trimmings_domain_local_info_pending = set()
    local = {}
    for record in records:
        for key in _domain_info_keys(record, session):
            _ensure_domain_info(session, key[1], key[0], key[2], priority=False)
        if record.get("origin") == "foldseek":
            info = record["local_info"]
            if info.get('_rendered_html_version') != _DOMAIN_INFO_RENDER_VERSION and id(info) not in pending:
                local[id(info)] = info
    if local:
        pending.update(local)
        def render():
            try:
                for info in local.values():
                    _domain_info_html(info)
            finally:
                pending.difference_update(local)
        Thread(target=render, name="ChimeraX local domain descriptions", daemon=True).start()


def _domain_exemplar_module():
    return _DOMAIN_EXEMPLAR_RUNTIME


def _domain_exemplar_state(session):
    from types import SimpleNamespace
    from threading import Semaphore
    state = getattr(session, '_trimmings_domain_exemplars', None)
    if state is None:
        state = session._trimmings_domain_exemplars = SimpleNamespace(pairs={}, jobs={}, locks={}, network=Semaphore(2))
    return state


def _domain_exemplar_key(record, spec):
    import json
    source = ({k: spec[k] for k in ('provider','id','database','accession','mapped')} if spec.get('mapped') else spec)
    return record['chain'], tuple(record['residues']), json.dumps(source, sort_keys=True)


def _domain_exemplar_pair(session, record, spec):
    pair = _domain_exemplar_state(session).pairs.get(_domain_exemplar_key(record, spec))
    if pair is not None and not any(m.deleted or m.id is None for m in (pair.group, pair.domain, pair.full) if m is not None):
        return pair
    return None


def _load_domain_exemplar(session, record, structure, *, context=False, notify=None):
    from pathlib import Path
    from threading import Thread, Lock
    from time import monotonic

    module = _domain_exemplar_module()
    spec = module.source_spec(record, structure)
    if spec is None or not _domain_record_is_current(record):
        raise ValueError('This domain or its structure mapping is no longer available.')
    key, state = _domain_exemplar_key(record, spec), _domain_exemplar_state(session)
    def report(message):
        session.logger.status(message, blank_after=0)
        if notify:
            notify(message)
    pair = _domain_exemplar_pair(session, record, spec)
    if pair is not None:
        if context and pair.full is not None:
            module.show_context(pair, pair.full is not None and not pair.full.display)
        elif not context:
            module.realign_pair(session, record, pair)
        pair.group.display = True
        module.focus_pair(session, pair)
        report(f"{spec['id'].upper()}: {pair.count} residues aligned; RMSD {pair.rmsd:.2f} Å."
               + (" Full context unavailable; showing the ECOD domain." if context and pair.full is None else ""))
        return pair
    if key in state.jobs:
        report(f"{spec['id'].upper()} is already loading.")
        return None
    token, started = object(), monotonic()
    state.jobs[key] = token
    report(f"Fetching {spec['id'].upper()} domain" + (" (0s)..." if spec.get('domain_only') else " and full context (0s)..."))
    def progress(phase):
        def show():
            if state.jobs.get(key) is token:
                report(f"{spec['id'].upper()}: {phase}; {monotonic()-started:.0f}s elapsed...")
        session.ui.thread_safe(show)
    def finish(path, resolved, error):
        if state.jobs.get(key) is not token:
            return
        state.jobs.pop(key)
        if error is None and not _domain_record_is_current(record):
            error = ValueError('The query domain changed or closed while the structure was loading.')
        if error is None:
            try:
                report(f"Opening and aligning {spec['id'].upper()}...")
                pair = module.open_aligned(session, record, resolved, path)
                state.pairs[key] = pair
                if context and pair.full is not None:
                    module.show_context(pair, True)
                module.focus_pair(session, pair)
                message = (f"{spec['id'].upper()} aligned: {pair.count} residues, RMSD {pair.rmsd:.2f} Å; "
                           + ('ECOD domain retrieved; full context unavailable.' if pair.full is None else
                              'full context shown.' if context else 'full context loaded and hidden.'))
                report(message)
                session.logger.info(message)
            except Exception as failure:
                error = failure
        if error is not None:
            message = f"Could not align {spec['id'].upper()}: {error}"
            report(message)
            session.logger.warning(message)
    def worker():
        path, resolved, error = None, spec, None
        try:
            with state.network, state.locks.setdefault((spec['provider'],spec['id']), Lock()):
                if spec.get('mapped') and not spec['regions']:
                    data = _interpro_request(f"/structure/pdb/{spec['id']}/entry/{spec['database']}/{spec['accession']}/")
                    checked = _domain_structure_match(data or {}, spec['id'], spec['database'], spec['accession'])
                    if checked is None:
                        raise ValueError('No current structure-region mapping was returned.')
                    resolved = module.source_spec(record, checked)
                path, resolved = module.fetch_source(resolved, Path(__file__).parent / '.domain_data/domain-structures', progress)
        except Exception as failure:
            error = failure
        session.ui.thread_safe(finish, path, resolved, error)
    Thread(target=worker, name='ChimeraX domain exemplar', daemon=True).start()
    return None


def _create_domain_browser(session):
    from Qt.QtCore import QTimer, Signal
    from Qt.QtGui import QColor, QIcon, QPainter, QPixmap
    from Qt.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QHeaderView, QHBoxLayout, QLabel, QPushButton, QTabBar, QTextBrowser, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget
    from chimerax.core.models import ADD_MODELS, MODEL_DISPLAY_CHANGED, MODEL_ID_CHANGED, MODEL_NAME_CHANGED, MODEL_POSITION_CHANGED, REMOVE_MODELS
    from chimerax.core.selection import SELECTION_CHANGED
    from chimerax.core.tools import ToolInstance
    from chimerax.std_commands.view import NamedView
    from chimerax.ui import MainToolWindow

    class OccurrenceItem(QTreeWidgetItem):
        def text(self, column=0):
            return super().text(column)

    class OccurrenceTree(QTreeWidget):
        # Flat occurrence navigation remains independent of group expansion.
        currentRowChanged = Signal(int)

        def __init__(self):
            super().__init__()
            self.setColumnCount(2)
            self.setHeaderLabels(("Domain / feature", "Copies"))
            self.header().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
            self.header().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
            self._leaves = []
            self.currentItemChanged.connect(lambda *_: self.currentRowChanged.emit(self.currentRow()))

        def clear(self):
            self._leaves = []
            super().clear()

        def count(self):
            return len(self._leaves)

        def item(self, index):
            return self._leaves[index]

        def currentRow(self):
            item = self.currentItem()
            return item.indices[0] if item is not None else -1

        def setCurrentRow(self, index):
            item = self._leaves[index] if 0 <= index < len(self._leaves) else None
            if item is not None:
                item.parent().setExpanded(True)
            self.setCurrentItem(item)

    class DomainBrowser(ToolInstance):
        SESSION_SAVE = False
        SESSION_ENDURING = False
        _domain_browser_version = _DOMAIN_BROWSER_VERSION
        _domain_browser_implementation = _DOMAIN_BROWSER_IMPLEMENTATION

        def __init__(self):
            super().__init__(session, "Protein Domains")
            self._closed, self._rows = False, []
            self._last_details_html = None
            self._last_details_record = None
            self._row_activation = None
            self._scope, self._nearby_index, self._last_cofr_key = 0, None, None
            self._kind = 0
            self._context_values = {}
            self._model_filter, self._remembered_domains = None, {}
            self._remembered_groups = {}
            view = session.main_view
            self._starting_view = NamedView(view, view.center_of_rotation.copy(), [])
            self._starting_cofr = view.center_of_rotation.copy()
            self._starting_cofr_method = view.center_of_rotation_method
            self.tool_window = MainToolWindow(self)
            root = QVBoxLayout(self.tool_window.ui_area)
            model_row = QHBoxLayout()
            model_row.addWidget(QLabel("Model:"))
            self.previous_model_button, self.next_model_button = QPushButton("<"), QPushButton(">")
            for button, label in ((self.previous_model_button, "Previous model"), (self.next_model_button, "Next model")):
                button.setFixedWidth(30)
                button.setToolTip(label + ": show only this model")
                button.setAccessibleName(label)
            self.previous_model_button.clicked.connect(lambda _checked=False: self.step_model(-1))
            self.next_model_button.clicked.connect(lambda _checked=False: self.step_model(1))
            model_row.addWidget(self.previous_model_button)
            self.model_menu = QComboBox()
            self.model_menu.currentIndexChanged.connect(self._model_changed)
            model_row.addWidget(self.model_menu, 1)
            model_row.addWidget(self.next_model_button)
            root.addLayout(model_row)
            self.show_only_model_button = QPushButton("Show only selected model")
            self.show_only_model_button.setToolTip(
                "Show this structure and hide other atomic models. In All analyzed models, use the current domain's model."
            )
            self.show_only_model_button.clicked.connect(lambda _checked=False: self.show_only_model())
            root.addWidget(self.show_only_model_button)
            self.kind_menu = QComboBox()
            self.kind_menu.addItems(("Domains", "Sequence features", "All annotations"))
            self.kind_menu.currentIndexChanged.connect(self._kind_changed)
            root.addWidget(self.kind_menu)
            self.scope_tabs = QTabBar()
            for label in ("All domains", "Nearby", "Selected"):
                self.scope_tabs.addTab(label)
            self.scope_tabs.currentChanged.connect(self._scope_changed)
            root.addWidget(self.scope_tabs)
            self.nearby_controls = QWidget()
            radius_row = QHBoxLayout(self.nearby_controls)
            radius_row.setContentsMargins(0, 0, 0, 0)
            radius_row.addWidget(QLabel("Within"))
            self.nearby_radius = QDoubleSpinBox()
            self.nearby_radius.setRange(1, 100)
            self.nearby_radius.setDecimals(1)
            self.nearby_radius.setSingleStep(.5)
            self.nearby_radius.setSuffix(" Å")
            self.nearby_radius.setValue(6)
            self.nearby_radius.setToolTip("Minimum distance from any annotated atom to the center of rotation, in visible models.")
            self.nearby_radius.valueChanged.connect(lambda *_args: self.refresh())
            radius_row.addWidget(self.nearby_radius)
            radius_row.addWidget(QLabel("of center of rotation"))
            radius_row.addStretch(1)
            self.nearby_controls.hide()
            root.addWidget(self.nearby_controls)
            self.summary = QLabel()
            root.addWidget(self.summary)
            self.domain_list = OccurrenceTree()
            self.domain_list.setAlternatingRowColors(True)
            self.domain_list.setUniformRowHeights(True)
            self.domain_list.setMinimumWidth(360)
            self.domain_list.currentRowChanged.connect(self._row_changed)
            self.domain_list.itemClicked.connect(self._item_clicked)
            root.addWidget(self.domain_list, 1)
            navigation = QHBoxLayout()
            self.previous_button, self.next_button = QPushButton("< Previous"), QPushButton("Next >")
            self.previous_button.clicked.connect(lambda _checked=False: self.step(-1))
            self.next_button.clicked.connect(lambda _checked=False: self.step(1))
            navigation.addWidget(self.previous_button)
            navigation.addWidget(self.next_button)
            root.addLayout(navigation)
            self.auto_center = QCheckBox("Auto-center on domain")
            self.auto_center.setChecked(True)
            root.addWidget(self.auto_center)
            self.details = QTextBrowser()
            self.details.setOpenExternalLinks(False)
            self.details.setOpenLinks(False)
            self.details.anchorClicked.connect(self._structure_link_clicked)
            self._structure_actions = {}
            self.details.setMinimumHeight(220)
            root.addWidget(self.details, 2)
            actions = QHBoxLayout()
            self.focus_button = QPushButton("Focus domain")
            self.focus_button.clicked.connect(lambda _checked=False: self.focus_current())
            restore = QPushButton("Restore view")
            restore.clicked.connect(lambda _checked=False: self.restore_view())
            self.reload_button = QPushButton("Reload details")
            self.reload_button.clicked.connect(lambda _checked=False: self.reload_details())
            for button in (self.focus_button, restore, self.reload_button):
                actions.addWidget(button)
            root.addLayout(actions)
            self.status = QLabel("Click a domain to select it. Details load from InterPro/Pfam.")
            self.status.setWordWrap(True)
            root.addWidget(self.status)
            self._context_timer = QTimer(self.tool_window.ui_area)
            self._context_timer.setSingleShot(True)
            self._context_timer.setInterval(100)
            self._context_timer.timeout.connect(self._refresh_context)
            from chimerax.atomic import get_triggers
            self._handlers = [
                session.triggers.add_handler(REMOVE_MODELS, self._models_removed),
                session.triggers.add_handler(ADD_MODELS, self._model_labels_changed),
                session.triggers.add_handler(MODEL_ID_CHANGED, self._model_labels_changed),
                session.triggers.add_handler(MODEL_NAME_CHANGED, self._model_labels_changed),
                session.triggers.add_handler(MODEL_POSITION_CHANGED, self._model_geometry_changed),
                session.triggers.add_handler(MODEL_DISPLAY_CHANGED, self._model_geometry_changed),
                session.triggers.add_handler(SELECTION_CHANGED, self._selection_changed),
                session.triggers.add_handler("new frame", self._new_frame),
                get_triggers().add_handler("changes", self._atomic_changes),
            ]
            self.tool_window.manage("side")
            self.refresh()

        def current(self):
            row = self.domain_list.currentRow()
            return self._rows[row] if 0 <= row < len(self._rows) else None

        def current_copies(self):
            item = self.domain_list.currentItem()
            return item.copies if item is not None else []

        def _current_is_live(self):
            checked = set()
            for record in self.current_copies():
                key = record["chain"], id(record["sequence_residues"])
                if key not in checked:
                    if not _domain_record_is_current(record):
                        return False
                    checked.add(key)
            return True

        def browser_model(self):
            record = self.current()
            if self._model_filter is None and len({r["structure"] for r in self.current_copies()}) > 1:
                return None
            model = self._model_filter if self._model_filter is not None else record["structure"] if record else None
            return model if model is not None and not model.deleted and model.id is not None else None

        def show_only_model(self):
            from chimerax.atomic import AtomicStructure
            from chimerax.core.undo import UndoState

            model = self.browser_model()
            if model is None:
                self.refresh()
                return
            displays = {other: other is model for other in session.models.list(type=AtomicStructure)}
            # A selected structure may be hidden by one or more parent groups.
            parent = model.parent
            while parent is not None and parent.id:
                displays[parent] = True
                parent = parent.parent
            undo = UndoState("show only browser model")
            for other, shown in displays.items():
                if other.display != shown:
                    undo.add(other, "display", other.display, shown)
                    other.display = shown
            if undo.state:
                session.undo.register(undo)
            self.status.setText(f"Showing only #{model.id_string} {model.name}. Undo restores prior model visibility.")

        def _record_key(self, record):
            return (record["chain"], record["accession"], record["start"], record["end"], record["name"])

        def _remember_current(self):
            current = self.current()
            if current is not None:
                key = self._model_filter, self._scope, self._kind
                self._remembered_domains[key] = self._record_key(current)
                self._remembered_groups[key] = self.domain_list.currentItem().parent() is None

        def _kind_changed(self, *_args):
            self._remember_current()
            self._kind = self.kind_menu.currentIndex()
            self.scope_tabs.setTabText(0, ("All domains", "All features", "All annotations")[self._kind])
            label = ("domain", "feature", "annotation")[self._kind]
            self.auto_center.setText("Auto-center on " + label)
            self.focus_button.setText("Focus " + label)
            self.refresh(remember=False)

        def _scope_changed(self, *_args):
            self._remember_current()
            self._scope = self.scope_tabs.currentIndex()
            self.nearby_controls.setVisible(self._scope == 1)
            self.refresh(remember=False)

        def _cofr(self):
            import numpy as np
            center = np.asarray(session.main_view.center_of_rotation, dtype=float)
            return center if center.shape == (3,) and np.isfinite(center).all() else None

        def _schedule_context(self):
            if not self._closed and self.tool_window.shown and not self._context_timer.isActive():
                self._context_timer.start()

        def _refresh_context(self):
            if not self._closed and self.tool_window.shown:
                self.refresh()

        def _selection_changed(self, *_args):
            if self._scope == 2:
                self._schedule_context()

        def _new_frame(self, *_args):
            if self._scope != 1 or self._closed or not self.tool_window.shown:
                return
            center = self._cofr()
            key = tuple(center.round(2)) if center is not None else None
            if key != self._last_cofr_key:
                self._schedule_context()

        def _model_labels_changed(self, *_args):
            self._schedule_context()

        def _model_geometry_changed(self, *_args):
            if self._scope == 1:
                self._schedule_context()

        def _atomic_changes(self, _trigger, changes):
            # Selection/color-only changes do not invalidate the geometry index.
            if (changes.num_deleted_atoms() or len(changes.created_atoms()) or changes.residue_reasons()
                    or changes.chain_reasons() or {"coord changed", "coordset changed"} & set(changes.atom_reasons())):
                self._nearby_index = None
                self._schedule_context()

        def _context_rows(self, records):
            from chimerax.atomic import concatenate, selected_atoms
            import numpy as np
            from types import SimpleNamespace

            self._context_values = {}
            if self._scope == 0:
                return records
            if self._scope == 2:
                selected = selected_atoms(session).unique_residues
                rows = []
                for record in records:
                    count = len(record["residues"].intersect(selected))
                    if count:
                        self._context_values[id(record)] = count
                        rows.append(record)
                return rows
            center = self._cofr()
            self._last_cofr_key = tuple(center.round(2)) if center is not None else None
            if center is None or not records:
                return []
            index = self._nearby_index
            models = list(dict.fromkeys(r["structure"] for r in records))
            counts = [(model.num_atoms, model.num_residues) for model in models]
            if (index is None or len(index.records) != len(records) or index.counts != counts
                    or any(a is not b for a, b in zip(index.records, records))):
                residues = concatenate([r["residues"] for r in records], remove_duplicates=True)
                atoms = residues.atoms
                index = self._nearby_index = SimpleNamespace(
                    records=tuple(records), counts=counts, atoms=atoms, atom_residues=residues.indices(atoms.residues),
                    domain_residues=[residues.indices(r["residues"]) for r in records], residue_count=len(residues))
            delta = index.atoms.scene_coords - center
            squared = np.einsum("ij,ij->i", delta, delta)
            distances = np.full(index.residue_count, np.inf)
            valid = np.isfinite(squared)
            np.minimum.at(distances, index.atom_residues[valid], squared[valid])
            rows = []
            for record, positions in zip(records, index.domain_residues):
                if not record["structure"].visible:
                    continue
                distance = float(np.sqrt(distances[positions].min()))
                if distance <= self.nearby_radius.value():
                    self._context_values[id(record)] = distance
                    rows.append(record)
            return sorted(rows, key=lambda record: self._context_values[id(record)])

        def _model_changed(self, *_args):
            self._remember_current()
            self._model_filter = self.model_menu.currentData()
            self.refresh(remember=False)

        def select_model(self, model):
            for index in range(self.model_menu.count()):
                if self.model_menu.itemData(index) is model:
                    self.model_menu.setCurrentIndex(index)
                    return

        def step_model(self, direction):
            self.refresh()
            count = self.model_menu.count() - 1  # Skip the combined-list entry.
            if not count:
                return
            current = self.model_menu.currentIndex()
            if current == 0:
                index = 1 if direction > 0 else count
            else:
                index = 1 + (current - 1 + direction) % count
            self.model_menu.setCurrentIndex(index)
            self.show_only_model()

        def refresh(self, remember=True):
            self._context_timer.stop()
            self._row_activation = None
            if remember:
                self._remember_current()
            models = _domain_analyzed_models(session)
            if self._model_filter is not None and self._model_filter not in models:
                self._model_filter = None
            self.model_menu.blockSignals(True)
            self.model_menu.clear()
            self.model_menu.addItem("All analyzed models", None)
            model_index = 0
            for index, model in enumerate(models, 1):
                self.model_menu.addItem(f"#{model.id_string} {model.name}", model)
                if model is self._model_filter:
                    model_index = index
            self.model_menu.setCurrentIndex(model_index)
            self.model_menu.blockSignals(False)
            self.previous_model_button.setEnabled(bool(models))
            self.next_model_button.setEnabled(bool(models))
            all_records = _live_domain_records(session)
            kind = self.kind_menu.currentIndex()
            self._rows = self._context_rows([record for record in all_records
                          if (self._model_filter is None or record["structure"] is self._model_filter)
                          and (kind == 2 or bool(_sequence_feature_kind(record)) == (kind == 1))])
            groups = _domain_groups(self._rows)
            if self._scope != 1:
                groups = dict(sorted(groups.items(), key=lambda item: (item[1][0]["name"].casefold(), item[0])))
            groups = {key: sorted(copies, key=lambda r: (r["structure"].id, r["chain"].chain_id, r["start"], r["end"]))
                      for key, copies in groups.items()}
            self._rows = [record for copies in groups.values() for record in copies]
            expanded = {self.domain_list.topLevelItem(i).group_key for i in range(self.domain_list.topLevelItemCount())
                        if self.domain_list.topLevelItem(i).isExpanded()}
            collapsed = {self.domain_list.topLevelItem(i).group_key for i in range(self.domain_list.topLevelItemCount())} - expanded
            self.domain_list.blockSignals(True)
            self.domain_list.clear()
            selected = 0 if self._rows else -1
            remembered = self._remembered_domains.get((self._model_filter, self._scope, kind))
            parents = {}
            for key, copies in groups.items():
                representative = min(copies, key=lambda r: bool(r.get("tentative")))
                label = f'{representative["name"]} ({representative["accession"]})'
                count_label = f'{len(copies)} ' + ("copy" if len(copies) == 1 else "copies")
                locations = f'{count_label} in {len({r["structure"] for r in copies})} model(s), {len({r["chain"] for r in copies})} chain(s)'
                parent = OccurrenceItem([label, count_label])
                parent.setToolTip(0, label + "; " + locations)
                parent.setToolTip(1, locations)
                parent.group_key, parent.copies, parent.indices = key, copies, []
                self.domain_list.addTopLevelItem(parent)
                parent.setExpanded(key in expanded)
                parents[key] = parent
            for index, record in enumerate(self._rows):
                if self._record_key(record) == remembered:
                    selected = index
                chain_spec = record["chain"].atomspec
                if not chain_spec.startswith("#"):
                    chain_spec = f'#{record["structure"].id_string}' + chain_spec
                text = (f'{chain_spec} | {_domain_display_name(record)} ({record["accession"]}) | '
                        f'{record["start"]}–{record["end"]}')
                if record.get("origin") in ("interpro", "foldseek", "prosite", "sequence", "tmsegments"):
                    text += f' | {record["source"]}'
                if record.get("alternatives"):
                    text += f' | +{len(record["alternatives"])} alternatives'
                if self._scope == 1:
                    text += f' | {self._context_values[id(record)]:.1f} Å'
                elif self._scope == 2:
                    count = self._context_values[id(record)]
                    text += f' | {count} selected residue' + ("s" if count != 1 else "")
                swatch = QPixmap(12, 12)
                rgba = _domain_record_color(session, record)
                if rgba is None:
                    text += ' | overlapping alternative (no separate color)'
                    swatch.fill(QColor('transparent'))
                    painter = QPainter(swatch)
                    painter.setPen(QColor('gray'))
                    painter.drawRect(1, 1, 9, 9)
                    painter.drawLine(1, 10, 10, 1)
                    painter.end()
                else:
                    swatch.fill(QColor.fromRgbF(*rgba))
                key = next(iter(_domain_groups([record])))
                parent = parents[key]
                item = OccurrenceItem([text])
                item.setToolTip(0, text)
                item.setIcon(0, QIcon(swatch))
                item.indices, item.copies = [index], [record]
                parent.addChild(item)
                parent.indices.append(index)
                if parent.childCount() == 1 or (rgba is not None and not getattr(parent, '_colored', False)):
                    parent.setIcon(0, QIcon(swatch))
                parent._colored = getattr(parent, '_colored', False) or rgba is not None
                self.domain_list._leaves.append(item)
            for parent in parents.values():
                if not getattr(parent, '_colored', False):
                    parent.setText(0, parent.text(0) + ' (overlapping alternative)')
            item = self.domain_list.item(selected) if selected >= 0 else None
            if item is not None and self._remembered_groups.get((self._model_filter, self._scope, kind)):
                item = item.parent()
            self.domain_list.setCurrentItem(item)
            if item is not None and item.parent() is not None and item.parent().group_key not in collapsed:
                item.parent().setExpanded(True)
            self.domain_list.blockSignals(False)
            scope = (f"{len(models)} analyzed model(s)" if self._model_filter is None
                     else f"#{self._model_filter.id_string} {self._model_filter.name}")
            tentative = sum(bool(r.get("tentative")) for r in self._rows)
            confidence = f" ({tentative} tentative)" if tentative else ""
            context = ("", f" within {self.nearby_radius.value():g} Å of center of rotation", " overlapping selected residues")[self._scope]
            self.summary.setText(f"{len(groups)} groups; {len(self._rows)} annotation matches{confidence} in {scope}{context}")
            self.summary.setWordWrap(True)
            for button in (self.previous_button, self.next_button, self.focus_button, self.reload_button):
                button.setEnabled(bool(self._rows))
            self.status.setText("Click an annotation to select it. Feature selection preserves domain colors."
                                if kind else "Click a domain to select it. Details load from InterPro/Pfam.")
            self.show_current()
            _prefetch_domain_info(session, all_records)

        def step(self, direction):
            item = self.domain_list.currentItem()
            if item is not None and item.parent() is None:
                index = (self.domain_list.indexOfTopLevelItem(item) + direction) % self.domain_list.topLevelItemCount()
                self.domain_list.setCurrentItem(self.domain_list.topLevelItem(index))
                self.domain_list.scrollToItem(self.domain_list.currentItem())
                return
            if self._rows:
                index = (self.domain_list.currentRow() + direction) % len(self._rows)
                if index == self.domain_list.currentRow():
                    self._activate_current()
                else:
                    self.domain_list.setCurrentRow(index)
                    self.domain_list.scrollToItem(self.domain_list.currentItem())

        def _row_changed(self, *_args):
            self._row_activation = self.current()
            self._activate_current()
            QTimer.singleShot(0, self._clear_row_activation)

        def _clear_row_activation(self):
            self._row_activation = None

        def _item_clicked(self, *_args):
            # A mouse click on a different row also emits currentRowChanged.
            # Handle it once; clicking the current row again still reselects/focuses.
            if self._row_activation is not None and self._row_activation is self.current():
                self._row_activation = None
                return
            self._activate_current()

        def show_current(self):
            self._structure_actions = {}
            self.show_only_model_button.setEnabled(self.browser_model() is not None)
            record = self.current()
            if record is not None and not self._current_is_live():
                self.refresh()
                return
            if record is None:
                self._last_details_html = None
                self._last_details_record = None
                if self._kind:
                    self.details.setHtml("<p>No sequence features match this view. Run "
                                         "<b>domains #1 source features</b> for local SEG and PROSITE analysis; "
                                         "Pfam/InterPro features also appear here when available.</p>")
                    return
                self.details.setHtml("<p>No annotated domains overlap the selected residues. Select residues in an analyzed model.</p>"
                                     if self._scope == 2 else
                                     f"<p>No annotated domains within {self.nearby_radius.value():g} Å of the center of rotation in visible models.</p>"
                                     if self._scope == 1 else
                                     "<p>No domain matches were found in the analyzed chains of this model.</p>"
                                     if self._model_filter is not None else
                                     "<p>No live domain matches. Run <b>domains #1</b> to annotate a structure.</p>")
                return
            info = _domain_info_for_record(session, record)
            html = _domain_details_html(session, record, info, loading=info is None)
            if _domain_record_color(session, record) is None:
                html = ('<p><b>Overlapping alternative:</b> this annotation has no separately colored residues. '
                        'Higher-priority matches control the colors of this region.</p>' + html)
            self._structure_actions = {}
            if info is not None and not _sequence_feature_kind(record):
                import hashlib, json
                for index, structure in enumerate(info.get('structures', [])):
                    try:
                        spec = _domain_exemplar_module().source_spec(record, structure)
                    except ValueError:
                        spec = None
                    if spec is None:
                        continue
                    token = str(id(record)) + '_' + hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:16]
                    self._structure_actions[token] = record, structure
                    pair = _domain_exemplar_pair(session, record, spec)
                    align_label = 'Realign domain' if pair else 'Fetch &amp; align domain'
                    context_label = 'Hide context' if pair and pair.full is not None and pair.full.display else 'Show full context'
                    context_control = (f'<td bgcolor="#dddddd"><a href="trimmings-domain:context:{token}">{context_label}</a></td>'
                                       if not spec.get('domain_only') and not (pair and pair.full is None)
                                       else '<td>Full context unavailable</td>')
                    controls = ('<br><table cellspacing="3"><tr>'
                        f'<td bgcolor="#dddddd"><a href="trimmings-domain:align:{token}">{align_label}</a></td>'
                        + context_control + '</tr></table>')
                    html = html.replace(f'<!--trimmings-structure-actions:{index}-->', controls)
            copies = self.current_copies()
            if len(copies) > 1:
                from html import escape
                from chimerax.atomic import concise_residue_spec
                html = (f'<h3>{len(copies)} copies</h3><p>Select this group to select every listed copy; '
                        'expand it to inspect individual occurrences.</p><ul>' + ''.join(
                            f'<li>#{r["structure"].id_string} {escape(r["structure"].name)}: '
                            f'{escape(concise_residue_spec(session, r["residues"]))}</li>' for r in copies)
                        + '</ul><p><b>Example copy — individual scores below apply to this copy:</b></p>' + html)
            if html != self._last_details_html:
                scroll = self.details.verticalScrollBar()
                position = scroll.value() if self._last_details_record is record else 0
                self.details.setHtml(html)
                scroll.setValue(position)
                self._last_details_html = html
                self._last_details_record = record
            for key in _domain_info_keys(record, session):
                _ensure_domain_info(session, key[1], key[0], key[2])
            return record

        def _structure_link_clicked(self, url):
            if url.scheme().lower() in ('http','https'):
                from Qt.QtGui import QDesktopServices
                QDesktopServices.openUrl(url)
                return
            if url.scheme() != 'trimmings-domain':
                return
            action, separator, token = url.path().partition(':')
            target = self._structure_actions.get(token)
            if not separator or action not in ('align','context') or target is None:
                return
            def notify(message):
                if not self._closed:
                    self.status.setText(message)
                    self.show_current()
            try:
                _load_domain_exemplar(session, *target, context=(action=='context'), notify=notify)
            except Exception as error:
                self.status.setText(str(error))
                session.logger.warning(str(error))

        def _activate_current(self, *_args):
            from chimerax.atomic import selected_atoms, concatenate
            from chimerax.core.objects import Objects
            from chimerax.std_commands.select import select

            record = self.show_current()
            if record is None:
                return
            atoms = concatenate([r["residues"].atoms for r in self.current_copies()], remove_duplicates=True)
            selected = selected_atoms(session)
            if len(selected) != len(atoms) or len(selected.intersect(atoms)) != len(atoms):
                select(session, Objects(atoms=atoms, bonds=atoms.intra_bonds))
            if self.auto_center.isChecked():
                self._focus_record(record, self.current_copies())

        def focus_current(self):
            record = self.current()
            if record is None:
                return
            if not self._current_is_live():
                self.refresh()
                return
            self._focus_record(record, self.current_copies())

        def _focus_record(self, record, copies=None):
            from chimerax.atomic import concatenate
            model = record["structure"]
            if not any(r["structure"].visible for r in copies or [record]):
                self.status.setText(f"Model #{model.id_string} is hidden; show it to inspect this annotation.")
                return
            atoms = [r["residues"].atoms for r in copies or [record] if r["structure"].visible]
            bounds = concatenate(atoms, remove_duplicates=True).scene_bounds if atoms else None
            if bounds is None:
                return
            view, center = session.main_view, bounds.center()
            view.view_all(bounds, pad=0.15)
            normal, depth = view.camera.view_direction(), max(bounds.radius(), 1.0)
            view.clip_planes.set_clip_position("near", center - depth * normal, view)
            view.clip_planes.set_clip_position("far", center + depth * normal, view)
            view.center_of_rotation = center
            label = "feature" if _sequence_feature_kind(record) else "domain"
            self.status.setText(f"Centered on {label}. Restore view returns to the starting camera and clipping.")

        def restore_view(self):
            view = session.main_view
            self._starting_view.set_view(view, [])
            view.center_of_rotation = self._starting_cofr.copy()
            view.center_of_rotation_method = self._starting_cofr_method

        def reload_details(self):
            record = self.current()
            if record is not None:
                if record.get("origin") == "foldseek":
                    record["local_info"].pop("_rendered_html", None)
                    record.pop("_enriched_info", None)
                for key in _domain_info_keys(record, session):
                    _ensure_domain_info(session, key[1], key[0], key[2], refresh=True)
                self.show_current()

        def _models_removed(self, *_args):
            self.refresh()

        def delete(self):
            if self._closed:
                return
            self._closed = True
            self._context_timer.stop()
            for handler in self._handlers:
                handler.remove()
            self._handlers = []
            if getattr(session, "_trimmings_domain_browser", None) is self:
                session._trimmings_domain_browser = None
            super().delete()

    return DomainBrowser()


def _domain_browser_is_current(tool):
    return (getattr(tool, '_domain_browser_version', 0) == _DOMAIN_BROWSER_VERSION
            and getattr(tool, '_domain_browser_implementation', None) is _DOMAIN_BROWSER_IMPLEMENTATION)


def domain_browser(session, models=None):
    """Open or refresh domain results, optionally choosing one analyzed model."""
    from chimerax.core.errors import UserError

    if not session.ui.is_gui:
        raise UserError("The domain browser requires the ChimeraX GUI.")
    if models is not None:
        models = list(models)
    if models is not None and (len(models) != 1 or models[0] not in _domain_analyzed_models(session)):
        raise UserError("Specify one analyzed atomic model, or omit the model to reopen the browser.")
    tool = getattr(session, "_trimmings_domain_browser", None)
    if tool is not None and not tool._closed and not _domain_browser_is_current(tool):
        # A script reload can leave an older browser instance alive.
        tool.delete()
        tool = None
    if tool is None or tool._closed:
        tool = _create_domain_browser(session)
        session._trimmings_domain_browser = tool
    else:
        tool.tool_window.shown = True
        tool.refresh()
    if models is not None:
        tool.select_model(models[0])
    return tool


def domains_cancel(session):
    job = getattr(session, "_trimmings_domain_job", None)
    if job is not None:
        job.set()
        session._trimmings_domain_job = None
        if getattr(session, "_trimmings_tm_job", None) is job:
            session._trimmings_tm_job = None
        session.logger.status("Domain lookup cancelled.")
        session.logger.info("Cancelled the domain lookup locally; submitted remote searches may still finish.")
    else:
        session.logger.info("No domain lookup is running.")


def _foldseek_module():
    return _DOMAIN_SEARCH_RUNTIME


def _foldseek_segments(positions):
    segments = []
    for position in sorted(set(positions)):
        if segments and position == segments[-1][1] + 1:
            segments[-1] = (segments[-1][0], position)
        else:
            segments.append((position, position))
    return segments


def _foldseek_display_hit(hit, entry):
    """Fill alignment gaps for display, preserving scoring and explicit query pieces."""
    hit = _foldseek_annotation_identity(hit)
    match = hit.get("structural_match", {})
    paired = (_foldseek_segments(match["aligned_positions"]) if "aligned_positions" in match
              else hit.get("aligned_segments", hit.get("segments", [(hit["start"], hit["end"])])))
    pieces = hit.get("query_segments")
    if pieces is None:
        # Upgrade cached annotations produced before display ranges were stored.
        occurrence = match.get("query", "").rsplit("_domain", 1)
        regions = entry.get("structural_query_regions", [])
        if len(occurrence) == 2 and occurrence[1].isdigit() and int(occurrence[1]) < len(regions):
            pieces = regions[int(occurrence[1])]
        else:
            pieces = [(hit["start"], hit["end"])]
    segments = [(max(a, hit["start"]), min(b, hit["end"])) for a, b in pieces
                if max(a, hit["start"]) <= min(b, hit["end"])]
    return {**hit, "segments": segments, "aligned_segments": paired, "display_query_segments": pieces}


def _domain_contact_boundary(positions, coordinates, residues, low, high, default):
    """Choose a well-supported loop cut in a local C-alpha contact graph."""
    import numpy as np

    first, last = int(positions[0]), int(positions[-1])
    delta = np.zeros(last - first + 2, dtype=np.int32)
    # A contact contributes to every possible sequence cut between its endpoints.
    # Ignore local backbone neighbors; use blocks to bound temporary memory.
    for start in range(0, len(positions), 128):
        vector = coordinates[start:start + 128, None, :] - coordinates[None, :, :]
        distance = np.einsum("ijk,ijk->ij", vector, vector)
        a, b = np.nonzero((distance <= 64.) & (positions[None, :] - positions[start:start + 128, None] > 3))
        np.add.at(delta, positions[start + a] - first, 1)
        np.add.at(delta, positions[b] - first, -1)
    costs = np.cumsum(delta)

    def splits_structure(cut):
        left, right = residues[cut - 1:cut + 1]
        return (left is not None and right is not None
                and getattr(left, "ss_type", 0) in (1, 2)
                and left.ss_type == getattr(right, "ss_type", 0)
                and getattr(left, "ss_id", 0) == getattr(right, "ss_id", -1))
    choices = [p for p in range(low, high + 1) if not splits_structure(p)]
    if not choices:
        return default, "overlap priority"
    best = min(choices, key=lambda p: (costs[p - first], abs(p - default), p))
    improvement = int(costs[default - first] - costs[best - first])
    if improvement >= max(3, .2 * costs[default - first]):
        return best, "local structural contacts"
    # Preserve the existing score-based decision when contacts are inconclusive.
    # A nearby loop may still prevent an otherwise arbitrary split through an SSE.
    if splits_structure(default):
        best = min(choices, key=lambda p: (abs(p - default), costs[p - first], p))
        return best, "secondary-structure continuity"
    return default, "overlap priority"


def _domain_structural_display_boundaries(entry, hits):
    """Reconcile adjacent structural footprints without altering search evidence.

    Only crossing/adjacent contiguous assignments compete. Inserted domains,
    explicit discontinuous pieces, sequence annotations and features stay intact.
    """
    import numpy as np

    selected = [dict(hit) for hit in hits]
    structural = sorted((h for h in selected if h.get("origin") == "foldseek" and len(h["segments"]) == 1),
                        key=lambda h: (h["start"], h["end"]))
    if len(structural) < 2:
        return hits
    residues = entry["residues"]
    atoms = {i: r.find_atom("CA") for i, r in enumerate(residues, 1) if r is not None}
    atoms = {i: a for i, a in atoms.items() if a is not None}

    def allowed(hit):
        return next(((a, b) for a, b in hit.get("display_query_segments", hit.get("query_segments", hit["segments"]))
                     if a <= hit["start"] <= hit["end"] <= b), hit["segments"][0])

    def update(hit, start, end, reason):
        if (start, end) == (hit["start"], hit["end"]):
            return
        hit.setdefault("display_boundary", {"original_range": [hit["start"], hit["end"]]})
        hit["display_boundary"] = {**hit["display_boundary"], "method": reason}
        hit.update(start=start, end=end, segments=[(start, end)])

    for left, right in zip(structural, structural[1:]):
        if not (left["start"] < right["start"] <= left["end"] + 1 and left["end"] < right["end"]):
            continue
        # Bound changes to a local boundary neighborhood, keeping substantial
        # seed cores on both sides and respecting the original query pieces.
        padding = min(30, max(6, int(.2 * min(left["end"] - left["start"] + 1,
                                             right["end"] - right["start"] + 1))))
        low = max(left["start"] + 19, min(left["end"], right["start"] - 1) - padding, allowed(right)[0] - 1)
        high = min(right["end"] - 20, max(left["end"], right["start"] - 1) + padding, allowed(left)[1])
        if low > high:
            continue
        positions = np.array([p for p in range(left["start"], right["end"] + 1) if p in atoms], dtype=np.int32)
        if len(positions) < 40:
            continue
        # Missing coordinates narrow the permitted window, rather than making
        # sequence-index counts stand in for modeled core residues.
        low, high = max(low, int(positions[19])), min(high, int(positions[-21]))
        if low > high:
            continue
        default = left["end"] if _domain_hit_sort_key(left) <= _domain_hit_sort_key(right) else right["start"] - 1
        default = max(low, min(high, default))
        coordinates = np.array([atoms[int(p)].coord for p in positions], dtype=np.float64)
        cut, reason = _domain_contact_boundary(positions, coordinates, residues, low, high, default)
        update(left, left["start"], cut, reason)
        update(right, cut + 1, right["end"], reason)
    return selected


_DOMAIN_LIBRARY_RELEASE = {
    "url": "https://github.com/olibclarke/chimerax-trimmings/releases/download/domains-2026-11-17/trimmings-ecod-cdd-v295.2.tar.gz",
    "sha256": "fe9d8d6a0cf780c9d9da2b7a5cb6f5d66c24a8f03afc6164bccc3a5b1a939bd1",
    "bytes": 675634414,
}


def _domain_library_user_path():
    from pathlib import Path
    import os
    import sys
    import chimerax

    directories = getattr(chimerax, "app_dirs_unversioned", None)
    if directories is not None:
        root = Path(directories.user_data_dir)
    elif sys.platform == "darwin":
        root = Path.home() / "Library" / "Application Support" / "ChimeraX"
    elif os.name == "nt":
        root = Path(os.environ.get("LOCALAPPDATA", str(Path.home() / "AppData" / "Local"))) / "ChimeraX"
    else:
        root = Path(os.environ.get("XDG_DATA_HOME", str(Path.home() / ".local" / "share"))) / "ChimeraX"
    return root / "trimmings" / "domain-libraries" / "ecod-cdd"


def _domain_library_permission_error(error):
    import errno
    from urllib.error import URLError
    import ssl

    return (not isinstance(error, (URLError, ssl.SSLError)) and
            (isinstance(error, PermissionError) or getattr(error, "errno", None) in
             (errno.EACCES, errno.EPERM, errno.EROFS)))


def _domain_library_download_error(error):
    import ssl
    from urllib.error import HTTPError, URLError

    reason = error.reason if isinstance(error, URLError) else error
    if isinstance(reason, ssl.SSLError):
        return "HTTPS certificate or TLS verification failed. Check your system/proxy certificate configuration, then retry."
    if isinstance(error, HTTPError):
        if error.code in (404, 410):
            return f"The domain library release was not found on GitHub (HTTP {error.code}). Update Trimmings, then retry."
        if error.code in (401, 403, 407):
            return f"The library download was blocked (HTTP {error.code}). Check GitHub access and proxy/firewall settings, then retry."
        if error.code == 429:
            return "GitHub is limiting downloads (HTTP 429). Wait and rerun the command to resume."
        return f"The library server returned HTTP {error.code}. Try again later; rerun the command to resume."
    if isinstance(reason, TimeoutError):
        return "The library download timed out. Check your internet connection, then rerun the command to resume."
    return "Cannot download the domain library. Check your internet connection and proxy/firewall access to GitHub, then rerun to resume."


def _domain_library_install_error(error):
    import errno
    import ssl
    from urllib.error import URLError

    if isinstance(error, (URLError, ssl.SSLError)):
        return _domain_library_download_error(error)
    if _domain_library_permission_error(error):
        return "Cannot write the domain library. Check folder permissions or move Trimmings to a writable folder, then retry."
    if getattr(error, "errno", None) in (errno.ENOSPC, getattr(errno, "EDQUOT", -1)):
        return "Not enough disk space for the domain library. Free at least 3 GB on the destination drive, then rerun to resume."
    if isinstance(error, OSError) and not isinstance(error, (InterruptedError, TimeoutError)):
        return "Cannot read or write the domain library files. Check the destination drive and folder, then retry."
    return str(error) or type(error).__name__


def _install_domain_library(library, cancel, progress):
    """Resume, verify and atomically install the pinned public library."""
    from hashlib import sha256
    import errno
    from http.client import IncompleteRead
    import json
    import os
    from pathlib import Path, PurePosixPath
    import shutil
    import tarfile
    import tempfile
    from time import monotonic
    import ssl
    from urllib.error import HTTPError, URLError
    from urllib.request import Request, urlopen

    library = Path(library)
    if (library / "manifest.json").is_file():
        return
    library.parent.mkdir(parents=True, exist_ok=True)
    release = _DOMAIN_LIBRARY_RELEASE
    archive = library.parent / (release["sha256"] + ".tar.gz.part")
    deadline = monotonic() + 7200

    def check():
        if cancel.is_set():
            raise InterruptedError("Domain library installation cancelled; rerun the command to resume.")
        if monotonic() > deadline:
            raise TimeoutError("Domain library download timed out; rerun the command to resume.")

    # OS locks are released on exit/crash; the small lock file may remain.
    with (library.parent / ".domain-library-install.lock").open("a+b") as lock:
        if os.name == "nt":
            import msvcrt
            lock.write(b"\0")
            lock.flush()
            lock.seek(0)
        try:
            if os.name == "nt":
                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            if error.errno in (errno.EACCES, errno.EAGAIN):
                raise ValueError("Another ChimeraX instance is installing the domain library. Retry when it finishes.") from None
            raise
        check()
        if (library / "manifest.json").is_file():
            return
        if library.exists():
            raise ValueError(f"An incomplete library already exists at {library}. Move it aside and rerun the command.")
        if shutil.disk_usage(library.parent).free < 3 * 1024 ** 3:
            raise ValueError("Installing the domain library needs at least 3 GB of free disk space. Free space on the destination drive, then retry.")
        if archive.is_file() and archive.stat().st_size > release["bytes"]:
            archive.unlink()
        for attempt in range(3):
            check()
            offset = archive.stat().st_size if archive.is_file() else 0
            if offset == release["bytes"]:
                break
            headers = {"User-Agent": "ChimeraX-Trimmings/1.0", "Accept-Encoding": "identity"}
            if offset:
                headers["Range"] = f"bytes={offset}-"
            progress(f"Downloading domain library: {offset / 1e6:.0f} / {release['bytes'] / 1e6:.0f} MB")
            response, download_error = None, None
            try:
                response = urlopen(Request(release["url"], headers=headers), timeout=30)
            except (OSError, URLError) as error:
                download_error = error
            if response is not None:
                with response:
                    append = response.status == 206 and offset > 0
                    if response.status == 206:
                        content_range = response.headers.get("Content-Range", "")
                        if not content_range.startswith(f"bytes {offset}-") or not content_range.endswith(f"/{release['bytes']}"):
                            raise ValueError("The library server returned an invalid download range. Rerun the command.")
                    elif response.status != 200:
                        raise ValueError(f"Library download returned HTTP {response.status}.")
                    if not append:
                        offset = 0
                    updated = 0
                    with archive.open("ab" if append else "wb") as stream:
                        while True:
                            check()
                            try:
                                chunk = response.read(1024 * 1024)
                            except (OSError, URLError, IncompleteRead) as error:
                                download_error = error
                                break
                            if not chunk:
                                break
                            if offset + len(chunk) > release["bytes"]:
                                raise ValueError("The library download is larger than the published release.")
                            stream.write(chunk)
                            offset += len(chunk)
                            now = monotonic()
                            if now - updated >= 0.25 or offset == release["bytes"]:
                                progress(f"Downloading domain library: {100 * offset / release['bytes']:.0f}% "
                                         f"({offset / 1e6:.0f} / {release['bytes'] / 1e6:.0f} MB)")
                                updated = now
                if download_error is None and offset != release["bytes"]:
                    download_error = OSError("The library download was interrupted.")
            if download_error is None:
                break
            check()
            reason = download_error.reason if isinstance(download_error, URLError) else download_error
            if (attempt == 2 or isinstance(reason, ssl.SSLError) or
                    isinstance(download_error, HTTPError) and 400 <= download_error.code < 500):
                raise ValueError(_domain_library_download_error(download_error)) from download_error
            progress("Domain library download interrupted; retrying…")
            cancel.wait(1)
        progress("Verifying domain library checksum…")
        digest = sha256()
        with archive.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                check()
                digest.update(chunk)
        if digest.hexdigest() != release["sha256"]:
            archive.unlink()
            raise ValueError("Domain library checksum failed. Rerun the command to download a fresh copy.")
        progress("Installing domain library…")
        with tempfile.TemporaryDirectory(prefix=".domain-library-", dir=library.parent) as temporary:
            destination = Path(temporary)
            extracted, updated = 0, 0
            with tarfile.open(archive, "r|gz") as bundle:
                for member in bundle:
                    check()
                    relative = PurePosixPath(member.name)
                    if (relative.is_absolute() or ".." in relative.parts or "\\" in member.name or ":" in member.name
                            or not relative.parts or relative.parts[0] != "ecod-cdd"
                            or not (member.isdir() or member.isfile())):
                        raise ValueError("Invalid file in the domain library archive.")
                    target = destination.joinpath(*relative.parts)
                    if member.isdir():
                        target.mkdir(parents=True, exist_ok=True)
                        continue
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with bundle.extractfile(member) as source, target.open("xb") as output:
                        while chunk := source.read(1024 * 1024):
                            check()
                            output.write(chunk)
                            extracted += len(chunk)
                            now = monotonic()
                            if now - updated >= 0.25:
                                progress(f"Installing domain library: {extracted / 1e6:.0f} MB extracted")
                                updated = now
            prepared = destination / "ecod-cdd"
            manifest = json.loads((prepared / "manifest.json").read_text())
            if manifest.get("schema") != 2 or manifest.get("reference_count", 0) <= 0:
                raise ValueError("The downloaded domain library has an invalid manifest.")
            database = _fs_library_file(prepared, manifest["database"])
            for path in [database, *[Path(str(database) + suffix) for suffix in
                         (".dbtype", ".index", "_ss", "_ss.dbtype", "_ss.index", "_ca", "_ca.dbtype", "_ca.index")],
                         _fs_library_file(prepared, manifest["reference_metadata"])]:
                if not path.is_file():
                    raise ValueError("The downloaded domain library is incomplete: " + path.name)
            check()
            prepared.rename(library)
        try:
            archive.unlink()
        except OSError:
            pass  # A leftover archive must not turn a completed installation into a failure.


def _install_domain_library_then_search(session, atoms, library, options):
    """Keep installation off the GUI thread and continue the original command."""
    from threading import Event, Thread
    from time import monotonic
    from chimerax.core.errors import UserError

    cancel, started = Event(), monotonic()
    session._trimmings_domain_job = cancel
    synchronous = options["wait"] or not session.ui.is_gui

    def status(message):
        if getattr(session, "_trimmings_domain_job", None) is cancel and not cancel.is_set():
            session.logger.status(f"{message} ({int(monotonic() - started)}s elapsed)", blank_after=0)

    def progress(message):
        if synchronous:
            status(message)
        else:
            session.ui.thread_safe(status, message)

    def finish(error):
        if getattr(session, "_trimmings_domain_job", None) is not cancel:
            return
        session._trimmings_domain_job = None
        if cancel.is_set():
            return
        if error is not None:
            from html import escape
            message = "Domain library installation failed: " + _domain_library_install_error(error)
            session.logger.status(message)
            session.logger.info(
                f'<p>Library destination: <code>{escape(str(library))}</code><br>'
                f'Details: {escape(str(error.__cause__ or error))}<br>'
                f'<a href="{escape(_DOMAIN_LIBRARY_RELEASE["url"], quote=True)}">Library download</a></p>',
                is_html=True)
            if synchronous:
                raise UserError(message)
            session.logger.error(message)
            return
        session.logger.info(f"Domain library installed at {library}. Starting Foldseek search.")
        session.logger.status("Domain library installed. Starting Foldseek search…")
        try:
            _annotate_foldseek(session, atoms, str(library), **options)
        except UserError as error:
            if synchronous:
                raise
            session.logger.status(str(error))
            session.logger.error(str(error))

    def install():
        nonlocal library
        try:
            try:
                _install_domain_library(library, cancel, progress)
            except OSError as problem:
                if not _domain_library_permission_error(problem) or cancel.is_set():
                    raise
                fallback = _domain_library_user_path()
                if fallback.resolve() == library.resolve():
                    raise
                progress("Library folder is not writable; trying your ChimeraX user-data folder…")
                library = fallback
                _install_domain_library(library, cancel, progress)
            error = None
        except Exception as error_value:
            error = error_value
        if synchronous:
            finish(error)
        else:
            session.ui.thread_safe(finish, error)

    session.logger.info("Downloading the domain library on first use (~676 MB). "
                        "The search will start automatically. Use 'domains cancel' to stop; rerun to resume.")
    status("Preparing domain library download…")
    if synchronous:
        install()
    else:
        Thread(target=install, name="ChimeraX domain library install", daemon=True).start()


def _foldseek_default_library():
    from pathlib import Path
    import json

    root = Path(__file__).parent / ".domain_data"
    settings = root / "default-library.json"
    if not settings.is_file():
        if (root / 'ecod-cdd' / 'manifest.json').is_file():
            return root / 'ecod-cdd'
        if (root / 'cdd-pilot' / 'manifest.json').is_file():
            return root / 'cdd-pilot'
        fallback = _domain_library_user_path()
        if (fallback / 'manifest.json').is_file():
            return fallback
        return root / "ecod-cdd"
    data = json.loads(settings.read_text())
    relative = data.get("library") if isinstance(data, dict) else None
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("The configured structural library must be a relative directory in .domain_data.")
    library = (root / relative).resolve()
    if not library.is_relative_to(root.resolve()) or not (library / "manifest.json").is_file():
        raise ValueError("The configured structural library is missing or outside .domain_data.")
    return library


def _foldseek_group_entries(entries, deduplicate):
    """Share one search across identical sequences with identical query index maps."""
    import json

    groups = {}
    for entry in entries:
        key = (entry["sequence"], json.dumps(entry["structural_query_regions"]),
               tuple(tuple(snapshot[0]) for snapshot in entry["structural_queries"]))
        representative = groups.setdefault(key if deduplicate else entry["lookup_key"], entry)
        entry["lookup_key"] = representative["lookup_key"]
        entry["structural_representative"] = representative["spec"]
    return {entry["lookup_key"]: entry for entry in groups.values()}


def _foldseek_merge_pfam(structural, pfam):
    """Prefer curated sequence assignments; retain overlapping structural evidence."""
    curated = [{**hit, "origin": "pfam", "source": hit.get("source", "Pfam"), "database": "pfam",
                "structural_support": []} for hit in pfam
               if not _sequence_feature_kind(hit) and not hit.get("tentative") and not hit["outcompeted"]]
    features = [{**hit, "origin": "pfam", "source": hit.get("source", "Pfam"), "database": "pfam"}
                for hit in pfam if _sequence_feature_kind(hit) and not hit.get("tentative")]
    remaining = []
    for hit in _domain_collapse_structural_hits(structural):
        positions = {p for a, b in hit.get("segments", [(hit["start"], hit["end"])]) for p in range(a, b + 1)}
        candidates = []
        for i, other in enumerate(curated):
            overlap = max(0, min(hit["end"], other["end"]) - max(hit["start"], other["start"]) + 1)
            span = overlap / max(hit["end"] - hit["start"] + 1, other["end"] - other["start"] + 1)
            paired = sum(other["start"] <= p <= other["end"] for p in positions) / max(
                len(positions), other["end"] - other["start"] + 1)
            if span >= .7 and paired >= .5:
                candidates.append((span, i))
        if candidates:
            _, i = max(candidates)
            curated[i]["structural_support"].extend([hit, *hit.get("alternatives", [])])
        else:
            remaining.append(hit)
    return curated + remaining + features


def _foldseek_reference_annotation(match):
    """Promote a supported Pfam identity while preserving the ECOD search hit."""
    import re

    reference = match["reference"]
    if reference.get("catalog") != "ECOD":
        return {}
    mappings = [m for m in reference.get("pfam_mappings", [
        {"accession": a, "method": "ecod"} for a in reference.get("pfam_accessions", [])])
                if isinstance(m.get("accession"), str) and re.fullmatch(r"PF\d+", m["accession"])]
    status = reference.get("pfam_mapping_status", "unresolved")
    validated = [m for m in mappings if m.get("method") == "hmmer" and .8 <= m.get("profile_coverage", 0) <= 1]
    candidates = []
    # ECOD family names and their Pfam links can refer to the source protein,
    # rather than this structural domain. Use the fold classification until a
    # localized sequence signature supports the part that actually aligned.
    classification = dict(piece.split(": ", 1) for piece in reference.get("description", "").split("; ") if ": " in piece)
    structural_name = (reference.get("ecod_topology_name") or classification.get("Topology")
                       or classification.get("Homology group") or classification.get("Architecture")
                       or reference["name"])
    fallback = {"name": structural_name + " (structural match)",
                "pfam_mapping_status": ("conflicting" if status == "conflicting" else "partial" if validated
                                        else "ambiguous" if len({m["accession"] for m in mappings}) > 1
                                        else "unverified" if mappings else "unresolved"),
                "pfam_accessions": [],
                "pfam_candidate_accessions": list(dict.fromkeys(m["accession"] for m in mappings))}
    if validated:
        paired, position = set(), match.get("tstart", 1)
        for query, target in zip(match.get("qaln", ""), match.get("taln", "")):
            if query != "-" and target != "-":
                paired.add(position)
            position += target != "-"
        families = {}
        for mapping in validated:
            families.setdefault(mapping["accession"], []).append(mapping)
        for family_mappings in families.values():
            footprint = {p for m in family_mappings for a, b in m.get("segments", []) for p in range(a, b + 1)}
            overlap = len(paired & footprint)
            # A complete family signature can identify a larger structural
            # domain containing insertions/extensions. It must account for
            # a majority of the actual alignment, rather than a small motif.
            if overlap > .5 * len(paired) and overlap >= .5 * len(footprint):
                candidates.extend(family_mappings)
    accessions = {m["accession"] for m in candidates}
    if len(accessions) > 1:
        # Competing signatures from the same curated Pfam clan can describe the
        # same fold. Resolve only mutually overlapping alternatives, not domains
        # at distinct positions or signatures from unrelated clans.
        groups = {a: [m for m in candidates if m['accession'] == a] for a in accessions}
        clans = {m.get('clan') for m in candidates}
        footprints = {a: {p for m in group for x, y in m['segments'] for p in range(x, y + 1)} for a, group in groups.items()}
        winner = min(accessions, key=lambda a: (-sum(m.get('bitscore', 0) for m in groups[a]),
                                               min(m.get('evalue', 1) for m in groups[a]), a))
        if len(clans) == 1 and None not in clans and '' not in clans and all(
                len(footprints[winner] & footprints[a]) >= .7 * min(len(footprints[winner]), len(footprints[a]))
                for a in accessions if a != winner):
            candidates, accessions = groups[winner], {winner}
    if len(accessions) != 1:
        if len(accessions) > 1:
            fallback["pfam_mapping_status"] = "ambiguous"
        return fallback
    mapping = min(candidates, key=lambda m: m.get("evalue", 0))
    accession = mapping["accession"]
    return {"database": "pfam", "accession": accession, "color_accession": accession,
            "name": mapping["name"] + " (structural match)", "source": "Foldseek / ECOD → Pfam",
            "entry_url": f"https://www.ebi.ac.uk/interpro/entry/pfam/{accession}/",
            "structures_url": f"https://www.ebi.ac.uk/interpro/structure/pdb/entry/pfam/{accession}/",
            "pfam_accessions": [accession], "pfam_mapping_status": status,
            "pfam_assignment": {"evidence": mapping["evidence"], "accession": reference["accession"],
                                "name": reference["name"], "entry_url": reference["entry_url"]},
            "mapped_description": {"source": accession, "database": "pfam", "text": mapping.get("description", "")}}


def _foldseek_annotation_identity(hit):
    """Derive labels and description identity from the actual matched reference.

    Cached display fields are never evidence. Rebuild them together so a stale
    Pfam accession or local description cannot survive a changed assignment.
    """
    match = hit.get('structural_match', {})
    reference = match.get('reference', {})
    if hit.get('origin') != 'foldseek' or not all(k in reference for k in ('accession', 'name', 'database', 'entry_url')):
        return hit
    result = {k: v for k, v in hit.items() if k not in ('pfam_assignment', 'pfam_accessions', 'pfam_candidate_accessions',
              'pfam_mapping_status', 'integrated', '_enriched_info', 'description')}
    result.update(database=reference['database'], accession=reference['accession'], name=reference['name'] + ' (structural match)',
                  color_accession=reference['accession'], source='Foldseek / ' + reference.get('catalog', 'CDD'),
                  entry_url=reference['entry_url'], structures_url=reference.get('structures_url'), description='', annotation_identity_version=3)
    local = hit.get('local_info', {})
    identifiers = {str(reference.get(k, '')).lower() for k in ('id', 'pdb', 'ecod_domain', 'uniprot') if reference.get(k)}
    structures = [s for s in local.get('structures', []) if not s.get('representative') and (
        (s.get('url') and s['url'] == (reference.get('structure_url') or reference.get('entry_url')))
        or str(s.get('accession', '')).lower() in identifiers)]
    result['local_info'] = {'structures': structures, 'descriptions': [
        {'source': reference.get('catalog_accession', reference.get('cdd_accession', reference['accession'])),
         'text': reference.get('description', ''), 'entry_url': reference['entry_url']}]}
    annotation = _foldseek_reference_annotation(match)
    description = annotation.pop('mapped_description', None)
    result.update(annotation)
    if description and description.get('text'):
        result['local_info']['descriptions'].insert(0, description)
    if hit.get('annotation_identity_version') == 3 and all(hit.get(k) == result.get(k) for k in (
            'database', 'accession', 'name', 'pfam_accessions', 'pfam_assignment', 'pfam_mapping_status', 'pfam_candidate_accessions',
            'description', 'color_accession', 'source', 'entry_url', 'structures_url', 'integrated')) \
            and local.get('descriptions') == result['local_info']['descriptions'] and local.get('structures') == structures:
        return hit
    return result


def _foldseek_search_entries(entries, library, module, timeout, cancel, executable, progress=None, pfam_results=None):
    from pathlib import Path
    import tempfile

    results = {entry["lookup_key"]: [] for entry in entries}
    with tempfile.TemporaryDirectory(prefix="trimmings-domains-") as directory:
        queries, owners, guidance = {}, {}, {}
        for index, entry in enumerate(entries):
            for occurrence, snapshot in enumerate(entry["structural_queries"]):
                identifier = f"chain{index}_domain{occurrence}"
                positions, names, coords = snapshot
                lines = [f'ATOM  {i:5d}  CA  {name:>3s} A{i:4d}    '
                         f'{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}           C\n'
                         for i, (name, (x, y, z)) in enumerate(zip(names, coords), 1)]
                (Path(directory) / (identifier + ".pdb")).write_text("".join(lines) + "TER\nEND\n")
                queries[identifier] = {"positions": positions, "label": entry["spec"] + f" query region {occurrence + 1}"}
                if "structural_query_regions" in entry:
                    queries[identifier]["segments"] = entry["structural_query_regions"][occurrence]
                owners[identifier] = entry["lookup_key"]
                if pfam_results is not None:
                    guidance[identifier] = [(h["start"], h["end"]) for h in pfam_results[entry["sequence"]]
                                            if not _sequence_feature_kind(h) and not h.get("tentative") and not h["outcompeted"]]
        options = {"guidance": guidance} if pfam_results is not None else {}
        result = module.run_search(directory, queries, library, timeout, cancel, executable, 4, progress, **options)
    if progress is not None:
        progress("preparing domain annotations")
    for match in module.select_matches(result["hits"]):
        reference = match["reference"]
        if reference.get("structure_source") == "alphafold":
            structure = {"accession": reference["uniprot"], "name": reference["structure_title"],
                         "representative": False, "chains": [reference["chain"]],
                         "url": reference["structure_url"], "label": "Matched AlphaFold reference"}
        elif reference.get("catalog") == "ECOD":
            structure = {"accession": reference["ecod_domain"] or match["target"], "name": reference["structure_title"],
                         "pdb": reference.get("pdb"),
                         "representative": False, "chains": [reference["chain"]],
                         "url": reference["structure_url"], "label": "Matched ECOD domain"}
        elif reference.get("structure_source") in ("prediction", "reused") and not reference.get("pdb"):
            structure = {"accession": reference["id"], "name": reference["structure_title"],
                         "representative": False, "chains": [reference.get("chain", "A")],
                         "url": reference.get("structure_url") or reference["entry_url"],
                         "label": "Matched predicted reference" if reference["structure_source"] == "prediction" else "Matched structural reference"}
        else:
            structure = {"accession": reference["pdb"], "name": reference["structure_title"],
                         "representative": False, "chains": [reference["chain"]]}
        segments = _foldseek_segments(match["aligned_positions"])
        accession = reference["accession"]
        hit = {
            "origin": "foldseek", "source": "Foldseek / " + reference.get("catalog", "CDD"), "database": reference["database"],
            "accession": accession, "name": reference["name"] + " (structural match)",
            "description": "", "type": "Structural match candidate", "outcompeted": False,
            "start": segments[0][0], "end": segments[-1][1], "segments": segments,
            **({"query_segments": queries[match["query"]]["segments"]}
               if "segments" in queries[match["query"]] else {}),
            "evalue": match["evalue"], "bitscore": match["bits"], "color_accession": accession,
            "entry_url": reference["entry_url"], "structural_match": match, "pilot": result["pilot"],
            **({"structures_url": reference["structures_url"]} if "structures_url" in reference else {}),
            "local_info": {"descriptions": [{"source": reference.get("catalog_accession", reference.get("cdd_accession", accession)),
                            "text": reference["description"], "entry_url": reference["entry_url"]}],
                           "structures": [structure]},
        }
        mapped = _foldseek_reference_annotation(match)
        description = mapped.pop("mapped_description", None)
        hit.update(mapped)
        if description is not None:
            hit["local_info"]["descriptions"].insert(0, description)
        results[owners[match["query"]]].append(hit)
    if pfam_results is not None:
        results = {entry["lookup_key"]: _foldseek_merge_pfam(results[entry["lookup_key"]], pfam_results[entry["sequence"]])
                   for entry in entries}
    return results, result


def _annotate_foldseek(session, atoms, library, regions, color, refresh, timeout, wait, gui, pfam_guided=False, deduplicate=True,
                       pfam_backend="auto", pfam_database=None, feature_options=None):
    from hashlib import sha256
    import json
    from pathlib import Path
    import shutil
    from threading import Event, Thread
    from time import monotonic
    import numpy as np
    from chimerax.core.errors import UserError

    executable = shutil.which("foldseek") or next((str(p) for p in (
        Path("/opt/homebrew/bin/foldseek"), Path("/usr/local/bin/foldseek")) if p.is_file()), None)
    if executable is None:
        message = "Foldseek not found. Install it, then restart ChimeraX (macOS: brew install foldseek). See Log for help."
        session.logger.info(
            '<p><b>Foldseek is required for structural domain searches.</b></p>'
            '<p>On macOS, run <code>brew install foldseek</code> in Terminal, then restart ChimeraX. '
            'For other systems, follow the '
            '<a href="https://github.com/steineggerlab/foldseek#installation">Foldseek installation instructions</a>.</p>'
            '<p>If already installed, make the <code>foldseek</code> executable available on PATH '
            'before starting ChimeraX. The standard Homebrew locations are also checked. '
            'Use <code>domains #1 source pfam</code> to search by sequence without Foldseek.</p>',
            is_html=True)
        session.logger.status(message)
        raise UserError(message)

    try:
        module = _foldseek_module()
        automatic_library = library is None
        library = Path(library).expanduser().resolve() if library else _foldseek_default_library()
        if automatic_library and not (library / "manifest.json").is_file():
            _domain_chain_inputs(atoms)
            return _install_domain_library_then_search(session, atoms, library, dict(
                regions=regions, color=color, refresh=refresh, timeout=timeout, wait=wait, gui=gui,
                pfam_guided=pfam_guided, deduplicate=deduplicate, pfam_backend=pfam_backend,
                pfam_database=pfam_database, feature_options=feature_options))
        manifest_text = (library / "manifest.json").read_text()
        manifest = json.loads(manifest_text)
        guidance_backend = _pfam_backend(pfam_backend, pfam_database, library) if pfam_guided else None
        pieces = module.parse_segments(regions) if regions is not None else None
        entries = _domain_chain_inputs(atoms)
        unique, geometry = {}, []
        for entry in entries:
            domains = pieces or [[(1, len(entry["sequence"]))]]
            if any(b > len(entry["sequence"]) for domain in domains for a, b in domain):
                raise ValueError(f"Domain regions exceed the sequence length of {entry['spec']}.")
            snapshots, targeted = [], set()
            for domain in domains:
                selected = [(i, residue, residue.find_atom("CA")) for i, residue in enumerate(entry["residues"], 1)
                            if residue is not None and residue in entry["targets"] and any(a <= i <= b for a, b in domain)]
                selected = [(i, residue, atom) for i, residue, atom in selected if atom is not None]
                if not 20 <= len(selected) <= 9999:
                    raise ValueError(f"Each query region in {entry['spec']} needs 20 to 9,999 modeled C-alpha atoms.")
                coords = np.array([atom.coord for _, _, atom in selected], dtype=np.float64)
                if not np.isfinite(coords).all() or (coords < -999).any() or (coords > 9999).any():
                    raise ValueError("Coordinates cannot be exported safely for Foldseek.")
                names = [r.name if len(r.name) <= 3 else "UNK" for _, r, _ in selected]
                snapshots.append(([i for i, _, _ in selected], names, coords.tolist()))
                geometry.extend((atom, np.array(atom.coord, copy=True)) for _, _, atom in selected)
                targeted.update(r for _, r, _ in selected)
            entry["targets"] = targeted
            entry["structural_queries"] = snapshots
            entry["structural_query_regions"] = domains
            entry["lookup_key"] = sha256((str(module.SEARCH_VERSION) + manifest_text + entry["sequence"]
                                           + json.dumps(snapshots) + json.dumps(domains) + "pfam-reference-labels-v3"
                                           + ("pfam-guided-v3" + guidance_backend["identity"] if pfam_guided else "")
                                           + ("separate:" + entry["spec"] if not deduplicate else "")).encode()).hexdigest()
        unique = _foldseek_group_entries(entries, deduplicate)
    except (OSError, ValueError, KeyError) as error:
        raise UserError(f"Cannot prepare Foldseek search: {error}") from None
    cache = getattr(session, "_trimmings_foldseek_cache", {})
    cached = {key: cache[key] for key in unique if not refresh and key in cache}
    missing = [entry for key, entry in unique.items() if key not in cached]
    model_count = len({entry["chain"].structure for entry in entries})
    if not missing and feature_options is None:
        records = _apply_domain_results(session, entries, cached, color, False, gui, source="foldseek")
        message = f"Foldseek complete (cached): {len(records)} matches across {model_count} model(s)."
        session.logger.status(message)
        session.logger.info(message)
        return
    cancel = Event()
    started = monotonic()
    session._trimmings_domain_job = cancel
    reference_count = manifest["reference_count"] if manifest.get("schema") == 2 else len(manifest["records"])
    session.logger.info((f"Searching {reference_count} local structural references with Foldseek. " if missing else
                         "Reusing cached Foldseek domain results; checking sequence features. ")
                        + ("Pilot library: coverage is limited. " if manifest.get("pilot") else "")
                        + f"{len(entries)} chain(s), {len(unique)} representative search(es). "
                        + "Use 'domains cancel' to stop waiting.")

    def update_status(phase):
        # Queued progress must not overwrite completion, cancellation, or a newer search.
        if getattr(session, "_trimmings_domain_job", None) is cancel and not cancel.is_set():
            seconds = int(monotonic() - started)
            session.logger.status(f"Foldseek running: {model_count} model(s), {seconds}s elapsed; {phase}...",
                                  blank_after=0)

    def progress(phase):
        if wait or not session.ui.is_gui:
            update_status(phase)
        else:
            session.ui.thread_safe(update_status, phase)

    update_status("preparing structure queries" if missing else "cached domains; checking sequence features")

    def finish(results, result, error):
        if getattr(session, "_trimmings_domain_job", None) is not cancel:
            return
        session._trimmings_domain_job = None
        if cancel.is_set():
            return
        if error is None and any(atom.deleted or not np.array_equal(atom.coord, xyz) for atom, xyz in geometry):
            error = ValueError("The query coordinates changed during the search; no annotations were applied.")
        if error is not None:
            state = "timed out" if isinstance(error, TimeoutError) else "failed"
            session.logger.status(f"Foldseek {state} after {monotonic() - started:.1f}s; see Log for details.")
            if wait or not session.ui.is_gui:
                raise UserError(f"Foldseek search failed: {error}")
            session.logger.error(f"Foldseek search failed: {error}")
            return
        cache = getattr(session, "_trimmings_foldseek_cache", {})
        if not result.get("guidance_error") and result.get("guidance_cacheable", True):
            cache.update(results)
        session._trimmings_foldseek_cache = cache
        if result.get("pfam_results") is not None:
            _pfam_cache_results(session, result["pfam_results"], result["pfam_backend"])
        if result.get("guidance_note"):
            session.logger.warning(result["guidance_note"])
        if result.get("guidance_error"):
            session.logger.warning("Pfam guidance unavailable; used the full local search: " + result["guidance_error"])
        for warning in result.get("feature_warnings", []):
            session.logger.warning(warning)
        records = _apply_domain_results(session, entries, {**cached, **results}, color, False, gui, source="foldseek",
                                        feature_results=result.get("feature_results"))
        duration = f"in {monotonic() - started:.1f}s" if missing else "(cached)"
        message = f"Foldseek complete {duration}: {_domain_completion_counts(records)} across {model_count} model(s)."
        session.logger.status(message)
        session.logger.info(message + " Names describe structurally matched reference families; inspect coverage and alignments.")

    def search():
        feature_task = None
        try:
            feature_task = _start_domain_features(entries, feature_options, timeout, cancel, progress)
            pfam_results, guidance_error, guidance_note = None, None, None
            used_backend = guidance_backend
            if pfam_guided and missing:
                sequences = list(dict.fromkeys(entry["sequence"] for entry in missing))
                cache = _pfam_cached(session, guidance_backend)
                pfam_results = {seq: cache[seq] for seq in sequences if not refresh and seq in cache}
                needed = [seq for seq in sequences if seq not in pfam_results]
                try:
                    if needed:
                        fresh, used_backend, guidance_note = _pfam_lookup(
                            needed, min(60, timeout / 3), cancel, progress, guidance_backend, all_sequences=sequences)
                        if used_backend["identity"] != guidance_backend["identity"]:
                            pfam_results = {}
                        pfam_results.update(fresh)
                except Exception as err:
                    if cancel.is_set():
                        raise
                    guidance_error, pfam_results = str(err) or type(err).__name__, None
                    progress("Pfam unavailable; starting full local search")
            if missing:
                remaining = timeout - (monotonic() - started)
                results, result = _foldseek_search_entries(missing, library, module, remaining, cancel, executable, progress,
                                                           pfam_results=pfam_results)
            else:
                results, result = {}, {}
            result.update(pfam_results=pfam_results, guidance_error=guidance_error, guidance_note=guidance_note,
                          pfam_backend=used_backend, guidance_cacheable=(not pfam_guided or
                          used_backend["identity"] == guidance_backend["identity"]))
            result["feature_results"], result["feature_warnings"] = _await_domain_features(feature_task, progress)
            error = None
        except Exception as err:
            if feature_task is not None:
                feature_task.stopped.set()
            results, result, error = None, None, err
        if wait or not session.ui.is_gui:
            finish(results, result, error)
        else:
            session.ui.thread_safe(finish, results, result, error)

    if wait or not session.ui.is_gui:
        search()
    else:
        Thread(target=search, name="ChimeraX Foldseek domains", daemon=True).start()


def _annotate_sequence_features(session, atoms, refresh, timeout, wait, gui,
                                prosite_profiles, prosite_common, prosite_database):
    from pathlib import Path
    from threading import Event, Thread
    from time import monotonic
    from types import SimpleNamespace
    from chimerax.core.errors import UserError

    entries = _domain_chain_inputs(atoms)
    sequences = list(dict.fromkeys(entry["sequence"] for entry in entries))
    try:
        module = _sequence_feature_backend()
        if prosite_database is not None and not Path(prosite_database).expanduser().is_file():
            raise ValueError("prosite_database must name a readable prosite.dat file.")
    except (OSError, ValueError) as error:
        raise UserError("Cannot prepare sequence features: " + str(error)) from None
    cancel, stopped = Event(), Event()
    started, phase = monotonic(), SimpleNamespace(text="preparing local scans")
    session._trimmings_domain_job = cancel
    session.logger.info(f"Sequence features: {len(sequences)} distinct sequence(s). "
                        + ("Online profile enrichment sends sequences to ScanProsite. " if prosite_profiles else
                           "SEG and PROSITE patterns run locally; only public definitions may be downloaded. ")
                        + "Use 'domains cancel' to stop waiting.")
    for entry in entries:
        if not entry["chain"].full_sequence_known:
            session.logger.warning(f"{entry['spec']}: features use the inferred chain sequence; "
                                   "missing sequence can affect motifs and low-complexity boundaries, "
                                   "and terminal matches may refer to fragment ends.")

    def update_status():
        if getattr(session, "_trimmings_domain_job", None) is cancel and not cancel.is_set():
            session.logger.status(f"Sequence features running: {int(monotonic() - started)}s elapsed; {phase.text}...",
                                  blank_after=0)

    def report(message):
        phase.text = message
        if wait or not session.ui.is_gui:
            update_status()
        else:
            session.ui.thread_safe(update_status)

    def heartbeat():
        while not stopped.wait(1) and not cancel.is_set():
            session.ui.thread_safe(update_status)

    def finish(results, warnings, error):
        if getattr(session, "_trimmings_domain_job", None) is not cancel:
            return
        session._trimmings_domain_job = None
        if cancel.is_set():
            return
        if error is not None:
            session.logger.status("Sequence-feature lookup failed.")
            if wait or not session.ui.is_gui:
                raise UserError("Sequence-feature lookup failed: " + str(error)) from None
            session.logger.error("Sequence-feature lookup failed: " + str(error))
            return
        for warning in warnings:
            session.logger.warning(warning)
        records = _apply_domain_results(session, entries, results, False, True, gui=gui, source="features",
                                        include_tentative=True)
        if gui and session.ui.is_gui:
            session._trimmings_domain_browser.kind_menu.setCurrentIndex(1)
        message = f"Sequence features complete in {monotonic() - started:.1f}s: {len(records)} matches across {len(entries)} chain(s)."
        session.logger.status(message)
        session.logger.info(message)

    def lookup():
        if session.ui.is_gui and not wait:
            Thread(target=heartbeat, name="ChimeraX feature progress", daemon=True).start()
        try:
            results, warnings = module.search(sequences, Path(__file__).parent / ".domain_data" / "sequence-features",
                                              timeout, cancel, report, database=prosite_database,
                                              common=prosite_common, profiles=prosite_profiles, refresh=refresh)
            error = None
        except Exception as err:
            results, warnings, error = None, [], err
        finally:
            stopped.set()
        if wait or not session.ui.is_gui:
            finish(results, warnings, error)
        else:
            session.ui.thread_safe(finish, results, warnings, error)

    update_status()
    if wait or not session.ui.is_gui:
        lookup()
    else:
        Thread(target=lookup, name="ChimeraX sequence features", daemon=True).start()


def annotate_domains(session, atoms, color=True, include_overlaps=False,
                     refresh=False, timeout=None, wait=False, gui=True, source="pfam", library=None, regions=None,
                     include_tentative=None, pfam_guided=False, deduplicate=True,
                     pfam_backend="auto", pfam_database=None,
                     prosite_profiles=False, prosite_common=False, prosite_database=None, features=True,
                     results=None, python=None):
    """Look up InterPro matches or search Pfam, and annotate only requested residues."""
    from math import isfinite
    from threading import Event, Thread
    from time import monotonic
    from chimerax.core.errors import UserError

    if timeout is None:
        timeout = 1200 if source == "tmsegments" else 600
    if not isfinite(timeout) or timeout <= 0:
        raise UserError("timeout must be finite and greater than zero (seconds).")
    if source not in ("pfam", "interpro", "foldseek", "features", "tmsegments"):
        raise UserError("source must be pfam, interpro, foldseek, features, or tmsegments.")
    if (results is not None or python is not None) and source != "tmsegments":
        raise UserError("results and python apply only to source tmsegments.")
    if not features and (source == "features" or prosite_profiles or prosite_common or prosite_database is not None):
        raise UserError("features false cannot be combined with source features or PROSITE scan options.")
    if prosite_database is not None:
        from pathlib import Path
        if not Path(prosite_database).expanduser().is_file():
            raise UserError("prosite_database must name a readable prosite.dat file.")
    if pfam_backend not in ("auto", "local", "online"):
        raise UserError("pfam_backend must be auto, local, or online.")
    if (pfam_backend != "auto" or pfam_database is not None) and not (
            source == "pfam" or source == "foldseek" and pfam_guided):
        raise UserError("pfam_backend and pfam_database apply to Pfam searches or Pfam-guided Foldseek searches.")
    if include_tentative is None:
        include_tentative = source == "pfam"
    if include_tentative and source != "pfam":
        raise UserError("include_tentative applies only to source pfam.")
    if getattr(session, "_trimmings_domain_job", None) is not None:
        raise UserError("A domain lookup is already running; use 'domains cancel' before starting another.")
    if source == "tmsegments":
        if library is not None or regions is not None or pfam_guided or not deduplicate:
            raise UserError("library, regions, pfam_guided and deduplicate apply only to source foldseek.")
        if prosite_profiles or prosite_common or prosite_database is not None:
            raise UserError("PROSITE options apply to local feature scans; use source features for those scans.")
        return annotate_tm_segments(session, atoms, color=False, results=results, python=python,
                                    refresh=refresh, timeout=timeout, wait=wait, gui=gui, _domain_source=True)
    if source == "features":
        if library is not None or regions is not None or pfam_guided or not deduplicate:
            raise UserError("library, regions, pfam_guided and deduplicate apply only to source foldseek.")
        return _annotate_sequence_features(session, atoms, refresh, timeout, wait, gui,
                                          prosite_profiles, prosite_common, prosite_database)
    if source == "foldseek":
        return _annotate_foldseek(session, atoms, library, regions, color, refresh, timeout, wait, gui, pfam_guided, deduplicate,
                                  pfam_backend, pfam_database,
                                  feature_options=({"database": prosite_database, "common": prosite_common,
                                                    "profiles": prosite_profiles, "refresh": False} if features else None))
    if pfam_guided or not deduplicate:
        raise UserError("pfam_guided and deduplicate apply only to source foldseek.")
    if library is not None or regions is not None:
        raise UserError("library and regions apply only to source foldseek.")
    entries = _domain_chain_inputs(atoms)
    sequences = list(dict.fromkeys(entry["sequence"] for entry in entries))
    try:
        backend = _pfam_backend(pfam_backend, pfam_database) if source == "pfam" else None
    except (OSError, ValueError, KeyError) as error:
        raise UserError(f"Cannot prepare Pfam lookup: {error}") from None
    cache_name = "_trimmings_interpro_cache" if source == "interpro" else "_trimmings_pfam_cache"
    cache = getattr(session, cache_name, {})
    if source == "interpro" and getattr(session, "_trimmings_interpro_cache_version", None) != _INTERPRO_CACHE_VERSION:
        cache = {}
    if source == "pfam":
        cache = _pfam_cached(session, backend)
    cached = {seq: cache[seq] for seq in sequences if not refresh and seq in cache}
    missing = [seq for seq in sequences if seq not in cached]
    if not missing and not features:
        _apply_domain_results(session, entries, cached, color, include_overlaps, gui=gui, source=source,
                              include_tentative=include_tentative)
        if source == "pfam":
            session.logger.info("Pfam results reused from cache (" + backend["kind"] + ").")
        return
    cancel = Event()
    started = monotonic()
    feature_options = ({"database": prosite_database, "common": prosite_common,
                        "profiles": prosite_profiles, "refresh": False} if features else None)
    session._trimmings_domain_job = cancel
    if not missing:
        description = f"Reusing cached {source} domains and checking sequence features"
    elif source == "interpro":
        description = "Looking up exact sequence hashes in InterPro Matches"
    elif backend["kind"] == "local":
        description = "Searching local Pfam " + backend["release"]
    else:
        description = "Searching Pfam at EMBL-EBI HMMER (curated gathering thresholds)"
    sequence_count = len(missing) if missing else len(sequences)
    session.logger.info(f"{description}: {sequence_count} distinct sequence(s). Use 'domains cancel' to stop waiting.")
    for entry in entries:
        if not entry["chain"].full_sequence_known:
            session.logger.warning(
                f"{entry['spec']}: full polymer sequence is unavailable; using ChimeraX's inferred chain sequence. "
                "Missing regions may prevent an exact InterPro lookup or reduce Pfam detection."
            )

    def finish(results, error, used_backend=None, note=None, feature_results=None, feature_warnings=()):
        # Cancellation/new requests invalidate queued callbacks as well as polling.
        if getattr(session, "_trimmings_domain_job", None) is not cancel:
            return
        session._trimmings_domain_job = None
        if cancel.is_set():
            return
        if error is not None:
            if wait or not session.ui.is_gui:
                raise UserError(f"Domain lookup failed: {error}")
            session.logger.error(f"Domain lookup failed: {error}")
            return
        if note:
            session.logger.warning(note)
        for warning in feature_warnings:
            session.logger.warning(warning)
        cache = getattr(session, cache_name, {})
        if source == "pfam":
            _pfam_cache_results(session, results, used_backend)
            effective_cache = cached if used_backend["identity"] == backend["identity"] else {}
        else:
            cache = dict(cache) if getattr(session, "_trimmings_interpro_cache_version", None) == _INTERPRO_CACHE_VERSION else {}
            cache.update(results)
            setattr(session, cache_name, cache)
            session._trimmings_interpro_cache_version = _INTERPRO_CACHE_VERSION
            effective_cache = cached
        records = _apply_domain_results(session, entries, {**effective_cache, **results}, color, include_overlaps, gui=gui, source=source,
                                        include_tentative=include_tentative, feature_results=feature_results)
        if source == "pfam":
            label = "local " + used_backend["release"] if used_backend["kind"] == "local" else "online"
            session.logger.info("Pfam complete (" + label + ").")
        name = "Pfam" if source == "pfam" else "InterPro"
        cached_note = " (cached domains)" if not missing else ""
        message = f"{name} complete{cached_note} in {monotonic() - started:.1f}s: {_domain_completion_counts(records)}."
        session.logger.status(message)
        session.logger.info(message)

    def lookup():
        used_backend, note, feature_task = backend, None, None
        feature_results, feature_warnings = None, []
        try:
            def update_status(message):
                if getattr(session, "_trimmings_domain_job", None) is cancel and not cancel.is_set():
                    name = "Pfam" if source == "pfam" else "InterPro"
                    session.logger.status(f"{name} running: {int(monotonic() - started)}s elapsed; {message}...", blank_after=0)
            def report(message):
                if wait or not session.ui.is_gui:
                    update_status(message)
                else:
                    session.ui.thread_safe(update_status, message)
            report("preparing lookups" if missing else "cached domains; checking sequence features")
            feature_task = _start_domain_features(entries, feature_options, timeout, cancel, report)
            if not missing:
                results = {}
            elif source == "pfam":
                results, used_backend, note = _pfam_lookup(missing, timeout, cancel, report, backend, all_sequences=sequences)
            else:
                results = _interpro_search(missing, timeout, cancel, report)
            feature_results, feature_warnings = _await_domain_features(feature_task, report)
            error = None
        except Exception as err:
            if feature_task is not None:
                feature_task.stopped.set()
            results, error = None, err
        if wait or not session.ui.is_gui:
            finish(results, error, used_backend, note, feature_results, feature_warnings)
        else:
            session.ui.thread_safe(finish, results, error, used_backend, note, feature_results, feature_warnings)

    if wait or not session.ui.is_gui:
        lookup()
    else:
        Thread(target=lookup, name="ChimeraX domains", daemon=True).start()


_TM_APP = "DTU/DeepTMHMM"
_TM_COLORS = {"M": (255, 145, 45, 255), "B": (40, 200, 230, 255), "": (190, 190, 190, 255)}
_TM_NAMES = {"M": "TM helix", "B": "TM beta strand"}


def _tm_parse_topologies(text):
    """Read DeepTMHMM 1.x three-line output, keyed by exact sequence.

    M and B are distinct predictions. S is a signal peptide, never a TM segment.
    Do not infer topology from the model's secondary-structure assignment.
    """
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or len(lines) % 3:
        raise ValueError("Expected DeepTMHMM predicted_topologies.3line output (header, sequence, topology).")
    predictions, identifiers = {}, set()
    for index in range(0, len(lines), 3):
        header, sequence, topology = lines[index:index + 3]
        sequence = sequence.upper()
        if not header.startswith(">") or not header[1:].strip():
            raise ValueError("Invalid DeepTMHMM sequence header.")
        identifier = header[1:].split()[0]
        if identifier in identifiers:
            raise ValueError(f"Duplicate DeepTMHMM sequence identifier: {identifier}")
        identifiers.add(identifier)
        if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
            raise ValueError(f"Invalid amino-acid sequence in DeepTMHMM result for {identifier}.")
        if len(sequence) != len(topology) or set(topology) - set("IOMBPS"):
            raise ValueError(f"Invalid DeepTMHMM topology or sequence/topology length mismatch for {identifier}.")
        if sequence in predictions and predictions[sequence] != topology:
            raise ValueError("Conflicting DeepTMHMM predictions for identical sequences.")
        predictions[sequence] = topology
    return predictions


def _tm_segments(topology):
    """Return inclusive, one-based ranges for contiguous M or B runs."""
    from itertools import groupby

    segments, offset = [], 1
    for label, group in groupby(topology):
        length = sum(1 for _ in group)
        if label in _TM_NAMES:
            segments.append({"label": label, "start": offset, "end": offset + length - 1,
                             "number": len(segments) + 1})
        offset += length
    return segments


def _tm_feature_hits(topology):
    """Expose the existing M/B topology segments as ordinary sequence features."""
    hits = []
    for segment in _tm_segments(topology):
        label = segment["label"]
        kind = "helix" if label == "M" else "strand"
        hits.append({"origin": "tmsegments", "feature_kind": "transmembrane_" + kind,
                     "database": "deeptmhmm", "source": "DeepTMHMM", "accession": "DeepTMHMM:" + label,
                     "name": _TM_NAMES[label], "type": "Predicted transmembrane " + kind,
                     "start": segment["start"], "end": segment["end"],
                     "tm_label": label, "tm_segment_number": segment["number"],
                     "evalue": None, "bitscore": None, "outcompeted": False,
                     "entry_url": "https://biolib.com/DTU/DeepTMHMM/", "structures_url": "",
                     "description": "DeepTMHMM sequence-based topology prediction. Inspect the segment in its structural context."})
    return hits


def _tm_biolib_worker():
    """Private child entry point; run in an external Python with pybiolib installed."""
    from pathlib import Path
    try:
        import biolib
    except ImportError:
        raise RuntimeError("Install pybiolib in the external Python environment: python -m pip install pybiolib") from None

    app = biolib.load(_TM_APP)
    job = app.cli(args="--fasta input.fasta", check=True)
    # Avoid downloading embeddings, plots, and other large outputs.
    job.save_files("output", path_filter=lambda path: Path(path).name == "predicted_topologies.3line")
    files = list(Path("output").rglob("predicted_topologies.3line"))
    if len(files) != 1:
        raise RuntimeError("DeepTMHMM did not return one predicted_topologies.3line file.")
    Path("topologies.3line").write_text(files[0].read_text(encoding="utf-8"), encoding="utf-8")


def _tm_python_executable(python=None):
    import os
    import shutil
    from pathlib import Path

    configured = python or os.environ.get("CHIMERAX_DEEPTMHMM_PYTHON")
    default = Path.home() / ".local/share/chimerax-trimmings/deeptmhmm/bin/python"
    candidate = str(Path(configured).expanduser()) if configured else str(default)
    resolved = shutil.which(candidate)
    if resolved is None:
        raise ValueError("BioLib Python was not found. Create the environment described in README.MD, "
                         "or supply 'python /absolute/path/to/venv/bin/python'. "
                         "Alternatively import a prediction with 'results /path/predicted_topologies.3line'.")
    # Do not resolve the symlink: a venv's Python must retain its venv path.
    return str(Path(resolved).absolute())


def _tm_predict(sequences, executable, timeout, cancel, status=None):
    """Only strings/files cross the process boundary; the child never sees models."""
    import os
    from pathlib import Path
    import subprocess
    import tempfile
    from time import monotonic

    environment = os.environ.copy()
    for key in ("PYTHONHOME", "PYTHONPATH", "PYTHONSTARTUP", "PYTHONEXECUTABLE", "__PYVENV_LAUNCHER__"):
        environment.pop(key, None)
    bootstrap = "import runpy,sys; runpy.run_path(sys.argv[1])['_tm_biolib_worker']()"
    script = str(Path(__file__).absolute())
    started = monotonic()
    deadline, next_status = started + timeout, started + 1
    with tempfile.TemporaryDirectory(prefix="chimerax-deeptmhmm-") as directory:
        directory = Path(directory)
        (directory / "input.fasta").write_text(
            "".join(f">query_{i}\n{seq}\n" for i, seq in enumerate(sequences)), encoding="utf-8")
        with (directory / "client.log").open("w+", encoding="utf-8") as log:
            if cancel.is_set():
                raise RuntimeError("DeepTMHMM prediction cancelled locally.")
            process = subprocess.Popen([executable, "-u", "-c", bootstrap, script], cwd=directory,
                                       env=environment, stdin=subprocess.DEVNULL, stdout=log, stderr=log)
            try:
                while process.poll() is None:
                    if cancel.wait(0.1):
                        raise RuntimeError("DeepTMHMM prediction cancelled locally; the remote job may continue.")
                    now = monotonic()
                    if now >= deadline:
                        raise TimeoutError("DeepTMHMM timed out; the remote job may continue. Increase timeout to wait longer.")
                    if status is not None and now >= next_status:
                        status(int(now - started))
                        next_status = now + 1
                if process.returncode:
                    log.seek(0)
                    detail = log.read()[-2000:].strip()
                    raise RuntimeError("BioLib client failed: " + (detail or f"exit code {process.returncode}"))
                text = (directory / "topologies.3line").read_text(encoding="utf-8")
            finally:
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
        predictions = _tm_parse_topologies(text)
        if set(predictions) != set(sequences):
            raise ValueError("DeepTMHMM returned sequences that do not exactly match the submitted sequences.")
        return predictions


def _apply_tm_results(session, entries, predictions, color=True, elapsed=None, gui=True):
    from html import escape
    from types import SimpleNamespace
    from urllib.parse import quote
    from chimerax.atomic import Residue, Residues, concise_residue_spec
    from chimerax.core.colors import Color
    from chimerax.core.objects import Objects
    from chimerax.core.undo import UndoState
    from chimerax.std_commands.color import color as color_objects

    for name, kind in (("tm_topology", str), ("tm_type", str), ("tm_segment", int)):
        Residue.register_attr(session, name, "Trimmings DeepTMHMM", attr_type=kind)
    previous = getattr(session, "_trimmings_tm_records", [])
    records, rows, updated, residue_count = [], [], set(), 0
    with session.undo.aggregate("annotate transmembrane segments"):
        undo = UndoState("DeepTMHMM annotations")
        for entry in entries:
            chain = entry["chain"]
            if not _domain_chain_is_current(chain, entry["sequence"], entry["residues"]):
                session.logger.warning(f"Skipped {entry['spec']}: chain changed during DeepTMHMM prediction.")
                continue
            updated.add(chain)
            topology = predictions[entry["sequence"]]
            assignments = {}
            for segment in _tm_segments(topology):
                members = [r for r in entry["residues"][segment["start"] - 1:segment["end"]] if r is not None]
                for residue in members:
                    assignments[residue] = segment
                modeled = Residues(members)
                records.append({**segment, "structure": chain.structure, "chain": chain,
                                "sequence": entry["sequence"], "sequence_residues": entry["residues"],
                                "residues": modeled})
                spec = concise_residue_spec(session, modeled) if members else ""
                if spec and not spec.startswith("#"):
                    spec = f"#{chain.structure.id_string}" + spec
                link = (f'<a href="cxcmd:{quote("select " + spec, safe="")}">{escape(spec)}</a>'
                        if spec else "No modeled residues")
                swatch = "#ff912d" if segment["label"] == "M" else "#28c8e6"
                rows.append(f'<tr><td>{escape(entry["spec"])}</td><td style="background:{swatch}">'
                            f'{_TM_NAMES[segment["label"]]}</td><td>{segment["number"]}</td>'
                            f'<td>{segment["start"]}&ndash;{segment["end"]}</td><td>{link}</td></tr>')
            groups = {}
            for residue, label in zip(entry["residues"], topology):
                if residue is None:
                    continue
                segment = assignments.get(residue)
                kind = {"M": "helix", "B": "strand"}.get(label, "")
                for name, value in (("tm_topology", label), ("tm_type", kind),
                                    ("tm_segment", segment["number"] if segment else 0)):
                    undo.add(residue, name, getattr(residue, name, None), value)
                    setattr(residue, name, value)
                groups.setdefault(label if label in _TM_NAMES else "", []).append(residue)
                residue_count += bool(segment)
            if color:
                for label, members in groups.items():
                    color_objects(session, Objects(atoms=Residues(members).atoms),
                                  Color(tuple(c / 255 for c in _TM_COLORS[label])), target="ac")
        combined = [record for record in previous if record["chain"] not in updated
                    and _domain_chain_is_current(record["chain"], record["sequence"], record["sequence_residues"])] + records
        state = SimpleNamespace(set_records=lambda value: setattr(session, "_trimmings_tm_records", value))
        undo.add(state, "set_records", previous, combined, option="M")
        state.set_records(combined)
        session.undo.register(undo)
        # tmsegments targets complete chains even when the atom spec touches
        # only one residue. Share that exact mapping with the feature browser.
        feature_entries = [{**entry, "targets": {r for r in entry["residues"] if r is not None}}
                           for entry in entries if entry["chain"] in updated]
        if feature_entries:
            feature_results = {entry["sequence"]: _tm_feature_hits(predictions[entry["sequence"]])
                               for entry in feature_entries}
            _apply_domain_results(session, feature_entries, feature_results, False, True, gui=gui,
                                  source="tmsegments", include_tentative=True)
            if gui and session.ui.is_gui:
                session._trimmings_domain_browser.kind_menu.setCurrentIndex(1)
    helices = sum(record["label"] == "M" for record in records)
    strands = sum(record["label"] == "B" for record in records)
    session.logger.info(f"DeepTMHMM: {helices} TM helices, {strands} TM beta strands; "
                        f"{residue_count} modeled TM residues across {len(updated)} chains. "
                        "Signal peptides (S) are annotated separately from membrane segments.")
    if rows:
        color_note = ("Orange = TM helix; cyan = TM beta strand; other residues are gray. " if color
                      else "Existing model colors were retained. Table colors identify segment types. ")
        session.logger.info('<b>DeepTMHMM predictions</b>. ' + color_note +
                            'Click structure residues to select a segment. Ranges refer to the submitted chain sequence.'
                            '<table border="1" cellpadding="3"><tr><th>Chain</th><th>Type</th><th>Segment</th>'
                            '<th>Sequence positions</th><th>Structure residues</th></tr>' + "".join(rows) + '</table>',
                            is_html=True)
    duration = f" in {elapsed:.0f}s" if elapsed is not None else ""
    session.logger.status(f"DeepTMHMM complete{duration}: {helices} TM helices, {strands} TM beta strands "
                          f"across {len(updated)} chains.")
    return records


def tm_segments_cancel(session):
    cancel = getattr(session, "_trimmings_tm_job", None)
    if cancel is not None:
        cancel.set()
        session._trimmings_tm_job = None
        if getattr(session, "_trimmings_domain_job", None) is cancel:
            session._trimmings_domain_job = None
        session.logger.status("DeepTMHMM cancelled locally; the remote job may continue.")
        session.logger.info("Cancelled DeepTMHMM locally; BioLib may finish a submitted remote job.")
    else:
        session.logger.info("No DeepTMHMM prediction is running.")


def annotate_tm_segments(session, atoms, color=True, results=None, python=None,
                         refresh=False, timeout=1200, wait=False, gui=True, _domain_source=False):
    """Predict or import topology for complete protein chains touched by atoms."""
    from math import isfinite
    from pathlib import Path
    from threading import Event, Thread
    from time import monotonic
    from chimerax.core.errors import UserError

    if not isfinite(timeout) or timeout <= 0:
        raise UserError("timeout must be finite and greater than zero (seconds).")
    if getattr(session, "_trimmings_tm_job", None) is not None:
        raise UserError("A DeepTMHMM prediction is running; use 'tmsegments cancel' before starting another.")
    entries = _domain_chain_inputs(atoms)
    for entry in entries:
        if set(entry["sequence"]) - set("ACDEFGHIKLMNPQRSTVWY"):
            raise UserError(f"{entry['spec']}: DeepTMHMM requires the 20 standard amino-acid letters; "
                            "this chain contains unknown or unsupported letters. Supply a complete, unambiguous sequence.")
        if not entry["chain"].full_sequence_known:
            session.logger.warning(f"{entry['spec']}: full polymer sequence is unavailable; "
                                   "using ChimeraX's inferred chain sequence. Missing sequence can affect topology prediction.")
    sequences = list(dict.fromkeys(entry["sequence"] for entry in entries))
    cache = getattr(session, "_trimmings_tm_cache", {})
    cached = {seq: cache[seq] for seq in sequences if not refresh and seq in cache}
    if results is not None:
        try:
            imported = _tm_parse_topologies(Path(results).expanduser().read_text(encoding="utf-8"))
            if any(seq not in imported for seq in sequences):
                raise ValueError("The result file must contain an exact sequence match for every requested chain.")
        except (OSError, ValueError) as error:
            raise UserError(f"Cannot import DeepTMHMM results: {error}") from None
        # Imported predictions take precedence but do not pollute the remote cache.
        return _apply_tm_results(session, entries, imported, color, gui=gui)
    missing = [seq for seq in sequences if seq not in cached]
    if not missing:
        return _apply_tm_results(session, entries, cached, color, gui=gui)
    try:
        executable = _tm_python_executable(python)
    except ValueError as error:
        raise UserError(str(error)) from None
    cancel = Event()
    started = monotonic()
    session._trimmings_tm_job = cancel
    if _domain_source:
        session._trimmings_domain_job = cancel
    session.logger.info(f"Submitting {len(missing)} distinct chain sequence(s) to {_TM_APP} on BioLib; "
                        "coordinates stay local. Use 'tmsegments cancel' to stop waiting.")

    def update_status(seconds):
        # A queued update must not replace completion/cancellation or a newer job's status.
        if getattr(session, "_trimmings_tm_job", None) is cancel and not cancel.is_set():
            session.logger.status(f"DeepTMHMM running: {len(missing)} distinct sequence(s), "
                                  f"{seconds}s elapsed (queued or processing on BioLib)...", blank_after=0)

    def progress(seconds):
        if wait or not session.ui.is_gui:
            update_status(seconds)
        else:
            session.ui.thread_safe(update_status, seconds)

    update_status(0)

    def finish(predictions, error):
        if getattr(session, "_trimmings_tm_job", None) is not cancel:
            return
        session._trimmings_tm_job = None
        if getattr(session, "_trimmings_domain_job", None) is cancel:
            session._trimmings_domain_job = None
        if cancel.is_set():
            return
        if error is not None:
            session.logger.status("DeepTMHMM timed out." if isinstance(error, TimeoutError)
                                  else "DeepTMHMM prediction failed; see Log for details.")
            if wait or not session.ui.is_gui:
                raise UserError(f"DeepTMHMM prediction failed: {error}")
            session.logger.error(f"DeepTMHMM prediction failed: {error}")
            return
        cache = getattr(session, "_trimmings_tm_cache", {})
        cache.update(predictions)
        session._trimmings_tm_cache = cache
        _apply_tm_results(session, entries, {**cached, **predictions}, color, elapsed=monotonic() - started, gui=gui)

    def predict():
        try:
            predictions = _tm_predict(missing, executable, timeout, cancel, status=progress)
            error = None
        except Exception as err:
            predictions, error = None, err
        if wait or not session.ui.is_gui:
            finish(predictions, error)
        else:
            session.ui.thread_safe(finish, predictions, error)

    if wait or not session.ui.is_gui:
        predict()
    else:
        Thread(target=predict, name="ChimeraX DeepTMHMM", daemon=True).start()


def _ligand_tool_module():
    """Load the optional browser beside this script, including runscript installs."""
    import importlib.util
    from pathlib import Path
    import sys
    from types import SimpleNamespace
    from chimerax.core.errors import UserError

    path = Path(__file__).resolve().with_name('ligand_search.py')
    if not path.is_file():
        raise UserError('Download ligand_search.py beside chimerax_trimmings.py to use ligandsearch.')
    name = '_trimmings_ligand_search'
    module = sys.modules.get(name)
    stamp = path.stat().st_mtime_ns
    if module is None or Path(module.__file__).resolve() != path or getattr(module, '_source_stamp', None) != stamp:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module._source_stamp = stamp
        sys.modules[name] = module
    module.SHAPE = SimpleNamespace(**globals())
    return module


def ligand_search(session, map_model=None, resolution=2.6, top=20,
                  library=None, workers=0, observed=True, codes=None, frequency=True, fast=False,
                  torsions=False, conformers=8, groups='all', group_match='any', sources=None,
                  refine=False, refine_conformers=64, smiles=None):
    return _ligand_tool_module().ligand_search(session, map_model, resolution,
                                              top, library, workers, observed, codes, frequency, fast,
                                              torsions=torsions, conformers=conformers, groups=groups,
                                              group_match=group_match, sources=sources,
                                              refine=refine, refine_conformers=refine_conformers,
                                              smiles=smiles)


def ligand_search_stop(session):
    return _ligand_tool_module().ligand_stop(session)


def _visible_atom_pick(view, x, y):
    from chimerax.atomic import PickedAtom
    from chimerax.graphics import PickedTriangle

    depth = None
    while True:
        # Inspect each hit so an opaque surface can stop the ray.
        pick = view.picked_object(x, y, beyond=depth, max_transparent_layers=0)
        if isinstance(pick, PickedAtom):
            return pick
        if pick is None:
            return None
        # Map and model picks wrap the triangle, which has the actual opacity.
        triangle = getattr(pick, "triangle_pick", getattr(pick, "picked_triangle", pick))
        if not (
            getattr(pick, "pick_through", False)
            or (isinstance(triangle, PickedTriangle) and triangle.is_transparent())
        ):
            return None
        if pick.distance is None or pick.distance >= 1 or (depth is not None and pick.distance <= depth):
            return None
        depth = pick.distance


def _visible_atom_box(view, start, end):
    from chimerax.atomic import Atoms
    from chimerax.atomic.structure import PickedAtoms
    from numpy import ones

    camera = view.camera
    width, height = view.window_size
    projection = camera.projection_matrix(view.near_far_distances(camera, None), None, (width, height))
    scene_to_clip = camera.position.inverse().opengl_matrix() @ projection
    left, right = sorted((start[0], end[0]))
    top, bottom = sorted((start[1], end[1]))
    visible = set()
    for group in view.rectangle_pick(*start, *end):
        if not isinstance(group, PickedAtoms):
            continue
        for structure, atoms in group.atoms.by_structure:
            coords = atoms.coords
            points = ones((len(atoms), 4))
            for position in structure.get_scene_positions(displayed_only=True):
                points[:, :3] = position * coords
                projected = points @ scene_to_clip
                for atom, (x, y, _, w) in zip(atoms, projected):
                    if atom in visible or w <= 0:
                        continue
                    px, py = (1 + x / w) * width / 2, (1 - y / w) * height / 2
                    if left <= px <= right and top <= py <= bottom:
                        pick = _visible_atom_pick(view, px, py)
                        if pick is not None and pick.atom == atom:
                            visible.add(atom)
    return [PickedAtoms(Atoms(list(visible)))] if visible else []


def _install_atom_selection(session):
    from types import MethodType

    from chimerax.mouse_modes import MouseMode
    from chimerax.mouse_modes.std_modes import SelectMouseMode, select_pick

    class VisibleAtomSelect(SelectMouseMode):
        def __init__(self, session):
            super().__init__(session)
            self.name = "select visible atoms"

        def mouse_up(self, event):
            if self.double_click:
                SelectMouseMode.mouse_up(self, event)
                return
            self._undraw_drag_rectangle()
            mode = "toggle" if event.shift_down() else "replace"
            if self._is_drag(event):
                pick = _visible_atom_box(self.view, self.mouse_down_position, event.position())
            else:
                pick = _visible_atom_pick(self.view, *event.position())
            MouseMode.mouse_up(self, event)
            select_pick(self.session, pick, mode)

    manager = session.ui.mouse_modes
    atom_select = manager.named_mode("select visible atoms")
    if atom_select is None:
        atom_select = VisibleAtomSelect(session)
        manager.add_mode(atom_select)
    else:
        # runscript uses fresh globals. Refresh the handler on the registered mode.
        atom_select.mouse_up = MethodType(VisibleAtomSelect.mouse_up, atom_select)
    for modifiers in (["control"], ["control", "shift"]):
        manager.bind_mouse_mode("left", modifiers, atom_select)
    manager.bind_mouse_mode("left", ["control", "alt", "shift"], manager.named_mode("select"))


def register_command(logger):
    from chimerax.atomic import AtomsArg
    from chimerax.core.commands import BoolArg, CmdDesc, EnumOf, FloatArg, IntArg, ModelsArg, OpenFileNameArg, StringArg, register

    _register_shape_align_command(logger)

    from chimerax.map import MapArg
    register('ligandsearch', CmdDesc(optional=[('map_model', MapArg)], keyword=[
        ('resolution', FloatArg), ('top', IntArg),
        ('library', StringArg), ('workers', IntArg), ('observed', BoolArg), ('codes', StringArg), ('frequency', BoolArg),
        ('fast', BoolArg), ('torsions', BoolArg), ('conformers', IntArg), ('groups', StringArg),
        ('group_match', EnumOf(('any', 'all'))), ('sources', EnumOf(('ccd', 'ccd_chebi'))),
        ('refine', BoolArg), ('refine_conformers', IntArg), ('smiles', StringArg),
    ], synopsis='Browse ligand fits from CCD and optional ChEBI metabolites against a density blob'),
             ligand_search, logger=logger)
    register('ligandsearch stop', CmdDesc(synopsis='Stop ligand screening or library update; retain completed top hits'),
             ligand_search_stop, logger=logger)

    register("tmsegments", CmdDesc(required=[("atoms", AtomsArg)], keyword=[
        ("color", BoolArg), ("results", OpenFileNameArg), ("python", OpenFileNameArg),
        ("refresh", BoolArg), ("timeout", FloatArg), ("wait", BoolArg), ("gui", BoolArg),
    ], synopsis="Annotate and color predicted transmembrane helices and beta strands with DeepTMHMM"),
        annotate_tm_segments, logger=logger)
    register("tmsegments cancel", CmdDesc(synopsis="Stop waiting for the current DeepTMHMM prediction"),
             tm_segments_cancel, logger=logger)

    register(
        "domains",
        CmdDesc(required=[("atoms", AtomsArg)], keyword=[
            ("color", BoolArg), ("include_overlaps", BoolArg), ("include_tentative", BoolArg),
            ("refresh", BoolArg), ("timeout", FloatArg), ("wait", BoolArg),
            ("gui", BoolArg),
            ("source", EnumOf(("pfam", "interpro", "foldseek", "features", "tmsegments"))),
            ("library", StringArg), ("regions", StringArg),
            ("pfam_guided", BoolArg), ("deduplicate", BoolArg),
            ("pfam_backend", EnumOf(("auto", "local", "online"))), ("pfam_database", StringArg),
            ("features", BoolArg), ("prosite_profiles", BoolArg), ("prosite_common", BoolArg), ("prosite_database", StringArg),
            ("results", OpenFileNameArg), ("python", OpenFileNameArg),
        ], synopsis="Annotate domains or sequence features using Pfam, InterPro, Foldseek, SEG, PROSITE and DeepTMHMM"),
        annotate_domains, logger=logger,
    )
    register("domains cancel", CmdDesc(synopsis="Cancel the current domain or sequence-feature lookup"),
             domains_cancel, logger=logger)
    register("domainbrowser", CmdDesc(optional=[("models", ModelsArg)],
                                      synopsis="Browse domains for all analyzed models or one specified model"),
             domain_browser, logger=logger)

    register(
        "hydrophobicpatches",
        CmdDesc(
            required=[("atoms", AtomsArg)],
            keyword=[
                ("binder_radius", FloatArg), ("binder_mass", FloatArg), ("context", AtomsArg),
                ("minimum_exposure", FloatArg), ("minimum_residues", IntArg),
                ("contact_distance", FloatArg), ("link_distance", FloatArg),
                ("probe_radius", FloatArg),
                ("samples", IntArg), ("max_patches", IntArg),
                ("approach", BoolArg), ("highlight", BoolArg), ("gui", BoolArg),
                ("residue_names", StringArg),
            ],
            synopsis="Highlight hydrophobic patches accessible to a specified globular binder sphere",
        ),
        hydrophobic_patches,
        logger=logger,
    )
    register(
        "patchbrowser",
        CmdDesc(synopsis="Browse the latest hydrophobic patches with previous/next navigation"),
        patch_browser,
        logger=logger,
    )

    register(
        "exposedhydrophobics",
        CmdDesc(
            required=[("atoms", AtomsArg)],
            keyword=[
                ("minimum_exposure", FloatArg),
                ("minimum_area", FloatArg),
                ("residue_names", StringArg),
                ("probe_radius", FloatArg),
            ],
            synopsis="Select hydrophobic residues with relative sidechain exposure > 25% by default",
        ),
        exposed_hydrophobics,
        logger=logger,
    )

    for name, function, synopsis in (
        ("active", model_motion_active, "Move only the specified models with mouse and trackpad"),
        ("freeze", model_motion_freeze, "Freeze the specified models during interactive movement"),
        ("unfreeze", model_motion_unfreeze, "Add the specified models to the active movement set"),
    ):
        register(
            f"modelmotion {name}",
            CmdDesc(required=[("models", ModelsArg)], synopsis=synopsis),
            function,
            logger=logger,
        )
    for name, function, synopsis in (
        ("all", model_motion_all, "Enable interactive movement for all current models"),
        ("off", model_motion_off, "Restore the mouse and trackpad bindings from before modelmotion"),
        ("status", model_motion_status, "List the models currently active for movement"),
    ):
        register(f"modelmotion {name}", CmdDesc(synopsis=synopsis), function, logger=logger)

    register(
        "nextmodel",
        CmdDesc(synopsis="Show the next atomic model"),
        next_model,
        logger=logger,
    )
    register(
        "prevmodel",
        CmdDesc(synopsis="Show the previous atomic model"),
        previous_model,
        logger=logger,
    )
    register(
        "nextitem",
        CmdDesc(synopsis="Show the next unskipped atomic model or volume map"),
        next_item,
        logger=logger,
    )
    register(
        "previtem",
        CmdDesc(synopsis="Show the previous unskipped atomic model or volume map"),
        previous_item,
        logger=logger,
    )
    register(
        "lowermapthreshold",
        CmdDesc(synopsis="Lower displayed map thresholds by 0.5 sigma without crossing zero"),
        lower_map_threshold,
        logger=logger,
    )
    register(
        "raisemapthreshold",
        CmdDesc(synopsis="Raise displayed map thresholds by 0.5 sigma without crossing zero"),
        raise_map_threshold,
        logger=logger,
    )
    register(
        "reorient",
        CmdDesc(
            synopsis="Recenter on visible inertia axes and cycle views along the principal axes"
        ),
        reorient,
        logger=logger,
    )
    register(
        "nextmap",
        CmdDesc(synopsis="Show the next map and hide the current map"),
        next_map,
        logger=logger,
    )
    register(
        "prevmap",
        CmdDesc(synopsis="Show the previous map and hide the current map"),
        previous_map,
        logger=logger,
    )
    register(
        "togglemodels",
        CmdDesc(synopsis="Hide shown atomic models, then restore the same set on the next use"),
        toggle_models,
        logger=logger,
    )
    register(
        "togglemaps",
        CmdDesc(synopsis="Hide shown volume maps, then restore the same set on the next use"),
        toggle_maps,
        logger=logger,
    )
    register(
        "cyclelighting",
        CmdDesc(synopsis="Cycle through lighting presets"),
        cycle_lighting,
        logger=logger,
    )
    register(
        "cyclemodeldisplay",
        CmdDesc(
            synopsis="Cycle shown atomic models through saved, CA-trace, ribbon, and related display presets"
        ),
        cycle_model_display,
        logger=logger,
    )


if "session" in globals():
    register_command(session.logger)
    if session.ui.is_gui:
        _install_atom_selection(session)

    # Reloading also upgrades an already-enabled controller in the session.
    if getattr(getattr(session, "_trimmings_model_motion", None), "enabled", False):
        _model_motion_controller(session)
        _watch_new_motion_models(session)
