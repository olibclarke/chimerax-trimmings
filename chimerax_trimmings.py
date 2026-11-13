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
_SHAPE_CACHE_BYTES = 256 * 1024 ** 2  # Retained preprocessing only, not fitter scratch space.
_SHAPE_AUTO_PARALLEL_MIN = 4
_MAX_GRID_VOXELS = 16_000_000  # 64 MB per float32 grid, before fitter scratch space.
_MAX_CORRELATION_VOXELS = 8_000_000  # Several moment arrays are needed at once.
_MAX_PARTIAL_FFT_VOXELS = 4_000_000
_PARTIAL_SIZE_RATIO = 0.75  # Atom-count proxy for unequal density mass/volume.
_CHAIN_BACKUP_CORRELATION = 0.4  # Weak partial fits may have unmatched parts on both sides.
_CHAIN_BACKUP_MIN_FRACTION = 0.2  # Ignore tiny chains and cap the extra searches below.
_MAX_CHAIN_BACKUP_STARTS = 4
_AUTO_FINE_LIMITS = (2.0, 8.0)  # Bound unnecessary atom-level detail and preserve helix-scale features.


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


def _shape_confidence_weights(atoms, mode):
    """Only interpret B-factors as pLDDT with provenance or explicit consent."""
    import numpy as np
    from chimerax.core.errors import UserError

    if mode not in ("auto", "off", "plddt"):
        raise UserError("confidence must be auto, off, or plddt.")
    database = getattr(atoms[0].structure, "database", None) if len(atoms) else None
    known = getattr(database, "database", None) == "AlphaFold"
    if mode == "off" or (mode == "auto" and not known):
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


def _shape_reference_fft(grid, field, fft_shape, objective, squared=False):
    def prepare():
        from scipy.fft import rfftn

        return grid, rfftn(field * field if squared else field, fft_shape)

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
                     confidence="auto", confidence_moving=None, confidence_reference=None):
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
                   partial=partial, fallback_correlation=fallback_correlation, objective=objective)
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
    try:
        if move:
            old = model.position
            parent_scene = model.parent.scene_position
            new = parent_scene.inverse() * transform * parent_scene * old
            undo.add(model, "position", old, new)
            model.position = new
        if scores is not None:
            _assign_shape_correlations(session, model, source_atoms, scores, undo, color_fit)
    except BaseException:
        undo.undo()
        raise
    if undo.state:
        session.undo.register(undo)
    if score_fit:
        stats["correlation_grid_seconds"] = stats["correlation_seconds"]
        stats["correlation_seconds"] = perf_counter() - scoring_started
    stats.update(moving_model_id=model.id_string, reference_model_id=reference.id_string)
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
                     confidence="auto", confidence_moving=None, confidence_reference=None):
    prepared = _prepare_shape_atomic_fit(
        session, atoms, to, resolution, max_points, extra_starts, move, radius, correlation,
        color_fit, correlation_resolution, correlation_window, sphere_center, partial,
        fallback_correlation, objective, confidence, confidence_moving, confidence_reference)
    result = _compute_shape_atomic_fit(_shape_atomic_snapshot(prepared), session.logger)
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
                if session.ui.is_gui:
                    from Qt.QtCore import QCoreApplication, QEventLoop

                    # Repaint progress without accepting reentrant user commands
                    # while this synchronous batch owns the pending snapshots.
                    QCoreApplication.processEvents(QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents)
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
                confidence="auto", confidence_moving=None, confidence_reference=None):
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
                   confidence=confidence, confidence_moving=confidence_moving, confidence_reference=confidence_reference)
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
        computed = _shape_atomic_workers(session, prepared, worker_count)
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
        indices = np.arange(start, min(start + 65536, matrix.size))
        k, rem = np.divmod(indices, nx * ny)
        j, i = np.divmod(rem, nx)
        points = np.column_stack((i, j, k)).astype(np.float32)
        points *= step
        points += origin
        matrix.ravel()[start:start + len(indices)] = volume_values(points, transform, grid.matrix())
    return matrix


class _ShapeDensityInput:
    """A selected atomic density or native full map, in the owner's local frame."""

    def __init__(self, model, atoms, center=None, radius=None, confidence="auto"):
        from itertools import product
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
            box = np.array(list(product(*[(v - radius, v + radius) for v in local_center])))
            ijk = data.xyz_to_ijk_transform.transform_points(box)
            lower = np.maximum(np.floor(ijk.min(axis=0)).astype(int), 0)
            upper = np.minimum(np.ceil(ijk.max(axis=0)).astype(int) + 1, data.size)
            if np.any(upper <= lower) or np.prod(upper - lower) > _MAX_GRID_VOXELS:
                raise UserError("The map CoFR region is empty or exceeds the working-grid limit.")
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
            raise UserError("The map density is collinear and cannot define a rigid orientation.")
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
            from scipy.fft import rfftn, irfftn, fftfreq, rfftfreq
            fz = fftfreq(work.shape[0])[:, None, None]
            fy = fftfreq(work.shape[1])[None, :, None]
            fx = rfftfreq(work.shape[2])[None, None, :]
            covariance = sigma * sigma * inverse @ inverse.T
            q = (covariance[0, 0] * fx * fx + covariance[1, 1] * fy * fy + covariance[2, 2] * fz * fz +
                 2 * covariance[0, 1] * fx * fy + 2 * covariance[0, 2] * fx * fz + 2 * covariance[1, 2] * fy * fz)
            spectrum = rfftn(work)
            spectrum *= np.exp(-2 * pi * pi * q)
            work = irfftn(spectrum, s=work.shape).astype(np.float32)
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
        confidence, getattr(getattr(model, "database", None), "database", None), _shape_array_key(atoms.bfactors))
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
    from chimerax.map_fit.fitmap import locate_maximum, volume_values

    points, weights = _shape_density_samples(source, maximum)
    field, second = _shape_reference_field(target, objective)
    if objective == "gradient":
        _, first = _gradient_density_field(source)
    else:
        first = np.einsum("ijk,ijk->", source.matrix(), source.matrix(), dtype=np.float64)
    denominator = sqrt(first * second)
    results = []
    for index, initial in enumerate(candidates):
        delta, _ = locate_maximum(initial.transform_points(points), weights, field, target.xyz_to_ijk_transform,
                                  max_steps=240, ijk_step_size_min=.01, ijk_step_size_max=.5, metric="sum product")
        pose = delta * initial
        values = volume_values(points, target.xyz_to_ijk_transform * pose, field)
        score = float(np.dot(weights.astype(np.float64), values)) / denominator
        results.append((score, pose, index))
    return sorted(results, key=lambda item: item[0], reverse=True), points


def _shape_map_fft_starts(source, target, source_landmarks, target_landmarks, objective, extra_starts):
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
        return [pose.inverse() for pose in _shape_map_fft_starts(
            target, source, target_landmarks, source_landmarks, objective, extra_starts)]
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
    target_fft = _shape_reference_fft(target, field, fft_shape, objective)
    peaks = []
    for initial in _shape_starts(source_landmarks, target_landmarks, 96 + extra_starts):
        rotation = Place(matrix=np.column_stack((initial.matrix[:, :3], np.zeros(3))))
        density = _resample_shape_density(source, origin, (side,) * 3, step, rotation)
        cross = irfftn(target_fft * rfftn(density[::-1, ::-1, ::-1], fft_shape), fft_shape) / sqrt(first * second)
        view = cross[half:half + field.shape[0], half:half + field.shape[1], half:half + field.shape[2]]
        for _ in range(2):
            k, j, i = np.unravel_index(np.argmax(view), view.shape)
            center = target.ijk_to_xyz_transform * np.array((i, j, k), dtype=float)
            peaks.append((float(view[k, j, i]), Place(matrix=np.column_stack((rotation.matrix[:, :3], center)))))
            distance = max(2, int(radius / (3 * step)))
            view[max(0, k-distance):k+distance+1, max(0, j-distance):j+distance+1,
                 max(0, i-distance):i+distance+1] = -np.inf
    peaks.sort(key=lambda item: item[0], reverse=True)
    kept = []
    probes = _shape_grid_corners(source)
    for _, pose in peaks:
        if any(np.sqrt(np.mean(np.sum((pose.transform_points(probes) - old.transform_points(probes)) ** 2, axis=1)))
               < step for old in kept):
            continue
        kept.append(pose)
        if len(kept) == 48:
            break
    return kept


def _fit_shape_density_inputs(source, target, options, center):
    import numpy as np
    from chimerax.geometry import translation

    objective, radius = options["objective"], options["radius"]
    levels, stats = _shape_resolution_schedule(source.frame_landmarks, target.frame_landmarks,
                                              options["resolution"], radius is not None)
    initial = (translation(-target.center) * target.model.scene_position.inverse() *
               source.model.scene_position * translation(source.center))
    starts = [initial] if radius is not None else [*_shape_starts(
        source.frame_landmarks, target.frame_landmarks, options["extra_starts"]), initial]
    search_radius = options.get("radius_search")
    if search_radius is not None:
        # A regional global search must explore translations even when PCA
        # already gives a good score. Only the reference density was cropped.
        source.model.session.logger.status("Shape alignment: global search of the reference CoFR sphere...")
        starts.extend(_shape_map_fft_starts(
            source.grid(levels[0]), target.grid(levels[0]), source.frame_landmarks,
            target.frame_landmarks, objective, options["extra_starts"]))

    def refine(candidates, first=0):
        for blur in levels[first:]:
            source.model.session.logger.status(f"Shape alignment: map density at {blur:.3g} A ({len(candidates)} starts)...")
            results, points = _refine_density_grids(source.grid(blur), target.grid(blur), candidates,
                                                   options["max_points"], objective)
            candidates = []
            probes = points[np.linspace(0, len(points) - 1, min(64, len(points)), dtype=int)]
            for _, pose, _ in results:
                if any(np.sqrt(np.mean(np.sum((pose.transform_points(probes) - old.transform_points(probes)) ** 2, axis=1)))
                       < blur / 3 for old in candidates):
                    continue
                candidates.append(pose)
                if len(candidates) == 48:
                    break
        return results[0][:2]

    score, pose = refine(starts)
    stats.update(objective=objective, partial=False, local=radius is not None, map_search=False,
                 moving_input_type="map" if source.is_map else "atoms",
                 reference_input_type="map" if target.is_map else "atoms", resolutions=levels)
    if search_radius is not None:
        stats.update(map_search=True, map_search_resolution=levels[0], radius_search=search_radius)
    elif radius is None and options["fallback_correlation"] > 0 and score < options["fallback_correlation"]:
        first = 1 if len(levels) > 2 else 0
        blur = levels[first]
        seeds = _shape_map_fft_starts(source.grid(blur), target.grid(blur), source.frame_landmarks,
                                      target.frame_landmarks, objective, options["extra_starts"])
        proposed_score, proposed = refine([pose, *seeds], first)
        stats.update(map_search=True, map_search_resolution=blur, map_initial_score=score)
        if proposed_score > score:
            score, pose = proposed_score, proposed
    stats["score"] = score
    if center is not None:
        stats["center_of_rotation"] = center.copy()
        if radius is not None:
            stats.update(radius=radius, radius_refine=radius)
    return translation(target.center) * pose * translation(-source.center), stats


def _shape_map_correlations(source, target, pose, resolution, window):
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
    _shape_log_confidence(session, source.model, source.weights, "moving")
    _shape_log_confidence(session, target.model, target.weights, "reference")
    pose, stats = _fit_shape_density_inputs(source, target, options, center)
    stats.update(moving_confidence="off" if source.weights is None else "plddt",
                 reference_confidence="off" if target.weights is None else "plddt")
    if options.get("radius_search") is not None:
        stats["radius_search_automatic"] = auto_search
    stats["seconds"] = perf_counter() - started
    scores = None
    if options["correlation"] or options["color_fit"]:
        scoring_started = perf_counter()
        scores, extra = _shape_map_correlations(source, target, pose, options["correlation_resolution"],
                                               options["correlation_window"])
        stats.update(extra, correlation_seconds=perf_counter() - scoring_started)
    model, reference = source.model, target.model
    transform = reference.scene_position * pose * model.scene_position.inverse()
    undo = UndoState("shape align")
    try:
        if options["move"]:
            parent = model.parent.scene_position
            old = model.position
            new = parent.inverse() * transform * parent * old
            undo.add(model, "position", old, new)
            model.position = new
        if scores is not None:
            _assign_shape_correlations(session, model, source.atoms, scores, undo, options["color_fit"])
    except BaseException:
        undo.undo()
        raise
    if undo.state:
        session.undo.register(undo)
    if scores is not None:
        stats["correlation_grid_seconds"] = stats["correlation_seconds"]
        stats["correlation_seconds"] = perf_counter() - scoring_started
    stats.update(moving_model_id=model.id_string, reference_model_id=reference.id_string)
    metric = "gradient correlation" if options["objective"] == "gradient" else "density overlap"
    region = (f"reference density in a {options['radius_search']:g} A CoFR search sphere."
              if options.get("radius_search") is not None else
              f"density in a {radius:g} A CoFR refinement sphere." if radius is not None else
              "full map density, independent of contour.")
    message = (f"Shape alignment #{model.id_string} ({stats['moving_input_type']}) to "
               f"#{reference.id_string} ({stats['reference_input_type']}): whole {metric} "
               f"score {stats['score']:.3f}; {stats['seconds']:.2f} s; {region}" +
               (" (not moved)" if not options["move"] else ""))
    session.logger.info(message)
    session.logger.status(message)
    session.logger.info("Fitting scales: " + " -> ".join(f"{v:.3g}" for v in stats["resolutions"]) +
                        " A; Gaussian blur is additional to the map's existing resolution.")
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
        ("objective density|gradient", "Scoring method; default: density."),
        ("resolution A", "Final fitting scale; default: automatic."),
        ("confidence auto|off|plddt", "Weight known AlphaFold models automatically; off disables; plddt forces B-factor confidence."),
        ("color_fit true", "Color moving atoms by local density agreement; default: off."),
        ("workers 0|N", "Atomic batches: 0 available CPUs minus one, 1 serial, or an explicit count."),
    )
    description = "Rigidly fit models/maps; only the moving inputs move. Supports Undo."
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
                confidence="auto", confidence_moving=None, confidence_reference=None):
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
    for model, _ in sources:
        if target in model.drawing_lineage or model in target.drawing_lineage:
            raise UserError("Moving and reference models must not contain one another.")
        if any(other in model.drawing_lineage[:-1] for other, _ in sources):
            raise UserError("Moving inputs must be independent models, not ancestor/descendant pairs.")
    options = dict(resolution=resolution, max_points=max_points, extra_starts=extra_starts, move=move,
                   radius=radius, correlation=correlation, color_fit=color_fit,
                   correlation_resolution=correlation_resolution, correlation_window=correlation_window,
                   partial=partial, fallback_correlation=fallback_correlation, objective=objective,
                   confidence=confidence, confidence_moving=confidence_moving, confidence_reference=confidence_reference)
    if all(atoms is not None for _, atoms in [*sources, *targets]):
        from chimerax.atomic import concatenate

        return _shape_align_atomic_models(session, concatenate([selected for _, selected in sources]),
                                          targets[0][1], workers=workers, **options)
    if workers > 1:
        raise UserError("Parallel workers currently support atomic model/model batches only; use workers 1 for maps.")
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
    with cache_context, undo_context as batch:
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
    from chimerax.core.commands import BoolArg, CmdDesc, EnumOf, FloatArg, IntArg, NoArg, ObjectsArg, Or, register

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
                 ("max_points", IntArg), ("extra_starts", IntArg), ("move", BoolArg),
                 ("radius_search", SearchRadiusArg(FloatArg, NoArg, name="radius in Angstroms or no value (automatic)")),
                 ("radius_refine", FloatArg), ("radius", FloatArg),
                 ("correlation", BoolArg), ("color_fit", BoolArg),
                 ("correlation_resolution", FloatArg), ("correlation_window", FloatArg),
                 ("partial", BoolArg), ("fallback_correlation", FloatArg),
                 ("workers", IntArg),
                 ("confidence", EnumOf(("auto", "off", "plddt"))),
                 ("confidence_moving", EnumOf(("auto", "off", "plddt"))),
                 ("confidence_reference", EnumOf(("auto", "off", "plddt"))),
                 ("objective", EnumOf(("density", "gradient")))],
        synopsis="Align models/maps by shape; radius_search searches a map region, radius_refine refines locally",
    ), shape_align, logger=logger)
    register("shapealign help", CmdDesc(synopsis="Show shapealign usage, options, and examples"),
             _shape_align_help, logger=logger)



_PFAM_API = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
_PFAM_POLL_SECONDS = 5
_DOMAIN_BROWSER_VERSION = 5


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
                if not domain["significant"]:
                    continue
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
                results.append({
                    "type": metadata["type"], "start": start, "end": end,
                    "evalue": evalue, "bitscore": score,
                    "accession": accession, "name": metadata["identifier"],
                    "description": metadata.get("description", ""),
                    "outcompeted": bool(domain["outcompeted"]),
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

    accepted_types = {"domain", "repeat", "homologous_superfamily", "family"}
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
            entry_type = entry.get("type") or signature.get("type") or ""
            if entry_type.lower() not in accepted_types:
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
                results.append({
                    "origin": "interpro", "source": provider, "database": database,
                    "accession": accession, "name": name, "type": entry_type,
                    "description": signature.get("description") or entry.get("description") or "",
                    "start": start, "end": end, "segments": segments,
                    "evalue": float(evalue) if evalue is not None else None,
                    "score": float(score) if score is not None else None, "bitscore": None,
                    "integrated": integrated, "color_accession": integrated or accession,
                    "entry_url": f"https://www.ebi.ac.uk/interpro/entry/{database}/{quote(accession, safe='')}/",
                    "outcompeted": False,
                })
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
    if hit.get("origin") == "interpro":
        # Scores from different member databases cannot be compared numerically.
        # Prefer structural superfamily coverage, then domain/repeat matches, then families.
        rank = {"ssf": 0, "cathgene3d": 1}.get(hit["database"], 2)
        if hit["type"].lower() == "family":
            rank = 3
        return (False, rank, -sum(b - a + 1 for a, b in hit["segments"]), hit["accession"], hit["start"], hit["end"])
    return (hit["outcompeted"], 0, -hit["bitscore"], hit["accession"], hit["start"], hit["end"])


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
    if tool is not None and not tool._closed and getattr(tool, "_domain_browser_version", 0) != _DOMAIN_BROWSER_VERSION:
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


def _live_domain_records(session):
    return [record for record in getattr(session, "_trimmings_domains", [])
            if len(record["residues"]) and _domain_chain_is_current(
                record["chain"], record["sequence"], record["sequence_residues"])]


def _apply_domain_results(session, entries, results, color, include_overlaps, gui=True, source="pfam"):
    from html import escape
    from types import SimpleNamespace
    from urllib.parse import quote
    from chimerax.atomic import Residue, Residues, concise_residue_spec
    from chimerax.core.colors import Color
    from chimerax.core.objects import Objects
    from chimerax.core.undo import UndoState
    from chimerax.std_commands.color import color as color_objects

    # One primary assignment controls display; retain all overlapping hits as text.
    attributes = ("domain_name", "domain_accession", "domain_names", "domain_accessions")
    for name in attributes:
        Residue.register_attr(session, name, "Trimmings domains", attr_type=str)
    rows, records, count, updated_chains = [], [], 0, set()
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
            hits = [h for h in results[entry["sequence"]]
                    if include_overlaps or not h["outcompeted"]]
            hits.sort(key=_domain_hit_sort_key)
            membership = {r: [] for r in entry["targets"]}
            for hit in hits:
                members = list(dict.fromkeys(
                    r for a, b in hit.get("segments", [(hit["start"], hit["end"])])
                    for r in entry["residues"][a - 1:b] if r in membership))
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
                        f'<td><a href="cxcmd:{quote("select " + spec, safe="")}">{escape(hit["name"])}</a></td>'
                        f'<td><a href="{_domain_entry_url(hit)}">'
                        f'{hit["accession"]}</a></td><td>{escape(hit["type"])}</td>'
                        f'<td>{hit["start"]}&ndash;{hit["end"]}</td><td>{escape(spec)}</td>'
                        f'<td>{_domain_number(hit["evalue"])}</td><td>{escape(hit.get("source", "Pfam"))}</td></tr>'
                    )
            groups = {}
            for residue, matches in membership.items():
                primary = matches[0] if matches else None
                values = (
                    primary["name"] if primary else "", primary["accession"] if primary else "",
                    " | ".join(dict.fromkeys(h["name"] for h in matches)),
                    " | ".join(dict.fromkeys(h["accession"] for h in matches)),
                )
                for name, value in zip(attributes, values):
                    # None is ChimeraX's unset attribute value and is undoable.
                    undo.add(residue, name, getattr(residue, name, None), value)
                    setattr(residue, name, value)
                groups.setdefault(primary.get("color_accession", values[1]) if primary else "", []).append(residue)
                count += bool(primary)
            if color:
                for accession, members in groups.items():
                    color_objects(session, Objects(atoms=Residues(members).atoms), Color(_domain_color(accession)), target="ac")
        records.sort(key=lambda r: (r["structure"].id, r["chain"].chain_id, r["start"], r["end"], r["accession"]))
        combined = [record for record in _live_domain_records(session) if record["chain"] not in updated_chains] + records
        combined.sort(key=lambda r: (r["structure"].id, r["chain"].chain_id, r["start"], r["end"], r["accession"]))
        # Keep the browser's results in step with annotation Undo/Redo.
        state = SimpleNamespace(set_records=lambda value: _set_domain_records(session, *value))
        undo.add(state, "set_records", (previous_records, previous_models), (combined, snapshots), option="M")
        state.set_records((combined, snapshots))
        session.undo.register(undo)
    if rows:
        session.logger.info(
            f'<b>{"InterPro" if source == "interpro" else "Pfam"} domain annotations</b> (click a name to select its footprint)'
            '<table border="1" cellpadding="3"><tr><th>Chain</th><th>Domain</th><th>Accession</th>'
            '<th>Hit type</th><th>Query positions</th><th>Structure residues</th><th>E-value</th><th>Source</th></tr>'
            + "".join(rows) + "</table>", is_html=True,
        )
    session.logger.info(f"{'InterPro' if source == 'interpro' else 'Pfam'}: {len(records)} annotation matches; "
                        f"annotated {count} target residues. Unassigned residues have empty domain_name/domain_accession.")
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
    return re.sub(r"\[+cite:[^\]]+\]+", "", text).strip()


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
    info["descriptions"] = [{**d, "source": accession, "database": database} for d in _domain_descriptions(metadata)]
    name = metadata.get("name") or {}
    info["name"] = name.get("name", "") if isinstance(name, dict) else str(name)
    integrated = metadata.get("integrated") or integrated_hint
    if isinstance(integrated, str) and re.fullmatch(r"IPR\d+", integrated):
        info["integrated"] = integrated
        entry = fetch(f"/entry/interpro/{integrated}/", "InterPro description")
        ipr = entry.get("metadata", {}) if isinstance(entry, dict) else {}
        info["descriptions"].extend({**d, "source": integrated} for d in _domain_descriptions(ipr))
        info["go_terms"] = [term for term in ipr.get("go_terms") or []
                            if isinstance(term, dict) and re.fullmatch(r"GO:\d+", term.get("identifier", ""))]

    representative = metadata.get("representative_structure") or {}
    pdb = representative.get("accession", "") if isinstance(representative, dict) else ""
    if re.fullmatch(r"[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+", pdb):
        structure = {"accession": pdb, "name": representative.get("name", ""), "representative": True}
        details = fetch(f"/structure/pdb/{pdb}/entry/{database}/{encoded_accession}/", "Representative structure details")
        if isinstance(details, dict):
            data = details.get("metadata", {})
            structure.update({"resolution": data.get("resolution"), "method": data.get("experiment_type"),
                              "chains": sorted({e["chain"] for e in details.get("entries", []) if e.get("chain")})})
        info["structures"].append(structure)
    else:
        examples = fetch(f"/structure/pdb/entry/{database}/{encoded_accession}/?page_size=3", "Example structures")
        if isinstance(examples, dict):
            for item in examples.get("results", [])[:3]:
                data = item.get("metadata", {})
                pdb = data.get("accession", "")
                if re.fullmatch(r"[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+", pdb):
                    info["structures"].append({
                        "accession": pdb, "name": data.get("name", ""), "representative": False,
                        "resolution": data.get("resolution"), "method": data.get("experiment_type"),
                        "chains": sorted({e["chain"] for e in item.get("entries", []) if e.get("chain")}),
                    })
    return info


def _domain_details_html(session, record, info=None, loading=False):
    from html import escape
    from chimerax.atomic import concise_residue_spec

    accession = record["accession"]
    database = record.get("database", "pfam")
    score_text = (f'{escape(record["source"])} score {_domain_number(record.get("score"))}'
                  if record.get("origin") == "interpro" else f'bit score {_domain_number(record.get("bitscore"))}')
    parts = [f'<h3>{escape(record["name"])} ({escape(accession)})</h3>',
             f'<p>{escape(concise_residue_spec(session, record["residues"]))}<br>'
             f'Query positions {record["start"]}&ndash;{record["end"]}; {len(record["residues"])} modeled residues<br>'
             f'{escape(record["type"])}; E-value {_domain_number(record.get("evalue"))}; {score_text}</p>',
             f'<p><a href="{_domain_entry_url(record)}">{escape(record.get("source", "Pfam"))} entry</a> | '
             f'<a href="https://www.ebi.ac.uk/interpro/structure/pdb/entry/{database}/{accession}/">All matching structures</a></p>']
    if len(record.get("segments", [])) > 1:
        parts.append('<p>Segments: ' + ', '.join(f'{a}–{b}' for a, b in record["segments"]) + '</p>')
    if record.get("description"):
        parts.append(f'<p>{escape(_domain_plain_text(record["description"]))}</p>')
    if info is None:
        parts.append("<p>Loading descriptions and structure links&hellip;</p>" if loading
                     else "<p>Details have not been fetched yet.</p>")
        return "".join(parts)
    if info.get("integrated"):
        ipr = info["integrated"]
        parts.append(f'<p>Integrated InterPro entry: '
                     f'<a href="https://www.ebi.ac.uk/interpro/entry/interpro/{ipr}/">{escape(ipr)}</a></p>')
    for description in info.get("descriptions", []):
        source = description["source"]
        description_database = description.get("database", "pfam" if source.startswith("PF") else "interpro")
        label = " (AI-generated database text)" if description.get("generated") else ""
        text = description["text"]
        if len(text) > 1400:
            text = text[:1400].rsplit(" ", 1)[0] + "…"
        parts.append(f'<p><b>{escape(source)} description{label}</b><br>{escape(text)} '
                     f'<a href="https://www.ebi.ac.uk/interpro/entry/{description_database}/{source}/">Full entry</a></p>')
    terms = info.get("go_terms", [])
    if terms:
        parts.append("<p><b>Associated GO annotations</b><br>" + "<br>".join(
            f'<a href="https://www.ebi.ac.uk/QuickGO/term/{term["identifier"]}">{term["identifier"]}</a>: '
            f'{escape(term.get("name", ""))}' for term in terms) + "</p>")
    structures = info.get("structures", [])
    if structures:
        parts.append("<p><b>Structures</b></p><ul>")
        for structure in structures:
            pdb = structure["accession"]
            label = "Database representative" if structure["representative"] else "Example match"
            details = [str(structure["method"])] if structure.get("method") else []
            if structure.get("resolution") is not None:
                details.append(f'{structure["resolution"]} Å')
            if structure.get("chains"):
                details.append("chains " + ", ".join(structure["chains"]))
            name = structure.get("name", "")
            if isinstance(name, dict):
                name = name.get("name", "")
            parts.append(f'<li>{label}: <a href="https://www.rcsb.org/structure/{pdb}">{pdb.upper()}</a> '
                         f'(<a href="https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb}">PDBe</a>)<br>'
                         f'{escape(name)}<br>{escape("; ".join(details))}</li>')
        parts.append("</ul>")
    elif not info.get("errors"):
        parts.append("<p>No experimental structure links are available for this entry.</p>")
    if not info.get("descriptions"):
        parts.append("<p>No additional functional description is available.</p>")
    if info.get("errors"):
        parts.append("<p>Some details could not be retrieved: " + escape("; ".join(info["errors"])) +
                     ". Use Reload details to retry.</p>")
    parts.append("<p><i>Annotations describe the domain family; exact function can differ between members.</i></p>")
    return "".join(parts)


def _ensure_domain_info(session, accession, database="pfam", integrated=None):
    from threading import Thread

    cache = getattr(session, "_trimmings_domain_info", {})
    if accession in cache or getattr(session, "_trimmings_domain_info_pending", None) is not None:
        return
    token = object()
    session._trimmings_domain_info_pending = token

    def finish(info):
        if getattr(session, "_trimmings_domain_info_pending", None) is not token:
            return
        session._trimmings_domain_info_pending = None
        cache = getattr(session, "_trimmings_domain_info", {})
        cache[accession] = info
        session._trimmings_domain_info = cache
        tool = getattr(session, "_trimmings_domain_browser", None)
        if tool is not None and not tool._closed:
            current = tool.current()
            if current is not None and (current["accession"] == accession or current["accession"] not in cache):
                tool.show_current()

    def lookup():
        try:
            info = _fetch_domain_info(accession, database, integrated)
        except Exception as error:
            info = {"errors": [str(error)]}
        session.ui.thread_safe(finish, info)

    Thread(target=lookup, name="ChimeraX domain details", daemon=True).start()


def _create_domain_browser(session):
    from Qt.QtGui import QColor, QIcon, QPixmap
    from Qt.QtWidgets import QCheckBox, QComboBox, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton, QTextBrowser, QVBoxLayout
    from chimerax.core.models import REMOVE_MODELS
    from chimerax.core.tools import ToolInstance
    from chimerax.std_commands.view import NamedView
    from chimerax.ui import MainToolWindow

    class DomainBrowser(ToolInstance):
        SESSION_SAVE = False
        SESSION_ENDURING = False
        _domain_browser_version = _DOMAIN_BROWSER_VERSION

        def __init__(self):
            super().__init__(session, "Protein Domains")
            self._closed, self._rows = False, []
            self._model_filter, self._remembered_domains = None, {}
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
                button.setToolTip(label)
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
            self.summary = QLabel()
            root.addWidget(self.summary)
            self.domain_list = QListWidget()
            self.domain_list.setAlternatingRowColors(True)
            self.domain_list.setMinimumWidth(360)
            self.domain_list.currentRowChanged.connect(self._activate_current)
            self.domain_list.itemClicked.connect(self._activate_current)
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
            self.details.setOpenExternalLinks(True)
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
            self._removed_handler = session.triggers.add_handler(REMOVE_MODELS, self._models_removed)
            self.tool_window.manage("side")
            self.refresh()

        def current(self):
            row = self.domain_list.currentRow()
            return self._rows[row] if 0 <= row < len(self._rows) else None

        def browser_model(self):
            record = self.current()
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
            return (record["chain"], record["accession"], record["start"], record["end"])

        def _remember_current(self):
            current = self.current()
            if current is not None:
                self._remembered_domains[self._model_filter] = self._record_key(current)

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

        def refresh(self, remember=True):
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
            self._rows = [record for record in _live_domain_records(session)
                          if self._model_filter is None or record["structure"] is self._model_filter]
            self.domain_list.blockSignals(True)
            self.domain_list.clear()
            selected = 0 if self._rows else -1
            remembered = self._remembered_domains.get(self._model_filter)
            for index, record in enumerate(self._rows):
                if self._record_key(record) == remembered:
                    selected = index
                chain_spec = record["chain"].atomspec
                if not chain_spec.startswith("#"):
                    chain_spec = f'#{record["structure"].id_string}' + chain_spec
                text = (f'{chain_spec} | {record["name"]} ({record["accession"]}) | '
                        f'{record["start"]}–{record["end"]}')
                if record.get("origin") == "interpro":
                    text += f' | {record["source"]}'
                swatch = QPixmap(12, 12)
                swatch.fill(QColor.fromRgbF(*_domain_color(record.get("color_accession", record["accession"]))))
                item = QListWidgetItem(QIcon(swatch), text)
                self.domain_list.addItem(item)
            self.domain_list.setCurrentRow(selected)
            self.domain_list.blockSignals(False)
            scope = (f"{len(models)} analyzed model(s)" if self._model_filter is None
                     else f"#{self._model_filter.id_string} {self._model_filter.name}")
            self.summary.setText(f"{len(self._rows)} annotation matches in {scope}")
            for button in (self.previous_button, self.next_button, self.focus_button, self.reload_button):
                button.setEnabled(bool(self._rows))
            self.status.setText("Click a domain to select it. Details load from InterPro/Pfam.")
            self.show_current()

        def step(self, direction):
            if self._rows:
                index = (self.domain_list.currentRow() + direction) % len(self._rows)
                if index == self.domain_list.currentRow():
                    self._activate_current()
                else:
                    self.domain_list.setCurrentRow(index)
                    self.domain_list.scrollToItem(self.domain_list.currentItem())

        def show_current(self):
            self.show_only_model_button.setEnabled(self.browser_model() is not None)
            record = self.current()
            if record is not None and not any(record is r for r in _live_domain_records(session)):
                self.refresh()
                return
            if record is None:
                self.details.setHtml("<p>No domain matches were found in the analyzed chains of this model.</p>"
                                     if self._model_filter is not None else
                                     "<p>No live domain matches. Run <b>domains #1</b> to annotate a structure.</p>")
                return
            cache = getattr(session, "_trimmings_domain_info", {})
            info = cache.get(record["accession"])
            self.details.setHtml(_domain_details_html(session, record, info, loading=info is None))
            self.details.verticalScrollBar().setValue(0)
            if info is None:
                _ensure_domain_info(session, record["accession"], record.get("database", "pfam"), record.get("integrated"))

        def _activate_current(self, *_args):
            from chimerax.atomic import selected_atoms
            from chimerax.core.objects import Objects
            from chimerax.std_commands.select import select

            record = self.current()
            if record is None:
                return
            if not any(record is r for r in _live_domain_records(session)):
                self.refresh()
                return
            atoms = record["residues"].atoms
            if set(selected_atoms(session)) != set(atoms):
                select(session, Objects(atoms=atoms, bonds=atoms.intra_bonds))
            self.show_current()
            if self.auto_center.isChecked():
                self.focus_current()

        def focus_current(self):
            record = self.current()
            if record is None:
                return
            if not any(record is r for r in _live_domain_records(session)):
                self.refresh()
                return
            model = record["structure"]
            if not model.visible:
                self.status.setText(f"Model #{model.id_string} is hidden; show it to inspect this domain.")
                return
            bounds = record["residues"].atoms.scene_bounds
            if bounds is None:
                return
            view, center = session.main_view, bounds.center()
            view.view_all(bounds, pad=0.15)
            normal, depth = view.camera.view_direction(), max(bounds.radius(), 1.0)
            view.clip_planes.set_clip_position("near", center - depth * normal, view)
            view.clip_planes.set_clip_position("far", center + depth * normal, view)
            view.center_of_rotation = center
            self.status.setText("Centered on domain. Restore view returns to the starting camera and clipping.")

        def restore_view(self):
            view = session.main_view
            self._starting_view.set_view(view, [])
            view.center_of_rotation = self._starting_cofr.copy()
            view.center_of_rotation_method = self._starting_cofr_method

        def reload_details(self):
            record = self.current()
            if record is not None:
                getattr(session, "_trimmings_domain_info", {}).pop(record["accession"], None)
                self.show_current()

        def _models_removed(self, *_args):
            self.refresh()

        def delete(self):
            if self._closed:
                return
            self._closed = True
            self._removed_handler.remove()
            if getattr(session, "_trimmings_domain_browser", None) is self:
                session._trimmings_domain_browser = None
            super().delete()

    return DomainBrowser()


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
    if tool is not None and not tool._closed and getattr(tool, "_domain_browser_version", 0) != _DOMAIN_BROWSER_VERSION:
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
        session.logger.info("Cancelled the Pfam lookup locally; EMBL-EBI may finish its submitted search.")
    else:
        session.logger.info("No Pfam lookup is running.")


def annotate_domains(session, atoms, color=True, include_overlaps=False,
                     refresh=False, timeout=600, wait=False, gui=True, source="pfam"):
    """Look up InterPro matches or search Pfam, and annotate only requested residues."""
    from math import isfinite
    from threading import Event, Thread
    from chimerax.core.errors import UserError

    if not isfinite(timeout) or timeout <= 0:
        raise UserError("timeout must be finite and greater than zero (seconds).")
    if source not in ("pfam", "interpro"):
        raise UserError("source must be pfam or interpro.")
    if getattr(session, "_trimmings_domain_job", None) is not None:
        raise UserError("A domain lookup is already running; use 'domains cancel' before starting another.")
    entries = _domain_chain_inputs(atoms)
    sequences = list(dict.fromkeys(entry["sequence"] for entry in entries))
    cache_name = "_trimmings_interpro_cache" if source == "interpro" else "_trimmings_pfam_cache"
    cache = getattr(session, cache_name, {})
    cached = {seq: cache[seq] for seq in sequences if not refresh and seq in cache}
    missing = [seq for seq in sequences if seq not in cached]
    if not missing:
        _apply_domain_results(session, entries, cached, color, include_overlaps, gui=gui, source=source)
        return
    cancel = Event()
    session._trimmings_domain_job = cancel
    description = ("Looking up exact sequence hashes in InterPro Matches" if source == "interpro"
                   else "Searching Pfam at EMBL-EBI HMMER (curated gathering thresholds)")
    session.logger.info(f"{description}: {len(missing)} distinct sequence(s). Use 'domains cancel' to stop waiting.")
    for entry in entries:
        if not entry["chain"].full_sequence_known:
            session.logger.warning(
                f"{entry['spec']}: full polymer sequence is unavailable; using ChimeraX's inferred chain sequence. "
                "Missing regions may prevent an exact InterPro lookup or reduce Pfam detection."
            )

    def finish(results, error):
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
        cache = getattr(session, cache_name, {})
        cache.update(results)
        setattr(session, cache_name, cache)
        _apply_domain_results(session, entries, {**cached, **results}, color, include_overlaps, gui=gui, source=source)

    def lookup():
        try:
            search = _interpro_search if source == "interpro" else _pfam_search
            results = search(missing, timeout, cancel,
                             lambda msg: session.ui.thread_safe(session.logger.status, msg))
            error = None
        except Exception as err:
            results, error = None, err
        if wait or not session.ui.is_gui:
            finish(results, error)
        else:
            session.ui.thread_safe(finish, results, error)

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


def _apply_tm_results(session, entries, predictions, color=True, elapsed=None):
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
        session.logger.status("DeepTMHMM cancelled locally; the remote job may continue.")
        session.logger.info("Cancelled DeepTMHMM locally; BioLib may finish a submitted remote job.")
    else:
        session.logger.info("No DeepTMHMM prediction is running.")


def annotate_tm_segments(session, atoms, color=True, results=None, python=None,
                         refresh=False, timeout=1200, wait=False):
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
        return _apply_tm_results(session, entries, imported, color)
    missing = [seq for seq in sequences if seq not in cached]
    if not missing:
        return _apply_tm_results(session, entries, cached, color)
    try:
        executable = _tm_python_executable(python)
    except ValueError as error:
        raise UserError(str(error)) from None
    cancel = Event()
    started = monotonic()
    session._trimmings_tm_job = cancel
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
        _apply_tm_results(session, entries, {**cached, **predictions}, color, elapsed=monotonic() - started)

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


def register_command(logger):
    from chimerax.atomic import AtomsArg
    from chimerax.core.commands import BoolArg, CmdDesc, EnumOf, FloatArg, IntArg, ModelsArg, OpenFileNameArg, StringArg, register

    _register_shape_align_command(logger)

    register("tmsegments", CmdDesc(required=[("atoms", AtomsArg)], keyword=[
        ("color", BoolArg), ("results", OpenFileNameArg), ("python", OpenFileNameArg),
        ("refresh", BoolArg), ("timeout", FloatArg), ("wait", BoolArg),
    ], synopsis="Annotate and color predicted transmembrane helices and beta strands with DeepTMHMM"),
        annotate_tm_segments, logger=logger)
    register("tmsegments cancel", CmdDesc(synopsis="Stop waiting for the current DeepTMHMM prediction"),
             tm_segments_cancel, logger=logger)

    register(
        "domains",
        CmdDesc(required=[("atoms", AtomsArg)], keyword=[
            ("color", BoolArg), ("include_overlaps", BoolArg),
            ("refresh", BoolArg), ("timeout", FloatArg), ("wait", BoolArg),
            ("gui", BoolArg),
            ("source", EnumOf(("pfam", "interpro"))),
        ], synopsis="Color and annotate domains using InterPro matches or a Pfam search"),
        annotate_domains, logger=logger,
    )
    register("domains cancel", CmdDesc(synopsis="Cancel the current Pfam domain lookup"),
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

    # Reloading also upgrades an already-enabled controller in the session.
    if getattr(getattr(session, "_trimmings_model_motion", None), "enabled", False):
        _model_motion_controller(session)
        _watch_new_motion_models(session)
