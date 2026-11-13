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
from math import pi, sqrt

_MAX_GRID_VOXELS = 16_000_000  # 64 MB per float32 grid, before fitter scratch space.
_MAX_CORRELATION_VOXELS = 8_000_000  # Several moment arrays are needed at once.
_MAX_PARTIAL_FFT_VOXELS = 4_000_000
_PARTIAL_SIZE_RATIO = 0.75  # Atom-count proxy for unequal density mass/volume.
_CHAIN_BACKUP_CORRELATION = 0.4  # Weak partial fits may have unmatched parts on both sides.
_CHAIN_BACKUP_MIN_FRACTION = 0.2  # Ignore tiny chains and cap the extra searches below.
_MAX_CHAIN_BACKUP_STARTS = 4
_AUTO_FINE_LIMITS = (2.0, 8.0)  # Bound unnecessary atom-level detail and preserve helix-scale features.


def _shape_spatial_size(xyz):
    """Equivalent diameter of a uniform PCA ellipsoid, independent of atom count."""
    import numpy as np

    centered = xyz - xyz.mean(axis=0)
    variances = np.linalg.eigvalsh(centered.T @ centered / len(centered))
    # A uniform solid ellipsoid has variance a^2/5 along semiaxis a. Its
    # equivalent-volume diameter uses all axes, not just the longest extent.
    diameters = 2 * np.sqrt(5 * np.maximum(variances, 0))
    # A planar atom cloud still has finite physical thickness at atomic scales.
    return float(np.exp(np.log(np.maximum(diameters, 2)).mean()))


def _shape_resolution_schedule(moving, reference, resolution=None, local=False):
    """Derive both endpoints and geometric intermediate levels from fitted size."""
    import numpy as np

    moving_size = _shape_spatial_size(moving)
    reference_size = _shape_spatial_size(reference)
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


def shape_correlations(moving, reference, resolution=8.0, window=12.0):
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
    for xyz, grid in ((moving, first_grid), (reference, second_grid)):
        points, weights = _shape_points(xyz, resolution / 6)
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


def _shape_points(xyz, spacing):
    """Collapse occupied cells to centroids, retaining the number of atoms."""
    import numpy as np

    cells = np.floor(xyz / spacing).astype(np.int64)
    cells -= cells.min(axis=0)
    keys = np.ravel_multi_index(cells.T, tuple(cells.max(axis=0) + 1))
    _, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
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


def _shape_starts(moving, reference, extra_starts=0):
    """Try every proper signed permutation of the principal-axis frames."""
    from itertools import permutations, product
    import numpy as np
    from chimerax.geometry import Place

    def frame(xyz):
        center = xyz.mean(axis=0, dtype=np.float64)
        centered = xyz - center
        _, axes = np.linalg.eigh(centered.T @ centered)
        return center, axes

    cm, am = frame(moving)
    cr, ar = frame(reference)
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
                          objective="density", max_starts=24):
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
    mx, mw = _shape_points(moving, step / 2)
    rx, rw = _shape_points(reference, step / 2)
    target_grid = _shape_grid(rx, rw, resolution)
    target = target_grid.matrix()
    gradient = objective == "gradient"
    if gradient:
        target, target_energy = _gradient_density_field(target_grid)
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
    target_fft = rfftn(target, fft_shape)
    target_squared_fft = None if gradient else rfftn(target * target, fft_shape)
    rotations = [pose.matrix[:, :3] for pose in _shape_starts(moving, reference)]
    def proper_axes(xyz):
        _, axes = np.linalg.eigh(xyz.T @ xyz)
        if np.linalg.det(axes) < 0:
            axes[:, 0] *= -1
        return axes

    moving_axes, reference_axes = proper_axes(moving), proper_axes(reference)
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
        anchor = reference.mean(axis=0)
        nearest = np.full(len(rx), np.inf)
        for _ in range(min(8, len(rx))):
            anchors.append(anchor)
            squared = np.sum((reference - anchor) ** 2, axis=1)
            indices = np.argpartition(squared, len(moving) - 1)[:len(moving)]
            rotations.extend(pose.matrix[:, :3] for pose in
                             _shape_starts(moving, reference[indices]))
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
                        partial=False, coarse=False, objective="density"):
    """Refine and rank poses against the same complete input densities."""
    import numpy as np
    from chimerax.map_fit.fitmap import locate_maximum, volume_values, overlap_and_correlation

    # Cell averaging preserves atom counts and hence relative density mass.
    mx, mw = _shape_points(moving, blur / 6)
    rx, rw = _shape_points(reference, blur / 6)
    mg, rg = _shape_grid(mx, mw, blur), _shape_grid(rx, rw, blur)
    points, weights, sampling, mnorm = _shape_samples(mg, max_points)
    rmatrix = rg.matrix()
    if objective == "gradient":
        _, moving_energy = _gradient_density_field(mg)
        rmatrix, reference_energy = _gradient_density_field(rg)
        denominator = sqrt(moving_energy * reference_energy)
    else:
        denominator = sqrt(mnorm * float(np.einsum("ijk,ijk->", rmatrix, rmatrix, dtype=np.float64)))
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
                          status, initial_transform, local, fallback_correlation):
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
    cm, cr = moving.mean(axis=0), reference.mean(axis=0)
    moving, reference = moving - cm, reference - cr
    levels, schedule_stats = _shape_resolution_schedule(moving, reference, resolution, local)
    if local:
        initial = Place() if initial_transform is None else initial_transform
        starts = [translation(-cr) * initial * translation(cm)]
    else:
        starts = _shape_starts(moving, reference, extra_starts)
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
        first_level = 1 if len(levels) > 2 else 0
        search_resolution = levels[first_level]
        seeds, search_stats = _partial_shape_starts(
            source, target, search_resolution, extra_starts, status,
            objective="gradient", max_starts=48,
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
                     partial=None, fallback_correlation=None, objective="density"):
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
    """
    from time import perf_counter
    import numpy as np
    from chimerax.geometry import Place, translation
    from chimerax.core.errors import UserError

    if objective not in ("density", "gradient"):
        raise UserError("objective must be density or gradient.")
    if fallback_correlation is None:
        fallback_correlation = 0.5 if objective == "gradient" else 0.7
    if objective == "gradient":
        if partial is True:
            raise UserError("objective gradient uses whole structures; do not combine it with partial true.")
        return _fit_gradient_points(moving, reference, resolution, max_points, extra_starts,
                                    status, initial_transform, local, fallback_correlation)

    started = perf_counter()
    moving = np.asarray(moving, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    use_partial = not local and (partial is True or (partial is None and
        min(len(moving), len(reference)) / max(len(moving), len(reference)) < _PARTIAL_SIZE_RATIO))
    if use_partial and len(moving) > len(reference):
        # Always fit the smaller template into the larger target. Invert the
        # result, not the user-visible choice of which model is allowed to move.
        reverse, stats = fit_shape_points(
            reference, moving, resolution, max_points, extra_starts, status,
            initial_transform.inverse() if initial_transform is not None else None,
            partial=True,
        )
        stats.update(moving_atoms=len(moving), reference_atoms=len(reference),
                     partial_template_reversed=True, seconds=perf_counter() - started,
                     partial_trigger="size" if partial is None else "requested")
        stats["moving_size_diameter"], stats["reference_size_diameter"] = (
            stats["reference_size_diameter"], stats["moving_size_diameter"])
        return reverse.inverse(), stats
    cm, cr = moving.mean(axis=0), reference.mean(axis=0)
    moving, reference = moving - cm, reference - cr
    levels, schedule_stats = _shape_resolution_schedule(moving, reference, resolution, local)
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
        starts, partial_stats = _partial_shape_starts(moving, reference, levels[0], extra_starts, status)
        if initial_transform is not None:
            starts.append(translation(-cr) * initial_transform * translation(cm))
    else:
        starts = _shape_starts(moving, reference, extra_starts)
        if initial_transform is not None:
            starts.append(translation(-cr) * initial_transform * translation(cm))
    candidates = starts
    counts = []
    for level, blur in enumerate(levels):
        if status:
            mode = "Local shape alignment" if local else "Partial shape alignment" if use_partial else "Shape alignment"
            status(f"{mode}: {blur:g} A density ({len(candidates)} starts)...")
        results, points = _refine_shape_level(moving, reference, blur, candidates,
                                              max_points, partial=use_partial, coarse=level == 0)
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
                        transform, stats, resolution, max_points, extra_starts, logger):
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
        try:
            pose, _ = fit_shape_points(source[indices], target, resolution, max_points, extra_starts,
                                       status=logger.status, partial=True, fallback_correlation=0)
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
    cm, cr = source.mean(axis=0), target.mean(axis=0)
    starts = [translation(-cr) * (pose.inverse() if reverse else pose) * translation(cm)
              for pose in poses]
    results, _ = _refine_shape_level(source - cm, target - cr, stats["resolutions"][-1],
                                      starts, max_points, partial=True)
    score, best, seed_index = results[0]
    if score <= stats["score"] + 1e-6:
        return transform
    result = translation(cr) * best * translation(-cm)
    result = result.inverse() if reverse else result
    stats.update(chain_backup_used=True, chain_backup_initial_score=stats["score"],
                 chain_backup_seed=labels[seed_index] or "original pose", score=score)
    return result


def _shape_align_one(session, atoms, to, resolution=None, max_points=12000,
                     extra_starts=0, move=True, radius=None, correlation=False,
                     color_fit=False, correlation_resolution=8.0, correlation_window=12.0,
                     sphere_center=None, partial=None, fallback_correlation=None, objective="density"):
    """Rigidly align one atomic model to another by shape, without sequences.

    Atom specifications restrict the fitting regions; the entire moving model
    moves. Hydrogens and water are excluded. The Gaussian search range is
    size-derived unless resolution overrides its final endpoint. max_points limits
    sampled density voxels per stage. extra_starts adds reproducible global
    orientations if the principal-axis starts are insufficient.
    move false reports a fit without applying it. The reference, coordinates,
    display, and camera are left unchanged.

    objective gradient is experimental whole-density gradient scoring, used in
    both optimization and ranking. It never switches to partial/chain fitting;
    its quality cutoff defaults to 0.5 and can trigger a broader gradient search.

    In density mode, weak partial fits try substantial selected polymer chains as
    pose seeds. Candidates are refined and compared on the same complete inputs
    at the final resolution; chains never move independently. partial false or
    fallback_correlation 0 disables this additional quality-triggered search.

    radius enables local refinement of an already aligned pair using only atom
    centers within that many Angstroms of the current center of rotation. The
    sphere is fixed in scene coordinates and each model's atoms are chosen once
    before fitting. Refine only the existing relative pose, with a size-derived
    range no broader than twice its final resolution. The whole model moves rigidly.

    correlation true optionally saves local density Pearson correlations as
    shape_correlation on fitted atoms and residue means. color_fit true implies
    correlation and colors fitted atoms/cartoons red (0) to blue (1); undefined
    values are gray. correlation_resolution sets the density resolution and
    correlation_window the Gaussian neighbourhood FWHM, both in Angstroms.
    These optional steps require move true, and their cost is reported separately.
    """
    from math import isfinite
    from time import perf_counter
    import numpy as np
    from chimerax.atomic import Residue
    from chimerax.core.errors import UserError
    from chimerax.core.undo import UndoState

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
    fitting_started = perf_counter()
    local_transform, stats = fit_shape_points(
        source_xyz, reference_xyz, resolution, max_points, extra_starts,
        status=session.logger.status,
        initial_transform=reference_scene.inverse() * source_scene,
        local=radius is not None,
        partial=partial, fallback_correlation=fallback_correlation, objective=objective,
    )
    if (objective == "density" and radius is None and partial is not False and stats["partial"] and
            stats["score"] < min(_CHAIN_BACKUP_CORRELATION, fallback_correlation)):
        local_transform = _chain_shape_backup(
            source_xyz, reference_xyz, _shape_chain_regions(source_atoms),
            _shape_chain_regions(reference_atoms), local_transform, stats,
            resolution, max_points, extra_starts, session.logger,
        )
    stats["seconds"] = perf_counter() - fitting_started
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


def _shape_align_atomic_models(session, atoms, to, resolution=None, max_points=12000,
                extra_starts=0, move=True, radius=None, correlation=False,
                color_fit=False, correlation_resolution=8.0, correlation_window=12.0,
                partial=None, fallback_correlation=None, objective="density"):
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
    from chimerax.core.errors import UserError

    targets = to.unique_structures
    moving = atoms.unique_structures
    if len(targets) != 1 or not len(moving):
        raise UserError("Specify moving atomic models and exactly one reference model.")
    reference = targets[0]
    if reference in moving:
        raise UserError("The reference must not be included in the moving selection.")
    for model in moving:
        if reference in model.drawing_lineage or model in reference.drawing_lineage:
            raise UserError("The moving and reference models must not contain one another.")
        if any(other in model.drawing_lineage[:-1] for other in moving):
            raise UserError("Moving models must not contain one another; specify independent models.")
    options = dict(resolution=resolution, max_points=max_points, extra_starts=extra_starts,
                   move=move, radius=radius, correlation=correlation, color_fit=color_fit,
                   correlation_resolution=correlation_resolution, correlation_window=correlation_window,
                   partial=partial, fallback_correlation=fallback_correlation, objective=objective)
    if len(moving) == 1:
        return _shape_align_one(session, atoms, to, **options)
    if radius is not None:
        # Read once: automatic CoFR modes must not change the sphere mid-batch.
        options["sphere_center"] = session.main_view.center_of_rotation.copy()
    results = []
    context = session.undo.aggregate("shape align models") if move else nullcontext()
    with context as batch:
        try:
            for _, selected in sorted(atoms.by_structure, key=lambda item: item[0].id):
                results.append(_shape_align_one(session, selected, to, **options))
        except BaseException:
            if batch is not None:
                for action in reversed(batch.actions):
                    action.undo()
                if batch.actions:
                    session.logger.warning("Shape alignment batch failed; earlier changes were undone.")
            raise
    return results


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
        return list(selection.by_structure)
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

    def __init__(self, model, atoms, center=None, radius=None):
        from itertools import product
        import numpy as np
        from chimerax.atomic import Residue
        from chimerax.core.errors import UserError
        from chimerax.map_data import ArrayGridData

        self.model, self.atoms = model, atoms
        self.is_map = atoms is None
        self._cache = {}
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
            self.xyz = xyz
            self.center = xyz.mean(axis=0)
            self.landmarks = xyz - self.center
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
        self.lower, self.upper = lower, upper

    def grid(self, resolution):
        import numpy as np
        from chimerax.core.errors import UserError
        from chimerax.map_data import ArrayGridData

        if resolution in self._cache:
            return self._cache[resolution]
        if not self.is_map:
            grid = _shape_grid(*_shape_points(self.landmarks, resolution / 6), resolution)
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
        self._cache.clear()  # One level per input, not a retained pyramid of large arrays.
        self._cache[resolution] = grid
        return grid


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
    field = target.matrix()
    if objective == "gradient":
        _, first = _gradient_density_field(source)
        field, second = _gradient_density_field(target)
    else:
        first = np.einsum("ijk,ijk->", source.matrix(), source.matrix(), dtype=np.float64)
        second = np.einsum("ijk,ijk->", field, field, dtype=np.float64)
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
    field = target.matrix()
    if objective == "gradient":
        _, first = _gradient_density_field(source)
        field, second = _gradient_density_field(target)
    else:
        first = np.sum(source.matrix().astype(float) ** 2)
        second = np.sum(field.astype(float) ** 2)
    target_fft = rfftn(field, fft_shape)
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
    levels, stats = _shape_resolution_schedule(source.landmarks, target.landmarks,
                                              options["resolution"], radius is not None)
    initial = (translation(-target.center) * target.model.scene_position.inverse() *
               source.model.scene_position * translation(source.center))
    starts = [initial] if radius is not None else [*_shape_starts(
        source.landmarks, target.landmarks, options["extra_starts"]), initial]
    search_radius = options.get("radius_search")
    if search_radius is not None:
        # A regional global search must explore translations even when PCA
        # already gives a good score. Only the reference density was cropped.
        source.model.session.logger.status("Shape alignment: global search of the reference CoFR sphere...")
        starts.extend(_shape_map_fft_starts(
            source.grid(levels[0]), target.grid(levels[0]), source.landmarks,
            target.landmarks, objective, options["extra_starts"]))

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
        seeds = _shape_map_fft_starts(source.grid(blur), target.grid(blur), source.landmarks,
                                      target.landmarks, objective, options["extra_starts"])
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
    source = _ShapeDensityInput(*source_item, center if radius is not None else None, radius)
    auto_search = options.get("radius_search") is True
    if auto_search:
        search_radius = _shape_search_radius(source)
        options = dict(options, radius_search=search_radius)
        session.logger.info(f"Automatic radius_search for #{source.model.id_string}: {search_radius:.3g} A "
                            "(1.5 times the moving input's enclosing radius).")
    target = target_input if target_input is not None else _ShapeDensityInput(
        *target_item, center, radius if radius is not None else options.get("radius_search"))
    pose, stats = _fit_shape_density_inputs(source, target, options, center)
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


def shape_align(session, moving, to, resolution=None, max_points=12000,
                extra_starts=0, move=True, radius=None, correlation=False,
                color_fit=False, correlation_resolution=8.0, correlation_window=12.0,
                partial=None, fallback_correlation=None, objective="density",
                radius_search=None, radius_refine=None):
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

    sources, targets = _shape_input_items(moving), _shape_input_items(to)
    if not sources or len(targets) != 1:
        raise UserError("Specify moving atomic models/maps and exactly one reference model or map.")
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
    for model, _ in sources:
        if target in model.drawing_lineage or model in target.drawing_lineage:
            raise UserError("The reference must be distinct from, and not contain, any moving model.")
        if any(other in model.drawing_lineage[:-1] for other, _ in sources):
            raise UserError("Moving inputs must be independent models, not ancestor/descendant pairs.")
    options = dict(resolution=resolution, max_points=max_points, extra_starts=extra_starts, move=move,
                   radius=radius, correlation=correlation, color_fit=color_fit,
                   correlation_resolution=correlation_resolution, correlation_window=correlation_window,
                   partial=partial, fallback_correlation=fallback_correlation, objective=objective)
    if all(atoms is not None for _, atoms in [*sources, *targets]):
        return _shape_align_atomic_models(session, moving.atoms if hasattr(moving, "atoms") else moving,
                                          to.atoms if hasattr(to, "atoms") else to, **options)
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
    with session.undo.aggregate("shape align models") if move and len(sources) > 1 else nullcontext() as batch:
        try:
            for item in sources:
                if item[1] is not None and targets[0][1] is not None:
                    result = _shape_align_one(session, item[1], targets[0][1], sphere_center=center, **options)
                else:
                    started = perf_counter()
                    # Automatic search sizes are per moving input, not inherited
                    # from the first member of a batch.
                    if target_input is None and radius_search is not True:
                        target_input = _ShapeDensityInput(*targets[0], center, region_radius)
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
        required=[("moving", ObjectsArg)],
        keyword=[("to", ObjectsArg), ("resolution", FloatArg),
                 ("max_points", IntArg), ("extra_starts", IntArg), ("move", BoolArg),
                 ("radius_search", SearchRadiusArg(FloatArg, NoArg, name="radius in Angstroms or no value (automatic)")),
                 ("radius_refine", FloatArg), ("radius", FloatArg),
                 ("correlation", BoolArg), ("color_fit", BoolArg),
                 ("correlation_resolution", FloatArg), ("correlation_window", FloatArg),
                 ("partial", BoolArg), ("fallback_correlation", FloatArg),
                 ("objective", EnumOf(("density", "gradient")))],
        required_arguments=["to"],
        synopsis="Align models/maps by shape; radius_search searches a map region, radius_refine refines locally",
    ), shape_align, logger=logger)



def register_command(logger):
    from chimerax.atomic import AtomsArg
    from chimerax.core.commands import BoolArg, CmdDesc, FloatArg, IntArg, ModelsArg, StringArg, register

    _register_shape_align_command(logger)

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
