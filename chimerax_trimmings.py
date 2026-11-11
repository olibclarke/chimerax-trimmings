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


def register_command(logger):
    from chimerax.atomic import AtomsArg
    from chimerax.core.commands import BoolArg, CmdDesc, FloatArg, IntArg, ModelsArg, StringArg, register

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


register_command(session.logger)

# Reloading the script also upgrades an already-enabled controller in the session.
if getattr(getattr(session, "_trimmings_model_motion", None), "enabled", False):
    _model_motion_controller(session)
    _watch_new_motion_models(session)
