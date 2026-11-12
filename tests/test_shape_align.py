"""Run with ChimeraX's bundled Python: python3.14 -m unittest discover -s tests."""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import chimerax_trimmings as shape

from chimerax.geometry import Place, Places, rotation, translation


def session():
    from chimerax.atomic import initialize_atomic
    from chimerax.core import core_settings
    from chimerax.core.attributes import RegAttrManager
    from chimerax.core.nogui import UI
    from chimerax.core.session import Session
    from chimerax.core.undo import Undo

    s = Session()
    s.ui = UI(s)
    core_settings.init(s)
    RegAttrManager(s)
    s.undo = Undo(s)
    initialize_atomic(s)
    return s


def cloud():
    rng = np.random.default_rng(7)
    centers = [(0, 0, 0), (35, 0, 5), (7, 20, -3), (-15, -8, 15)]
    return np.concatenate([rng.normal(center, 3, size=(180, 3)) for center in centers])


def model(s, name, xyz):
    from chimerax.atomic import AtomicStructure, Element

    m = AtomicStructure(s, name=name, auto_style=False, log_info=False)
    residue = m.new_residue("LIG", "A", 1)
    for index, coord in enumerate(xyz):
        atom = m.new_atom(f"C{index}", Element.get_element("C"))
        residue.add_atom(atom)
        atom.coord = coord
    s.models.add([m])
    return m


class ShapeMathTests(unittest.TestCase):
    def test_voxel_centroids_preserve_atom_count_and_centroid(self):
        xyz = cloud()
        points, weights = shape._shape_points(xyz, 5)
        self.assertEqual(weights.sum(), len(xyz))
        np.testing.assert_allclose(np.average(points, weights=weights, axis=0), xyz.mean(axis=0), atol=1e-5)

    def test_starts_are_proper_and_reproducible(self):
        xyz = cloud()
        reference = xyz + (100, 200, -300)
        starts = shape._shape_starts(xyz, reference, 12)
        self.assertEqual(len(starts), 38)
        repeated = shape._shape_starts(xyz, reference, 12)
        for a, b in zip(starts, repeated):
            self.assertAlmostEqual(np.linalg.det(a.matrix[:, :3]), 1)
            np.testing.assert_array_equal(a.matrix, b.matrix)

    def test_rigid_fit_does_not_modify_arrays(self):
        xyz = cloud()
        reference = xyz.copy()
        displaced = (translation((110, -50, 80)) * rotation((2, -1, 3), 137)).transform_points(xyz)
        original = displaced.copy()
        transform, stats = shape.fit_shape_points(displaced, reference)
        rmsd = np.sqrt(np.mean(np.sum((transform.transform_points(displaced) - reference) ** 2, axis=1)))
        self.assertLess(rmsd, 0.5)
        self.assertGreater(stats["score"], 0.9)
        np.testing.assert_array_equal(displaced, original)
        np.testing.assert_array_equal(xyz, reference)

    def test_fourfold_shape(self):
        from scipy.spatial import cKDTree

        rng = np.random.default_rng(5)
        arm = rng.normal((60, 4, 0), (18, 5, 7), (500, 3))
        reference = np.concatenate([rotation((0, 0, 1), angle).transform_points(arm)
                                    for angle in (0, 90, 180, 270)])
        moving = (translation((300, 50, -20)) * rotation((1, 2, 3), 123)).transform_points(reference)
        transform, _ = shape.fit_shape_points(moving, reference)
        distances = cKDTree(reference).query(transform.transform_points(moving))[0]
        # Any of the four symmetry-equivalent orientations is legitimate.
        self.assertLess(np.quantile(distances, 0.95), 0.7)

    def test_grid_memory_guard(self):
        from chimerax.core.errors import UserError

        with self.assertRaises(UserError):
            shape._shape_grid(np.array([[0., 0., 0.], [1000., 1000., 1000.]]), np.ones(2), 2)

    def test_local_refinement_does_not_recenter_or_search(self):
        from unittest.mock import patch

        moving = cloud()
        reference = cloud()[:400] + (10, 20, 30)
        initial = translation((2, -1, 3)) * rotation((1, 2, 3), 7)
        with patch.object(shape, "_shape_starts", side_effect=AssertionError("global search")), \
             patch("chimerax.map_fit.fitmap.locate_maximum", return_value=(Place(), {})):
            transform, stats = shape.fit_shape_points(
                moving, reference, local=True, initial_transform=initial)
        np.testing.assert_allclose(transform.matrix, initial.matrix, atol=1e-12)
        self.assertEqual(stats["starts"], 1)
        self.assertEqual(stats["resolution_mode"], "automatic")
        self.assertLessEqual(stats["resolutions"][0], 2 * stats["resolutions"][-1])
        self.assertTrue(stats["local"])


class ShapeCommandTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = session()
        shape.register_command(cls.s.logger)

    def setUp(self):
        self.s.models.close(self.s.models.list())
        self.s.undo.clear()
        xyz = cloud()[::5]
        self.reference = model(self.s, "reference", xyz)
        self.moving = model(self.s, "moving", xyz)
        self.moving.position = translation((100, -70, 40)) * rotation((1, 3, -2), 101)

    def test_command_parse_apply_undo_redo_and_display_preservation(self):
        from chimerax.core.commands import run

        source, target = self.moving, self.reference
        old = source.position
        original_coords = source.atoms.coords.copy()
        original_display = source.atoms.displays.copy()
        original_colors = source.atoms.colors.copy()
        source.display = False
        target_position = target.position
        camera = self.s.main_view.camera.position
        ids = [m.id for m in self.s.models.list()]
        run(self.s, f"shapealign #{source.id_string} to #{target.id_string}")
        np.testing.assert_allclose(source.atoms.scene_coords, target.atoms.scene_coords, atol=0.5)
        np.testing.assert_array_equal(original_coords, source.atoms.coords)
        np.testing.assert_array_equal(original_display, source.atoms.displays)
        np.testing.assert_array_equal(original_colors, source.atoms.colors)
        np.testing.assert_array_equal(target_position.matrix, target.position.matrix)
        np.testing.assert_array_equal(camera.matrix, self.s.main_view.camera.position.matrix)
        self.assertFalse(source.display)
        self.assertEqual(ids, [m.id for m in self.s.models.list()])
        fitted = source.position
        self.s.undo.undo()
        np.testing.assert_array_equal(old.matrix, source.position.matrix)
        self.s.undo.redo()
        np.testing.assert_array_equal(fitted.matrix, source.position.matrix)

    def test_dry_run(self):
        old = self.moving.position
        shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, move=False)
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        self.assertFalse(self.s.undo.undo_stack)

    def test_selected_region_moves_entire_model(self):
        shape.shape_align(self.s, self.moving.atoms[:120], self.reference.atoms[:120], extra_starts=8)
        np.testing.assert_allclose(self.moving.atoms.scene_coords, self.reference.atoms.scene_coords, atol=0.5)

    def test_water_and_hydrogen_are_excluded(self):
        from chimerax.atomic import Element

        source = self.moving
        residue = source.residues[0]
        hydrogen = source.new_atom("H", Element.get_element("H"))
        residue.add_atom(hydrogen)
        hydrogen.coord = (100000, 100000, 100000)
        water = source.new_residue("HOH", "W", 1)
        oxygen = source.new_atom("O", Element.get_element("O"))
        water.add_atom(oxygen)
        oxygen.coord = (-100000, -100000, -100000)
        _, stats = shape.shape_align(self.s, source.atoms, self.reference.atoms)
        self.assertEqual(stats["moving_atoms"], len(self.reference.atoms))
        np.testing.assert_allclose(source.atoms.scene_coords[:-2], self.reference.atoms.scene_coords, atol=0.5)

    def test_scene_transform_does_not_change_fit(self):
        _, first = shape.shape_align(self.s, self.moving.atoms, self.reference.atoms)
        fitted = self.moving.atoms.scene_coords.copy()
        self.moving.position = rotation((3, 1, -5), 158) * translation((90, 70, -40))
        _, second = shape.shape_align(self.s, self.moving.atoms, self.reference.atoms)
        np.testing.assert_allclose(self.moving.atoms.scene_coords, fitted, atol=1e-5)
        self.assertAlmostEqual(first["score"], second["score"], places=8)

    def test_parent_transform(self):
        from chimerax.core.models import Model

        group = Model("group", self.s)
        group.position = translation((50, -30, 10)) * rotation((0, 1, 1), 43)
        self.s.models.add([group])
        group.add([self.moving])
        old = self.moving.position
        group_old = group.position
        shape.shape_align(self.s, self.moving.atoms, self.reference.atoms)
        np.testing.assert_allclose(self.moving.atoms.scene_coords, self.reference.atoms.scene_coords, atol=0.5)
        np.testing.assert_array_equal(group_old.matrix, group.position.matrix)
        self.s.undo.undo()
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)

    def test_reject_invalid_input(self):
        from chimerax.core.errors import UserError

        for kwargs in ({"resolution": 0}, {"resolution": float("nan")},
                       {"max_points": 3}, {"extra_starts": -1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(UserError):
                shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, **kwargs)
        with self.assertRaises(UserError):
            shape.shape_align(self.s, self.moving.atoms, self.moving.atoms)
        self.moving.positions = Places([Place(), translation((2, 0, 0))])
        with self.assertRaises(UserError):
            shape.shape_align(self.s, self.moving.atoms, self.reference.atoms)


class LocalShapeCommandTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = session()
        shape.register_command(cls.s.logger)

    def setUp(self):
        self.s.models.close(self.s.models.list())
        self.s.undo.clear()
        core = cloud()[::4] * 0.5
        outside = np.random.default_rng(5).normal((150, 0, 0), 10, (500, 3))
        self.core_count = len(core)
        self.reference = model(self.s, "reference", np.concatenate((core, outside)))
        # Make the surrounding structure disagree strongly with the local core.
        self.moving = model(self.s, "moving", np.concatenate((core, outside + (0, 70, 0))))
        self.common = translation((250, -130, 40)) * rotation((1, 2, -3), 39)
        self.reference.position = self.common
        self.moving.position = self.common * translation((0.8, -0.6, 0.4)) * rotation((1, 2, -1), 5)
        self.s.main_view.center_of_rotation = np.array((250, -130, 40), dtype=float)
        self.radius = 30

    def test_local_command_ignores_outside_and_preserves_state(self):
        from chimerax.core.commands import run
        from unittest.mock import patch

        old = self.moving.position
        coords = self.moving.atoms.coords.copy()
        target = self.reference.position
        center = self.s.main_view.center_of_rotation.copy()
        cofr_method = self.s.main_view.center_of_rotation_method
        camera = self.s.main_view.camera.position
        display = self.moving.atoms.displays.copy()
        colors = self.moving.atoms.colors.copy()
        self.moving.display = False
        ids = [m.id for m in self.s.models.list()]
        with patch.object(shape, "_shape_starts", side_effect=AssertionError("global search")):
            transform, stats = run(self.s,
                f"shapealign #{self.moving.id_string} to #{self.reference.id_string} radius 30 resolution 4")
        self.assertEqual(stats["moving_atoms"], self.core_count)
        self.assertEqual(stats["reference_atoms"], self.core_count)
        self.assertEqual(stats["starts"], 1)
        self.assertEqual(stats["resolutions"][-1], 4)
        self.assertLessEqual(stats["resolutions"][0], 8)
        self.assertEqual(stats["resolution_mode"], "explicit")
        self.assertEqual(stats["radius"], 30)
        np.testing.assert_allclose(
            self.moving.atoms.scene_coords[:self.core_count],
            self.reference.atoms.scene_coords[:self.core_count], atol=0.3)
        # The full model, including the excluded surroundings, gets one transform.
        expected = transform.transform_points(old.transform_points(coords))
        np.testing.assert_allclose(self.moving.atoms.scene_coords, expected, atol=1e-10)
        np.testing.assert_array_equal(coords, self.moving.atoms.coords)
        np.testing.assert_array_equal(target.matrix, self.reference.position.matrix)
        np.testing.assert_array_equal(center, self.s.main_view.center_of_rotation)
        self.assertEqual(cofr_method, self.s.main_view.center_of_rotation_method)
        np.testing.assert_array_equal(camera.matrix, self.s.main_view.camera.position.matrix)
        np.testing.assert_array_equal(display, self.moving.atoms.displays)
        np.testing.assert_array_equal(colors, self.moving.atoms.colors)
        np.testing.assert_array_equal(center, stats["center_of_rotation"])
        self.assertFalse(self.moving.display)
        self.assertEqual(ids, [m.id for m in self.s.models.list()])
        fitted = self.moving.position
        self.s.undo.undo()
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        self.s.undo.redo()
        np.testing.assert_array_equal(fitted.matrix, self.moving.position.matrix)

    def test_local_dry_run(self):
        old = self.moving.position
        _, stats = shape.shape_align(self.s, self.moving.atoms, self.reference.atoms,
                                     radius=self.radius, move=False)
        self.assertEqual(stats["moving_atoms"], self.core_count)
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        self.assertFalse(self.s.undo.undo_stack)

    def test_local_sphere_uses_model_hierarchy(self):
        from chimerax.core.models import Model

        group = Model("group", self.s)
        self.s.models.add([group])
        group.add([self.moving])
        self.moving.position = self.common.inverse() * self.moving.position
        group.position = self.common
        group_old = group.position
        _, stats = shape.shape_align(self.s, self.moving.atoms, self.reference.atoms,
                                     radius=self.radius, resolution=4)
        self.assertEqual(stats["moving_atoms"], self.core_count)
        np.testing.assert_allclose(
            self.moving.scene_position.transform_points(self.moving.atoms.coords[:self.core_count]),
            self.reference.atoms.scene_coords[:self.core_count], atol=0.3)
        np.testing.assert_array_equal(group_old.matrix, group.position.matrix)

    def test_local_sphere_intersects_atom_specification(self):
        _, stats = shape.shape_align(self.s, self.moving.atoms[:80], self.reference.atoms[:80],
                                     radius=self.radius, resolution=4, move=False)
        self.assertEqual(stats["moving_atoms"], 80)
        self.assertEqual(stats["reference_atoms"], 80)

    def test_local_rejects_invalid_radius_or_empty_sphere(self):
        from chimerax.core.errors import UserError

        old = self.moving.position
        for radius in (0, -5, float("inf"), float("nan")):
            with self.subTest(radius=radius), self.assertRaises(UserError):
                shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, radius=radius)
        with self.assertRaises(UserError):
            shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, radius=30, extra_starts=8)
        self.s.main_view.center_of_rotation = np.array((10000, 10000, 10000), dtype=float)
        with self.assertRaisesRegex(UserError, "within 30 A"):
            shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, radius=30)
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        self.assertFalse(self.s.undo.undo_stack)


if __name__ == "__main__":
    unittest.main()
