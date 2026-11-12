"""Optional density correlation and independent multi-model alignment tests."""

import unittest
from unittest.mock import patch

import numpy as np

from test_shape_align import cloud, model, session, shape, rotation, translation


class CorrelationMathTests(unittest.TestCase):
    def test_field_matches_direct_weighted_pearson(self):
        rng = np.random.default_rng(3)
        first = rng.random((21, 21, 21)).astype(np.float32)
        second = (first * .7 + rng.random(first.shape) * .3).astype(np.float32)
        sigma = 1.3
        field, valid = shape._local_correlation_field(first.copy(), second.copy(), sigma)
        radius = int(3 * sigma + .5)
        axis = np.arange(-radius, radius + 1)
        kernel = np.exp(-axis ** 2 / (2 * sigma ** 2))
        weights = kernel[:, None, None] * kernel[None, :, None] * kernel[None, None, :]
        weights /= weights.sum()
        cut = slice(10 - radius, 11 + radius)
        a, b = first[cut, cut, cut], second[cut, cut, cut]
        a = a - np.sum(weights * a)
        b = b - np.sum(weights * b)
        expected = np.sum(weights * a * b) / np.sqrt(np.sum(weights * a * a) * np.sum(weights * b * b))
        self.assertEqual(valid[10, 10, 10], 1)
        self.assertAlmostEqual(field[10, 10, 10], expected, places=5)

    def test_scale_offset_and_anticorrelation(self):
        a = np.random.default_rng(5).random((25, 25, 25)).astype(np.float32)
        field, valid = shape._local_correlation_field(a.copy(), 3 * a + 7, 1.5)
        np.testing.assert_allclose(field[8:-8, 8:-8, 8:-8], 1, atol=1e-4)
        self.assertTrue(valid[8:-8, 8:-8, 8:-8].all())
        field, _ = shape._local_correlation_field(a.copy(), 1 - a, 1.5)
        np.testing.assert_allclose(field[8:-8, 8:-8, 8:-8], -1, atol=1e-4)

    def test_undefined_variance(self):
        zero = np.zeros((25, 25, 25), dtype=np.float32)
        field, valid = shape._local_correlation_field(zero.copy(), zero.copy(), 1)
        self.assertFalse(valid.any())
        self.assertTrue(np.isfinite(field).all())
        field, valid = shape._local_correlation_field(np.ones_like(zero), np.ones_like(zero), 1)
        self.assertEqual(valid[12, 12, 12], 0)

    def test_atom_sampling_self_fit_and_no_input_mutation(self):
        points = cloud()
        before = points.copy()
        scores, stats = shape.shape_correlations(points, points)
        self.assertTrue(np.isfinite(scores).all())
        np.testing.assert_allclose(scores, 1, atol=1e-5)
        np.testing.assert_array_equal(before, points)
        self.assertEqual(stats["correlation_defined_atoms"], len(points))

    def test_mismatch_is_not_reported_as_a_good_fit(self):
        points = cloud()
        scores, _ = shape.shape_correlations(points, points[:540])
        self.assertGreater(np.nanmedian(scores[:180]), .9)
        absent = scores[540:]
        self.assertGreater(np.mean(~np.isfinite(absent) | (absent < .5)), .8)

    def test_correlation_memory_guard(self):
        from chimerax.core.errors import UserError

        with patch.object(shape, "_MAX_CORRELATION_VOXELS", 10), self.assertRaises(UserError):
            shape.shape_correlations(cloud(), cloud())


class _CorrelationFixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = session()
        shape.register_command(cls.s.logger)

    def setUp(self):
        self.s.models.close(self.s.models.list())
        self.s.undo.clear()
        self.xyz = cloud()[::5]
        self.reference = model(self.s, "reference", self.xyz)
        self.moving = model(self.s, "moving", self.xyz)
        self.moving.position = translation((80, -90, 100)) * rotation((2, 3, 1), 110)


class CorrelationCommandTests(_CorrelationFixture):
    def test_default_skips_correlation_entirely(self):
        with patch.object(shape, "shape_correlations", side_effect=AssertionError("scoring enabled")):
            _, stats = shape.shape_align(self.s, self.moving.atoms, self.reference.atoms)
        self.assertNotIn("correlation_seconds", stats)
        self.assertFalse(any(hasattr(a, "shape_correlation") for a in self.moving.atoms))

    def test_color_attributes_and_undo_redo(self):
        from chimerax.core.commands import run

        source = self.moving
        source.atoms.colors = np.tile((17, 33, 52, 123), (len(source.atoms), 1)).astype(np.uint8)
        old_colors = source.atoms.colors.copy()
        old_ribbon = source.residues.ribbon_colors.copy()
        old = source.position
        reference_colors = self.reference.atoms.colors.copy()
        source.atoms[0].shape_correlation = .123
        source.residues[0].shape_correlation = .456
        _, stats = run(self.s, f"shapealign #{source.id_string} to #{self.reference.id_string} color_fit true")
        scores = np.array([a.shape_correlation for a in source.atoms])
        np.testing.assert_allclose(scores, 1, atol=1e-4)
        self.assertAlmostEqual(source.residues[0].shape_correlation, scores.mean(), places=6)
        self.assertFalse(np.array_equal(old_colors, source.atoms.colors))
        np.testing.assert_array_equal(source.atoms.colors[:, 3], old_colors[:, 3])
        np.testing.assert_array_equal(reference_colors, self.reference.atoms.colors)
        self.assertEqual(len(self.s.undo.undo_stack), 1)
        self.assertIn("correlation_seconds", stats)
        new_colors = source.atoms.colors.copy()
        self.s.undo.undo()
        np.testing.assert_array_equal(old.matrix, source.position.matrix)
        np.testing.assert_array_equal(old_colors, source.atoms.colors)
        np.testing.assert_array_equal(old_ribbon, source.residues.ribbon_colors)
        self.assertEqual(source.atoms[0].shape_correlation, .123)
        self.assertFalse(hasattr(source.atoms[1], "shape_correlation"))
        self.assertEqual(source.residues[0].shape_correlation, .456)
        self.s.undo.redo()
        np.testing.assert_array_equal(new_colors, source.atoms.colors)
        np.testing.assert_allclose([a.shape_correlation for a in source.atoms], scores)

    def test_assign_only_does_not_color(self):
        old_colors = self.moving.atoms.colors.copy()
        old_ribbon = self.moving.residues.ribbon_colors.copy()
        shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, correlation=True)
        self.assertTrue(hasattr(self.moving.atoms[0], "shape_correlation"))
        np.testing.assert_array_equal(old_colors, self.moving.atoms.colors)
        np.testing.assert_array_equal(old_ribbon, self.moving.residues.ribbon_colors)

    def test_grid_failure_does_not_move_model(self):
        from chimerax.core.errors import UserError

        old = self.moving.position
        with patch.object(shape, "_MAX_CORRELATION_VOXELS", 10), self.assertRaises(UserError):
            shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, color_fit=True)
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        self.assertFalse(self.s.undo.undo_stack)

    def test_attribute_failure_rolls_back_position(self):
        old = self.moving.position
        with patch.object(shape, "_assign_shape_correlations", side_effect=RuntimeError("test")), \
             self.assertRaises(RuntimeError):
            shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, color_fit=True)
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        self.assertFalse(self.s.undo.undo_stack)

    def test_late_coloring_failure_restores_attributes_and_colors(self):
        old = self.moving.position
        colors = self.moving.atoms.colors.copy()
        assign = shape._assign_shape_correlations

        def fail_after_assignment(*args, **kwargs):
            assign(*args, **kwargs)
            raise RuntimeError("test late failure")

        with patch.object(shape, "_assign_shape_correlations", side_effect=fail_after_assignment), \
             self.assertRaises(RuntimeError):
            shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, color_fit=True)
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        np.testing.assert_array_equal(colors, self.moving.atoms.colors)
        self.assertFalse(hasattr(self.moving.atoms[0], "shape_correlation"))
        self.assertFalse(hasattr(self.moving.residues[0], "shape_correlation"))
        self.assertFalse(self.s.undo.undo_stack)

    def test_local_scoring_scope_and_stale_attribute_cleanup(self):
        source = self.moving
        source.position = translation((.5, -.3, .2))
        self.s.main_view.center_of_rotation = np.zeros(3)
        source.atoms[-1].shape_correlation = .5
        before = source.atoms.colors.copy()
        mask = np.linalg.norm(source.position.transform_points(source.atoms.coords), axis=1) <= 10
        shape.shape_align(self.s, source.atoms, self.reference.atoms, radius=10, color_fit=True)
        self.assertFalse(hasattr(source.atoms[-1], "shape_correlation"))
        np.testing.assert_array_equal(before[~mask], source.atoms.colors[~mask])
        self.s.undo.undo()
        self.assertEqual(source.atoms[-1].shape_correlation, .5)

    def test_undefined_is_gray_and_has_no_attribute(self):
        from chimerax.core.undo import UndoState

        undo = UndoState("test")
        scores = np.full(len(self.moving.atoms), np.nan)
        shape._assign_shape_correlations(self.s, self.moving, self.moving.atoms, scores, undo, True)
        self.assertFalse(hasattr(self.moving.atoms[0], "shape_correlation"))
        np.testing.assert_array_equal(self.moving.atoms.colors[:, :3], np.full((len(scores), 3), 160))

    def test_invalid_optional_arguments(self):
        from chimerax.core.errors import UserError

        for kwargs in ({"correlation_resolution": 0}, {"correlation_window": float("nan")},
                       {"correlation": True, "move": False}, {"color_fit": True, "move": False}):
            with self.subTest(kwargs=kwargs), self.assertRaises(UserError):
                shape.shape_align(self.s, self.moving.atoms, self.reference.atoms, **kwargs)


class BatchShapeTests(_CorrelationFixture):
    def test_multiple_models_align_independently_and_undo_together(self):
        from chimerax.core.commands import run

        second = model(self.s, "second", self.xyz)
        second.position = rotation((1, -2, 3), 81) * translation((40, 80, -60))
        originals = [m.position for m in (self.moving, second)]
        result = run(self.s, f"shapealign #{self.moving.id_string}-{second.id_string} "
                     f"to #{self.reference.id_string} color_fit true")
        self.assertEqual(len(result), 2)
        for m in (self.moving, second):
            np.testing.assert_allclose(m.atoms.scene_coords, self.reference.atoms.scene_coords, atol=.5)
            self.assertTrue(hasattr(m.atoms[0], "shape_correlation"))
        self.assertEqual(len(self.s.undo.undo_stack), 1)
        self.s.undo.undo()
        for m, original in zip((self.moving, second), originals):
            np.testing.assert_array_equal(m.position.matrix, original.matrix)
            self.assertFalse(hasattr(m.atoms[0], "shape_correlation"))
        self.s.undo.redo()
        for m in (self.moving, second):
            np.testing.assert_allclose(m.atoms.scene_coords, self.reference.atoms.scene_coords, atol=.5)

    def test_batch_failure_rolls_back_earlier_models(self):
        from chimerax.core.errors import UserError

        invalid = model(self.s, "too few atoms", self.xyz[:2])
        old = self.moving.position
        colors = self.moving.atoms.colors.copy()
        with self.assertRaises(UserError):
            shape.shape_align(self.s, self.moving.atoms | invalid.atoms, self.reference.atoms, color_fit=True)
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        np.testing.assert_array_equal(colors, self.moving.atoms.colors)
        self.assertFalse(hasattr(self.moving.atoms[0], "shape_correlation"))
        self.assertFalse(self.s.undo.undo_stack)

    def test_batch_dry_run_and_reference_rejection(self):
        from chimerax.core.errors import UserError

        second = model(self.s, "second", self.xyz)
        old = self.moving.position
        result = shape.shape_align(self.s, self.moving.atoms | second.atoms, self.reference.atoms, move=False)
        self.assertEqual(len(result), 2)
        np.testing.assert_array_equal(old.matrix, self.moving.position.matrix)
        self.assertFalse(self.s.undo.undo_stack)
        with self.assertRaises(UserError):
            shape.shape_align(self.s, self.moving.atoms | self.reference.atoms, self.reference.atoms)

    def test_local_batch_keeps_one_sphere(self):
        second = model(self.s, "second", self.xyz)
        self.moving.position = translation((.5, -.3, .2))
        second.position = translation((-.4, .2, -.1))
        self.s.main_view.center_of_rotation = np.array((4, 0, 0), dtype=float)
        result = shape.shape_align(self.s, self.moving.atoms | second.atoms, self.reference.atoms, radius=20)
        np.testing.assert_array_equal(result[0][1]["center_of_rotation"], result[1][1]["center_of_rotation"])
        self.assertEqual(len(self.s.undo.undo_stack), 1)


if __name__ == "__main__":
    unittest.main()
