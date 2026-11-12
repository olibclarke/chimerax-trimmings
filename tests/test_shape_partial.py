"""Partial-template placement and quality-triggered fallback regressions."""

import unittest
from unittest.mock import patch

import numpy as np

from test_shape_align import cloud, model, session, shape, rotation, translation


class PartialShapeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.small = cloud()
        background = np.random.default_rng(12).normal((110, 15, -10), (15, 23, 12), (1500, 3))
        cls.large = np.concatenate((cls.small, background))
        cls.displacement = translation((80, -60, 120)) * rotation((1, -2, 3), 137)

    def test_subcomplex_is_placed_on_known_corresponding_points(self):
        moving = self.displacement.transform_points(self.small)
        transform, stats = shape.fit_shape_points(moving, self.large)
        error = np.sqrt(np.mean(np.sum((transform.transform_points(moving) - self.small) ** 2, axis=1)))
        self.assertLess(error, 0.7)
        self.assertTrue(stats["partial"])
        self.assertEqual(stats["partial_trigger"], "size")
        self.assertFalse(stats["partial_template_reversed"])
        self.assertGreater(stats["partial_orientations"], 24)

    def test_larger_moving_model_uses_inverse_template_placement(self):
        moving = self.displacement.transform_points(self.large)
        transform, stats = shape.fit_shape_points(moving, self.small)
        error = np.sqrt(np.mean(np.sum((transform.transform_points(moving[:len(self.small)]) - self.small) ** 2, axis=1)))
        self.assertLess(error, 0.7)
        self.assertTrue(stats["partial_template_reversed"])
        self.assertEqual(stats["moving_atoms"], len(self.large))
        self.assertEqual(stats["reference_atoms"], len(self.small))

    def test_partial_search_is_not_tied_to_file_axes(self):
        target_pose = translation((-70, 20, 10)) * rotation((.3, -2, 4), 117)
        moving = self.displacement.transform_points(self.small)
        reference = target_pose.transform_points(self.large)
        transform, _ = shape.fit_shape_points(moving, reference)
        desired = target_pose.transform_points(self.small)
        error = np.sqrt(np.mean(np.sum((transform.transform_points(moving) - desired) ** 2, axis=1)))
        self.assertLess(error, 0.7)

    def test_low_correlation_triggers_backup(self):
        moving = self.displacement.transform_points(self.small)
        transform, stats = shape.fit_shape_points(moving, self.small, fallback_correlation=1)
        self.assertTrue(stats["partial"])
        self.assertEqual(stats["partial_trigger"], "low_correlation")
        self.assertLess(stats["fallback_initial_score"], 1)
        error = np.sqrt(np.mean(np.sum((transform.transform_points(moving) - self.small) ** 2, axis=1)))
        self.assertLess(error, .7)

    def test_good_whole_fit_skips_partial_search(self):
        with patch.object(shape, "_partial_shape_starts", side_effect=AssertionError("partial search")):
            _, stats = shape.fit_shape_points(self.small, self.small)
        self.assertFalse(stats["partial"])

    def test_explicit_false_disables_size_and_quality_triggers(self):
        with patch.object(shape, "_partial_shape_starts", side_effect=AssertionError("partial search")):
            _, stats = shape.fit_shape_points(self.small, self.large, partial=False, fallback_correlation=1)
        self.assertFalse(stats["partial"])

    def test_local_mode_never_restarts_globally(self):
        with patch.object(shape, "_partial_shape_starts", side_effect=AssertionError("partial search")):
            _, stats = shape.fit_shape_points(self.small, self.small, local=True, fallback_correlation=1)
        self.assertFalse(stats["partial"])
        self.assertEqual(stats["starts"], 1)

    def test_partial_fft_memory_guard(self):
        from chimerax.core.errors import UserError

        with patch.object(shape, "_MAX_PARTIAL_FFT_VOXELS", 10), self.assertRaises(UserError):
            shape._partial_shape_starts(self.small, self.large, 40)


class PartialCommandTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = session()
        shape.register_command(cls.s.logger)

    def setUp(self):
        self.s.models.close(self.s.models.list())
        self.s.undo.clear()
        self.source = model(self.s, "source", cloud()[::5])
        self.target = model(self.s, "target", cloud()[::5])
        self.source.position = translation((50, 70, -80)) * rotation((1, 2, -3), 95)

    def test_partial_keyword_and_grouped_undo(self):
        from chimerax.core.commands import run

        old = self.source.position
        _, stats = run(self.s, f"shapealign #{self.source.id_string} to #{self.target.id_string} partial true")
        self.assertTrue(stats["partial"])
        self.assertEqual(stats["partial_trigger"], "requested")
        self.s.undo.undo()
        np.testing.assert_array_equal(old.matrix, self.source.position.matrix)

    def test_invalid_cutoff_and_local_partial_request(self):
        from chimerax.core.errors import UserError

        for value in (-.1, 1.1, float("nan")):
            with self.subTest(value=value), self.assertRaises(UserError):
                shape.shape_align(self.s, self.source.atoms, self.target.atoms, fallback_correlation=value)
        with self.assertRaises(UserError):
            shape.shape_align(self.s, self.source.atoms, self.target.atoms, radius=30, partial=True)


if __name__ == "__main__":
    unittest.main()
