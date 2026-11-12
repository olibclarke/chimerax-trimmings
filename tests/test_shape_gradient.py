"""Whole-density gradient scoring, search isolation, and native command checks."""

import unittest
from unittest.mock import patch

import numpy as np

from test_shape_align import cloud, model, session, shape, Place, rotation, translation


class GradientMathTests(unittest.TestCase):
    def test_discrete_integration_by_parts(self):
        from chimerax.map_data import ArrayGridData

        rng = np.random.default_rng(4)
        density = np.pad(rng.random((7, 6, 5)).astype(np.float32), 1)
        grid = ArrayGridData(density, (0, 0, 0), (1.5,) * 3)
        field, energy = shape._gradient_density_field(grid)
        padded = np.pad(density.astype(np.float64), 1)
        expected = sum(np.sum(np.diff(padded, axis=axis) ** 2) for axis in range(3)) / 1.5 ** 2
        self.assertAlmostEqual(energy / expected, 1, places=6)
        self.assertTrue(np.any(field < 0))
        self.assertTrue(np.any(field > 0))
        np.testing.assert_array_equal(grid.matrix(), density)

    def test_self_score_and_native_refinement(self):
        xyz = cloud()
        results, _ = shape._refine_shape_level(xyz, xyz, 4, [Place()], 12000, objective="gradient")
        self.assertAlmostEqual(results[0][0], 1, delta=.01)
        initial = translation((1, -1, .5)) * rotation((1, 2, -3), 5)
        fitted, stats = shape.fit_shape_points(xyz, xyz, local=True, initial_transform=initial,
                                               objective="gradient")
        error = np.sqrt(np.mean(np.sum((fitted.transform_points(xyz) - xyz) ** 2, axis=1)))
        self.assertLess(error, .3)
        self.assertFalse(stats["partial"])
        self.assertFalse(stats["gradient_search"])

    def test_unmatched_density_contributes_to_global_normalization(self):
        xyz = cloud()
        enlarged = np.concatenate((xyz, xyz + (120, 0, 0)))
        with patch("chimerax.map_fit.fitmap.locate_maximum", return_value=(Place(), {})):
            original, _ = shape._refine_shape_level(xyz, xyz, 4, [Place()], 12000, objective="gradient")
            extra, _ = shape._refine_shape_level(xyz, enlarged, 4, [Place()], 12000, objective="gradient")
        self.assertGreater(original[0][0], .99)
        self.assertLess(extra[0][0], .8)
        self.assertGreater(extra[0][0], .6)

    def test_known_rigid_transform_without_array_changes(self):
        xyz = cloud()
        moving = (translation((80, -60, 40)) * rotation((1, -2, 3), 137)).transform_points(xyz)
        saved = moving.copy()
        fitted, stats = shape.fit_shape_points(moving, xyz, objective="gradient", fallback_correlation=0)
        error = np.sqrt(np.mean(np.sum((fitted.transform_points(moving) - xyz) ** 2, axis=1)))
        self.assertLess(error, .7)
        self.assertEqual(stats["objective"], "gradient")
        self.assertFalse(stats["partial"])
        np.testing.assert_array_equal(moving, saved)

    def test_arbitrary_reference_coordinate_frame(self):
        xyz = cloud()
        moving = (translation((80, -60, 40)) * rotation((1, -2, 3), 137)).transform_points(xyz)
        reference = (translation((-60, 40, 25)) * rotation((-3, 1, 2), 76)).transform_points(xyz)
        fitted, _ = shape.fit_shape_points(moving, reference, objective="gradient", fallback_correlation=0)
        error = np.sqrt(np.mean(np.sum((fitted.transform_points(moving) - reference) ** 2, axis=1)))
        self.assertLess(error, .7)

    def test_local_mode_preserves_initial_pose_without_any_global_search(self):
        initial = translation((2, -1, 3)) * rotation((1, 2, 3), 7)
        with patch.object(shape, "_shape_starts", side_effect=AssertionError("PCA search")), \
                patch.object(shape, "_partial_shape_starts", side_effect=AssertionError("FFT search")), \
                patch("chimerax.map_fit.fitmap.locate_maximum", return_value=(Place(), {})):
            fitted, stats = shape.fit_shape_points(cloud(), cloud()[:400] + (10, 20, 30), local=True,
                                                   initial_transform=initial, objective="gradient")
        np.testing.assert_allclose(fitted.matrix, initial.matrix, atol=1e-12)
        self.assertTrue(stats["local"])

    def test_fft_backup_keeps_gradient_objective_and_cannot_degrade_score(self):
        xyz = cloud()
        search_stats = {"partial_orientations": 1, "partial_fft_voxels": 100}
        count = 0

        def refine(*args, **kwargs):
            nonlocal count
            count += 1
            self.assertEqual(kwargs["objective"], "gradient")
            self.assertFalse(kwargs.get("partial", False))
            # Three initial levels, then two FFT-seeded levels that score worse.
            return [(.4 if count <= 3 else .1, Place(), 0)], xyz

        with patch.object(shape, "_refine_shape_level", side_effect=refine), \
                patch.object(shape, "_partial_shape_starts", return_value=([Place()], search_stats)) as fft:
            _, stats = shape.fit_shape_points(xyz, xyz, objective="gradient")
        self.assertTrue(stats["gradient_search"])
        self.assertFalse(stats["gradient_search_improved"])
        self.assertEqual(stats["score"], .4)
        self.assertEqual(fft.call_args.kwargs["objective"], "gradient")

    def test_zero_cutoff_disables_search_even_for_negative_scores(self):
        xyz = cloud()
        with patch.object(shape, "_refine_shape_level", return_value=([(-.1, Place(), 0)], xyz)), \
                patch.object(shape, "_partial_shape_starts", side_effect=AssertionError("FFT search")):
            _, stats = shape.fit_shape_points(xyz, xyz, objective="gradient", fallback_correlation=0)
        self.assertFalse(stats["gradient_search"])

    def test_actual_fft_search_places_the_whole_template(self):
        xyz = cloud()
        xyz -= xyz.mean(axis=0)
        shift = np.array((20, -15, 10))
        starts, _ = shape._partial_shape_starts(xyz, xyz + shift, 12,
                                               objective="gradient", max_starts=4)
        errors = [np.sqrt(np.mean(np.sum((pose.transform_points(xyz) - (xyz + shift)) ** 2, axis=1)))
                  for pose in starts]
        self.assertLess(min(errors), 3)


class GradientCommandTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = session()
        shape.register_command(cls.s.logger)

    def setUp(self):
        self.s.models.close(self.s.models.list())
        self.s.undo.clear()
        self.source = model(self.s, "source", cloud()[::5])
        self.target = model(self.s, "target", cloud()[::5])

    def test_command_and_undo_without_partial_or_chain_fallback(self):
        from chimerax.core.commands import run

        self.source.position = translation((80, -60, 40)) * rotation((1, -2, 3), 137)
        old = self.source.position
        fixed = self.target.position
        coords = self.source.atoms.coords.copy()
        ids = [m.id for m in self.s.models.list()]
        with patch.object(shape, "_chain_shape_backup", side_effect=AssertionError("chain fit")):
            _, stats = run(self.s, f"shapealign #{self.source.id_string} to #{self.target.id_string} "
                          "objective gradient fallback_correlation 0")
        self.assertFalse(stats["partial"])
        self.assertEqual(stats["objective"], "gradient")
        np.testing.assert_array_equal(coords, self.source.atoms.coords)
        np.testing.assert_array_equal(fixed.matrix, self.target.position.matrix)
        self.assertEqual(ids, [m.id for m in self.s.models.list()])
        fitted = self.source.position
        self.s.undo.undo()
        np.testing.assert_array_equal(old.matrix, self.source.position.matrix)
        self.s.undo.redo()
        np.testing.assert_array_equal(fitted.matrix, self.source.position.matrix)

    def test_partial_true_is_rejected_not_silently_applied(self):
        from chimerax.core.errors import UserError

        with self.assertRaises(UserError):
            shape.shape_align(self.s, self.source.atoms, self.target.atoms, objective="gradient", partial=True)
        with self.assertRaises(UserError):
            shape.fit_shape_points(cloud(), cloud(), objective="gradient", partial=True)

    def test_invalid_objective(self):
        from chimerax.core.errors import UserError

        with self.assertRaises(UserError):
            shape.shape_align(self.s, self.source.atoms, self.target.atoms, objective="unknown")

    def test_correlation_coloring_remains_independent_and_dry_run_preserves_pose(self):
        before = self.source.position
        _, stats = shape.shape_align(self.s, self.source.atoms, self.target.atoms,
                                      objective="gradient", move=False)
        np.testing.assert_array_equal(self.source.position.matrix, before.matrix)
        self.assertFalse(self.s.undo.undo_stack)
        _, stats = shape.shape_align(self.s, self.source.atoms, self.target.atoms,
                                      objective="gradient", color_fit=True)
        self.assertEqual(stats["objective"], "gradient")
        self.assertGreater(stats["correlation_defined_atoms"], 0)
        self.assertTrue(hasattr(self.source.atoms[0], "shape_correlation"))


if __name__ == "__main__":
    unittest.main()
