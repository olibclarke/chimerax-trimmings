"""Size-driven resolution ranges, without molecule-class-specific presets."""

import unittest
import numpy as np

from test_shape_align import shape, rotation, translation


class AdaptiveResolutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.points = np.random.default_rng(15).normal(size=(1000, 3)) * (7, 9, 13)

    def test_small_and_large_objects_get_different_ranges(self):
        small, a = shape._shape_resolution_schedule(self.points, self.points)
        large, b = shape._shape_resolution_schedule(self.points * 5, self.points * 5)
        self.assertGreater(large[0], small[0])
        self.assertGreater(large[-1], small[-1])
        self.assertLess(small[0], 20)
        self.assertEqual(a["resolution_mode"], "automatic")
        self.assertGreater(b["object_size"], a["object_size"])

    def test_both_endpoints_follow_size_continuously(self):
        first, _ = shape._shape_resolution_schedule(self.points, self.points)
        second, _ = shape._shape_resolution_schedule(self.points * 1.05, self.points * 1.05)
        self.assertAlmostEqual(second[0] / first[0], 1.05)
        self.assertAlmostEqual(second[-1] / first[-1], 1.05)

    def test_duplicate_atoms_do_not_change_size_or_range(self):
        duplicated = np.tile(self.points, (10, 1))
        first, _ = shape._shape_resolution_schedule(self.points, self.points)
        second, _ = shape._shape_resolution_schedule(duplicated, duplicated)
        np.testing.assert_allclose(first, second, atol=1e-12)

    def test_rotation_and_translation_do_not_change_size(self):
        moved = (translation((100, -90, 40)) * rotation((1, 2, 3), 127)).transform_points(self.points)
        np.testing.assert_allclose(shape._shape_spatial_size(self.points),
                                   shape._shape_spatial_size(moved), rtol=1e-8)

    def test_smaller_object_sets_range_in_either_argument_order(self):
        first, a = shape._shape_resolution_schedule(self.points, self.points * 5)
        second, b = shape._shape_resolution_schedule(self.points * 5, self.points)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(a["object_size"], b["object_size"])

    def test_long_thin_objects_use_more_than_longest_dimension(self):
        thin = np.random.default_rng(12).normal(size=(1000, 3)) * (100, 3, 3)
        levels, _ = shape._shape_resolution_schedule(thin, thin)
        self.assertLess(levels[0], 25)

    def test_explicit_final_override_keeps_size_based_coarse_scale(self):
        automatic, _ = shape._shape_resolution_schedule(self.points, self.points)
        explicit, info = shape._shape_resolution_schedule(self.points, self.points, resolution=3)
        self.assertEqual(explicit[-1], 3)
        self.assertEqual(automatic[0], explicit[0])
        self.assertEqual(info["resolution_mode"], "explicit")

    def test_geometric_levels_have_no_large_jumps(self):
        for scale in (.05, 1, 5, 50):
            levels, _ = shape._shape_resolution_schedule(self.points * scale, self.points * scale)
            self.assertTrue(2 <= levels[-1] <= 8)
            ratios = np.array(levels[:-1]) / np.array(levels[1:])
            self.assertTrue((ratios <= 2 + 1e-10).all())
            self.assertTrue((ratios >= 1).all())

    def test_planar_and_tiny_clouds_have_finite_ranges(self):
        for points in (self.points * .001, np.column_stack((self.points[:, :2], np.zeros(len(self.points))))):
            levels, _ = shape._shape_resolution_schedule(points, points)
            self.assertTrue(np.isfinite(levels).all())
            self.assertGreaterEqual(levels[-1], 2)

    def test_local_range_scales_with_roi_and_stays_narrow(self):
        small, _ = shape._shape_resolution_schedule(self.points, self.points, local=True)
        large, _ = shape._shape_resolution_schedule(self.points * 5, self.points * 5, local=True)
        self.assertGreater(large[-1], small[-1])
        self.assertLessEqual(small[0], 2 * small[-1])
        self.assertLessEqual(large[0], 2 * large[-1])


if __name__ == "__main__":
    unittest.main()
