"""Chain-derived seeds for complexes with unmatched components on both sides."""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from test_shape_align import cloud, model, session, shape, translation


class ChainRegionTests(unittest.TestCase):
    def regions(self, chains, polymer):
        # Special methods are looked up on the type, not the instance.
        class Atoms:
            coords = cloud()[:len(chains)]
            residues = SimpleNamespace(chain_ids=np.array(chains), polymer_types=np.array(polymer))

            def __len__(self):
                return len(chains)

        return shape._shape_chain_regions(Atoms())

    def test_only_substantial_selected_polymer_chains(self):
        regions = self.regions(['A'] * 40 + ['B'] * 40 + ['C'] * 5 + ['L'] * 15,
                               [1] * 85 + [0] * 15)
        self.assertEqual([label for label, _ in regions], ['A', 'B'])
        np.testing.assert_array_equal(regions[1][1], np.arange(40, 80))

    def test_one_chain_is_not_a_smaller_seed(self):
        self.assertEqual(self.regions(['A'] * 40, [1] * 40), [])

    def test_nonpolymers_do_not_create_chain_seeds(self):
        self.assertEqual(self.regions(['A'] * 40 + ['B'] * 40, [0] * 80), [])


class ChainBackupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = session()

    def test_non_nested_shapes_keep_one_rigid_transform(self):
        common = cloud()
        moving = np.concatenate((common, common + (90, 0, 0)))
        rng = np.random.default_rng(8)
        reference = np.concatenate((common, rng.normal((0, -80, 15), (8, 9, 13), (500, 3))))
        originals = moving.copy(), reference.copy()
        regions = [('A', np.arange(len(common))), ('B', np.arange(len(common), len(moving)))]
        stats = {'score': .01, 'resolutions': (4,), 'partial_template_reversed': True}
        pose = shape._chain_shape_backup(moving, reference, regions, [], translation((50, 0, 0)),
                                          stats, None, 12000, 0, self.s.logger)
        errors = [np.sqrt(np.mean(np.sum((pose.transform_points(moving[indices]) - common) ** 2, axis=1)))
                  for _, indices in regions]
        self.assertLess(min(errors), .7)
        self.assertTrue(stats['chain_backup_used'])
        self.assertEqual(stats['chain_backup_attempts'], 2)
        np.testing.assert_array_equal(moving, originals[0])
        np.testing.assert_array_equal(reference, originals[1])

    def test_lower_whole_input_score_cannot_replace_baseline(self):
        xyz = cloud()
        stats = {'score': .8, 'resolutions': (4,), 'partial_template_reversed': False}
        baseline, candidate = translation((10, 20, 30)), translation((4, 5, 6))
        with patch.object(shape, 'fit_shape_points', return_value=(candidate, {})), \
                patch.object(shape, '_refine_shape_level', return_value=([(.7, candidate, 1)], xyz)):
            result = shape._chain_shape_backup(xyz, xyz, [('A', np.arange(400))], [], baseline,
                                               stats, None, 12000, 0, self.s.logger)
        self.assertIs(result, baseline)
        self.assertFalse(stats['chain_backup_used'])

    def test_extra_searches_are_bounded_and_use_whole_other_input(self):
        xyz = cloud()
        regions = [(str(i), np.arange(150)) for i in range(10)]
        stats = {'score': .8, 'resolutions': (4,), 'partial_template_reversed': False}
        with patch.object(shape, 'fit_shape_points', return_value=(translation((0, 0, 0)), {})) as fit, \
                patch.object(shape, '_refine_shape_level', return_value=([(.7, translation((0, 0, 0)), 0)], xyz)):
            shape._chain_shape_backup(xyz, xyz, regions, [], translation((0, 0, 0)),
                                      stats, None, 12000, 0, self.s.logger)
        self.assertEqual(fit.call_count, shape._MAX_CHAIN_BACKUP_STARTS)
        for call in fit.call_args_list:
            self.assertEqual(len(call.args[0]), 150)
            np.testing.assert_array_equal(call.args[1], xyz)

    def test_optional_seed_allocation_failure_preserves_original_fit(self):
        from chimerax.core.errors import UserError

        xyz = cloud()
        stats = {'score': .2, 'resolutions': (4,), 'partial_template_reversed': False}
        baseline = translation((0, 0, 0))
        with patch.object(shape, 'fit_shape_points', side_effect=UserError('grid too large')):
            result = shape._chain_shape_backup(xyz, xyz, [('A', np.arange(400))], [], baseline,
                                               stats, None, 12000, 0, self.s.logger)
        self.assertIs(result, baseline)
        self.assertFalse(stats['chain_backup_used'])

    def test_good_partial_fit_does_not_generate_seeds(self):
        self.s.models.close(self.s.models.list())
        source, target = model(self.s, 'source', cloud()), model(self.s, 'target', cloud())
        with patch.object(shape, '_chain_shape_backup', side_effect=AssertionError('backup')):
            shape.shape_align(self.s, source.atoms, target.atoms, partial=True)

    def test_disabled_quality_search_and_local_mode_do_not_generate_seeds(self):
        self.s.models.close(self.s.models.list())
        source, target = model(self.s, 'source', cloud()), model(self.s, 'target', cloud())
        stats = dict(score=.2, seconds=0, resolutions=(4,), moving_atoms=720,
                     reference_atoms=720, resolution_mode='explicit', object_size=40,
                     partial=True, partial_template_reversed=False)
        for options in ({'partial': False}, {'fallback_correlation': 0}, {'radius': 30}):
            with self.subTest(options=options), \
                    patch.object(shape, 'fit_shape_points', return_value=(translation((0, 0, 0)), stats.copy())), \
                    patch.object(shape, '_chain_shape_backup', side_effect=AssertionError('backup')):
                shape.shape_align(self.s, source.atoms, target.atoms, move=False, **options)

    def test_backup_pose_is_applied_once_and_is_undoable(self):
        self.s.models.close(self.s.models.list())
        self.s.undo.clear()
        source, target = model(self.s, 'source', cloud()), model(self.s, 'target', cloud())
        coords = source.atoms.coords.copy()
        old, fixed = source.position, target.position
        ids = [m.id for m in self.s.models.list()]
        pose = translation((10, 20, 30))
        stats = dict(score=.2, seconds=0, resolutions=(4,), moving_atoms=720,
                     reference_atoms=720, resolution_mode='explicit', object_size=40,
                     partial=True, partial_template_reversed=False)
        with patch.object(shape, 'fit_shape_points', return_value=(translation((0, 0, 0)), stats)), \
                patch.object(shape, '_chain_shape_backup', return_value=pose) as backup:
            shape.shape_align(self.s, source.atoms, target.atoms)
        backup.assert_called_once()
        np.testing.assert_allclose(source.position.matrix, pose.matrix)
        np.testing.assert_array_equal(source.atoms.coords, coords)
        np.testing.assert_array_equal(target.position.matrix, fixed.matrix)
        self.assertEqual([m.id for m in self.s.models.list()], ids)
        self.s.undo.undo()
        np.testing.assert_array_equal(source.position.matrix, old.matrix)
        self.s.undo.redo()
        np.testing.assert_allclose(source.position.matrix, pose.matrix)


if __name__ == '__main__':
    unittest.main()
