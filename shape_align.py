"""Sequence-independent rigid shape fitting using ChimeraX's native map fitter.

Load with ``runscript /path/to/shape_align.py``, then ``shapealign #1 to #2``.
The numerical helpers also work without a GUI or an open ChimeraX session.
"""

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
            raise UserError("extra_starts cannot be used with radius; local fitting refines only the current pose.")
        if partial:
            raise UserError("partial true cannot be used with radius; local fitting must keep the current pose.")
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
        stats.update(radius=radius, center_of_rotation=center.copy())
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


def shape_align(session, atoms, to, resolution=None, max_points=12000,
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


def register_command(logger):
    from chimerax.atomic import AtomsArg
    from chimerax.core.commands import BoolArg, CmdDesc, EnumOf, FloatArg, IntArg, register

    register("shapealign", CmdDesc(
        required=[("atoms", AtomsArg)],
        keyword=[("to", AtomsArg), ("resolution", FloatArg),
                 ("max_points", IntArg), ("extra_starts", IntArg), ("move", BoolArg),
                 ("radius", FloatArg), ("correlation", BoolArg), ("color_fit", BoolArg),
                 ("correlation_resolution", FloatArg), ("correlation_window", FloatArg),
                 ("partial", BoolArg), ("fallback_correlation", FloatArg),
                 ("objective", EnumOf(("density", "gradient")))],
        required_arguments=["to"],
        synopsis="Align models by density or experimental gradient shape scoring; optionally refine at CoFR",
    ), shape_align, logger=logger)


if "session" in globals():
    register_command(session.logger)
