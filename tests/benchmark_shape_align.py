"""Benchmark local mmCIF files with ChimeraX's Python; never download anything."""

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shape_align import shape_align
from test_shape_align import session


def main():
    from chimerax.geometry import quaternion_rotation, translation
    from chimerax.mmcif import open_mmcif

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("moving", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--runs", type=int, default=6)
    parser.add_argument("--objective", choices=("density", "gradient"), default="density")
    parser.add_argument("--color-fit", action="store_true")
    parser.add_argument("--correlation", action="store_true")
    parser.add_argument("--correlation-resolution", type=float, default=8)
    parser.add_argument("--correlation-window", type=float, default=12)
    args = parser.parse_args()
    if args.runs < 1:
        parser.error("--runs must be positive")
    for path in (args.moving, args.reference):
        if not path.is_file():
            parser.error(f"No such mmCIF file: {path}")

    s = session()
    models = []
    for path in (args.moving, args.reference):
        opened, _ = open_mmcif(s, str(path), auto_style=False, log_info=False, slider=False)
        if len(opened) != 1:
            parser.error("Each mmCIF file must contain exactly one model")
        s.models.add(opened)
        models.extend(opened)
    moving, reference = models
    source_coords = moving.atoms.coords.copy()
    reference_position = reference.position
    model_ids = [model.id for model in s.models.list()]
    print(json.dumps({"moving_atoms": moving.num_atoms, "reference_atoms": reference.num_atoms}), flush=True)

    rng = np.random.default_rng(11)
    times = []
    for index in range(args.runs):
        quaternion = rng.normal(size=4)
        quaternion /= np.linalg.norm(quaternion)
        moving.position = translation(rng.uniform(-500, 500, size=3)) * quaternion_rotation(quaternion)
        started = perf_counter()
        _, stats = shape_align(s, moving.atoms, reference.atoms,
                              objective=args.objective,
                              color_fit=args.color_fit, correlation=args.correlation,
                              correlation_resolution=args.correlation_resolution,
                              correlation_window=args.correlation_window)
        total_seconds = perf_counter() - started
        times.append(stats["seconds"])
        print(json.dumps({
            "run": index + 1, "seconds": stats["seconds"], "overlap": stats["score"],
            "total_seconds": total_seconds,
            "correlation_seconds": stats.get("correlation_seconds"),
            "correlation_defined_atoms": stats.get("correlation_defined_atoms"),
            "partial": stats["partial"], "partial_trigger": stats.get("partial_trigger"),
            "objective": args.objective, "gradient_search": stats.get("gradient_search", False),
        }), flush=True)
        np.testing.assert_array_equal(source_coords, moving.atoms.coords)
        np.testing.assert_array_equal(reference_position.matrix, reference.position.matrix)
        assert model_ids == [model.id for model in s.models.list()]
    print(json.dumps({"median_seconds": float(np.median(times)), "runs": args.runs}), flush=True)


if __name__ == "__main__":
    main()
