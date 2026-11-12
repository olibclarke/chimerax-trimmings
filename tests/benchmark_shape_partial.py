"""Compare shape placement with Matchmaker on explicitly selected chains.

Use ChimeraX's bundled Python and local mmCIF files. No downloads are performed.
"""

import argparse
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from chimerax_trimmings import shape_align
from test_shape_align import session


def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def main():
    import chimerax
    from chimerax.geometry import quaternion_rotation, translation
    from chimerax.match_maker.match import cmd_match
    from chimerax.mmcif import open_mmcif

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("moving", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--moving-chain", default="S2")
    parser.add_argument("--equivalent-moving-chain", action="append", default=[],
                        help="Additional symmetry-equivalent chain for the Matchmaker baseline; repeatable")
    parser.add_argument("--reference-chain", default="S2")
    parser.add_argument("--fit-chains", action="store_true", help="Restrict shape fitting to the same selected chains")
    parser.add_argument("--resolution", type=float, help="Override the automatic fine endpoint")
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--partial", choices=("auto", "true", "false"), default="auto")
    parser.add_argument("--objective", choices=("density", "gradient"), default="density")
    parser.add_argument("--fallback-correlation", type=float, help="Override the objective-specific search cutoff")
    parser.add_argument("--max-placement-error", type=float)
    args = parser.parse_args()
    if args.runs < 1:
        parser.error("--runs must be positive")
    s = session()
    models = []
    for path in (args.moving, args.reference):
        if not path.is_file():
            parser.error(f"No such mmCIF file: {path}")
        opened, _ = open_mmcif(s, str(path), auto_style=False, log_info=False, slider=False)
        if len(opened) != 1:
            parser.error("Each file must contain one model")
        s.models.add(opened)
        models.extend(opened)
    moving, reference = models
    source_chain = moving.residues.filter(moving.residues.chain_ids == args.moving_chain).atoms
    target_chain = reference.residues.filter(reference.residues.chain_ids == args.reference_chain).atoms
    if not len(source_chain) or not len(target_chain):
        parser.error("The requested chains were not found")

    # A standalone Python session has no application directories. Matrix lookup
    # still uses ChimeraX's stock matrices, with an empty user-matrix directory.
    with TemporaryDirectory(prefix="shapealign-mmaker-") as data_dir:
        had_dirs = hasattr(chimerax, "app_dirs_unversioned")
        if not had_dirs:
            chimerax.app_dirs_unversioned = SimpleNamespace(user_data_dir=data_dir, user_config_dir=data_dir)
        baselines = []
        try:
            for chain_id in dict.fromkeys([args.moving_chain, *args.equivalent_moving_chain]):
                selected = moving.residues.filter(moving.residues.chain_ids == chain_id).atoms
                if not len(selected):
                    parser.error(f"Equivalent moving chain {chain_id!r} was not found")
                started = perf_counter()
                matched = cmd_match(s, selected, to=target_chain, pairing="ss",
                                    show_alignment=False, verbose=False, log_parameters=False)[0]
                baselines.append({
                    "chain": chain_id, "truth": moving.scene_position,
                    "source": matched["final match atoms"].coords.copy(),
                    "target": reference.scene_position.transform_points(matched["final ref atoms"].coords),
                    "full_source": matched["full match atoms"].coords.copy(),
                    "full_target": reference.scene_position.transform_points(matched["full ref atoms"].coords),
                })
                print(json.dumps({"mmaker_chain": chain_id, "mmaker_seconds": perf_counter() - started,
                                  "retained_pairs": len(baselines[-1]["source"]),
                                  "full_pairs": len(baselines[-1]["full_source"]),
                                  "mmaker_retained_rmsd_A": matched["final RMSD"],
                                  "mmaker_full_rmsd_A": matched["full RMSD"]}), flush=True)
        finally:
            if not had_dirs:
                del chimerax.app_dirs_unversioned
    source_coords = moving.atoms.coords.copy()
    reference_position = reference.position
    ids = [m.id for m in s.models.list()]
    rng = np.random.default_rng(11)
    for index in range(args.runs):
        quaternion = rng.normal(size=4)
        quaternion /= np.linalg.norm(quaternion)
        moving.position = translation(rng.uniform(-500, 500, size=3)) * quaternion_rotation(quaternion)
        _, stats = shape_align(s, source_chain if args.fit_chains else moving.atoms,
                              target_chain if args.fit_chains else reference.atoms,
                              resolution=args.resolution,
                              objective=args.objective,
                              partial={"auto": None, "true": True, "false": False}[args.partial],
                              fallback_correlation=args.fallback_correlation)
        baseline = min(baselines, key=lambda item: rmsd(
            moving.scene_position.transform_points(item["source"]), item["truth"].transform_points(item["source"])))
        truth, paired_source, paired_target = baseline["truth"], baseline["source"], baseline["target"]
        full_source, full_target = baseline["full_source"], baseline["full_target"]
        truth_points = truth.transform_points(paired_source)
        fitted = moving.scene_position.transform_points(paired_source)
        placement_error = rmsd(fitted, truth_points)
        print(json.dumps({"run": index + 1, "seconds": stats["seconds"], "score": stats["score"],
                          "partial": stats["partial"], "trigger": stats.get("partial_trigger"),
                          "matched_chain": baseline["chain"],
                          "chain_backup_used": stats.get("chain_backup_used", False),
                          "objective": args.objective, "gradient_search": stats.get("gradient_search", False),
                          "resolutions": stats["resolutions"], "size_A": stats["object_size"],
                          "same_retained_pair_rmsd_A": rmsd(fitted, paired_target),
                          "full_matched_pair_rmsd_A": rmsd(moving.scene_position.transform_points(full_source), full_target),
                          "placement_rmsd_from_mmaker_A": placement_error,
                          "angle_from_mmaker_degrees": (moving.scene_position * truth.inverse()).axis_center_angle_shift()[2]}), flush=True)
        np.testing.assert_array_equal(source_coords, moving.atoms.coords)
        np.testing.assert_array_equal(reference_position.matrix, reference.position.matrix)
        assert ids == [m.id for m in s.models.list()]
        if args.objective == "gradient":
            assert not stats["partial"] and not stats.get("chain_backup_used", False)
        if args.max_placement_error is not None:
            assert placement_error < args.max_placement_error, "Shape placement differs too far from the Matchmaker baseline"


if __name__ == "__main__":
    main()
