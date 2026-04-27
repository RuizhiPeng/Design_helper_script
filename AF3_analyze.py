#!/usr/bin/env python3
"""
AlphaFold3 Score Aggregator

Walks a directory of AlphaFold3 prediction folders, finds every
*_confidences.json (top-level rank-0 plus all seed-*_sample-*/ samples), and
writes a single CSV summarizing average pLDDT and pae_interaction per sample.

pae_interaction:
    - If the prediction has >1 chain, mean of PAE entries between residues
      belonging to different chains.
    - Otherwise, mean of the full PAE matrix excluding the diagonal.

Usage:
    python AF3_analyze.py -i /path/to/af3_outputs
    python AF3_analyze.py -i /path/to/af3_outputs -o results.csv
"""

import argparse
import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob

import numpy as np


def detect_cpus() -> int:
    """Detect usable CPUs, preferring SLURM allocation when available."""
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE"):
        val = os.environ.get(var)
        if val:
            try:
                return max(1, int(val.split("(")[0]))
            except ValueError:
                pass
    if hasattr(os, "sched_getaffinity"):
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except OSError:
            pass
    return max(1, os.cpu_count() or 1)


def parse_sample_id(json_path: str, design_dir: str) -> str:
    """Return 'rank0' for the top-level JSON, else the seed-*_sample-* dir name."""
    parent = os.path.basename(os.path.dirname(json_path))
    if os.path.abspath(os.path.dirname(json_path)) == os.path.abspath(design_dir):
        return "rank0"
    return parent


def find_confidences_jsons(design_dir: str):
    """Yield (json_path, sample_id) for every confidences JSON in a design folder."""
    for fname in os.listdir(design_dir):
        if (
            fname.endswith("_confidences.json")
            and not fname.endswith("_summary_confidences.json")
        ):
            full = os.path.join(design_dir, fname)
            if os.path.isfile(full):
                yield full, "rank0"

    for sub in sorted(glob(os.path.join(design_dir, "seed-*_sample-*"))):
        if not os.path.isdir(sub):
            continue
        for fname in os.listdir(sub):
            if (
                fname.endswith("_confidences.json")
                and not fname.endswith("_summary_confidences.json")
            ):
                yield os.path.join(sub, fname), os.path.basename(sub)


def compute_pae_interaction(pae: np.ndarray, chain_ids) -> float:
    """Mean inter-chain PAE; falls back to mean off-diagonal if single chain."""
    chains = np.asarray(chain_ids)
    if chains.shape[0] != pae.shape[0]:
        raise ValueError(
            f"token_chain_ids length {chains.shape[0]} != pae dim {pae.shape[0]}"
        )

    inter_mask = chains[:, None] != chains[None, :]
    if inter_mask.any():
        return float(pae[inter_mask].mean())

    n = pae.shape[0]
    if n < 2:
        return float("nan")
    off_diag_sum = pae.sum() - np.trace(pae)
    return float(off_diag_sum / (n * n - n))


def score_json(json_path: str):
    """Return dict of metrics for a single confidences JSON, or None on error."""
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except (IOError, json.JSONDecodeError) as e:
        return {"error": f"read/parse: {e}"}

    atom_plddts = data.get("atom_plddts")
    pae = data.get("pae")
    chain_ids = data.get("token_chain_ids")

    if atom_plddts is None or pae is None or chain_ids is None:
        return {"error": "missing required keys"}

    plddt_arr = np.asarray(atom_plddts, dtype=float)
    pae_arr = np.asarray(pae, dtype=float)

    return {
        "avg_plddt": float(plddt_arr.mean()),
        "pae_interaction": compute_pae_interaction(pae_arr, chain_ids),
        "n_tokens": int(pae_arr.shape[0]),
        "n_atoms": int(plddt_arr.shape[0]),
    }


def _worker(task):
    """Pool worker: task = (design_name, sample_id, json_path) -> result row."""
    design_name, sample_id, json_path = task
    metrics = score_json(json_path)
    if "error" in metrics:
        return {
            "design_folder": design_name,
            "sample_id": sample_id,
            "json_path": json_path,
            "error": metrics["error"],
        }
    return {
        "design_folder": design_name,
        "sample_id": sample_id,
        "avg_plddt": round(metrics["avg_plddt"], 3),
        "pae_interaction": round(metrics["pae_interaction"], 3),
        "n_tokens": metrics["n_tokens"],
        "n_atoms": metrics["n_atoms"],
        "json_path": json_path,
    }


def load_already_processed(csv_path: str):
    """Return set of (design_folder, sample_id) pairs already in the CSV."""
    done = set()
    if not os.path.exists(csv_path):
        return done
    try:
        with open(csv_path, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                done.add((row["design_folder"], row["sample_id"]))
    except (IOError, csv.Error, KeyError) as e:
        print(f"Warning: could not read existing CSV {csv_path}: {e}", file=sys.stderr)
    return done


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Aggregate AlphaFold3 confidence scores into a single CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python %(prog)s -i /path/to/af3_outputs
  python %(prog)s -i /path/to/af3_outputs -o /path/to/scores.csv
  python %(prog)s -i /path/to/af3_outputs --force-reprocess
        """,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Directory containing AlphaFold3 design subfolders.",
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output CSV path (default: <input>/af3_scores.csv).",
    )
    parser.add_argument(
        "--force-reprocess", action="store_true",
        help="Reprocess all samples, ignoring existing rows in the output CSV.",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true",
        help="Reduce output verbosity.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="List the JSONs that would be processed without scoring.",
    )
    parser.add_argument(
        "-j", "--workers", type=int, default=0,
        help="Number of parallel worker processes (0 = auto-detect from SLURM/affinity).",
    )
    parser.add_argument(
        "--chunksize", type=int, default=8,
        help="Pool chunksize for task dispatch (default: 8).",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    input_dir = os.path.abspath(args.input)

    if not os.path.isdir(input_dir):
        print(f"Error: input directory does not exist: {input_dir}", file=sys.stderr)
        sys.exit(1)

    output_csv = args.output or os.path.join(input_dir, "af3_scores.csv")

    if not args.quiet:
        print(f"Input:  {input_dir}")
        print(f"Output: {output_csv}")

    design_dirs = sorted(
        os.path.join(input_dir, d) for d in os.listdir(input_dir)
        if os.path.isdir(os.path.join(input_dir, d))
    )
    if not args.quiet:
        print(f"Found {len(design_dirs)} design folders.")

    already_done = set()
    if not args.force_reprocess:
        already_done = load_already_processed(output_csv)
        if already_done and not args.quiet:
            print(f"Resuming: {len(already_done)} (folder, sample) pairs already in CSV.")

    fieldnames = [
        "design_folder", "sample_id", "avg_plddt", "pae_interaction",
        "n_tokens", "n_atoms", "json_path",
    ]
    file_missing_or_empty = (
        not os.path.exists(output_csv) or os.path.getsize(output_csv) == 0
    )
    write_header = file_missing_or_empty or args.force_reprocess
    mode = "w" if write_header else "a"

    expected_header = ",".join(fieldnames)
    if not write_header:
        try:
            with open(output_csv, "r") as f_check:
                first_line = f_check.readline().rstrip("\r\n")
            if first_line != expected_header:
                print(
                    f"Warning: existing CSV {output_csv} does not start with the expected header.\n"
                    f"  Expected: {expected_header}\n"
                    f"  Found:    {first_line}\n"
                    f"  Consider deleting it and re-running, or pass --force-reprocess.",
                    file=sys.stderr,
                )
        except IOError as e:
            print(f"Warning: could not peek at existing CSV: {e}", file=sys.stderr)

    if args.dry_run:
        total = 0
        for d in design_dirs:
            for jp, sid in find_confidences_jsons(d):
                key = (os.path.basename(d), sid)
                if not args.force_reprocess and key in already_done:
                    continue
                total += 1
                if not args.quiet:
                    print(f"  would process: {os.path.basename(d)}  {sid}")
        print(f"Dry run total: {total} JSONs")
        return

    out_dir = os.path.dirname(output_csv)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    tasks = []
    n_skipped = 0
    for d in design_dirs:
        design_name = os.path.basename(d)
        for jp, sid in find_confidences_jsons(d):
            key = (design_name, sid)
            if not args.force_reprocess and key in already_done:
                n_skipped += 1
                continue
            tasks.append((design_name, sid, jp))

    n_workers = args.workers if args.workers > 0 else detect_cpus()
    n_workers = max(1, min(n_workers, len(tasks))) if tasks else 1
    if not args.quiet:
        print(f"Tasks to process: {len(tasks)} (skipping {n_skipped} already done)")
        print(f"Using {n_workers} worker process(es).")

    n_written = 0
    n_failed = 0
    with open(output_csv, mode, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
            f.flush()

        if not tasks:
            pass
        elif n_workers == 1:
            for task in tasks:
                row = _worker(task)
                if "error" in row:
                    n_failed += 1
                    print(f"  Skip {row['json_path']}: {row['error']}", file=sys.stderr)
                    continue
                writer.writerow(row)
                n_written += 1
                if n_written % 200 == 0:
                    f.flush()
                    if not args.quiet:
                        print(f"  ...{n_written}/{len(tasks)} rows written")
        else:
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                for row in pool.map(_worker, tasks, chunksize=args.chunksize):
                    if "error" in row:
                        n_failed += 1
                        print(f"  Skip {row['json_path']}: {row['error']}", file=sys.stderr)
                        continue
                    writer.writerow(row)
                    n_written += 1
                    if n_written % 200 == 0:
                        f.flush()
                        if not args.quiet:
                            print(f"  ...{n_written}/{len(tasks)} rows written")
        f.flush()

    if not args.quiet:
        print(f"\nDone. Wrote {n_written} rows, skipped {n_skipped} already-done, {n_failed} failed.")
        print(f"CSV: {output_csv}")


if __name__ == "__main__":
    main()
