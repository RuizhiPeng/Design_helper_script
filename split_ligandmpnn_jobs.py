#!/usr/bin/env python
"""
Generate split LigandMPNN sbatch scripts from a template, with checkpointing.

Usage:
    python split_mpnn_jobs.py \
        --job_dir ./TNT_helical/4_mpnn/IL9 \
        --template ./TNT_helical/4_mpnn/IL9/IL9_ligandmpnn.sh \
        --n_splits 8 \
        [--job_name IL9]   # override base name for scripts/jobs

What it does:
1. Reads pdb_ids.json + fixed_residues_multi.json from --job_dir.
2. Checkpoint: each generated sbatch script touches {job_name}_split_i.done
   on clean exit. Re-running this script reads those sentinel files (O(K) stat
   calls, not an O(N) directory listing) and excludes already-finished PDBs.
3. Splits the remaining PDBs evenly across --n_splits jobs.
4. For each split i, writes:
     pdb_ids_split_i.json
     fixed_residues_multi_split_i.json
     {job_name}_split_i.sh   (patched copy of the template)
5. Prints ready-to-paste sbatch commands.

Checkpoint note:
  A split is considered done only when its sbatch script exits cleanly (exit
  code 0) — a killed/failed job leaves no sentinel and will be re-queued.
  On re-run the remaining structures are re-split evenly into --n_splits new
  jobs; completed splits are excluded entirely.
"""

import argparse
import json
import math
import os
import re
import stat
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_completed_pdbs(job_dir: str, job_name: str) -> set[str]:
    """Return the set of PDB paths already processed, using sentinel files.

    Each split sbatch script touches {job_name}_split_i.done on clean exit.
    We stat those files (O(K), K = number of splits) and read the corresponding
    pdb_ids_split_i.json to learn which PDBs they covered — no seqs/ scan needed.
    """
    completed: set[str] = set()
    sentinel_re = re.compile(rf"^{re.escape(job_name)}_split_(\d+)\.done$")
    try:
        entries = os.listdir(job_dir)
    except OSError:
        return completed
    for fname in entries:
        m = sentinel_re.match(fname)
        if not m:
            continue
        split_json = os.path.join(job_dir, f"pdb_ids_split_{m.group(1)}.json")
        if os.path.exists(split_json):
            with open(split_json) as f:
                completed.update(json.load(f).keys())
    return completed


def patch_template(template: str, split_idx: int, job_name: str,
                   pdb_json_path: str, fixed_json_path: str) -> str:
    """Return a copy of template with SBATCH headers, JSON paths, and sentinel updated."""
    text = template

    # SBATCH directives
    text = re.sub(
        r'(#SBATCH\s+--job-name=)\S+',
        lambda m: m.group(1) + f"{job_name}_split_{split_idx}",
        text,
    )
    text = re.sub(
        r'(#SBATCH\s+--output=)\S+',
        lambda m: m.group(1) + f"{job_name}_split_{split_idx}_out",
        text,
    )
    text = re.sub(
        r'(#SBATCH\s+--error=)\S+',
        lambda m: m.group(1) + f"{job_name}_split_{split_idx}_err",
        text,
    )

    # INPUT_PDB_JSON — double-quoted, single-quoted, or bare
    for pat, repl in [
        (r'INPUT_PDB_JSON="[^"]*"', f'INPUT_PDB_JSON="{pdb_json_path}"'),
        (r"INPUT_PDB_JSON='[^']*'", f"INPUT_PDB_JSON='{pdb_json_path}'"),
        (r'INPUT_PDB_JSON=\S+',     f'INPUT_PDB_JSON="{pdb_json_path}"'),
    ]:
        new_text = re.sub(pat, repl, text)
        if new_text != text:
            text = new_text
            break

    # FIXED_RESIDUES_JSON
    for pat, repl in [
        (r'FIXED_RESIDUES_JSON="[^"]*"', f'FIXED_RESIDUES_JSON="{fixed_json_path}"'),
        (r"FIXED_RESIDUES_JSON='[^']*'", f"FIXED_RESIDUES_JSON='{fixed_json_path}'"),
        (r'FIXED_RESIDUES_JSON=\S+',     f'FIXED_RESIDUES_JSON="{fixed_json_path}"'),
    ]:
        new_text = re.sub(pat, repl, text)
        if new_text != text:
            text = new_text
            break

    return text


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate split LigandMPNN sbatch scripts with checkpointing"
    )
    parser.add_argument(
        "--job_dir", required=True,
        help="Directory containing pdb_ids.json and fixed_residues_multi.json",
    )
    parser.add_argument(
        "--template", required=True,
        help="Template sbatch script to base each split on",
    )
    parser.add_argument(
        "--n_splits", type=int, required=True,
        help="Number of parallel job splits to create",
    )
    parser.add_argument(
        "--job_name", default=None,
        help="Base name used for job scripts and SBATCH --job-name. "
             "Defaults to the basename of --job_dir.",
    )
    args = parser.parse_args()

    job_dir = os.path.abspath(args.job_dir)
    pdb_ids_path = os.path.join(job_dir, "pdb_ids.json")
    fixed_res_path = os.path.join(job_dir, "fixed_residues_multi.json")
    template_path = os.path.abspath(args.template)

    # Validate
    for p, label in [(pdb_ids_path, "pdb_ids.json"),
                     (fixed_res_path, "fixed_residues_multi.json"),
                     (template_path, "template")]:
        if not os.path.exists(p):
            print(f"Error: {label} not found at {p}", file=sys.stderr)
            sys.exit(1)

    with open(template_path) as f:
        template_text = f.read()
    with open(pdb_ids_path) as f:
        pdb_ids: dict = json.load(f)
    with open(fixed_res_path) as f:
        fixed_residues: dict = json.load(f)

    job_name = args.job_name or os.path.basename(job_dir)

    # Checkpoint: exclude PDBs from splits whose sentinel file exists
    all_pdbs = list(pdb_ids.keys())
    completed = get_completed_pdbs(job_dir, job_name)
    pending = [p for p in all_pdbs if p not in completed]
    n_done = len(all_pdbs) - len(pending)
    if n_done:
        print(f"Checkpoint: {n_done}/{len(all_pdbs)} structures already done, "
              f"{len(pending)} remaining.")
    else:
        print(f"No sentinel files found; processing all {len(all_pdbs)} structures.")

    if not pending:
        print("All structures already processed. Nothing to do.")
        sys.exit(0)

    # Cap splits to number of pending structures
    n_splits = min(args.n_splits, len(pending))
    if n_splits < args.n_splits:
        print(
            f"Note: only {len(pending)} pending structures; "
            f"reducing to {n_splits} split(s)."
        )

    chunk_size = math.ceil(len(pending) / n_splits)
    chunks = [pending[i:i + chunk_size] for i in range(0, len(pending), chunk_size)]

    print(f"\nSplitting {len(pending)} structures into {len(chunks)} job(s):")
    for i, chunk in enumerate(chunks):
        print(f"  split_{i}: {len(chunk)} structures")

    # Write per-split JSONs and sbatch scripts
    sbatch_scripts = []
    for i, chunk in enumerate(chunks):
        pdb_split_path   = os.path.join(job_dir, f"pdb_ids_split_{i}.json")
        fixed_split_path = os.path.join(job_dir, f"fixed_residues_multi_split_{i}.json")
        sbatch_path      = os.path.join(job_dir, f"{job_name}_split_{i}.sh")
        sentinel_path    = os.path.join(job_dir, f"{job_name}_split_{i}.done")

        with open(pdb_split_path, "w") as f:
            json.dump({p: pdb_ids[p] for p in chunk}, f, indent=2)

        with open(fixed_split_path, "w") as f:
            json.dump({p: fixed_residues[p] for p in chunk}, f, indent=2)

        patched = patch_template(
            template_text, i, job_name, pdb_split_path, fixed_split_path,
        )
        with open(sbatch_path, "w") as f:
            f.write(patched)
        os.chmod(sbatch_path, os.stat(sbatch_path).st_mode | stat.S_IXUSR | stat.S_IXGRP)

        sbatch_scripts.append(sbatch_path)

    print("\nFiles written:")
    for s in sbatch_scripts:
        print(f"  {s}")

    print("\nTo submit all splits:")
    for s in sbatch_scripts:
        print(f"  sbatch {s}")


if __name__ == "__main__":
    main()
