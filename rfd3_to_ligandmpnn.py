#!/usr/bin/env python
"""
Prepare LigandMPNN input JSONs (pdb_ids.json and fixed_residues_multi.json)
from RFdiffusion3 output.

Usage:
    python prepare_mpnn_inputs.py \
        --rfd3_output_dir /path/to/rfd3_output \
        --mpnn_output_dir /path/to/mpnn_dir

This script:
1. Reads each JSON from rfd3 output to determine which residues to fix.
2. Converts .cif.gz files to .pdb files (BioPython).
3. Writes pdb_ids.json and fixed_residues_multi.json for LigandMPNN.

Logic for determining fixed residues:
- The output structure has all residues on chain A, numbered 1..N.
- `diffused_index_map` maps original-chain residues (e.g. C1, B5) to output
  A-chain numbering (e.g. A1, A142).
- `sampled_contig` lists all residues. Entries that are just a number (no chain
  prefix) are de novo generated backbone and should be DESIGNED (not fixed).
- `select_unfixed_sequence` lists original-chain residues whose sequence should
  be redesigned by MPNN (not fixed).
- Fixed residues = all mapped original-chain residues MINUS those in
  select_unfixed_sequence. De novo residues are never fixed.
"""

import argparse
import glob
import gzip
import json
import os
import re
import sys
from multiprocessing import Pool, cpu_count

from Bio.PDB import MMCIFParser, PDBIO


def parse_range_spec(spec: str) -> list[str]:
    """Parse a range spec like 'B1-26,B28-30,C4-6' into individual residue IDs.

    Each segment is like 'B1-26' meaning chain B residues 1 through 26,
    or 'C4-6' meaning chain C residues 4 through 6.
    """
    residues = []
    for segment in spec.split(","):
        segment = segment.strip()
        if not segment:
            continue
        # Match chain letter + start number, optionally dash + end number
        m = re.match(r"([A-Z])(\d+)(?:-(\d+))?$", segment)
        if not m:
            raise ValueError(f"Cannot parse range segment: {segment}")
        chain = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3)) if m.group(3) else start
        for i in range(start, end + 1):
            residues.append(f"{chain}{i}")
    return residues


def get_fixed_residues(json_data: dict) -> str:
    """Determine fixed residues in output A-chain numbering.

    Returns a space-separated string like 'A1 A2 A3 A81 A82'.
    """
    index_map = json_data["diffused_index_map"]
    spec = json_data["specification"]
    sampled_contig = spec["extra"]["sampled_contig"]

    # Parse select_unfixed_sequence to get original-chain residues to design
    unfixed_spec = spec.get("select_unfixed_sequence", "")
    if unfixed_spec:
        unfixed_original = set(parse_range_spec(unfixed_spec))
    else:
        unfixed_original = set()

    # All mapped original-chain residues that are NOT in the unfixed set
    # should be fixed. Map them to output A-chain numbering.
    fixed = []
    for orig_res, out_res in index_map.items():
        if orig_res not in unfixed_original:
            fixed.append(out_res)

    # Sort by residue number
    fixed.sort(key=lambda x: int(re.search(r"\d+", x).group()))

    return " ".join(fixed)


def cif_gz_to_pdb(cif_gz_path: str, pdb_path: str):
    """Convert a .cif.gz file to .pdb using BioPython."""
    parser = MMCIFParser(QUIET=True)
    with gzip.open(cif_gz_path, "rt") as f:
        structure = parser.get_structure("s", f)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)


def process_one(args):
    """Process a single JSON/CIF pair. Returns (pdb_path, fixed) or None if skipped."""
    jf, rfd3_dir, pdb_dir = args
    basename = os.path.splitext(os.path.basename(jf))[0]
    cif_gz = os.path.join(rfd3_dir, basename + ".cif.gz")

    if not os.path.exists(cif_gz):
        print(f"Warning: no CIF.gz for {basename}, skipping", flush=True)
        return None

    with open(jf) as f:
        data = json.load(f)

    pdb_path = os.path.join(pdb_dir, basename + ".pdb")
    cif_gz_to_pdb(cif_gz, pdb_path)

    fixed = get_fixed_residues(data)
    print(f"  {basename}: {len(fixed.split())} fixed residues", flush=True)
    return (pdb_path, fixed)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare LigandMPNN inputs from RFdiffusion3 output"
    )
    parser.add_argument(
        "--rfd3_output_dir", required=True,
        help="Directory containing rfd3 output (.json and .cif.gz files)"
    )
    parser.add_argument(
        "--mpnn_output_dir", required=True,
        help="Directory where PDBs and JSON inputs will be written"
    )
    parser.add_argument(
        "--workers", type=int, default=cpu_count(),
        help="Number of parallel worker processes (default: all CPUs)"
    )
    args = parser.parse_args()

    rfd3_dir = args.rfd3_output_dir
    mpnn_dir = args.mpnn_output_dir
    pdb_dir = os.path.join(mpnn_dir, "pdbs")
    os.makedirs(pdb_dir, exist_ok=True)

    json_files = sorted(glob.glob(os.path.join(rfd3_dir, "*.json")))
    if not json_files:
        print(f"No JSON files found in {rfd3_dir}", file=sys.stderr)
        sys.exit(1)

    tasks = [(jf, rfd3_dir, pdb_dir) for jf in json_files]

    pdb_ids = {}
    fixed_residues_multi = {}
    skipped = 0

    with Pool(processes=args.workers) as pool:
        for result in pool.imap_unordered(process_one, tasks):
            if result is None:
                skipped += 1
            else:
                pdb_path, fixed = result
                pdb_ids[pdb_path] = ""
                fixed_residues_multi[pdb_path] = fixed

    # Write output JSONs
    pdb_ids_path = os.path.join(mpnn_dir, "pdb_ids.json")
    fixed_res_path = os.path.join(mpnn_dir, "fixed_residues_multi.json")

    with open(pdb_ids_path, "w") as f:
        json.dump(pdb_ids, f, indent=2)

    with open(fixed_res_path, "w") as f:
        json.dump(fixed_residues_multi, f, indent=2)

    print(f"\nDone! Processed {len(pdb_ids)} structures ({skipped} skipped)")
    print(f"  PDBs:            {pdb_dir}/")
    print(f"  pdb_ids.json:    {pdb_ids_path}")
    print(f"  fixed_residues:  {fixed_res_path}")


if __name__ == "__main__":
    main()
