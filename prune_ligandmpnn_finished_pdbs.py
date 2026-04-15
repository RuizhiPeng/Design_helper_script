#!/usr/bin/env python3
"""
Pre-step for LigandMPNN: move already-finished PDB entries out of the
input pdb_ids JSON so that run.py only processes what remains.

Checks for existing .fa files in <out_folder>/seqs/ and moves matching
entries from pdb_ids.json into a companion *_done.json file.

Usage (inside sbatch, before run.py):
    python prune_ligandmpnn_finished_pdbs.py \
        --pdb_path_multi  /path/to/pdb_ids.json \
        --out_folder      /path/to/mpnn_output \
        [--file_ending ""]
"""

import argparse
import json
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_path_multi", required=True)
    parser.add_argument("--out_folder", required=True)
    parser.add_argument("--file_ending", default="")
    args = parser.parse_args()

    seqs_dir = os.path.join(args.out_folder, "seqs")
    pdb_json = args.pdb_path_multi
    pdb_done_path = pdb_json.replace(".json", "_done.json")

    with open(pdb_json) as f:
        pdb_ids = json.load(f)

    if os.path.exists(pdb_done_path):
        with open(pdb_done_path) as f:
            pdb_done = json.load(f)
    else:
        pdb_done = {}

    # Find finished PDBs by checking for .fa output
    newly_done = []
    for pdb_path in list(pdb_ids.keys()):
        name = pdb_path[pdb_path.rfind("/") + 1:]
        if name.endswith(".pdb"):
            name = name[:-4]
        fa_path = os.path.join(seqs_dir, name + args.file_ending + ".fa")
        if os.path.exists(fa_path):
            newly_done.append(pdb_path)

    if not newly_done:
        print(f"[prune] 0 newly finished, {len(pdb_ids)} remaining")
        return

    for p in newly_done:
        pdb_done[p] = pdb_ids.pop(p)

    # Write back (tmp + rename for safety)
    for path, data in [(pdb_json, pdb_ids), (pdb_done_path, pdb_done)]:
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(data, f, indent=2)
        os.replace(tmp, path)

    print(f"[prune] moved {len(newly_done)} finished entries "
          f"({len(pdb_ids)} remaining, {len(pdb_done)} total done)")


if __name__ == "__main__":
    main()
