#!/usr/bin/env python3
"""Convert LigandMPNN FASTA output to AlphaFold3 JSON input files."""

import argparse
import glob
import json
import os
import string
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert LigandMPNN FASTA output to AlphaFold3 JSON input files."
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing .fa files from LigandMPNN",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Base path for output directories (splits appended as _0, _1, ...)",
    )
    parser.add_argument(
        "--num_splits",
        type=int,
        default=1,
        help="Number of output subdirectories for batch processing (default: 1)",
    )
    parser.add_argument(
        "--model_seeds",
        type=int,
        nargs="+",
        default=[1],
        help="Model seeds for AF3 JSON (default: 1)",
    )
    return parser.parse_args()


def parse_header(header_line):
    """Parse '>name, key=val, key=val, ...' into a dict."""
    line = header_line.lstrip(">").strip()
    tokens = line.split(", ")
    metadata = {"name": tokens[0]}
    for token in tokens[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            metadata[key] = value
    return metadata


def parse_fasta(filepath):
    """Parse .fa file, return list of (metadata_dict, sequence_str) for designed sequences only."""
    entries = []
    current_header = None
    current_seq_parts = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    entries.append((current_header, "".join(current_seq_parts)))
                current_header = parse_header(line)
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

    if current_header is not None:
        entries.append((current_header, "".join(current_seq_parts)))

    # Only return designed sequences (those with 'id' key in header)
    return [(meta, seq) for meta, seq in entries if "id" in meta]


def sequence_to_chains(sequence):
    """Split multi-chain sequence by '/', return list of chain sequences."""
    if "/" in sequence:
        return [c for c in sequence.split("/") if c]
    return [sequence]


def build_af3_json(name, chains, model_seeds):
    """Build an AlphaFold3-format JSON dict."""
    return {
        "dialect": "alphafold3",
        "version": 1,
        "name": name,
        "modelSeeds": model_seeds,
        "sequences": [
            {"protein": {"id": chain_id, "sequence": chain_seq}}
            for chain_id, chain_seq in zip(string.ascii_uppercase, chains)
        ],
    }


def process_file(fa_path, model_seeds):
    """Process one .fa file, return list of (json_filename, json_dict)."""
    stem = os.path.splitext(os.path.basename(fa_path))[0]
    results = []

    designed_sequences = parse_fasta(fa_path)
    if not designed_sequences:
        print(f"  Warning: no designed sequences in {fa_path}", file=sys.stderr)
        return results

    for meta, seq in designed_sequences:
        seq_id = meta.get("id", "unknown")
        json_name = f"{stem}_id{seq_id}"
        json_filename = f"{json_name}.json"

        chains = sequence_to_chains(seq)
        json_dict = build_af3_json(json_name, chains, model_seeds)
        results.append((json_filename, json_dict))

    return results


def main():
    args = parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: input directory does not exist: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    # Create output directories
    split_dirs = []
    for i in range(args.num_splits):
        d = f"{args.output_dir}_{i}"
        os.makedirs(d, exist_ok=True)
        split_dirs.append(d)

    # Collect and sort .fa files
    fa_files = sorted(glob.glob(os.path.join(args.input_dir, "*.fa")))
    if not fa_files:
        print(f"Warning: no .fa files found in {args.input_dir}", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(fa_files)} .fa files in {args.input_dir}")

    # Process and distribute
    total_jsons = 0
    split_counts = [0] * args.num_splits

    for fa_path in fa_files:
        results = process_file(fa_path, args.model_seeds)
        for json_filename, json_dict in results:
            split_idx = total_jsons % args.num_splits
            out_path = os.path.join(split_dirs[split_idx], json_filename)
            with open(out_path, "w") as f:
                json.dump(json_dict, f, indent=2)
                f.write("\n")
            total_jsons += 1
            split_counts[split_idx] += 1

    print(f"\nDone. Wrote {total_jsons} JSON files across {args.num_splits} split(s):")
    for i, count in enumerate(split_counts):
        print(f"  {split_dirs[i]}: {count} files")


if __name__ == "__main__":
    main()
