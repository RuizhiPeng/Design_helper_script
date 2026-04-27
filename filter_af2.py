#!/usr/bin/env python3
"""
Script to copy PDB files based on CSV criteria.

Filters directories where:
- plddt > 90
- iptm >= 0.8

Copies only .pdb files containing 'rank_001' in filename.
"""

import csv
import shutil
from pathlib import Path
import argparse
import sys


def read_csv_and_filter(csv_path, plddt_threshold=90, iptm_threshold=0.8):
    """
    Read CSV file and return directories that meet the criteria.

    Args:
        csv_path: Path to the CSV file
        plddt_threshold: Minimum plddt value (exclusive)
        iptm_threshold: Minimum iptm value (inclusive)

    Returns:
        List of directory names that pass the criteria
    """
    passing_dirs = []

    try:
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)

            # Check if required columns exist
            if not reader.fieldnames:
                print("Error: CSV file appears to be empty")
                return []

            # Check for both 'plddt' and 'avg_plddt' column names
            plddt_col = None
            if 'avg_plddt' in reader.fieldnames:
                plddt_col = 'avg_plddt'
            elif 'plddt' in reader.fieldnames:
                plddt_col = 'plddt'

            if not plddt_col or 'iptm' not in reader.fieldnames:
                print(f"Error: CSV must contain 'iptm' and either 'plddt' or 'avg_plddt' columns")
                print(f"Found columns: {reader.fieldnames}")
                return []

            # Get the first column name (directory name)
            dir_col = reader.fieldnames[0]

            for row_num, row in enumerate(reader, start=2):
                try:
                    dir_name = row[dir_col].strip()
                    plddt = float(row[plddt_col])
                    iptm = float(row['iptm'])

                    if plddt > plddt_threshold and iptm >= iptm_threshold:
                        passing_dirs.append(dir_name)
                        print(f"[PASS] {dir_name}: plddt={plddt}, iptm={iptm}")
                    else:
                        print(f"[SKIP] {dir_name}: plddt={plddt}, iptm={iptm} (filtered out)")

                except (ValueError, KeyError) as e:
                    print(f"Warning: Skipping row {row_num} due to error: {e}")
                    continue

    except FileNotFoundError:
        print(f"Error: CSV file not found: {csv_path}")
        return []
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return []

    return passing_dirs


def copy_rank_001_pdb_files(source_base, dir_names, output_path, prefer_v1=False):
    """
    Copy rank_001 PDB files from specified directories to output path.

    Args:
        source_base: Base directory containing the subdirectories
        dir_names: List of directory names to search
        output_path: Destination path for copied files
        prefer_v1: If True, when directory names with same prefix exist (differing only in _v1 vs _v2),
                   only process the _v1 directories
    """
    source_base = Path(source_base)
    output_path = Path(output_path)

    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)

    # If prefer_v1 is enabled, filter out v2 directories when v1 exists
    if prefer_v1:
        filtered_dirs = []
        v1_prefixes = set()

        # First pass: collect all v1 prefixes
        for dir_name in dir_names:
            if '_v1' in dir_name:
                # Extract prefix up to and including _v1
                prefix = dir_name.split('_v1')[0]
                v1_prefixes.add(prefix)

        # Second pass: only keep v1 directories, and v2 only if no v1 exists
        for dir_name in dir_names:
            if '_v2' in dir_name:
                prefix = dir_name.split('_v2')[0]
                if prefix in v1_prefixes:
                    print(f"[SKIP] {dir_name} (corresponding v1 version exists)")
                    continue
            filtered_dirs.append(dir_name)

        dir_names = filtered_dirs
        print(f"\nAfter v1 preference filtering: {len(dir_names)} directories to process")
        print(f"{'='*60}")

    copied_count = 0
    not_found_dirs = []

    for dir_name in dir_names:
        dir_path = source_base / dir_name

        if not dir_path.exists():
            print(f"Warning: Directory not found: {dir_path}")
            not_found_dirs.append(dir_name)
            continue

        if not dir_path.is_dir():
            print(f"Warning: Not a directory: {dir_path}")
            continue

        # Find all .pdb files with 'rank_001' in the filename
        pdb_files = list(dir_path.glob("*rank_001*.pdb"))

        if not pdb_files:
            print(f"  No rank_001 PDB files found in {dir_name}")
            continue

        for pdb_file in pdb_files:
            try:
                # Keep the same filename
                output_filename = pdb_file.name
                dest_path = output_path / output_filename

                shutil.copy2(pdb_file, dest_path)
                print(f"  Copied: {pdb_file.name}")
                copied_count += 1

            except Exception as e:
                print(f"  Error copying {pdb_file}: {e}")

    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Total files copied: {copied_count}")
    print(f"  Directories processed: {len(dir_names)}")
    if not_found_dirs:
        print(f"  Directories not found: {len(not_found_dirs)}")
        print(f"  Missing: {', '.join(not_found_dirs[:5])}")
        if len(not_found_dirs) > 5:
            print(f"           ... and {len(not_found_dirs) - 5} more")
    print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(
        description='Copy PDB files based on CSV criteria',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python copy_pdb_files.py data.csv /path/to/structures /path/to/output
  python copy_pdb_files.py data.csv . ./selected_pdbs --plddt 85 --iptm 0.75
        """
    )

    parser.add_argument('csv_file', help='Path to CSV file with criteria')
    parser.add_argument('source_base', help='Base directory containing structure subdirectories')
    parser.add_argument('output_path', help='Output directory for selected PDB files')
    parser.add_argument('--plddt', type=float, default=90,
                        help='Minimum pLDDT threshold (exclusive, default: 90)')
    parser.add_argument('--iptm', type=float, default=0.8,
                        help='Minimum ipTM threshold (inclusive, default: 0.8)')
    parser.add_argument('--prefer-v1', action='store_true',
                        help='Only process v1 directories when both v1 and v2 exist with same prefix')

    args = parser.parse_args()

    print(f"Reading CSV file: {args.csv_file}")
    print(f"Criteria: pLDDT > {args.plddt}, ipTM >= {args.iptm}")
    print(f"{'='*60}")

    # Read CSV and filter directories
    passing_dirs = read_csv_and_filter(args.csv_file, args.plddt, args.iptm)

    if not passing_dirs:
        print("\nNo directories passed the filtering criteria.")
        sys.exit(1)

    print(f"\n{len(passing_dirs)} directories passed the criteria")
    print(f"{'='*60}")

    # Copy PDB files
    print(f"\nCopying rank_001 PDB files to: {args.output_path}")
    if args.prefer_v1:
        print("Preference: Only v1 directories (when both v1 and v2 exist)")
    print(f"{'='*60}")
    copy_rank_001_pdb_files(args.source_base, passing_dirs, args.output_path, prefer_v1=args.prefer_v1)


if __name__ == "__main__":
    main()
