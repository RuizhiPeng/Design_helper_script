#!/usr/bin/env python3
"""
Script to sort unprocessed files from a source split directory into multiple target split directories.
Reads .completion_cache.json to identify already processed files and redistributes unprocessed ones.

Usage:
    python sort_files.py --base <base_dir> --source split3 --targets split1 split2 split3
    python sort_files.py --base <base_dir> --source split3 --targets split1 split2 split3 --dry-run
"""

import json
import os
import shutil
import random
import argparse
from pathlib import Path
from typing import Set, List


def load_completion_cache(cache_path: str) -> Set[str]:
    """
    Load the completion cache JSON file and extract processed file identifiers.

    Args:
        cache_path: Path to the .completion_cache.json file

    Returns:
        Set of processed file identifiers (without extensions)
    """
    with open(cache_path, 'r') as f:
        cache_data = json.load(f)

    # Extract the keys (file identifiers) from the cache
    # These are in format like "conf_10_cirbp_86_108_1_12_10_21_0_v1_be280"
    # We need to match them with input files which have format like "fill_3_cirbp_86_108_1_10_1_10_0_v1.fa"
    processed_files = set(cache_data.keys())

    return processed_files


def is_file_processed(filename: str, processed_set: Set[str]) -> bool:
    """
    Check if a file has been processed based on the completion cache.

    Args:
        filename: Input filename to check (e.g., "fill_3_cirbp_86_108_1_10_1_10_0_v1.fa")
        processed_set: Set of processed file identifiers from cache

    Returns:
        True if the file has been processed, False otherwise
    """
    # Remove .fa extension to get base name
    # Example: fill_3_cirbp_86_108_1_10_1_10_0_v1.fa -> fill_3_cirbp_86_108_1_10_1_10_0_v1
    base = os.path.splitext(filename)[0]

    # Check if any cache entry starts with this base name followed by underscore and hash
    # Example cache entry: fill_3_cirbp_86_108_1_10_1_10_0_v1_ed6a9
    for processed in processed_set:
        if processed.startswith(base + '_'):
            return True

    return False


def get_files_in_directory(directory: str) -> List[str]:
    """
    Get all files in a directory (excluding subdirectories).

    Args:
        directory: Path to the directory

    Returns:
        List of filenames
    """
    path = Path(directory)
    if not path.exists():
        return []

    files = [f.name for f in path.iterdir() if f.is_file() and not f.name.startswith('.') and f.name != 'desktop.ini']
    return files


def redistribute_files(unprocessed_files: List[str], source_dir: str, target_dirs: List[str], dry_run: bool = False):
    """
    Randomly and evenly redistribute unprocessed files into target directories.

    Args:
        unprocessed_files: List of filenames to redistribute
        source_dir: Source directory path
        target_dirs: List of target directory paths
        dry_run: If True, only print what would be done without actually moving files
    """
    # Shuffle files randomly
    random.shuffle(unprocessed_files)

    # Calculate how many files per directory
    num_targets = len(target_dirs)
    files_per_dir = len(unprocessed_files) // num_targets
    remainder = len(unprocessed_files) % num_targets

    print(f"\nRedistributing {len(unprocessed_files)} files into {num_targets} directories")
    print(f"Base files per directory: {files_per_dir}")
    if remainder > 0:
        print(f"First {remainder} directories will get 1 extra file")

    # Distribute files
    file_index = 0
    for dir_index, target_dir in enumerate(target_dirs):
        # First 'remainder' directories get one extra file
        num_files = files_per_dir + (1 if dir_index < remainder else 0)

        print(f"\n{os.path.basename(target_dir)}: {num_files} files")

        for i in range(num_files):
            if file_index >= len(unprocessed_files):
                break

            filename = unprocessed_files[file_index]
            source_path = os.path.join(source_dir, filename)
            target_path = os.path.join(target_dir, filename)

            if dry_run:
                print(f"  [DRY RUN] Would move: {filename}")
            else:
                shutil.move(source_path, target_path)
                if i < 3 or i == num_files - 1:  # Print first 3 and last file
                    print(f"  Moved: {filename}")
                elif i == 3:
                    print(f"  ... ({num_files - 4} more files)")

            file_index += 1


def main():
    parser = argparse.ArgumentParser(
        description='Redistribute unprocessed AF2 input files evenly across split directories.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sort_files.py --base my_project --source af2_input_split3 --targets af2_input_split1 af2_input_split2 af2_input_split3
  python sort_files.py --base my_project --source af2_input_split3 --targets af2_input_split1 af2_input_split2 af2_input_split3 --dry-run
        """
    )
    parser.add_argument('--base', required=True,
                        help='Base project directory containing af2_output and split directories')
    parser.add_argument('--source', required=True,
                        help='Source split directory name (relative to base), e.g. af2_input_split3')
    parser.add_argument('--targets', nargs='+', required=True,
                        help='Target split directory names (relative to base) to redistribute files into')
    parser.add_argument('--cache', default=None,
                        help='Path to .completion_cache.json (default: <base>/af2_output/.completion_cache.json)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without moving any files')

    args = parser.parse_args()

    BASE_DIR = args.base
    SOURCE_DIR = os.path.join(BASE_DIR, args.source)
    TARGET_DIRS = [os.path.join(BASE_DIR, t) for t in args.targets]
    CACHE_FILE = args.cache if args.cache else os.path.join(BASE_DIR, "af2_output", ".completion_cache.json")
    DRY_RUN = args.dry_run

    print("=" * 70)
    print("File Redistribution Script")
    print("=" * 70)

    # Load completion cache
    print(f"\nLoading completion cache from: {CACHE_FILE}")
    if not os.path.exists(CACHE_FILE):
        print(f"ERROR: Cache file not found: {CACHE_FILE}")
        return

    processed_files = load_completion_cache(CACHE_FILE)
    print(f"Found {len(processed_files)} processed files in cache")

    # Get all files in source directory
    print(f"\nScanning source directory: {SOURCE_DIR}")
    all_files = get_files_in_directory(SOURCE_DIR)
    print(f"Found {len(all_files)} total files in {os.path.basename(SOURCE_DIR)}")

    # Filter out processed files
    print("\nChecking which files are unprocessed...")
    unprocessed_files = []
    for filename in all_files:
        if not is_file_processed(filename, processed_files):
            unprocessed_files.append(filename)

    print(f"Found {len(unprocessed_files)} unprocessed files")
    print(f"Found {len(all_files) - len(unprocessed_files)} already processed files")

    if len(unprocessed_files) == 0:
        print("\nNo unprocessed files to redistribute. Exiting.")
        return

    # Show sample of unprocessed files
    print("\nSample unprocessed files:")
    for f in unprocessed_files[:5]:
        print(f"  - {f}")
    if len(unprocessed_files) > 5:
        print(f"  ... and {len(unprocessed_files) - 5} more")

    # Verify target directories exist
    print("\nVerifying target directories...")
    for target_dir in TARGET_DIRS:
        if not os.path.exists(target_dir):
            print(f"Creating directory: {target_dir}")
            os.makedirs(target_dir)
        else:
            print(f"  [OK] {os.path.basename(target_dir)}")

    # Redistribute files
    if DRY_RUN:
        print("\n" + "=" * 70)
        print("DRY RUN MODE - No files will be moved")
        print("=" * 70)

    redistribute_files(unprocessed_files, SOURCE_DIR, TARGET_DIRS, dry_run=DRY_RUN)

    print("\n" + "=" * 70)
    if DRY_RUN:
        print("DRY RUN COMPLETE")
        print("Run without --dry-run to actually move files")
    else:
        print("REDISTRIBUTION COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
