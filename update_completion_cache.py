#!/usr/bin/env python3
"""
Script to detect .done.txt files in subdirectories and update .completion_cache.json
to mark corresponding items as completed (true).

Usage: python update_completion_cache.py  af2_output/  af2_output/.completion_cache.json
"""

import os
import json
import argparse
from pathlib import Path
from typing import Set, Dict, Any


def find_done_files(base_path: str) -> Set[str]:
    """
    Find all .done.txt files in subdirectories.

    Args:
        base_path: Root directory to search

    Returns:
        Set of directory names containing .done.txt files
    """
    done_dirs = set()
    base_path_obj = Path(base_path)

    if not base_path_obj.exists():
        print(f"Warning: Path '{base_path}' does not exist")
        return done_dirs

    # Search for all .done.txt files
    for done_file in base_path_obj.rglob("*.done.txt"):
        # Get the parent directory name
        parent_dir = done_file.parent.name
        done_dirs.add(parent_dir)
        print(f"Found: {done_file.relative_to(base_path_obj)}")

    return done_dirs


def update_completion_cache(cache_path: str, completed_items: Set[str]) -> None:
    """
    Update the .completion_cache.json file to mark items as completed.

    Args:
        cache_path: Path to the .completion_cache.json file
        completed_items: Set of item names that should be marked as true
    """
    cache_file = Path(cache_path)

    # Load existing cache or create new one
    if cache_file.exists():
        with open(cache_file, 'r', encoding='utf-8') as f:
            try:
                cache_data = json.load(f)
            except json.JSONDecodeError:
                print(f"Warning: Could not parse {cache_path}, creating new cache")
                cache_data = {}
    else:
        print(f"Creating new cache file: {cache_path}")
        cache_data = {}

    # Update cache with completed items
    updated_count = 0
    for item in completed_items:
        if item not in cache_data or cache_data[item] != True:
            cache_data[item] = True
            updated_count += 1
            print(f"Marked '{item}' as completed")

    # Write updated cache back to file
    with open(cache_file, 'w', encoding='utf-8') as f:
        json.dump(cache_data, f, indent=2, ensure_ascii=False)

    print(f"\nUpdated {updated_count} items in cache")
    print(f"Total completed items: {len([v for v in cache_data.values() if v == True])}")


def main():
    """Main function to process completion status."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Detect .done.txt files and update completion cache"
    )
    parser.add_argument(
        "search_path",
        help="Path to search for .done.txt files (e.g., af2_output)"
    )
    parser.add_argument(
        "cache_file",
        help="Path to .completion_cache.json file"
    )

    args = parser.parse_args()

    print(f"Scanning directory: {args.search_path}")
    print(f"Cache file: {args.cache_file}")
    print("=" * 60)

    # Find all directories with .done.txt files
    completed_dirs = find_done_files(args.search_path)

    if not completed_dirs:
        print("\nNo .done.txt files found")
        return

    print(f"\nFound {len(completed_dirs)} completed items")
    print("=" * 60)

    # Update the completion cache
    update_completion_cache(args.cache_file, completed_dirs)

    print("=" * 60)
    print("Done!")


if __name__ == "__main__":
    main()
