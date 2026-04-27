#!/usr/bin/env python3
"""
AlphaFold2 Score Processor

This script processes AlphaFold2 result folders to extract and analyze scores from JSON files.
It finds the model with the highest average pLDDT for each result folder and outputs
the metrics to a CSV file.

Usage:
    python AF2_analyze.py --input /path/to/af2_results --output results.csv

"""

import os
import json
import csv
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import statistics


def calculate_average_plddt(plddt_scores: List[float]) -> float:
    """Calculate average pLDDT from a list of scores."""
    return statistics.mean(plddt_scores)


def parse_score_file(file_path: str) -> Optional[Dict]:
    """
    Parse a JSON score file and extract metrics.

    Returns:
        Dictionary with plddt, ptm, iptm values, or None if parsing fails
    """
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Extract metrics
        plddt_scores = data.get('plddt', [])
        ptm = data.get('ptm')
        iptm = data.get('iptm')

        if not plddt_scores:
            print(f"Warning: No pLDDT scores found in {file_path}")
            return None

        avg_plddt = calculate_average_plddt(plddt_scores)

        return {
            'avg_plddt': avg_plddt,
            'ptm': ptm,
            'iptm': iptm,
            'file_path': file_path
        }

    except (json.JSONDecodeError, IOError) as e:
        print(f"Error parsing {file_path}: {e}")
        return None


def find_score_files(result_folder: str) -> List[str]:
    """
    Find all score files in a result folder.

    Returns:
        List of score file paths
    """
    score_files = []

    try:
        for file in os.listdir(result_folder):
            if file.endswith('.json') and 'scores_rank_' in file and 'alphafold2_multimer_v3_model_' in file:
                score_files.append(os.path.join(result_folder, file))
    except OSError as e:
        print(f"Error accessing folder {result_folder}: {e}")

    return score_files


def process_result_folder(result_folder: str) -> Optional[Dict]:
    """
    Process a single result folder and return the best model's metrics.

    Returns:
        Dictionary with best model's metrics, or None if no valid scores found
    """
    score_files = find_score_files(result_folder)

    if not score_files:
        return None

    best_model = None
    best_avg_plddt = -1

    for score_file in score_files:
        metrics = parse_score_file(score_file)
        if metrics and metrics['avg_plddt'] > best_avg_plddt:
            best_avg_plddt = metrics['avg_plddt']
            best_model = metrics

    if best_model:
        # Add folder name
        folder_name = os.path.basename(result_folder)
        best_model['folder_name'] = folder_name

        return best_model

    return None


def load_processed_folders(csv_file: str) -> set:
    """
    Load already processed folders from existing CSV file.

    Returns:
        Set of folder names that have been processed
    """
    processed = set()

    if os.path.exists(csv_file):
        try:
            with open(csv_file, 'r', newline='') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    processed.add(row['result_folder'])
        except (IOError, csv.Error) as e:
            print(f"Warning: Could not read existing CSV file {csv_file}: {e}")

    return processed


def write_results_to_csv(results: List[Dict], csv_file: str, append_mode: bool = False):
    """
    Write results to CSV file.

    Args:
        results: List of result dictionaries
        csv_file: Output CSV file path
        append_mode: If True, append to existing file
    """
    if not results:
        print("No results to write.")
        return

    # Ensure output directory exists
    output_dir = os.path.dirname(csv_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    mode = 'a' if append_mode and os.path.exists(csv_file) else 'w'
    write_header = mode == 'w' or not os.path.exists(csv_file)

    try:
        with open(csv_file, mode, newline='') as f:
            fieldnames = ['result_folder', 'avg_plddt', 'ptm', 'iptm']
            writer = csv.DictWriter(f, fieldnames=fieldnames)

            if write_header:
                writer.writeheader()

            for result in results:
                writer.writerow({
                    'result_folder': result['folder_name'],
                    'avg_plddt': round(result['avg_plddt'], 2),
                    'ptm': result['ptm'],
                    'iptm': result['iptm']
                })

        print(f"Results written to {csv_file}")

    except IOError as e:
        print(f"Error writing to CSV file {csv_file}: {e}")


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Process AlphaFold2 result folders to extract and analyze scores.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process default directory and output to default file
  python %(prog)s

  # Specify custom input directory
  python %(prog)s --input /path/to/alphafold/results

  # Specify both input and output
  python %(prog)s --input /path/to/af2_results --output /path/to/results.csv

  # Use short options
  python %(prog)s -i /path/to/af2_results -o results.csv

  # Force reprocessing of all folders (ignore existing CSV)
  python %(prog)s --input /path/to/af2_results --force-reprocess
        """
    )

    parser.add_argument(
        '-i', '--input',
        default='af2_output',
        help='Directory containing AlphaFold2 result folders (default: af2_output)'
    )

    parser.add_argument(
        '-o', '--output',
        default='af2_results_summary.csv',
        help='Output CSV file path (default: af2_results_summary.csv)'
    )

    parser.add_argument(
        '--force-reprocess',
        action='store_true',
        help='Reprocess all folders, ignoring existing CSV file'
    )

    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Reduce output verbosity'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be processed without actually processing'
    )

    return parser.parse_args()


def main():
    """Main function to process AlphaFold2 results."""

    # Parse command line arguments
    args = parse_arguments()

    input_dir = args.input
    output_csv = args.output

    if not os.path.exists(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.")
        sys.exit(1)

    if not args.quiet:
        print(f"Processing AlphaFold2 results from: {input_dir}")
        print(f"Output CSV file: {output_csv}")

    # Load already processed folders (unless force reprocessing)
    processed_folders = set()
    if not args.force_reprocess:
        processed_folders = load_processed_folders(output_csv)
        if processed_folders and not args.quiet:
            print(f"Found {len(processed_folders)} already processed folders. Skipping them.")

    # Get all result folders
    try:
        all_folders = [f for f in os.listdir(input_dir)
                      if os.path.isdir(os.path.join(input_dir, f))]
    except OSError as e:
        print(f"Error accessing input directory {input_dir}: {e}")
        sys.exit(1)

    # Filter out already processed folders (unless force reprocessing)
    if args.force_reprocess:
        folders_to_process = all_folders
    else:
        folders_to_process = [f for f in all_folders if f not in processed_folders]

    if not folders_to_process:
        if not args.quiet:
            print("No new folders to process.")
        return

    if not args.quiet:
        print(f"Found {len(folders_to_process)} {'total' if args.force_reprocess else 'new'} result folders to process.")

    # Dry run mode - just show what would be processed
    if args.dry_run:
        print("\nDry run mode - showing folders that would be processed:")
        for i, folder_name in enumerate(sorted(folders_to_process), 1):
            folder_path = os.path.join(input_dir, folder_name)
            score_files = find_score_files(folder_path)
            status = f"OK ({len(score_files)} score files)" if score_files else "SKIP (no score files)"
            print(f"  {i:3d}. {folder_name} - {status}")
        print(f"\nTotal folders to process: {len(folders_to_process)}")
        print("Run without --dry-run to actually process these folders.")
        return

    # Process folders
    results = []
    processed_count = 0
    skipped_count = 0

    for folder_name in sorted(folders_to_process):
        folder_path = os.path.join(input_dir, folder_name)

        if not args.quiet:
            print(f"Processing {folder_name}...", end=" ")

        result = process_result_folder(folder_path)
        if result:
            results.append(result)
            processed_count += 1
            if not args.quiet:
                print(f"OK (avg_pLDDT: {result['avg_plddt']:.2f})")
        else:
            skipped_count += 1
            if not args.quiet:
                print("SKIP (no valid score files)")

    # Write results to CSV
    if results:
        append_mode = len(processed_folders) > 0 and not args.force_reprocess
        write_results_to_csv(results, output_csv, append_mode)

        if not args.quiet:
            print(f"\nSummary:")
            print(f"  Processed: {processed_count} folders")
            print(f"  Skipped: {skipped_count} folders (no score files)")
            print(f"  Results written to: {output_csv}")

            if args.force_reprocess and append_mode:
                print(f"  Note: Results were appended to existing file.")
    else:
        if not args.quiet:
            print("\nNo valid results found.")


if __name__ == "__main__":
    main()