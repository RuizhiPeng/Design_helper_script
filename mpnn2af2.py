#!/usr/bin/env python3
"""
Usage: mpnn2af2.py --directory <fasta_dir> <jsonl_file> -o <output_dir>
"""

"""
Script to process FASTA and JSONL files and generate output FASTA files.
Combines sequences from FASTA (designed chains) with sequences from JSONL (fixed chains).
"""

import json
import os
import sys
from typing import Dict, List, Tuple
import argparse


def parse_fasta(fasta_file: str) -> Dict[str, Tuple[str, str]]:
    """
    Parse FASTA file and extract sequences with their headers.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence index to (header, sequence) tuple
    """
    sequences = {}
    current_header = None
    current_seq = []
    seq_index = 0
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences[seq_index] = (current_header, ''.join(current_seq))
                    seq_index += 1
                
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_header is not None:
            sequences[seq_index] = (current_header, ''.join(current_seq))
    
    return sequences


def parse_jsonl(jsonl_file: str) -> List[Dict]:
    """
    Parse JSONL file and extract structural data.
    
    Args:
        jsonl_file: Path to JSONL file
        
    Returns:
        List of dictionaries containing structural data
    """
    data = []
    with open(jsonl_file, 'r') as f:
        for line in f:
            if line.strip():
                data.append(json.loads(line.strip()))
    return data


def parse_header_chain_info(header: str) -> Tuple[List[str], List[str]]:
    """
    Parse FASTA header to extract fixed_chains and designed_chains information.

    Args:
        header: FASTA header string

    Returns:
        Tuple of (fixed_chains, designed_chains) lists
    """
    import re

    fixed_chains = []
    designed_chains = []

    # Extract fixed_chains
    fixed_match = re.search(r"fixed_chains=\[([^\]]*)\]", header)
    if fixed_match:
        chains_str = fixed_match.group(1)
        # Parse chain IDs from string like "'A', 'B'" or "'A'"
        chains = re.findall(r"'([A-Z])'", chains_str)
        fixed_chains = chains

    # Extract designed_chains
    designed_match = re.search(r"designed_chains=\[([^\]]*)\]", header)
    if designed_match:
        chains_str = designed_match.group(1)
        # Parse chain IDs from string like "'A', 'B'" or "'B'"
        chains = re.findall(r"'([A-Z])'", chains_str)
        designed_chains = chains

    return fixed_chains, designed_chains


def extract_chain_sequences(jsonl_data: Dict) -> Dict[str, str]:
    """
    Extract chain sequences from JSONL data.

    Args:
        jsonl_data: Dictionary containing structural data

    Returns:
        Dictionary mapping chain ID to sequence
    """
    chain_sequences = {}

    # Extract chain A sequence
    if 'seq_chain_A' in jsonl_data:
        chain_sequences['A'] = jsonl_data['seq_chain_A']

    # Extract chain B sequence
    if 'seq_chain_B' in jsonl_data:
        chain_sequences['B'] = jsonl_data['seq_chain_B']

    # Handle additional chains if present
    for key in jsonl_data:
        if key.startswith('seq_chain_') and key not in ['seq_chain_A', 'seq_chain_B']:
            chain_id = key.split('_')[-1]
            chain_sequences[chain_id] = jsonl_data[key]

    return chain_sequences


def get_base_name(fasta_file: str) -> str:
    """
    Extract base name from FASTA file path.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Base name without extension
    """
    return os.path.splitext(os.path.basename(fasta_file))[0]


def process_single_file(fasta_file: str, jsonl_file: str, output_dir: str = "."):
    """
    Process a single FASTA and JSONL file pair to generate output FASTA files.

    Args:
        fasta_file: Path to input FASTA file
        jsonl_file: Path to input JSONL file
        output_dir: Directory to save output files
    """
    # Parse input files
    jsonl_data_list = parse_jsonl(jsonl_file)

    # Use the first (and typically only) JSONL entry
    if not jsonl_data_list:
        print(f"Warning: No data found in {jsonl_file}")
        return

    jsonl_data = jsonl_data_list[0]  # Assuming one entry per JSONL file

    # Use the shared processing function
    process_single_file_with_data(fasta_file, jsonl_data, output_dir)


def find_fasta_files(input_dir: str) -> List[str]:
    """
    Find all FASTA files in a directory.
    
    Args:
        input_dir: Directory containing input files
        
    Returns:
        List of FASTA file paths
    """
    fasta_files = []
    
    # Get all files in directory
    files = os.listdir(input_dir)
    
    # Find FASTA files
    for file in files:
        if file.endswith('.fa') or file.endswith('.fasta'):
            fasta_path = os.path.join(input_dir, file)
            fasta_files.append(fasta_path)
    
    return fasta_files


def find_jsonl_data_by_name(jsonl_data_list: List[Dict], target_name: str) -> Dict:
    """
    Find JSONL data entry that matches the target name.
    
    Args:
        jsonl_data_list: List of all JSONL entries
        target_name: Name to search for (e.g., 'sst14_001')
        
    Returns:
        Matching JSONL data dictionary, or None if not found
    """
    for data in jsonl_data_list:
        # Check if 'name' field matches
        if data.get('name') == target_name:
            return data
        
        # Also check if the name is part of any identifier
        for key, value in data.items():
            if isinstance(value, str) and target_name in value:
                return data
    
    return None


def process_directory_with_single_jsonl(input_dir: str, jsonl_file: str, output_dir: str = "."):
    """
    Process all FASTA files in a directory using a single JSONL file.
    
    Args:
        input_dir: Directory containing FASTA files
        jsonl_file: Path to the single JSONL file containing all data
        output_dir: Directory to save output files
    """
    if not os.path.exists(input_dir):
        print(f"Error: Input directory '{input_dir}' not found.")
        return
    
    if not os.path.exists(jsonl_file):
        print(f"Error: JSONL file '{jsonl_file}' not found.")
        return
    
    # Parse the single JSONL file
    print(f"Loading JSONL data from: {jsonl_file}")
    jsonl_data_list = parse_jsonl(jsonl_file)
    print(f"Found {len(jsonl_data_list)} entries in JSONL file")
    
    # Find all FASTA files
    fasta_files = find_fasta_files(input_dir)
    
    if not fasta_files:
        print("No FASTA files found in directory.")
        return
    
    print(f"Processing {len(fasta_files)} FASTA files...")
    
    # Process each FASTA file
    for fasta_file in fasta_files:
        base_name = get_base_name(fasta_file)
        print(f"\nProcessing: {os.path.basename(fasta_file)} (looking for '{base_name}' in JSONL)")
        
        # Find corresponding JSONL data
        jsonl_data = find_jsonl_data_by_name(jsonl_data_list, base_name)
        
        if jsonl_data is None:
            print(f"Warning: No matching JSONL data found for {base_name}")
            continue
        
        try:
            process_single_file_with_data(fasta_file, jsonl_data, output_dir)
        except Exception as e:
            print(f"Error processing {fasta_file}: {e}")
            continue
    
    print(f"\nProcessing completed! Output files saved to: {output_dir}")


def process_single_file_with_data(fasta_file: str, jsonl_data: Dict, output_dir: str = "."):
    """
    Process a single FASTA file with provided JSONL data.

    Args:
        fasta_file: Path to input FASTA file
        jsonl_data: Dictionary containing structural data
        output_dir: Directory to save output files
    """
    # Parse FASTA file
    fasta_sequences = parse_fasta(fasta_file)

    # Get base name for output files
    base_name = get_base_name(fasta_file)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Extract chain sequences from JSONL data
    chain_sequences = extract_chain_sequences(jsonl_data)

    # Get the reference sequence (first sequence) to extract chain information
    if 0 not in fasta_sequences:
        print(f"Warning: No reference sequence found in {fasta_file}")
        return

    reference_header, _ = fasta_sequences[0]
    fixed_chains, designed_chains = parse_header_chain_info(reference_header)

    # Skip the first sequence (reference) and process generated sequences
    generated_sequences = {k: v for k, v in fasta_sequences.items() if k > 0}

    for seq_idx, (header, designed_seq) in generated_sequences.items():
        # Create output filename
        output_filename = f"{base_name}_v{seq_idx}.fa"
        output_path = os.path.join(output_dir, output_filename)

        with open(output_path, 'w') as f:
            # Write all chains in order, handling fixed vs designed chains
            all_chains = set(fixed_chains + designed_chains)

            for chain_id in sorted(all_chains):
                if chain_id in fixed_chains:
                    # Fixed chain: get sequence from JSONL data
                    if chain_id in chain_sequences:
                        f.write(f">{base_name}_{seq_idx}_{chain_id}\n")
                        f.write(f"{chain_sequences[chain_id]}\n")
                    else:
                        print(f"Warning: Fixed chain {chain_id} not found in JSONL data for {base_name}")

                elif chain_id in designed_chains:
                    # Designed chain: use the sequence from FASTA
                    f.write(f">{base_name}_{seq_idx}_{chain_id}\n")
                    f.write(f"{designed_seq}\n")

        print(f"Created: {output_path}")


def process_files(fasta_file: str, jsonl_file: str, output_dir: str = "."):
    """
    Process FASTA and JSONL files to generate output FASTA files.
    Wrapper function for backward compatibility.
    
    Args:
        fasta_file: Path to input FASTA file
        jsonl_file: Path to input JSONL file
        output_dir: Directory to save output files
    """
    process_single_file(fasta_file, jsonl_file, output_dir)


def main():
    """Main function to handle command line arguments and process files."""
    parser = argparse.ArgumentParser(
        description="Process FASTA and JSONL files to generate output FASTA files."
    )
    
    # Add mutually exclusive group for different processing modes
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--single", 
        nargs=2, 
        metavar=("FASTA", "JSONL"),
        help="Process a single FASTA and JSONL file pair"
    )
    group.add_argument(
        "--directory", 
        nargs=2,
        metavar=("FASTA_DIR", "JSONL_FILE"),
        help="Process all FASTA files in directory with a single JSONL file"
    )
    
    parser.add_argument(
        "-o", "--output", 
        default=".", 
        help="Output directory (default: current directory)"
    )
    
    args = parser.parse_args()
    
    try:
        if args.single:
            fasta_file, jsonl_file = args.single
            
            # Validate input files
            if not os.path.exists(fasta_file):
                print(f"Error: FASTA file '{fasta_file}' not found.")
                sys.exit(1)
            
            if not os.path.exists(jsonl_file):
                print(f"Error: JSONL file '{jsonl_file}' not found.")
                sys.exit(1)
            
            process_single_file(fasta_file, jsonl_file, args.output)
            print("Processing completed successfully!")
            
        elif args.directory:
            fasta_dir, jsonl_file = args.directory
            process_directory_with_single_jsonl(fasta_dir, jsonl_file, args.output)
            
    except Exception as e:
        print(f"Error processing files: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
