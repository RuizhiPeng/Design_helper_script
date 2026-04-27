#!/usr/bin/env python3
"""
Renumber PDB residues to avoid duplicates across chains.
Chain B (and subsequent chains) residues will be renumbered to follow chain A.

Usage:
    python renumber_pdb_chains.py input.pdb output.pdb
    python renumber_pdb_chains.py input.pdb  # overwrites input.pdb
"""

import sys
import argparse
from collections import defaultdict


def get_chain_info(pdb_file):
    """
    Extract chain information from PDB file.
    Returns a dict mapping chain ID to list of residue numbers.
    """
    chain_residues = defaultdict(set)

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                res_num = int(line[22:26].strip())
                chain_residues[chain_id].add(res_num)

    # Convert sets to sorted lists
    for chain in chain_residues:
        chain_residues[chain] = sorted(chain_residues[chain])

    return chain_residues


def calculate_offsets(chain_info):
    """
    Calculate offset for each chain so residues are sequential.
    Returns dict mapping chain ID to offset value.
    """
    offsets = {}
    next_residue = 1  # Start numbering from 1

    # Sort chains alphabetically to ensure consistent processing
    sorted_chains = sorted(chain_info.keys())

    for chain_id in sorted_chains:
        min_res = min(chain_info[chain_id])
        max_res = max(chain_info[chain_id])

        # Offset to make this chain start at next_residue
        offsets[chain_id] = next_residue - min_res

        # Calculate where the next chain should start
        # Number of residues in current chain
        num_residues = max_res - min_res + 1
        next_residue = next_residue + num_residues

    return offsets


def renumber_pdb(input_file, output_file, offsets):
    """
    Renumber PDB file residues based on calculated offsets.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract current residue number and chain ID
                chain_id = line[21]
                res_num = int(line[22:26].strip())

                # Calculate new residue number
                new_res_num = res_num + offsets.get(chain_id, 0)

                # Reconstruct line with new residue number
                # PDB format: residue number is in columns 23-26 (right-justified)
                new_line = line[:22] + f"{new_res_num:4d}" + line[26:]
                outfile.write(new_line)
            else:
                # Write non-ATOM/HETATM lines as-is
                outfile.write(line)


def main():
    parser = argparse.ArgumentParser(
        description='Renumber PDB residues to avoid duplicates across chains.',
        epilog='Example: python renumber_pdb_chains.py input.pdb output.pdb'
    )
    parser.add_argument('input_pdb', help='Input PDB file')
    parser.add_argument('output_pdb', nargs='?', help='Output PDB file (default: overwrites input)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print detailed information')

    args = parser.parse_args()

    input_file = args.input_pdb
    output_file = args.output_pdb if args.output_pdb else input_file

    # Get chain information
    chain_info = get_chain_info(input_file)

    if args.verbose:
        print(f"Processing {input_file}...")
        print(f"Found {len(chain_info)} chain(s):")
        for chain_id in sorted(chain_info.keys()):
            res_list = chain_info[chain_id]
            print(f"  Chain {chain_id}: {len(res_list)} residues (range: {min(res_list)}-{max(res_list)})")

    # Calculate offsets
    offsets = calculate_offsets(chain_info)

    if args.verbose:
        print("\nRenumbering scheme:")
        for chain_id in sorted(offsets.keys()):
            res_list = chain_info[chain_id]
            old_range = f"{min(res_list)}-{max(res_list)}"
            new_min = min(res_list) + offsets[chain_id]
            new_max = max(res_list) + offsets[chain_id]
            new_range = f"{new_min}-{new_max}"
            print(f"  Chain {chain_id}: {old_range} -> {new_range} (offset: +{offsets[chain_id]})")

    # If overwriting, use temporary file
    if input_file == output_file:
        temp_file = output_file + '.tmp'
        renumber_pdb(input_file, temp_file, offsets)

        # Replace original with renumbered version
        import os
        os.replace(temp_file, output_file)
    else:
        renumber_pdb(input_file, output_file, offsets)

    if args.verbose:
        print(f"\nRenumbered PDB written to: {output_file}")
    else:
        print(f"Successfully renumbered {output_file}")


if __name__ == "__main__":
    main()
