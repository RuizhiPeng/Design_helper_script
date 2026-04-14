#!/usr/bin/env python
"""
Standalone script to thread sequences from a FASTA file onto PDB structures.
Replaces amino acid sequences while preserving backbone coordinates.

Usage:
    Single file mode:
        python thread_sequence_to_pdb.py --pdb input.pdb --fasta sequences.fasta --out output.pdb

    Batch mode (multiple files):
        python thread_sequence_to_pdb.py --pdb_dir ./pdbs --fasta_dir ./fastas --out_dir ./output

    In batch mode, PDB files have the shortest prefix (e.g., "protein.pdb")
    FASTA files have suffixes like "_v1", "_v2" (e.g., "protein_v1.fasta", "protein_v2.fasta")
    Output PDBs will be named after the FASTA files (e.g., "protein_v1.pdb", "protein_v2.pdb")
"""

import os
import sys
import argparse
import pyrosetta
from pyrosetta.rosetta.std import ostringstream
from pathlib import Path

# Initialize PyRosetta
PYROSETTA_OPTIONS = "-mute all -beta_nov16 -in:file:silent_struct_type binary" \
    " -use_terminal_residues true -precompute_ig" \
    " -dunbrack_prob_buried 0.8 -dunbrack_prob_nonburied 0.8" \
    " -dunbrack_prob_buried_semi 0.8 -dunbrack_prob_nonburied_semi 0.8" \
    " -optimization:default_max_cycles 200"

pyrosetta.init(PYROSETTA_OPTIONS)

# Amino acid conversion dictionaries
alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','GAP']

aa_1_3 = {a:b for a,b in zip(alpha_1,alpha_3)}
aa_3_1 = {b:a for a,b in zip(alpha_1,alpha_3)}

def read_fasta(fasta_file):
    """
    Read sequences from a FASTA file.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        List of tuples (header, sequence)
    """
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))

                # Start new sequence
                current_header = line[1:].strip()  # Remove '>' and whitespace
                current_seq = []
            else:
                # Add to current sequence
                current_seq.append(line)

        # Don't forget the last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences

def thread_mpnn_seq(pose, binder_seq):
    """
    Thread a sequence onto a pose, replacing residue identities while preserving backbone.

    Args:
        pose: PyRosetta pose object
        binder_seq: Sequence string (1-letter amino acid codes)

    Returns:
        Modified pose with threaded sequence
    """
    rsd_set = pose.residue_type_set_for_pose(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)

    for resi, mut_to in enumerate(binder_seq):
        resi += 1  # PyRosetta uses 1-based indexing

        # Skip non-standard amino acids
        if mut_to not in aa_1_3:
            print(f"Warning: Skipping non-standard amino acid '{mut_to}' at position {resi}")
            continue

        name3 = aa_1_3[mut_to]
        new_res = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(
            rsd_set.name_map(name3)
        )
        pose.replace_residue(resi, new_res, True)  # True = orient backbone

    return pose

def remove_clashes(pose):
    """
    Remove major clashes using simple repacking.

    Args:
        pose: PyRosetta pose object

    Returns:
        Pose with clashes removed
    """
    xml = """<ROSETTASCRIPTS>
        <SCOREFXNS>
            <ScoreFunction name="sfxn" weights="beta_nov16" />
        </SCOREFXNS>
        <RESIDUE_SELECTORS>
            <Chain name="chainA" chains="1"/>
            <Chain name="chainB" chains="2"/>
            <Neighborhood name="interface_chA" selector="chainB" distance="10.0" />
            <Neighborhood name="interface_chB" selector="chainA" distance="10.0" />
            <And name="AB_interface" selectors="interface_chA,interface_chB" />
            <Not name="Not_interface" selector="AB_interface" />
            <And name="actual_interface_chB" selectors="AB_interface,chainB" />
            <And name="not_interface_chB" selectors="Not_interface,chainB" />
            <True name="all" />
        </RESIDUE_SELECTORS>
        <TASKOPERATIONS>
            <IncludeCurrent name="current" />
            <LimitAromaChi2 name="limitchi2" chi2max="110" chi2min="70" include_trp="True" />
            <ExtraRotamersGeneric name="ex1_ex2" ex1="1" ex2="1" />
            <OperateOnResidueSubset name="restrict2repacking" selector="all">
                <RestrictToRepackingRLT/>
            </OperateOnResidueSubset>
            <OperateOnResidueSubset name="restrict_target_not_interface" selector="not_interface_chB">
                <PreventRepackingRLT/>
            </OperateOnResidueSubset>
        </TASKOPERATIONS>
        <MOVERS>
            <PackRotamersMover name="remove_massive_clashes" scorefxn="sfxn"
                task_operations="current,restrict2repacking,restrict_target_not_interface"/>
        </MOVERS>
    </ROSETTASCRIPTS>
    """

    objs = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(xml)
    remove_massive_clashes = objs.get_mover('remove_massive_clashes')
    remove_massive_clashes.apply(pose)

    return pose

def validate_sequence_length(pose, sequence, chain_num=1):
    """
    Validate that the sequence length matches the chain length in the pose.

    Args:
        pose: PyRosetta pose object
        sequence: Sequence string
        chain_num: Chain number to validate against (1-indexed)

    Returns:
        True if valid, False otherwise
    """
    chain_length = pose.chain_end(chain_num) - pose.chain_begin(chain_num) + 1
    seq_length = len(sequence)

    if chain_length != seq_length:
        print(f"ERROR: Chain {chain_num} has {chain_length} residues but sequence has {seq_length} residues")
        return False

    return True

def find_matching_files(pdb_dir, fasta_dir):
    """
    Find and match PDB and FASTA files based on filename prefixes.
    PDB files have the shortest prefix, FASTA files have suffixes like _v1, _v2.

    Args:
        pdb_dir: Directory containing PDB files
        fasta_dir: Directory containing FASTA files

    Returns:
        List of tuples (pdb_file, fasta_file, output_name)
    """
    matches = []

    # Get all PDB files
    pdb_files = list(Path(pdb_dir).glob('*.pdb'))

    # Get all FASTA files (support common extensions)
    fasta_files = []
    for ext in ['*.fasta', '*.fa', '*.faa']:
        fasta_files.extend(Path(fasta_dir).glob(ext))

    if not pdb_files:
        print(f"WARNING: No PDB files found in {pdb_dir}")
        return matches

    if not fasta_files:
        print(f"WARNING: No FASTA files found in {fasta_dir}")
        return matches

    print(f"Found {len(pdb_files)} PDB file(s) and {len(fasta_files)} FASTA file(s)")

    # For each FASTA file, find matching PDB
    for fasta_file in fasta_files:
        fasta_stem = fasta_file.stem  # Filename without extension

        # Try to find the PDB file by removing common suffixes (_v1, _v2, etc.)
        # Look for patterns like: protein_v1 -> protein
        matched_pdb = None

        for pdb_file in pdb_files:
            pdb_stem = pdb_file.stem

            # Check if FASTA filename starts with PDB filename
            # This handles: protein.pdb matches protein_v1.fasta, protein_v2.fasta, etc.
            if fasta_stem.startswith(pdb_stem):
                # Additional check: ensure it's followed by underscore or end of string
                # to avoid partial matches (e.g., "prot" shouldn't match "protein")
                if len(fasta_stem) == len(pdb_stem) or fasta_stem[len(pdb_stem)] == '_':
                    matched_pdb = pdb_file
                    break

        if matched_pdb:
            # Output will have the same name as the FASTA file but with .pdb extension
            output_name = fasta_stem + '.pdb'
            matches.append((str(matched_pdb), str(fasta_file), output_name))
            print(f"  Matched: {matched_pdb.name} + {fasta_file.name} -> {output_name}")
        else:
            print(f"  WARNING: No matching PDB found for {fasta_file.name}")

    return matches

def process_single_pair(pdb_file, fasta_file, output_file, args):
    """
    Process a single PDB-FASTA pair and create output PDB.

    If FASTA contains multiple sequence entries, they are mapped to sequential chains.
    If --split_by_slash is used, sequences are split by '/' and mapped to chains.

    Args:
        pdb_file: Path to input PDB file
        fasta_file: Path to input FASTA file
        output_file: Path to output PDB file
        args: Argument namespace with processing options

    Returns:
        True if successful, False otherwise
    """
    # Read sequences from FASTA
    print(f"\nReading sequences from {fasta_file}")
    sequences = read_fasta(fasta_file)

    if len(sequences) == 0:
        print(f"ERROR: No sequences found in {fasta_file}")
        return False

    # Load PDB structure
    print(f"Loading PDB structure from {pdb_file}")
    pose = pyrosetta.pose_from_pdb(pdb_file)
    print(f"  Structure has {pose.total_residue()} residues in {pose.num_chains()} chain(s)")

    # Check if we have multiple sequences in FASTA
    if len(sequences) > 1:
        print(f"  FASTA file contains {len(sequences)} sequences")
        print(f"  Mapping to chains 1-{len(sequences)}")

        # Map each FASTA sequence to sequential chains
        for seq_idx, (header, sequence) in enumerate(sequences):
            chain_num = seq_idx + 1

            if chain_num > pose.num_chains():
                print(f"  WARNING: Chain {chain_num} not found in PDB (only {pose.num_chains()} chains available)")
                print(f"  Skipping sequence: {header}")
                continue

            print(f"  Chain {chain_num}: {header}")
            print(f"    Sequence length: {len(sequence)}")

            # Validate sequence length
            if not validate_sequence_length(pose, sequence, chain_num):
                print(f"    Skipping due to length mismatch")
                continue

            # Get pose indices for this chain
            chain_start = pose.chain_begin(chain_num)

            # Thread sequence onto this chain
            rsd_set = pose.residue_type_set_for_pose(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)
            for res_idx, aa in enumerate(sequence):
                pose_idx = chain_start + res_idx
                if aa in aa_1_3:
                    name3 = aa_1_3[aa]
                    new_res = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(
                        rsd_set.name_map(name3)
                    )
                    pose.replace_residue(pose_idx, new_res, True)
                elif aa not in ['-', 'X']:
                    print(f"    WARNING: Skipping non-standard amino acid '{aa}' at position {res_idx + 1}")

            print(f"    Successfully threaded sequence onto chain {chain_num}")

    else:
        # Single sequence in FASTA
        header, sequence = sequences[0]
        print(f"  Sequence: {header}")
        print(f"  Length: {len(sequence)}")

        # Handle multi-chain sequences (split by '/')
        if args.split_by_slash and '/' in sequence:
            chain_sequences = sequence.split('/')
            print(f"  Detected multi-chain sequence: {len(chain_sequences)} chains")

            for chain_idx, chain_seq in enumerate(chain_sequences):
                chain_num = chain_idx + 1
                print(f"    Threading chain {chain_num}: {chain_seq[:50]}..." if len(chain_seq) > 50 else f"    Threading chain {chain_num}: {chain_seq}")

                # Validate length
                if not validate_sequence_length(pose, chain_seq, chain_num):
                    print(f"    Skipping due to length mismatch")
                    continue

                # Get pose indices for this chain
                chain_start = pose.chain_begin(chain_num)

                # Thread sequence onto this chain
                rsd_set = pose.residue_type_set_for_pose(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)
                for res_idx, aa in enumerate(chain_seq):
                    pose_idx = chain_start + res_idx
                    if aa in aa_1_3:
                        name3 = aa_1_3[aa]
                        new_res = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(
                            rsd_set.name_map(name3)
                        )
                        pose.replace_residue(pose_idx, new_res, True)
        else:
            # Single chain threading
            print(f"  Threading sequence onto chain {args.chain}")

            # Validate sequence length
            if not validate_sequence_length(pose, sequence, args.chain):
                print(f"ERROR: Sequence length mismatch for {fasta_file}")
                return False

            # Thread the sequence
            pose = thread_mpnn_seq(pose, sequence)

    # Remove clashes if requested
    if args.remove_clashes:
        print("  Removing clashes...")
        pose = remove_clashes(pose)

    # Save the threaded structure
    print(f"  Saving to {output_file}")
    pose.dump_pdb(output_file)
    print(f"  Successfully saved threaded structure")

    return True

def main():
    parser = argparse.ArgumentParser(
        description='Thread sequences from FASTA file onto PDB structure (single or batch mode)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Single file mode arguments
    parser.add_argument('--pdb', type=str,
                        help='Input PDB file (single file mode)')
    parser.add_argument('--fasta', type=str,
                        help='Input FASTA file with sequences to thread (single file mode)')
    parser.add_argument('--out', type=str, default='',
                        help='Output PDB file (single file mode)')
    parser.add_argument('--out_prefix', type=str, default='threaded',
                        help='Prefix for output PDB files (used when multiple sequences in FASTA)')

    # Batch mode arguments
    parser.add_argument('--pdb_dir', type=str,
                        help='Directory containing PDB files (batch mode)')
    parser.add_argument('--fasta_dir', type=str,
                        help='Directory containing FASTA files (batch mode)')
    parser.add_argument('--out_dir', type=str,
                        help='Output directory for generated PDB files (batch mode)')

    # Processing options
    parser.add_argument('--chain', type=int, default=1,
                        help='Chain number to thread sequence onto (1-indexed)')
    parser.add_argument('--remove_clashes', action='store_true',
                        help='Remove clashes after threading (recommended)')
    parser.add_argument('--split_by_slash', action='store_true',
                        help='Split sequences by "/" character for multi-chain designs')

    args = parser.parse_args()

    # Determine mode: batch or single file
    batch_mode = args.pdb_dir is not None and args.fasta_dir is not None
    single_mode = args.pdb is not None and args.fasta is not None

    if not batch_mode and not single_mode:
        print("ERROR: Must specify either:")
        print("  Single file mode: --pdb and --fasta")
        print("  Batch mode: --pdb_dir and --fasta_dir")
        sys.exit(1)

    if batch_mode and single_mode:
        print("ERROR: Cannot use both single file mode and batch mode simultaneously")
        sys.exit(1)

    # ============= BATCH MODE =============
    if batch_mode:
        print("="*60)
        print("BATCH MODE: Processing multiple PDB-FASTA pairs")
        print("="*60)

        # Check directories exist
        if not os.path.exists(args.pdb_dir):
            print(f"ERROR: PDB directory not found: {args.pdb_dir}")
            sys.exit(1)

        if not os.path.exists(args.fasta_dir):
            print(f"ERROR: FASTA directory not found: {args.fasta_dir}")
            sys.exit(1)

        # Create output directory if needed
        out_dir = args.out_dir if args.out_dir else '.'
        if not os.path.exists(out_dir):
            print(f"Creating output directory: {out_dir}")
            os.makedirs(out_dir)

        # Find matching files
        print(f"\nMatching files:")
        print(f"  PDB directory: {args.pdb_dir}")
        print(f"  FASTA directory: {args.fasta_dir}")
        print(f"  Output directory: {out_dir}")
        print()

        matches = find_matching_files(args.pdb_dir, args.fasta_dir)

        if not matches:
            print("ERROR: No matching PDB-FASTA pairs found")
            sys.exit(1)

        print(f"\nFound {len(matches)} matching pair(s)")
        print("="*60)

        # Process each pair
        successful = 0
        failed = 0

        for idx, (pdb_file, fasta_file, output_name) in enumerate(matches):
            print(f"\n{'='*60}")
            print(f"Processing pair {idx+1}/{len(matches)}")
            print(f"  PDB:    {os.path.basename(pdb_file)}")
            print(f"  FASTA:  {os.path.basename(fasta_file)}")
            print(f"  Output: {output_name}")
            print('='*60)

            output_path = os.path.join(out_dir, output_name)

            try:
                success = process_single_pair(pdb_file, fasta_file, output_path, args)
                if success:
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                print(f"ERROR processing pair: {e}")
                failed += 1

        # Summary
        print(f"\n{'='*60}")
        print("BATCH PROCESSING COMPLETE")
        print(f"  Successful: {successful}/{len(matches)}")
        print(f"  Failed:     {failed}/{len(matches)}")
        print("="*60)

    # ============= SINGLE FILE MODE =============
    else:
        print("="*60)
        print("SINGLE FILE MODE")
        print("="*60)

        # Check files exist
        if not os.path.exists(args.pdb):
            print(f"ERROR: PDB file not found: {args.pdb}")
            sys.exit(1)

        if not os.path.exists(args.fasta):
            print(f"ERROR: FASTA file not found: {args.fasta}")
            sys.exit(1)

        # Read sequences from FASTA
        print(f"Reading sequences from {args.fasta}")
        sequences = read_fasta(args.fasta)
        print(f"Found {len(sequences)} sequence(s)")

        if len(sequences) == 0:
            print("ERROR: No sequences found in FASTA file")
            sys.exit(1)

        # Load original PDB
        print(f"Loading PDB structure from {args.pdb}")
        pose = pyrosetta.pose_from_pdb(args.pdb)
        print(f"Structure has {pose.total_residue()} residues in {pose.num_chains()} chain(s)")

        # Check if we have multiple sequences in FASTA
        if len(sequences) > 1:
            print(f"\nFASTA file contains {len(sequences)} sequences")
            print(f"Mapping to chains 1-{len(sequences)}")

            # Map each FASTA sequence to sequential chains
            for seq_idx, (header, sequence) in enumerate(sequences):
                chain_num = seq_idx + 1

                if chain_num > pose.num_chains():
                    print(f"  WARNING: Chain {chain_num} not found in PDB (only {pose.num_chains()} chains available)")
                    print(f"  Skipping sequence: {header}")
                    continue

                print(f"\nChain {chain_num}: {header}")
                print(f"  Sequence length: {len(sequence)}")

                # Validate sequence length
                if not validate_sequence_length(pose, sequence, chain_num):
                    print(f"  Skipping due to length mismatch")
                    continue

                # Get pose indices for this chain
                chain_start = pose.chain_begin(chain_num)

                # Thread sequence onto this chain
                rsd_set = pose.residue_type_set_for_pose(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)
                for res_idx, aa in enumerate(sequence):
                    pose_idx = chain_start + res_idx
                    if aa in aa_1_3:
                        name3 = aa_1_3[aa]
                        new_res = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(
                            rsd_set.name_map(name3)
                        )
                        pose.replace_residue(pose_idx, new_res, True)
                    elif aa not in ['-', 'X']:
                        print(f"  WARNING: Skipping non-standard amino acid '{aa}' at position {res_idx + 1}")

                print(f"  Successfully threaded sequence onto chain {chain_num}")

            # Remove clashes if requested
            if args.remove_clashes:
                print("\nRemoving clashes...")
                pose = remove_clashes(pose)

            # Determine output filename
            if args.out:
                out_file = args.out
            else:
                out_file = f"{args.out_prefix}_multichain.pdb"

            # Save the threaded structure
            print(f"\nSaving to {out_file}")
            pose.dump_pdb(out_file)
            print(f"Successfully saved threaded structure with {len(sequences)} chain(s)")

        else:
            # Single sequence in FASTA
            header, sequence = sequences[0]
            print(f"\nSequence: {header}")
            print(f"Length: {len(sequence)}")

            # Handle multi-chain sequences (split by '/')
            if args.split_by_slash and '/' in sequence:
                chain_sequences = sequence.split('/')
                print(f"Detected multi-chain sequence: {len(chain_sequences)} chains")

                for chain_idx, chain_seq in enumerate(chain_sequences):
                    chain_num = chain_idx + 1
                    print(f"  Threading chain {chain_num}: {chain_seq[:50]}..." if len(chain_seq) > 50 else f"  Threading chain {chain_num}: {chain_seq}")

                    # Validate length
                    if not validate_sequence_length(pose, chain_seq, chain_num):
                        print(f"  Skipping due to length mismatch")
                        continue

                    # Get pose indices for this chain
                    chain_start = pose.chain_begin(chain_num)

                    # Thread sequence onto this chain
                    rsd_set = pose.residue_type_set_for_pose(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)
                    for res_idx, aa in enumerate(chain_seq):
                        pose_idx = chain_start + res_idx
                        if aa in aa_1_3:
                            name3 = aa_1_3[aa]
                            new_res = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(
                                rsd_set.name_map(name3)
                            )
                            pose.replace_residue(pose_idx, new_res, True)
            else:
                # Single chain threading
                print(f"Threading sequence onto chain {args.chain}")

                # Validate sequence length
                if not validate_sequence_length(pose, sequence, args.chain):
                    print("ERROR: Sequence length mismatch")
                    sys.exit(1)

                # Thread the sequence
                pose = thread_mpnn_seq(pose, sequence)

            # Remove clashes if requested
            if args.remove_clashes:
                print("Removing clashes...")
                pose = remove_clashes(pose)

            # Determine output filename
            if args.out:
                out_file = args.out
            else:
                # Use header for filename, sanitize it
                safe_header = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in header)
                out_file = f"{args.out_prefix}_{safe_header}.pdb"

            # Save the threaded structure
            print(f"Saving to {out_file}")
            pose.dump_pdb(out_file)
            print(f"Successfully saved threaded structure")

        print(f"\n{'='*60}")
        print("COMPLETED")
        print("="*60)

if __name__ == "__main__":
    main()
