#!/usr/bin/env python3
"""
Extract nucleotide and amino acid reference sequences from GenBank files.

USAGE:
  python extract_gbk_references.py reference.gbk -o output_dir

  Creates:
    - output_dir/nuc_ref.fasta (nucleotide sequence from ORIGIN)
    - output_dir/aa_ref.fasta (amino acid sequences from CDS translations)

EXAMPLE:
  python scripts/extract_gbk_references.py \\
      resources/references/rsva/nextclade/reference.gbk \\
      -o resources/references/rsva

OPTIONS:
  -o, --output-dir DIR     Output directory (default: same as input)
  --nuc-output FILE        Custom nucleotide FASTA output path
  --aa-output FILE         Custom amino acid FASTA output path
  --wrap INT               Wrap sequences at N characters (0=no wrap, default=70)

OUTPUT FORMAT:

  nuc_ref.fasta:
    >LOCUS_NAME
    FULL_GENOME_SEQUENCE (70 chars per line)

  aa_ref.fasta:
    >GENE_NAME
    PROTEIN_SEQUENCE (70 chars per line, includes stop codons *)
    >GENE_NAME2
    PROTEIN_SEQUENCE2
    ...

IMPLEMENTATION:
  - Parses GenBank LOCUS name
  - Extracts sequence from ORIGIN section
  - Identifies all CDS features with /translation qualifiers
  - Preserves gene names and stop codons
  - Outputs standard FASTA format
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple


def parse_genbank_file(gbk_path: Path) -> Tuple[str, str, List[Dict[str, str]]]:
    """
    Parse a GenBank file to extract sequences.

    Args:
        gbk_path: Path to the GenBank file

    Returns:
        Tuple of (locus_name, nucleotide_sequence, list_of_cds_features)
    """
    with open(gbk_path, 'r') as f:
        content = f.read()

    # Extract LOCUS name
    locus_match = re.search(r'^LOCUS\s+(\S+)', content, re.MULTILINE)
    locus_name = locus_match.group(1) if locus_match else "Unknown"

    # Extract nucleotide sequence from ORIGIN section
    origin_match = re.search(r'ORIGIN\s*\n(.*?)^//', content, re.MULTILINE | re.DOTALL)
    if not origin_match:
        raise ValueError("No ORIGIN section found in GenBank file")

    origin_text = origin_match.group(1)
    # Remove line numbers and whitespace, keep only nucleotides
    nuc_sequence = re.sub(r'[\s\d]', '', origin_text).upper()

    # Extract CDS features with their translations
    cds_features = []

    # Find all CDS features
    cds_pattern = r'^\s+CDS\s+(\S+)\n(.*?)(?=^\s+(?:gene|CDS|mRNA|5\'UTR|3\'UTR|source|ORIGIN)|^//)'
    cds_matches = re.finditer(cds_pattern, content, re.MULTILINE | re.DOTALL)

    for match in cds_matches:
        location = match.group(1)
        feature_text = match.group(2)

        # Extract gene name
        gene_match = re.search(r'/gene="([^"]+)"', feature_text)
        gene_name = gene_match.group(1) if gene_match else None

        # Extract product name
        product_match = re.search(r'/product="([^"]+)"', feature_text)
        product = product_match.group(1) if product_match else None

        # Extract translation
        translation_match = re.search(r'/translation="([^"]+)"', feature_text)
        if translation_match:
            translation = translation_match.group(1).replace('\n', '').replace(' ', '')

            # Create feature entry
            # Use gene name if available, otherwise use product or location
            feature_id = gene_name if gene_name else (product if product else location)

            cds_features.append({
                'id': feature_id,
                'location': location,
                'gene': gene_name,
                'product': product,
                'translation': translation
            })

    return locus_name, nuc_sequence, cds_features


def write_fasta(sequences: Dict[str, str], output_path: Path, wrap_length: int = 70):
    """
    Write sequences to a FASTA file.

    Args:
        sequences: Dictionary mapping sequence IDs to sequences
        output_path: Path to output FASTA file
        wrap_length: Number of characters per line (0 for no wrapping)
    """
    with open(output_path, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n")

            if wrap_length > 0:
                # Wrap sequence to specified length
                for i in range(0, len(sequence), wrap_length):
                    f.write(sequence[i:i+wrap_length] + '\n')
            else:
                f.write(sequence + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Extract reference sequences from GenBank files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example:
  python extract_gbk_references.py reference.gbk -o output_dir

  This will create:
    output_dir/nuc_ref.fasta - Nucleotide reference sequence
    output_dir/aa_ref.fasta  - Amino acid reference sequences (CDS translations)
        '''
    )

    parser.add_argument('gbk_file', type=Path, help='Input GenBank (.gbk) file')
    parser.add_argument('-o', '--output-dir', type=Path,
                        help='Output directory (default: same as input file)',
                        default=None)
    parser.add_argument('--nuc-output', type=Path,
                        help='Custom output path for nucleotide FASTA',
                        default=None)
    parser.add_argument('--aa-output', type=Path,
                        help='Custom output path for amino acid FASTA',
                        default=None)
    parser.add_argument('--wrap', type=int, default=70,
                        help='Wrap sequences at this length (0 for no wrapping, default: 70)')

    args = parser.parse_args()

    # Validate input file
    if not args.gbk_file.exists():
        raise FileNotFoundError(f"Input file not found: {args.gbk_file}")

    # Determine output directory
    if args.output_dir:
        output_dir = args.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = args.gbk_file.parent

    # Determine output file paths
    nuc_output = args.nuc_output if args.nuc_output else output_dir / 'nuc_ref.fasta'
    aa_output = args.aa_output if args.aa_output else output_dir / 'aa_ref.fasta'

    print(f"Reading GenBank file: {args.gbk_file}")
    locus_name, nuc_sequence, cds_features = parse_genbank_file(args.gbk_file)

    print(f"Found locus: {locus_name}")
    print(f"Nucleotide sequence length: {len(nuc_sequence)} bp")
    print(f"Found {len(cds_features)} CDS features with translations")

    # Write nucleotide reference
    nuc_sequences = {locus_name: nuc_sequence}
    write_fasta(nuc_sequences, nuc_output, wrap_length=args.wrap)
    print(f"Wrote nucleotide reference to: {nuc_output}")

    # Write amino acid references
    aa_sequences = {}
    for feature in cds_features:
        # Use gene name as ID, fall back to product
        seq_id = feature['gene'] if feature['gene'] else feature['id']
        aa_sequences[seq_id] = feature['translation']

    write_fasta(aa_sequences, aa_output, wrap_length=args.wrap)
    print(f"Wrote {len(aa_sequences)} amino acid sequences to: {aa_output}")

    # Print summary
    print("\nAmino acid sequences:")
    for seq_id, seq in aa_sequences.items():
        print(f"  {seq_id}: {len(seq)} aa")


if __name__ == '__main__':
    main()
