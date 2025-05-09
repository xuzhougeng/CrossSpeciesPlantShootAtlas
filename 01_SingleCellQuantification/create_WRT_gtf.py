#!/usr/bin/env python3

import sys

def parse_fasta(fasta_path):
    """
    Parse a FASTA file and yield sequence records
    
    Args:
        fasta_path (str): Path to input FASTA file
        
    Yields:
        tuple: (sequence_id, sequence)
    """
    seq_id = None
    sequence = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # If we already have a sequence, yield it
                if seq_id:
                    yield seq_id, ''.join(sequence)
                
                # Start a new sequence
                seq_id = line[1:].split()[0]  # Get ID (first part after '>')
                sequence = []
            else:
                sequence.append(line)
        
        # Yield the last sequence if there is one
        if seq_id:
            yield seq_id, ''.join(sequence)

def write_gtf(fasta_path, gtf_path):
    """
    Create a simplified GTF file from transcriptome FASTA for WRT pipeline.
    Each transcript is treated as a single exon with gene and transcript IDs.
    
    Args:
        fasta_path (str): Path to input FASTA file
        gtf_path (str): Path to output GTF file
    """
    try:
        with open(gtf_path, 'w') as gtf_out:
            for seq_id, sequence in parse_fasta(fasta_path):
                # Extract transcript information
                length = len(sequence)
                
                # GTF fields
                chrom = seq_id  # Each transcript as a reference sequence
                source = "WRT"
                feature = "exon"
                start = 1
                end = length
                score = "."
                strand = "+"  # Default to positive strand
                frame = "."
                
                # Generate gene and transcript IDs
                gene_id = f"gene_{seq_id}"
                transcript_id = f"transcript_{seq_id}"
                
                # Format GTF attributes
                attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
                
                # Create GTF line
                gtf_line = "\t".join([
                    chrom, source, feature, str(start), str(end),
                    score, strand, frame, attributes
                ])
                gtf_out.write(gtf_line + "\n")
                
        print(f"Successfully created GTF file: {gtf_path}")
        
    except FileNotFoundError:
        print(f"Error: Could not find input FASTA file: {fasta_path}")
        sys.exit(1)
    except Exception as e:
        print(f"Error creating GTF file: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_wrt_gtf.py input.fasta output.gtf")
        sys.exit(1)

    fasta_input = sys.argv[1]
    gtf_output = sys.argv[2]
    write_gtf(fasta_input, gtf_output)
