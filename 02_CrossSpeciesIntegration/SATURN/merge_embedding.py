import json
import pickle
from pathlib import Path
from typing_extensions import Literal
import glob
import torch
from tap import Tap
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor



# constants
BLAST_HEADER = [
    'qseqid',
    'sseqid',
    'pident',
    'length',
    'mismatch',
    'gapopen',
    'qstart',
    'qend',
    'sstart',
    'send',
    'evalue',
    'bitscore'
]

# Vocabulary of amino acid characters
# Note: 0 reserved for pad token; '-' for class token where full sequence embedding will be extracted
# Note: '*' is for stop codon; 'X' is for any amino acid; 'U' is for Selenocysteine
VOCAB = {
    character: index
    for index, character in enumerate([
        '-', '*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y'
    ])
}

# Last layer of pretrained transformer
LAST_LAYER = 33 # ESM1b
MSA_LAST_LAYER = 12 # MSA
LAST_LAYER_2 = 48 # ESM2


def load_proteins_embedding(embedding_dir: str, last_layer: int, max_workers: int = 10) -> dict:
    """Load protein embedding from file"""

    # Get protein embedding paths
    protein_embedding_paths = [Path(x) for x in glob.glob(embedding_dir + "*.pt")]

    # Define a function to load embedding for a single protein
    def load_single_embedding(protein_embedding_path: Path):
        return protein_embedding_path.stem, torch.load(protein_embedding_path)['mean_representations'][last_layer]

    # Use ThreadPoolExecutor to parallelize loading
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all loading jobs to the executor
        future_to_path = {executor.submit(load_single_embedding, path): path for path in protein_embedding_paths if path.exists()}

        # Collect results as they complete
        protein_id_to_embedding = {}
        for future in tqdm(as_completed(future_to_path), total=len(future_to_path)):
            path = future_to_path[future]
            try:
                protein_id, embedding = future.result()
                protein_id_to_embedding[protein_id] = embedding
            except Exception as exc:
                print(f"{path.name} generated an exception: {exc}")

    return protein_id_to_embedding



def load_from_multiple_directories(directories: list, last_layer: int, max_workers: int = 4):
    """Load embeddings from multiple directories in parallel"""
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Map each directory to a separate task
        futures = [executor.submit(load_proteins_embedding, d, last_layer) for d in directories]
        
        # Wait for all futures to complete and collect results
        results = []
        for future in as_completed(futures):
            results.append(future.result())

    # Merge dictionaries from each directory

    return results


def merge_embedding(protein_id_to_embedding):

    # key是基因，value是一组转录本
    gene_symbol_to_protein_ids = {}
    for transcript_id in protein_id_to_embedding.keys():
        gene_id = transcript_id.split(".")[0]
        if gene_id in gene_symbol_to_protein_ids:
            gene_symbol_to_protein_ids[gene_id].append(transcript_id)
        else:
            gene_symbol_to_protein_ids[gene_id] = [transcript_id]


    # Average protein embeddings for each gene to get gene embeddings

    protein_id_to_embedding_avg = {
        gene_symbol: torch.mean(
            torch.stack([
                protein_id_to_embedding[protein_id] for protein_id in protein_ids
            ]),
            dim=0
        )
        for gene_symbol, protein_ids in tqdm(gene_symbol_to_protein_ids.items())
    }

    return protein_id_to_embedding_avg


def main(input_dir, output_dir, prefix, threads=20):
    input_dir = f"{input_dir}/"
    
    pt_dict = load_proteins_embedding(input_dir, LAST_LAYER_2, threads)
    pt_avg = merge_embedding(pt_dict)


    with open(f'{output_dir}/{prefix}_pt.pkl', 'wb') as f:
        pickle.dump(pt_dict, f)

    # Save gene symbol to embedding map
    torch.save(pt_avg, f"{output_dir}/{prefix}_protein_id_to_embedding.pt")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process FASTA file and remove duplicate headers, keeping the longest sequence.')
    parser.add_argument('input_dir', type=str, help='Path to input Protein Embedding file')
    parser.add_argument('output_dir', type=str, help='Path to output')
    parser.add_argument('prefix', type=str, help='prefix of output filename')

    args = parser.parse_args()

    main(args.input_dir, args.output_dir, args.prefix)
