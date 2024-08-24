import pandas as pd
import json
import warnings
import torch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.utils import set_seed
import os
import argparse

warnings.filterwarnings('ignore')
set_seed(42)

def initialize_model(args_file, model_file, vocab_file):
    if not all(os.path.exists(f) for f in [args_file, model_file, vocab_file]):
        raise FileNotFoundError(f"Required model files not found: {args_file}, {model_file}, {vocab_file}")

    vocab = GeneVocab.from_file(vocab_file)
    special_tokens = ["<pad>", "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    
    with open(args_file, "r") as f:
        model_configs = json.load(f)
    
    ntokens = len(vocab)
    model = TransformerModel(
        ntokens,
        model_configs["embsize"],
        model_configs["nheads"],
        model_configs["d_hid"],
        model_configs["nlayers"],
        vocab=vocab,
        pad_value=model_configs.get("pad_value", -2),
        n_input_bins=model_configs.get("n_bins", 51),
    )
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    state_dict = torch.load(model_file, map_location=device)
    
    model.load_state_dict(state_dict, strict=False)
    model.to(device)
    model.eval()
    
    return model, vocab, device

def generate_and_save_embeddings(file_path, model, vocab, device, output_file):
    expression_data = pd.read_csv(file_path, index_col=0)
    gene_names = expression_data.index.tolist()
    tokenized_genes = []
    for gene in gene_names:
        upper_gene = gene.upper()
        if upper_gene in vocab:
            tokenized_genes.append(vocab[upper_gene])
        else:
            tokenized_genes.append(vocab["<pad>"])  # Use <pad> token for unknown genes

    gene_ids = torch.tensor(tokenized_genes, dtype=torch.long).to(device)
    with torch.no_grad():
        gene_embeddings = model.encoder(gene_ids)
    gene_embeddings = gene_embeddings.detach().cpu().numpy()
    
    gene_embeddings_dict = {gene: gene_embeddings[i] for i, gene in enumerate(gene_names)}
    embeddings_df = pd.DataFrame(gene_embeddings_dict).T
    embeddings_df.to_csv(output_file)
    print(f'Saved gene embeddings for {len(gene_embeddings_dict)} genes to {output_file}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate gene embeddings')
    parser.add_argument('--input', required=True, help='Path to input expression data CSV file')
    parser.add_argument('--model_dir', default='/app', help='Directory containing model files')
    args = parser.parse_args()

    # Generate output path based on input path
    input_dir = os.path.dirname(args.input)
    input_filename = os.path.basename(args.input)
    output_filename = "EmbeddingsData.csv"
    output_path = os.path.join(input_dir, output_filename)
    
    print("Debug Information:")
    print(f"Input file path: {args.input}")
    print(f"Output file path: {output_path}")
    print(f"Model directory: {args.model_dir}")

    args_file = os.path.join(args.model_dir, "args.json")
    model_file = os.path.join(args.model_dir, "best_model.pt")
    vocab_file = os.path.join(args.model_dir, "vocab.json")

    try:
        model, vocab, device = initialize_model(args_file, model_file, vocab_file)
        generate_and_save_embeddings(args.input, model, vocab, device, output_path)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure that the model files (args.json, best_model.pt, vocab.json) are in the correct directory.")
        print(f"Current model directory: {args.model_dir}")
        exit(1)