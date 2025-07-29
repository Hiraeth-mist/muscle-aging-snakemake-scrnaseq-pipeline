# scripts/01_initial_qc.py

import matplotlib
matplotlib.use('Agg') # Use the 'Agg' backend for non-interactive plotting


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse # Import the argparse module

# --- Add Argument Parsing ---
parser = argparse.ArgumentParser(description='Initial QC for single-cell RNA-seq data.')
parser.add_argument('--input_barcodes', required=True, help='Path to barcodes.tsv.gz file.')
parser.add_argument('--input_features', required=True, help='Path to features.tsv.gz file.')
parser.add_argument('--input_matrix', required=True, help='Path to matrix.mtx.gz file.')
parser.add_argument('--output_adata', required=True, help='Path to output AnnData raw file (.h5ad).')
parser.add_argument('--output_qc_plot', required=True, help='Path to output QC plot (.png).')
parser.add_argument('--qc_plot_dir', required=True, help='Directory for QC plots.')
parser.add_argument('--sample_name', required=True, help='Name of the sample.')

args = parser.parse_args()

# Get inputs and outputs from parsed arguments - Keep these lines, they correctly use 'args'
input_barcodes = args.input_barcodes
input_features = args.input_features
input_matrix = args.input_matrix
output_adata = args.output_adata
output_qc_plot = args.output_qc_plot
qc_plot_dir = args.qc_plot_dir
sample_name = args.sample_name

# This line now uses 'input_barcodes' which came from 'args'
input_sample_dir = os.path.dirname(input_barcodes)

# Create output directories if they don't exist
os.makedirs(os.path.dirname(output_adata), exist_ok=True)
os.makedirs(qc_plot_dir, exist_ok=True)

# Read 10x Genomics data from the sample directory
adata = sc.read_10x_mtx(
    input_sample_dir,                # Path to directory containing barcodes.tsv, features.tsv, matrix.mtx
    var_names='gene_symbols',        # Use gene symbols from features.tsv
    cache=True                       # Cache for faster loading 
)

# Add sample name to adata object (important for distinguishing samples later)
adata.obs['sample'] = sample_name

# Basic QC
# Identify mitochondrial genes (adjust pattern if needed, e.g., for human 'MT-')
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Plot QC metrics
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, show=False, ax=axes[0])
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False, ax=axes[1]) 
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False, ax=axes[2])
plt.suptitle(f"QC Metrics for {sample_name}") # Add a title to the plot
plt.tight_layout()
plt.savefig(output_qc_plot)
plt.close(fig)

# Save the raw AnnData object
adata.write(output_adata)
print(f"Initial QC and raw adata saved for {sample_name}")