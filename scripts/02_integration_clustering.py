# scripts/02_integration_clustering.py

import matplotlib
matplotlib.use('Agg') # Use the 'Agg' backend for non-interactive plotting

import scanpy as sc
import anndata as ad # Keep import as 'ad' for general AnnData usage if needed
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Integration and Clustering of scRNA-seq data.')
parser.add_argument('--input_adatas', nargs='+', required=True, help='Paths to input AnnData raw files (.h5ad), space-separated.')
parser.add_argument('--output_integrated_adata', required=True, help='Path to output integrated AnnData file (.h5ad).')
parser.add_argument('--umap_plot', required=True, help='Path to output UMAP cluster plot (.png).')
parser.add_argument('--age_group_umap_plot', required=True, help='Path to output UMAP age group plot (.png).')
parser.add_argument('--fig_dir', required=True, help='Directory to save plots.')

args = parser.parse_args()

input_adatas = args.input_adatas # This will be a list of paths
output_integrated_adata = args.output_integrated_adata
output_umap_cluster_plot = args.umap_plot
output_umap_age_group_plot = args.age_group_umap_plot
fig_dir = args.fig_dir

# Create output directories if they don't exist
os.makedirs(fig_dir, exist_ok=True)
os.makedirs(os.path.dirname(output_integrated_adata), exist_ok=True)

# Load individual AnnData objects
adatas = [sc.read_h5ad(f) for f in input_adatas]

# Concatenate AnnData objects without directly adding batch_key during concat
adata = ad.concat(adatas, axis=0, join='outer', merge='unique')
adata.obs_names_make_unique() # Ensure observation names (cell IDs) are unique across combined samples

# Add Age Group metadata based on sample names (this logic remains fine)
# The 'sample' column now exists thanks to the concatenate operation with batch_key="sample"
adata.obs['Age_Group'] = ['Young' if 'Young' in s else 'Old' for s in adata.obs['sample'].astype(str)] # Convert to string for 'in' check

# Filter cells and genes (example thresholds, adjust as needed based on your QC plots)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Filter out cells with high mitochondrial content (example threshold, adjust based on QC)
adata = adata[adata.obs['pct_counts_mt'] < 10, :] # Use 'pct_counts_mt' as it's the more common Scanpy output

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True)

# Scale the data (important for PCA)
sc.pp.scale(adata, max_value=10)

# PCA (Dimensionality Reduction)
sc.tl.pca(adata, svd_solver='arpack')

# Harmony for batch correction (recommended for multiple samples/batches to remove technical variation)
# Harmony usually works on the PCA output (adata.obsm['X_pca']).
# Harmony for batch correction
print("Performing Harmony integration...")
# Remove 'inplace=False'. The function will implicitly return the result.
harmony_result = sc.external.pp.harmony_integrate(adata, key='sample', basis='X_pca')

# Check if Harmony returned a valid result.
if harmony_result is not None:
    adata.obsm['X_harmony'] = harmony_result
    print("Harmony integration successful, X_harmony created.")
else:
    # Fallback: if Harmony returns None, it means it couldn't integrate.
    print("WARNING: Harmony integration returned None. Falling back to X_pca for neighbors.")
    adata.obsm['X_harmony'] = adata.obsm['X_pca'] # Use PCA if Harmony failed

# Neighbors graph (used for UMAP and clustering)
# Use the Harmony embedding for a corrected embedding space
# We can be flexible here if X_harmony was successfully created, otherwise use X_pca
if 'X_harmony' in adata.obsm:
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20, use_rep='X_harmony')
    print("Neighbors graph computed using X_harmony.")
else:
    # This block should ideally not be reached if the above fallback worked
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20, use_rep='X_pca')
    print("Neighbors graph computed using X_pca (Harmony was skipped or failed).")

# UMAP (Another dimensionality reduction for visualization)
sc.tl.umap(adata)

# Clustering
sc.tl.leiden(adata, resolution=0.5) # Adjust resolution for desired number of clusters

# Plot UMAP by clusters
fig_clusters, ax_clusters = plt.subplots(figsize=(8, 7)) # Create figure and axes
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='UMAP by Clusters', show=False, ax=ax_clusters)
plt.savefig(output_umap_cluster_plot, bbox_inches='tight') # Save with tight bounding box
plt.close(fig_clusters)

# Plot UMAP by age group
fig_age_group, ax_age_group = plt.subplots(figsize=(8, 7))
sc.pl.umap(adata, color='Age_Group', title='UMAP by Age Group', show=False, ax=ax_age_group)
plt.savefig(output_umap_age_group_plot, bbox_inches='tight')
plt.close(fig_age_group)

# Save the integrated AnnData object for the next step
adata.write(output_integrated_adata)
print("Integration, clustering, and UMAP completed and data saved.")