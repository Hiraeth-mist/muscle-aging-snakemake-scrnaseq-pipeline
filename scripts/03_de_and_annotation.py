# scripts/03_de_and_annotation.py

import matplotlib
matplotlib.use('Agg') # Use the 'Agg' backend for non-interactive plotting

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Differential Expression and Cell Type Annotation.')
parser.add_argument('--input_integrated_adata', required=True, help='Path to input integrated AnnData file (.h5ad).')
parser.add_argument('--output_diff_exp_csv', required=True, help='Path to output differential expression CSV.')
parser.add_argument('--output_marker_genes_csv', required=True, help='Path to output marker genes CSV.')
parser.add_argument('--output_annotated_umap_plot', required=True, help='Path to output annotated UMAP plot (.png).')
parser.add_argument('--output_final_annotated_adata', required=True, help='Path to output final annotated AnnData file (.h5ad).')
# Add these new arguments, needed because they are used by os.makedirs
parser.add_argument('--output_dir', required=True, help='Base output directory for results (e.g., diff_exp_results).')
parser.add_argument('--fig_dir', required=True, help='Base directory for figures.')

args = parser.parse_args()

# --- THESE ARE THE ONLY LINES THAT SHOULD ASSIGN INPUT/OUTPUTS/PARAMS ---
input_adata = args.input_integrated_adata
output_diff_exp_csv = args.output_diff_exp_csv
output_marker_genes_csv = args.output_marker_genes_csv
output_annotated_umap_plot = args.output_annotated_umap_plot
output_final_annotated_adata = args.output_final_annotated_adata
output_dir = args.output_dir
fig_dir = args.fig_dir
# --- END ASSIGNMENT ---

# Create output directories if they don't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

# Load the integrated AnnData object
adata = sc.read_h5ad(input_adata)

# --- Differential Expression Analysis (Young vs Old across all cells) ---
print("Performing differential expression analysis between Young and Old groups...")
# Ensure 'Age_Group' is a categorical variable needed for rank_genes_groups
if not pd.api.types.is_categorical_dtype(adata.obs['Age_Group']):
    adata.obs['Age_Group'] = adata.obs['Age_Group'].astype('category')

sc.tl.rank_genes_groups(adata, 'Age_Group', method='wilcoxon')
de_results_young_vs_old = pd.DataFrame({
    group + '_' + key: adata.uns['rank_genes_groups'][key][group]
    for group in adata.obs['Age_Group'].cat.categories
    for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']
})
de_results_young_vs_old.to_csv(output_diff_exp_csv, index=False)
print(f"Differential expression results saved to {output_diff_exp_csv}")

# --- Find marker genes for each cluster ---
print("Finding marker genes for each cluster...")
if not pd.api.types.is_categorical_dtype(adata.obs['leiden']):
    adata.obs['leiden'] = adata.obs['leiden'].astype('category')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
marker_genes_df = pd.DataFrame({
    cluster + '_' + key: adata.uns['rank_genes_groups'][key][cluster]
    for cluster in adata.obs['leiden'].cat.categories
    for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']
})
marker_genes_df.to_csv(output_marker_genes_csv, index=False)
print(f"Cluster marker genes saved to {output_marker_genes_csv}")


# --- Cell Type Annotation (Manual or using reference datasets) ---
print("Applying placeholder cell type annotation. Please review 'cluster_marker_genes.csv' for biological annotation.")
cluster_names_map = {
    "0": "Cluster 0 (Review Markers)",
    "1": "Cluster 1 (Review Markers)",
    "2": "Cluster 2 (Review Markers)",
    "3": "Cluster 3 (Review Markers)",
    "4": "Cluster 4 (Review Markers)",
    "5": "Cluster 5 (Review Markers)"
}
# Ensure 'leiden' is a categorical type for renaming
# It should already be categorical from sc.tl.leiden, but explicit cast is safe.
adata.obs['leiden'] = adata.obs['leiden'].astype('category')
adata.rename_categories('leiden', [cluster_names_map.get(str(c), str(c)) for c in adata.obs['leiden'].cat.categories])
adata.obs['cell_type'] = adata.obs['leiden'] # Create a new column for annotated cell types

# Plot UMAP with annotated cell types
fig_annotated_umap, ax_annotated_umap = plt.subplots(figsize=(8, 7))
sc.pl.umap(adata, color='cell_type', legend_loc='on data', title='UMAP with Annotated Cell Types', show=False, ax=ax_annotated_umap)
plt.savefig(output_annotated_umap_plot, bbox_inches='tight')
plt.close(fig_annotated_umap)
print(f"Annotated UMAP plot saved to {output_annotated_umap_plot}")

# Save the final AnnData object with annotations
adata.write(output_final_annotated_adata)
print("Final analysis complete and annotated adata saved.")