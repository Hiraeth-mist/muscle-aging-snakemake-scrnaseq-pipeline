# Snakefile

import pandas as pd
import os

# Load configuration
configfile: "config.yaml"

SAMPLES = config["samples"]
DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]

# --- ALL rule (what you want to achieve at the end) ---
rule all:
    input:
        # Expected outputs from the last rule
        f"{OUTPUT_DIR}/integrated_adata_annotated.h5ad", # The final annotated data object
        f"{OUTPUT_DIR}/figures/umap_clusters.png",
        f"{OUTPUT_DIR}/figures/umap_age_group.png",
        f"{OUTPUT_DIR}/figures/umap_annotated.png",
        f"{OUTPUT_DIR}/diff_exp_results/young_vs_old_DEGs.csv",
        f"{OUTPUT_DIR}/diff_exp_results/cluster_marker_genes.csv"

# --- Rule to load and preprocess individual samples ---
rule load_and_initial_qc:
    input:
        barcodes=f"{DATA_DIR}/{{sample}}/barcodes.tsv.gz",
        features=f"{DATA_DIR}/{{sample}}/features.tsv.gz",
        matrix=f"{DATA_DIR}/{{sample}}/matrix.mtx.gz"
    output:
        adata_raw=f"{OUTPUT_DIR}/processed_data/{{sample}}_raw.h5ad",
        qc_plot=f"{OUTPUT_DIR}/qc_plots/{{sample}}_qc_metrics.png"
    params:
        output_dir=f"{OUTPUT_DIR}/processed_data",
        qc_plot_dir=f"{OUTPUT_DIR}/qc_plots",
        script_path="scripts/01_initial_qc.py"
    conda: "envs/scanpy_env.yaml"
    shell:
        """
        

        env -u MPLBACKEND python {params.script_path} \
               --input_barcodes {input.barcodes} \
               --input_features {input.features} \
               --input_matrix {input.matrix} \
               --output_adata {output.adata_raw} \
               --output_qc_plot {output.qc_plot} \
               --qc_plot_dir {params.qc_plot_dir} \
               --sample_name {wildcards.sample}
        """

# --- Rule to integrate and cluster all samples ---
rule integrate_and_cluster:
    input:
        adatas=expand(f"{OUTPUT_DIR}/processed_data/{{sample}}_raw.h5ad", sample=SAMPLES)
    output:
        integrated_adata=f"{OUTPUT_DIR}/integrated_adata.h5ad",
        umap_plot=f"{OUTPUT_DIR}/figures/umap_clusters.png",
        age_group_umap_plot=f"{OUTPUT_DIR}/figures/umap_age_group.png"
    params:
        output_dir=OUTPUT_DIR,
        fig_dir=f"{OUTPUT_DIR}/figures",
        script_path="scripts/02_integration_clustering.py"
    conda: "envs/scanpy_env.yaml"
    shell: # <-- This line is correct
        """ # <-- This opening triple quote MUST be indented 8 spaces (4 from 'shell:')
        
        env -u MPLBACKEND python {params.script_path} \
               --input_adatas {input.adatas} \
               --output_integrated_adata {output.integrated_adata} \
               --umap_plot {output.umap_plot} \
               --age_group_umap_plot {output.age_group_umap_plot} \
               --fig_dir {params.fig_dir}
        """ # <-- This closing triple quote MUST be indented 8 spaces

# --- Rule for Differential Expression Analysis and Cell Type Annotation ---
rule differential_expression_and_annotation:
    input:
        integrated_adata=f"{OUTPUT_DIR}/integrated_adata.h5ad"
    output:
        diff_exp_csv=f"{OUTPUT_DIR}/diff_exp_results/young_vs_old_DEGs.csv",
        marker_genes_csv=f"{OUTPUT_DIR}/diff_exp_results/cluster_marker_genes.csv",
        annotated_umap_plot=f"{OUTPUT_DIR}/figures/umap_annotated.png",
        final_annotated_adata=f"{OUTPUT_DIR}/integrated_adata_annotated.h5ad"
    params:
        output_dir=f"{OUTPUT_DIR}/diff_exp_results",
        fig_dir=f"{OUTPUT_DIR}/figures",
        script_path="scripts/03_de_and_annotation.py"
    conda: "envs/scanpy_env.yaml"
    shell: # <-- This line is correct
        """ # <-- This opening triple quote MUST be indented 8 spaces (4 from 'shell:')
        
        env -u MPLBACKEND python {params.script_path} \
               --input_integrated_adata {input.integrated_adata} \
               --output_diff_exp_csv {output.diff_exp_csv} \
               --output_marker_genes_csv {output.marker_genes_csv} \
               --output_annotated_umap_plot {output.annotated_umap_plot} \
               --output_final_annotated_adata {output.final_annotated_adata} \
               --output_dir {OUTPUT_DIR}/diff_exp_results \
               --fig_dir {OUTPUT_DIR}/figures
        """ # <-- This closing triple quote MUST be indented 8 spaces