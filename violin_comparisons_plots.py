import os
import gzip
import argparse
import numpy as np
import scanpy as sc
import pandas as pd
import scipy.stats as stat
import matplotlib.pyplot as plt

# Global
FOLD_THRESHOLD   = 5
PVALUE_THRESHOLD = 15 

def run(args):

    rna_count_path      = args.rna_count_path
    metadata_path       = args.metadata_path
    day                 = args.day
    celltype            = args.celltype
    enriched_genes_path = args.enriched_genes_path
    out_dir             = args.out_dir

    label = f"{celltype}_day{day}"

    print(label,flush=True)

    metadata_df = pd.read_csv(metadata_path)
    metadata_df = metadata_df[metadata_df["technology"] == "citeseq"]

    enriched_genes_df        = pd.read_csv(enriched_genes_path,compression="gzip")
    enriched_genes_subset_df = enriched_genes_df[(enriched_genes_df["cell_type"] == celltype) & (enriched_genes_df["day"] == day)]
    gene_list                = [row["gene"] for _,row in enriched_genes_subset_df.iterrows()]

    cite_rna_df    = pd.read_hdf(rna_count_path)
    cite_rna_adata = sc.AnnData(cite_rna_df)
    
    day_dict      = {}
    donor_dict    = {}
    celltype_dict = {}
    for idx,row in metadata_df.iterrows():
        day_dict[     row["cell_id"]] = row["day"]
        donor_dict[   row["cell_id"]] = row["donor"]
        celltype_dict[row["cell_id"]] = row["cell_type"]

    cite_rna_adata.obs[      "day"] = [day_dict[     cell_id] for cell_id in cite_rna_adata.obs_names]
    cite_rna_adata.obs[    "donor"] = [donor_dict[   cell_id] for cell_id in cite_rna_adata.obs_names]
    cite_rna_adata.obs["cell_type"] = [celltype_dict[cell_id] for cell_id in cite_rna_adata.obs_names]

    del metadata_df, enriched_genes_df, cite_rna_df, day_dict, donor_dict, celltype_dict

    print("RNA expression count data and metadata processed...",flush=True)

    # filter all-zero columns (genes) with respect to celltype (across all days)
    # print(f"Filtering all-zero columns (genes) with respect to {celltype}...",flush=True)
    celltype_adata  = cite_rna_adata[cite_rna_adata.obs.cell_type == celltype]
    # nonzero_indices = [sum(celltype_adata.X[:,j]) > 0 for j in range(celltype_adata.X.shape[1])]
    # filtered_adata  = cite_rna_adata[:,nonzero_indices]
    subset_adata    = cite_rna_adata[cite_rna_adata.obs.day == day]

    sc.tl.rank_genes_groups(subset_adata,"cell_type",method='t-test_overestim_var')

    with plt.rc_context({'figure.figsize': (2, 3)}):
        sc.pl.rank_genes_groups_violin(
            subset_adata,
            n_genes=3,
            jitter=True,
            groups = [celltype],
            save=f"rna_rank_genes_groups_violin_{label}_naive.pdf"
        )
    
    with plt.rc_context({'figure.figsize': (2, 3)}):
        sc.pl.rank_genes_groups_violin(
            subset_adata,
            gene_names = gene_list,
            jitter=True,
            groups = [celltype],
            save=f"rna_rank_genes_groups_violin_{label}_volcano.pdf"
        )

    print("Done!",flush=True)

def main():
    parser = argparse.ArgumentParser(description="",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r','--rna_count_path',type=str,required=True,help='')
    parser.add_argument('-m','--metadata_path',type=str,required=True,help='')
    parser.add_argument('-d','--day',type=int,required=True,help='')
    parser.add_argument('-c','--celltype',type=str,required=True,help='')
    parser.add_argument('-e','--enriched_genes_path',type=str,required=True,help='')
    parser.add_argument('-o','--out_dir',type=str,required=True,help='')
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

