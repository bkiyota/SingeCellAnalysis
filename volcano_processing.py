import os
import gzip
import argparse
import numpy as np
import scanpy as sc
import pandas as pd
import scipy.stats as stat

# Global
PSEUDOCOUNT=0.0001

def run(args):

    rna_count_path = args.rna_count_path
    metadata_path  = args.metadata_path
    day            = args.day
    celltype       = args.celltype
    out_dir        = args.out_dir

    print(f"Volcano processing: day {day}, {celltype}", flush=True)

    metadata       = pd.read_csv(metadata_path)
    metadata       = metadata[metadata["technology"] == "citeseq"]
    cite_rna_df    = pd.read_hdf(rna_count_path)
    cite_rna_adata = sc.AnnData(cite_rna_df)
    
    day_dict       = {}
    donor_dict     = {}
    celltype_dict  = {}
    
    for i,row in metadata.iterrows():
        day_dict[row["cell_id"]] = row["day"]
        donor_dict[row["cell_id"]] = row["donor"]
        celltype_dict[row["cell_id"]] = row["cell_type"]
    
    cite_rna_adata.obs["day"] = [day_dict[cell_id] for cell_id in cite_rna_adata.obs_names]
    cite_rna_adata.obs["donor"] = [donor_dict[cell_id] for cell_id in cite_rna_adata.obs_names]
    cite_rna_adata.obs["cell_type"] = [celltype_dict[cell_id] for cell_id in cite_rna_adata.obs_names]

    del metadata, cite_rna_df, day_dict, donor_dict, celltype_dict

    day_subset_adata = cite_rna_adata[cite_rna_adata.obs.day == day]
    celltype_adata   = day_subset_adata[day_subset_adata.obs.cell_type == celltype]
    rest_adata       = day_subset_adata[day_subset_adata.obs.cell_type != celltype]

    n,m = celltype_adata.X.shape
    nonzero_indices = [sum(celltype_adata.X[:,j]) > 0 for j in range(m)]    
    celltype_subset_adata = celltype_adata[:,nonzero_indices]
    rest_subset_adata     = rest_adata[:,nonzero_indices]

    celltype_subset_adata.X += PSEUDOCOUNT
    rest_subset_adata.X     += PSEUDOCOUNT

    del celltype_adata,rest_adata

    foldchanges = np.log2(
        np.divide(
            np.mean(celltype_subset_adata.X,axis=0),
            np.mean(rest_subset_adata.X,axis=0)
        )
    )
    gene_list = celltype_subset_adata.var_names
    num_subset_genes = len(foldchanges)

    ttest_results = stat.ttest_ind(
        celltype_subset_adata.X,
        rest_subset_adata.X,
        axis=0,
        equal_var=True
    )
    pvalues = ttest_results[1]
    transformed_pvalues = -1*np.log10(np.array(pvalues))
    bonferroni_pvalues  = -1*np.log10(num_subset_genes*np.array(pvalues)) # Bonferri correction

    out_path = f"cite_rna_day{day}_{celltype}_volcano_processing.txt.gz"
    with gzip.open(os.path.join(out_dir,out_path),"wt") as handle:
        handle.write("gene,foldchange,pval,transformed_pval,bonferroni_pval\n")
        for gene,foldchange,pval,tpval,bpval in zip(gene_list,foldchanges,pvalues,transformed_pvalues,bonferroni_pvalues):
            handle.write(f"{gene},{foldchange},{pval},{tpval},{bpval}\n")
    
    print("Done!",flush=True)

def main():
    parser = argparse.ArgumentParser(description="",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r','--rna_count_path',type=str,required=True,help='')
    parser.add_argument('-m','--metadata_path',type=str,required=True,help='')
    parser.add_argument('-d','--day',type=int,required=True,help='')
    parser.add_argument('-c','--celltype',type=str,required=True,help='')
    parser.add_argument('-o','--out_dir',type=str,required=True,help='')
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

