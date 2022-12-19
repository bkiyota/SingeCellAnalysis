import scanpy as sc
import pandas as pd

print("Libraries loaded successfully.",flush=True)

cite_rna_df = pd.read_hdf("/path/to/raw/train_cite_inputs.h5")
outpath = "/path/to/processed/train_cite_inputs.csv.gz"
cite_rna_df.to_csv(outpath,compression="gzip")
del cite_rna_df


cite_protein_df = pd.read_hdf("/path/to/raw/train_cite_targets.h5")
outpath = "/path/to/processed/train_cite_targets.csv.gz"
cite_protein_df.to_csv(outpath,compression="gzip")

print("Done!",flush=True)