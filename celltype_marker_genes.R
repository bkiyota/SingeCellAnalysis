library(Seurat)
library(tidyverse)

args     <- commandArgs(trailingOnly=T)
celltype <- args[1]

rna_counts_path <- "path/to/processed/train_cite_inputs.csv.gz"
metadata_path   <- "path/to/raw/metadata.csv"
out_dir         <- "path/to/processed/seurat_marker_genes"

rna_df <- read.csv(rna_counts_path,row.names=1)
rna_df <- t(rna_df)

metadata_df <- read.csv(metadata_path,row.names=1) 
metadata_df <- metadata_df %>%
    filter(rownames(metadata_df) %in% colnames(rna_df))

rna <- CreateSeuratObject(
    counts=rna_df,
    project="cs547",
    assay="RNA",
    meta.data=metadata_df
)
Idents(rna) <- "cell_type"
rna <- SetAssayData(rna, slot="scale.data",new.data=rna_df)

markers_df <- FindMarkers(rna,ident.1=celltype)
outpath    <- paste0(out_dir,"/",celltype,"_marker_genes.csv")
write.csv(markers_df,outpath)