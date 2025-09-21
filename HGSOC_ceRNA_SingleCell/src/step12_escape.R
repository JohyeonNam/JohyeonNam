#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(dplyr)
library(Seurat)
library(escape)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)
library(ggplot2)
library(msigdb)
library(BiocParallel)
library(GSVA)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory

geneset_enrichment <- function(seurat_obj, dir_out) {
    print("Function geneset_enrichment")
    
    seurat_obj <- DietSeurat(
        seurat_obj,
        assays = "RNA",
        dimreducs = c("integrated.harmony.full", "umap.full")
    )
    
    print("Diet finished")

    hallmark <- getGeneSets(library = "H")
    kegg <- getGeneSets(library = "C2", subcategory = "CP:KEGG")
    reactome <- getGeneSets(library = "C2", subcategory = "CP:REACTOME")

    reflist <- list(
        hallmark = hallmark, 
        kegg = kegg, 
        reactome = reactome
    )

    print("Pathway download complete")

    for (name in names(reflist)) {
        ref <- reflist[[name]]

        print(paste0("Processing ", name))

        escape_mat <- escape.matrix(
            seurat_obj,
            gene.sets = ref,
            groups = 1000,
            min.size = 5,
            BPPARAM = MulticoreParam(workers = 8)
        )       

        seurat_obj <- performNormalization(
            seurat_obj,
            enrichment.data = escape_mat,
            gene.sets = ref
        )

        seurat_obj <- RenameAssays(
            seurat_obj,
            assay.name = "escape_normalized",
            new.assay.name = paste0(name, "_normalised")
        )
    }

    return(seurat_obj)
}

seurat_obj <- geneset_enrichment(seurat_obj, dir_out)

saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_9.rds"))