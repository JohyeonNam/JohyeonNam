#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(harmony)
library(clustree)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory

seurat_harmony <- function(seurat_obj) {
    print("Function seurat_harmony")

    options(future.globals.maxSize = Inf)

    DefaultAssay(seurat_obj) <- "sketch"

    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)

    seurat_obj <- IntegrateLayers(
        object = seurat_obj, 
        method = HarmonyIntegration,
        orig = "pca",
        new.reduction = "integrated.harmony",
        verbose = TRUE
    )
    
    return(seurat_obj)
}

seurat_obj <- seurat_harmony(seurat_obj)

saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_5.rds"))