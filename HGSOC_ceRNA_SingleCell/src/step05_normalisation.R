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

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory

seurat_normalisation <- function(seurat_obj) {
        print("Function seurat_normalisation")

        options(future.globals.maxSize = Inf)

        seurat_obj <- NormalizeData(seurat_obj)
        seurat_obj <- FindVariableFeatures(seurat_obj)
        
        return(seurat_obj)
}

seurat_obj <- seurat_normalisation(seurat_obj)

saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_3.rds"))