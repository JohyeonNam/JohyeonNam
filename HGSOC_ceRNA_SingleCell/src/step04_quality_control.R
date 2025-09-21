#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(dplyr)
library(tibble)
library(data.table)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory

seurat_quality_control <- function(seurat_obj, dir_out) {
    print("Function seurat_quality_control")

    seurat_obj <- seurat_obj %>%
        subset(nFeature_RNA > 250 & nFeature_RNA < 5000) %>%
        subset(nCount_RNA > 1000 & nCount_RNA < 20000) %>%
        subset(percent.mt < 15) %>%
        subset(percent.ribo < 50)

    ribo_genes <- grep("^RPS|^RPL", rownames(seurat_obj), value = T)
    mito_genes <- grep("^MT-", rownames(seurat_obj), value = T)

    seurat_obj <- seurat_obj %>%
        subset(features = setdiff(rownames(.), mito_genes)) %>%
        subset(features = setdiff(rownames(.), ribo_genes))
    
    pdf_path <- file.path(dir_out, "feature_vlnplot_filt.pdf")
    cairo_pdf(pdf_path, width = 15, height = 15)
    print(VlnPlot(
        seurat_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        ncol = 2,
        pt.size = 0,
        group.by = "dataset"
    ))
    dev.off()

    return(seurat_obj)
}

seurat_obj <- seurat_quality_control(seurat_obj, dir_out)

saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_2.rds"))