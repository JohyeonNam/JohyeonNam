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
library(SingleR)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory
dim <- args[3] # Dimensions to use
res <- args[4] # Resolution

seurat_projection <- function(seurat_obj, dir_out, dim, res) {
    print("Function seurat_projection")

    seurat_obj <- ProjectIntegration(
        object = seurat_obj,
        sketched.assay = "sketch",
        assay = "RNA",
        reduction = "integrated.harmony"
    )

    seurat_obj <- ProjectData(
        object = seurat_obj,
        sketched.assay = "sketch",
        assay = "RNA",
        sketched.reduction = "integrated.harmony",
        full.reduction = "integrated.harmony.full",
        dims = 1:dim,
        refdata = list(seurat_clusters_full = "seurat_clusters")
    )

    DefaultAssay(seurat_obj) <- "RNA"

    seurat_obj <- JoinLayers(
        seurat_obj,
        assay = "RNA"
    )

    seurat_obj <- RunUMAP(
        seurat_obj,
        reduction = "integrated.harmony.full",
        dims = 1:dim,
        reduction.name = "umap.full"
    )

    plot1 <- DimPlot(
        object = seurat_obj,
        reduction = "umap.full",
        raster = FALSE,
        label = TRUE,
        group.by = "seurat_clusters_full",
        label.size = 2
    ) +
    coord_fixed() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(text = element_text(family = "Arial", size = 14),
        plot.title = element_blank()
    ) +
    NoLegend()

    plot2 <- DimPlot(
        object = seurat_obj,
        reduction = "umap.full",
        raster = FALSE,
        label = FALSE,
        group.by = "dataset" 
    ) +
    coord_fixed() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(text = element_text(family = "Arial", size = 36),
        plot.title = element_blank()
    )

    plot3 <- DimPlot(
        object = seurat_obj,
        reduction = "umap.full",
        raster = FALSE,
        label = FALSE,
        group.by = "orig.ident" 
    ) +
    coord_fixed() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(text = element_text(family = "Arial", size = 36),
        plot.title = element_blank()
    )

    pdf_path1 <- file.path(dir_out, "UMAP_clusters.pdf")
    cairo_pdf(pdf_path1, width = 4, height = 4)
    print(plot1)
    dev.off()

    pdf_path2 <- file.path(dir_out, "UMAP_patients_full.pdf")
    cairo_pdf(pdf_path2, width = 4, height = 4)
    print(plot2)
    dev.off()

    pdf_path3 <- file.path(dir_out, "UMAP_datasets_full.pdf")
    cairo_pdf(pdf_path3, width = 4, height = 4)
    print(plot3)
    dev.off()

    return(seurat_obj)
}

seurat_obj <- seurat_projection(seurat_obj, dir_out, dim, res)

saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_7.rds"))