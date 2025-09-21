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
library(tibble)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory

singler_annotation <- function(seurat_obj, dir_out) {
    print("Function singler_annotation")

    seurat_RNA_counts <- GetAssayData(
        seurat_obj,
        layer = "data",
        assay = "RNA"
    )

    hpca_se <- celldex::HumanPrimaryCellAtlasData()

    pred_counts <- SingleR(
        test = seurat_RNA_counts,
        ref = hpca_se,
        assay.type.test = 1,
        labels = hpca_se$label.main,
        clusters = seurat_obj$seurat_clusters_full,
        num.threads = 32
    )

    pred_counts <- pred_counts %>%
        as.data.frame() %>%
        rownames_to_column(var = "clusters")

    print(table(pred_counts$labels))

    seurat_obj[["Celltype.label"]] <- pred_counts$labels[match(seurat_obj$seurat_clusters_full, pred_counts$clusters)]

    csv_path <- file.path(dir_out, "Celltype_annotation_results.csv")
    write.csv(pred_counts, csv_path, quote = FALSE, row.names = FALSE)

    seurat_obj$Celltype.label.plot <- seurat_obj$Celltype.label %>% 
        gsub("T_cells", "T cells", .) %>%
        gsub("NK_cell", "NK cells", .) %>%
        gsub("Epithelial_cells", "Epithelial\ncells", .) %>%
        gsub("Macrophage", "Macrophages", .) %>%
        gsub("B_cell", "B cells", .) %>%
        gsub("Monocyte", "Monocytes", .) %>%
        gsub("Endothelial_cells", "Endothelial\ncells", .) %>%
        gsub("Smooth_muscle_cells", "SMCs", .) %>%
        gsub("iPS_cells", "iPS cells", .)

    plot1 <- DimPlot(
        object = seurat_obj,
        reduction = "umap.full",
        raster = FALSE,
        label = TRUE,
        group.by = "Celltype.label.plot",
        label.size = 2.4
    ) +
    coord_fixed() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(text = element_text(
        family = "Arial",
        size = 14),
        plot.title = element_blank()
    ) + 
    NoLegend()

    plot2 <- FeaturePlot(
        seurat_obj,
        features = "EPCAM",
        raster = FALSE
    ) +
    coord_fixed() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        text = element_text(family = "Arial", size = 24), plot.title = element_text(family = "Arial", size = 28),
        legend.text = element_text(family = "Arial", size = 12)
    )

    plot3 <- FeaturePlot(
        seurat_obj,
        features = "PAX8",
        raster = FALSE
    ) +
    coord_fixed() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        text = element_text(family = "Arial", size = 24), plot.title = element_text(family = "Arial", size = 28),
        legend.text = element_text(family = "Arial", size = 12)
    )

    pdf_path1 <- file.path(dir_out, "UMAP_celltypes.pdf")
    cairo_pdf(pdf_path1, width = 4, height = 4)
    print(plot1)
    dev.off()

    pdf_path2 <- file.path(dir_out, "UMAP_EPCAM.pdf")
    cairo_pdf(pdf_path2, width = 4, height = 4)
    print(plot2)
    dev.off()

    pdf_path3 <- file.path(dir_out, "UMAP_PAX8.pdf")
    cairo_pdf(pdf_path3, width = 4, height = 4)
    print(plot3)
    dev.off()

    return(seurat_obj)
}

seurat_obj <- singler_annotation(seurat_obj, dir_out)

saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_8.rds"))