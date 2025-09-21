#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: Making Seurat object

set.seed(231116)

library(dplyr)
library(tibble)
library(data.table)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
dir_in <- args[1]  # CellRanger directory
dir_out <- args[2]  # Output directory
file_sampleinfo <- args[3]  # Samplesheet file

seurat_make_object <- function(dir_in, dir_out, file_sampleinfo) {
    message("Function seurat_make_object")

    seurat_objects <- list()

    for (sample in file_sampleinfo$Patient_ID) {
        message(paste0("Working on Sample ", sample))

        data_dir <- file.path(dir_in, sample, "outs", "filtered_feature_bc_matrix")
        data <- Read10X(data.dir = data_dir)

        seurat_objects[[sample]] <- CreateSeuratObject(
            counts = data,
            project = sample,
            min.cells = 10,
            min.features = 250
        )
    }

    data_merged <- merge(
        seurat_objects[[1]],
        y = seurat_objects[-1],
        add.cell.ids = names(seurat_objects)
    )

    data_merged[["percent.mt"]] <- PercentageFeatureSet(data_merged, pattern = "^MT-")
    data_merged[["percent.ribo"]] <- PercentageFeatureSet(data_merged, pattern = "^RP[SL]")

    data_merged[["dataset"]] <- file_sampleinfo$Source[match(data_merged$orig.ident, file_sampleinfo$Patient_ID)]
    data_merged[["state"]] <- file_sampleinfo$State[match(data_merged$orig.ident, file_sampleinfo$Patient_ID)]
    data_merged[["orig.ident"]] <- substr(data_merged$orig.ident, 1, 10)

    pdf_path <- file.path(dir_out, "feature_vlnplot.pdf")
    cairo_pdf(pdf_path, width = 15, height = 15)
    print(VlnPlot(
        data_merged,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        ncol = 2,
        pt.size = 0,
        group.by = "dataset"
    ))
    dev.off()

    return(data_merged)
}

seurat_obj <- seurat_make_object(dir_in, dir_out, file_sampleinfo)

saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_1.rds"))