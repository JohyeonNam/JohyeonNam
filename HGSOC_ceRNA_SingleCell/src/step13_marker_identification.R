#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(dplyr)
library(Seurat)
library(presto)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory
idnametype <- read.csv(args[3], header = TRUE) # Custom GTF file

marker_identification <- function(seurat_obj, dir_out, idnametype) {
    print("Function marker_identification")
    
    Idents(seurat_obj) <- "Celltype.label"
    
    all_genes <- rownames(seurat_obj)
    pcg <- all_genes[all_genes %in% idnametype$GeneName[idnametype$GeneType == "protein_coding"]]
    lncRNA <- all_genes[all_genes %in% idnametype$GeneName[idnametype$GeneType == "lncRNA"]]

    pcg_markers <- FindAllMarkers(
        object = seurat_obj, 
        assay = "RNA",
        only.pos = TRUE, 
        logfc.threshold = 1, 
        return.thresh = 0.01,
        features = pcg
    )

    lnc_markers <- FindAllMarkers(
        object = seurat_obj, 
        assay = "RNA",
        only.pos = TRUE, 
        logfc.threshold = 1, 
        return.thresh = 0.01,
        features = lncRNA
    )

    csv_path1 <- file.path(dir_out, "pcg_markers.csv")
    write.csv(data.frame(pcg_markers), csv_path1, quote = F)
    csv_path2 <- file.path(dir_out, "lnc_markers.csv")
    write.csv(data.frame(lnc_markers), csv_path2, quote = F)

    return(list(pcg_markers, lnc_markers, seurat_obj))
}

result2 <- marker_identification(seurat_obj, dir_out, idnametype)
pcg_markers <- result2[[1]]
lnc_markers <- result2[[2]]
seurat_obj <- result2[[3]]

saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_10.rds"))
saveRDS(pcg_markers, file.path(dir_out, "pcg_markers.rds"))
saveRDS(lnc_markers, file.path(dir_out, "lnc_markers.rds"))