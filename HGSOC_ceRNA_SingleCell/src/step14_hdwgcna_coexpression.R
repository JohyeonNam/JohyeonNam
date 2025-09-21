#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(Seurat)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory
pcg_markers <- readRDS(args[3])
lnc_markers <- readRDS(args[4])

hdwgcna_coexp <- function(seurat_obj, dir_out, pcg_markers, lnc_markers) {
    print("Function hdwgcna_coexp")
    
    set.seed(231116)
    
    marker_genes <- c(pcg_markers$gene[pcg_markers$cluster == "Epithelial_cells"], lnc_markers$gene[lnc_markers$cluster == "Epithelial_cells"])

    wgcna_obj_cancer <- SetupForWGCNA(
        seurat_obj,
        wgcna_name = "wgcna_cancer",
        gene_select = "custom",
        gene_list = marker_genes
    )

    wgcna_obj_cancer <- MetacellsByGroups(
        seurat_obj = wgcna_obj_cancer,
        group.by = c("Celltype.label", "orig.ident"),
        reduction = "integrated.harmony.full",
        ident.group = "orig.ident"
    )

    wgcna_obj_cancer <- NormalizeMetacells(wgcna_obj_cancer)

    wgcna_obj_cancer <- SetDatExpr(
        wgcna_obj_cancer,
        group_name = "Epithelial_cells",
        group.by = "Celltype.label"
    )

    wgcna_obj_cancer <- TestSoftPowers(
        wgcna_obj_cancer,
        networkType = 'signed'
    )

    plot_list <- PlotSoftPowers(wgcna_obj_cancer)

    plot1 <- wrap_plots(plot_list, ncol = 2)

    wgcna_obj_cancer <- ConstructNetwork(
        wgcna_obj_cancer,
        tom_name = "Cancer Cells",
        overwrite_tom = TRUE
    )

    wgcna_obj_cancer <- ScaleData(
        wgcna_obj_cancer, 
        features = VariableFeatures(wgcna_obj_cancer)
    )

    wgcna_obj_cancer <- ModuleEigengenes(
        wgcna_obj_cancer,
        group.by.vars = "orig.ident"
    )

    wgcna_obj_cancer <- ModuleConnectivity(
        wgcna_obj_cancer,
        group.by = "Celltype.label", 
        group_name = "Epithelial_cells"
    )

    wgcna_obj_cancer <- ResetModuleNames(
        wgcna_obj_cancer,
        new_name = "Cancer"
    )

    plot2 <- PlotKMEs(
        wgcna_obj_cancer, 
        ncol = 4,
        text_size = 2.4,
        n_hubs = 15
        ) 

    MEs <- GetMEs(wgcna_obj_cancer, harmonized = TRUE)
    modules <- GetModules(wgcna_obj_cancer)
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
    wgcna_obj_cancer@meta.data <- cbind(wgcna_obj_cancer@meta.data, MEs)

    plot3 <- DotPlot(
        wgcna_obj_cancer,
        features = mods,
        group.by = "Celltype.label",
        dot.scale = 12
    ) +
    coord_flip() +
    scale_color_gradient2(high='red', mid='grey95', low='blue') +
    theme(
        axis.text.x = element_text(size = 24, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 24),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 24),
        text = element_text(family = "Arial", size = 36),
        panel.background = element_rect(fill = "transparent", colour = "Black")
    )

    pdf_path1 <- file.path(dir_out, "Softpower.pdf")
    cairo_pdf(pdf_path1, width = 15, height = 15)
    print(plot1)
    dev.off()

    pdf_path2 <- file.path(dir_out, "Dendrogram.pdf")
    cairo_pdf(pdf_path2, width = 12, height = 5)
    PlotDendrogram(wgcna_obj_cancer, main='Cancer Cells hdWGCNA Dendrogram')
    dev.off()

    pdf_path3 <- file.path(dir_out, "KMEs.pdf")
    cairo_pdf(pdf_path3, width = 10, height = 15)
    print(plot2)
    dev.off()

    pdf_path4 <- file.path(dir_out, "ExpDot.pdf")
    cairo_pdf(pdf_path4, width = 10, height = 20)
    print(plot3)
    dev.off()

    hubgenes_table <- wgcna_obj_cancer %>% 
        GetHubGenes(n_hubs = 15) %>%
        mutate(GeneType = ifelse(gene_name %in% pcg_markers$gene[pcg_markers$cluster == "Epithelial_cells"], "pcg", "lncRNA"))

    allgenes_table <- wgcna_obj_cancer %>% 
        GetHubGenes(n_hubs = Inf) %>%
        mutate(GeneType = ifelse(gene_name %in% pcg_markers$gene[pcg_markers$cluster == "Epithelial_cells"], "pcg", "lncRNA"))

    csv_path1 <- file.path(dir_out, "hubgenes_table.csv")
    write.csv(hubgenes_table, csv_path1, quote = FALSE, row.names = FALSE)
    csv_path2 <- file.path(dir_out, "allgenes_table.csv")
    write.csv(allgenes_table, csv_path2, quote = FALSE, row.names = FALSE)

    return(list(wgcna_obj_cancer, hubgenes_table, allgenes_table))
}

result <- hdwgcna_coexp(seurat_obj, dir_out, pcg_markers, lnc_markers)
wgcna_obj_cancer <- result[[1]]
hubgenes_table <- result[[2]]
allgenes_table <- result[[3]]

saveRDS(wgcna_obj_cancer, file.path(dir_out, "wgcna_obj_cancer.rds"))
saveRDS(hubgenes_table, file.path(dir_out, "hubgenes_table.rds"))
saveRDS(allgenes_table, file.path(dir_out, "allgenes_table.rds"))