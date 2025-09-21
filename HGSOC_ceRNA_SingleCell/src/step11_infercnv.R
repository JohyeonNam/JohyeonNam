#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(dplyr)
library(Seurat)
library(infercnv)
library(tibble)
library(devtools)
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj <- readRDS(args[1])  # Seurat object
dir_out <- args[2]  # Output directory
gene_order <- read.table(args[3], row.names = 1) # Gene order file for inferCNV
celltype_anno <- read.csv(args[4], header = TRUE) # Celltype annotation result

cnv_inference <- function(seurat_obj, dir_out, gene_order, celltype_anno) {
    print("Function cnv_inference")
    
    options(scipen = 100)

    genes <- intersect(rownames(seurat_obj), rownames(gene_order))

    counts <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")[genes, ]
    idents <- as.data.frame(seurat_obj$seurat_clusters_full)
    ref_groups <- as.vector(unique(seurat_obj$seurat_clusters_full[seurat_obj$Celltype.label != "Epithelial_cells"]))

    print(all_of(rownames(counts) %in% rownames(gene_order)))

    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = counts,
        annotations_file = idents,
        gene_order_file = gene_order,
        ref_group_names = ref_groups,
        max_cells_per_group = 1000
    )

    infercnv_result <- infercnv::run(
        infercnv_obj,
        cutoff = 0.1,
        out_dir = dir_out,
        cluster_by_groups = TRUE,
        HMM = TRUE,
        denoise = TRUE,
        write_expr_matrix = TRUE,
        num_threads = 32,
        output_format = "pdf"
    )

    Idents(seurat_obj) <- "seurat_clusters_full"

    expr_data <- infercnv_result@expr.data
    normal_means <- rowMeans(expr_data[, unlist(infercnv_result@reference_grouped_cell_indices)])
    cnv_scores <- expr_data %>% 
        apply(2, function(cell_expr) {sum(abs(cell_expr - normal_means))}) %>%
        merge(y = as.data.frame(seurat_obj$seurat_clusters_full), all.x = T, by = 0) %>%
        column_to_rownames(var = 'Row.names') %>%
        merge(y = as.data.frame(seurat_obj$Celltype.label), all.x = T, by = 0) %>%
        setNames(c('ident', 'log2_CNVscore', 'cluster', 'celltype')) %>%
        mutate(log2_CNVscore = log2(log2_CNVscore))
    
    csv_path <- file.path(dir_out, "CNV_score.csv")
    write.csv(cnv_scores, csv_path, row.names = F)

    cnv_scores_cancer <- cnv_scores %>%
        mutate(cluster = ifelse(celltype == "Epithelial_cells", cluster, "Ref"))

    cnv_scores$celltype.plot <- cnv_scores$celltype %>% 
        gsub("T_cells", "T cells", .) %>%
        gsub("NK_cell", "NK cells", .) %>%
        gsub("Epithelial_cells", "Epithelial cells", .) %>%
        gsub("Macrophage", "Macrophages", .) %>%
        gsub("B_cell", "B cells", .) %>%
        gsub("Monocyte", "Monocytes", .) %>%
        gsub("Endothelial_cells", "Endothelial cells", .) %>%
        gsub("Smooth_muscle_cells", "SMCs", .) %>%
        gsub("iPS_cells", "iPS cells", .)

    plot1 <- ggplot(cnv_scores, 
            aes(x = celltype.plot, y = log2_CNVscore, fill = celltype.plot)
        ) +
        geom_violin(trim = TRUE) +
        labs(x = "", y = "log2 CNV score") +
        theme(
            text = element_text(
                family = "Arial",
                size = 10
            ),
            plot.background = element_rect(fill = 'transparent'), 
            panel.background = element_rect(fill = 'transparent'), 
            legend.background = element_rect(fill = 'transparent'), 
            axis.line = element_line(), 
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 10, hjust = 1, angle = 45), 
            legend.position = "none"
        )

    plot2 <- ggplot(cnv_scores_cancer, 
            aes(x = cluster, y = log2_CNVscore, fill = cluster)
        ) +
        geom_violin(trim = TRUE) +
        labs(title = "CNV score of cancer clusters", x = "Cancer clusters", y = "log2 CNV score") +
        theme(
            text = element_text(
                family = "Arial",
                size = 20
            ),
            plot.background = element_rect(fill = 'transparent'), 
            panel.background = element_rect(fill = 'transparent'), 
            legend.background = element_rect(fill = 'transparent'), 
            plot.margin = margin(0.5,0.25,0.5,0.5,'cm'), 
            axis.line = element_line(), 
            axis.text.x = element_text(size = 20, hjust = 1, angle = 45), 
            plot.title = element_text(size = 24, hjust = 0.5, face='bold'),
            legend.position = "none"
        )

    pdf_path1 <- file.path(dir_out, "cnv_celltypes.pdf")
    cairo_pdf(pdf_path1, width = 4, height = 2)
    print(plot1)
    dev.off()

    pdf_path2 <- file.path(dir_out, "cnv_cancerclusters.pdf")
    cairo_pdf(pdf_path2, width = 8, height = 4)
    print(plot2)
    dev.off()

    return(infercnv_result)
}

infercnv_result <- cnv_inference(seurat_obj, dir_out, gene_order, celltype_anno)