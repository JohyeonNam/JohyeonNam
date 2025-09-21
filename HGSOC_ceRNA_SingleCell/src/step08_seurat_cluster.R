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

jnam_selectelbow <- function(seurat_obj, dir_out) {
    print("Function jnam_selectelbow")

    seurat_obj@reductions$integrated.harmony@stdev = as.numeric(apply(seurat_obj@reductions$integrated.harmony@cell.embeddings, 2, stats::sd))

    stdev_harmony <- seurat_obj[["integrated.harmony"]]@stdev

    print(stdev_harmony)

    n <- length(stdev_harmony)
    x <- 1:n
    y <- stdev_harmony
    p1 <- c(1, y[1])
    p2 <- c(n, y[n])
    distances <- abs((p2[2] - p1[2]) * x - (p2[1] - p1[1]) * y + p2[1]*p1[2] - p2[2]*p1[1]) / sqrt((p2[2] - p1[2])^2 + (p2[1] - p1[1])^2)
    elbow <- which.max(distances)

    df1 <- data.frame(
        Component = x,
        Stdev = y,
        Distance = distances
    )

    elbow_plot <- ggplot(df1, aes(x = Component, y = Stdev)) +
        geom_line(color = "steelblue") +
        geom_point(color = "steelblue") +
        geom_point(
            data = df1[elbow, ], 
            aes(x = Component, y = Stdev), 
            color = "red", 
            size = 3
        ) +
        geom_text(data = df1[elbow, ],
            aes(x = Component, y = Stdev, label = paste("Elbow = ", Component)),
            vjust = -1,
            color = "red"
        ) +
        labs(title = "Elbow Plot of Harmony Components",
            x = "Harmony Component",
            y = "Standard Deviation"
        ) +
        theme(plot.background = element_rect(fill = "transparent"),
            text = element_text(family = "Helvetica", size = 30),
            panel.grid = element_line(color = "#EBEBEB"),
            panel.background = element_rect(fill = "transparent", colour = "Black"),
            axis.title = element_text(size = 30, family = "Helvetica"),
            legend.background = element_rect(colour = "white", fill = "transparent"),
            legend.title = element_text(size = 24),
            legend.key = element_rect(fill = "transparent", color = "white"),
            legend.key.size = unit(2, units = "cm"),
            legend.text = element_text(size = 18),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
        )

    pdf_path <- file.path(dir_out, "elbow_plot.pdf")
    cairo_pdf(pdf_path, width = 15, height = 15)
    print(elbow_plot)
    dev.off()

    var_explained <- stdev_harmony^2
    prop_var <- var_explained / sum(var_explained)
    cum_var <- cumsum(prop_var)

    df2 <- data.frame(
        PC = 1:length(stdev_harmony),
        CumulativeVariance = cum_var
    )

    pc_90 <- which(cum_var >= 0.9)[1]

    cum_plot <- ggplot(df2, aes(x = PC, y = CumulativeVariance)) +
        geom_line(color = "#0072B2", size = 1.2) +
        geom_point(color = "#D55E00", size = 2) +
        geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray50") +
        labs(title = "Cumulative Variance Explained by Harmony Components",
            x = "Harmony Component",
            y = "Cumulative Proportion of Variance Explained") +
        theme(plot.background = element_rect(fill = "transparent"),
            text = element_text(family = "Helvetica", size = 30),
            panel.grid = element_line(color = "#EBEBEB"),
            panel.background = element_rect(fill = "transparent", colour = "Black"),
            axis.title = element_text(size = 30, family = "Helvetica"),
            legend.background = element_rect(colour = "white", fill = "transparent"),
            legend.title = element_text(size = 24),
            legend.key = element_rect(fill = "transparent", color = "white"),
            legend.key.size = unit(2, units = "cm"),
            legend.text = element_text(size = 18),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
        ) 

    pdf_path <- file.path(dir_out, "cum_plot.pdf")
    cairo_pdf(pdf_path, width = 15, height = 15)
    print(cum_plot)
    dev.off()

    print(paste0("Elbow cut: ", elbow))
    print(paste0("Cumulative cut: ", pc_90))

    dim <- floor(median(elbow:pc_90))

    print(paste0("Dims for use: ", dim))
    return(dim)
}

seurat_clustercells <- function(seurat_obj, dir_out, dim) {
    print("Function seurat_clustercells")

    seurat_obj <- FindNeighbors(
        seurat_obj,
        reduction = "integrated.harmony",
        dims = 1:dim
    )

    seurat_obj <- FindClusters(
        seurat_obj,
        resolution = c(0.5, 0.8, 1.1, 1.4, 1.7, 2.0)
    )

    seurat_obj$seurat_clusters <- seurat_obj$sketch_snn_res.1.1

    return(seurat_obj)
}

clustree_drawclustree <- function(seurat_obj, dir_out) {
    print("Function clustree_drawclustree")

    sketch_cells <- WhichCells(seurat_obj, colnames(seurat_obj[["sketch"]]))
    seurat_sketch_only <- subset(seurat_obj, cells = sketch_cells)

    clustree_plot <- clustree(seurat_sketch_only, 
        prefix = "sketch_snn_res."
    )

    pdf_path <- file.path(dir_out, "clustree_plot.pdf")
    cairo_pdf(pdf_path, width = 30, height = 15)
    print(clustree_plot)
    dev.off()
}

dim <- jnam_selectelbow(seurat_obj, dir_out)
seurat_obj <- seurat_clustercells(seurat_obj, dir_out, dim)
clustree_drawclustree(seurat_obj, dir_out)

print(dim)
saveRDS(seurat_obj, file.path(dir_out, "seurat_obj_6.rds"))