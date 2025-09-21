#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
dir_out <- args[1]  # Output directory
pred_df <- readRDS(args[2])

visualise_prediction <- function(pred_df, dir_out) {
    print("Function visualise_prediction")

    set.seed(981214)

    edges_long <- pred_df %>%
        pivot_longer(cols = c("lncRNA", "miRNA", "pcg"), names_to = "type", values_to = "gene") %>%
        group_by(module) %>%
        mutate(next_gene = lead(gene),
            next_type = lead(type)) %>%
        filter(!is.na(next_gene)) %>%
        ungroup()

    nodes <- edges_long %>%
        select(name = gene, type) %>%
        bind_rows(
            edges_long %>% select(name = next_gene, type = next_type)
        ) %>%
        distinct()

    g <- graph_from_data_frame(
        d = edges_long %>% select(gene, next_gene),
        vertices = nodes,
        directed = FALSE
    )

    plot1 <- ggraph(g, layout = 'fr') +
        geom_edge_link(color = "gray80", width = 1) +
        geom_node_point(aes(fill = type), color = "gray40", shape = 21, size = 10, stroke = 1.2) +
        geom_node_text(aes(label = name), repel = TRUE, size = 7) +
        theme_void() +
        scale_fill_manual(values = c(
            lncRNA = "#E41A1C",   
            miRNA = "#377EB8",    
            pcg   = "#4DAF4A"     
        ),
        name = "Gene Type"
    ) +
    theme(
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 28)
    )

    pdf_path <- file.path(dir_out, "Network.pdf")
    cairo_pdf(pdf_path, width = 15, height = 15)
    print(plot1)
    dev.off()
}

visualise_prediction(pred_df, dir_out)