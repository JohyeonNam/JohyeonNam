#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(scales)
library(data.table)
library(ReactomePA)

args <- commandArgs(trailingOnly = TRUE)
dir_out <- args[1]  # Output directory
go_df <- readRDS(args[2])
kegg_df <- readRDS(args[3])
reactome_df <- readRDS(args[4])
wiki_df <- readRDS(args[5])

plot_ora <- function(dir_out, go_df, kegg_df, reactome_df, wiki_df) {
    print("Function plot_ora")

    dfs <- list(
        GO = go_df, 
        KEGG = kegg_df, 
        Reactome = reactome_df, 
        WikiPathways = wiki_df
    )
    cols <- c(
        GO = "steelblue",
        KEGG = "forestgreen",
        WikiPathways = "firebrick",
        Reactome = "darkorchid"
    )
    num <- 1

    for (name in names(dfs)) {
        col <- cols[[name]]
        df <- dfs[[name]]

        df <- df %>%
            arrange(desc(p.adjust))
        df$Description <- factor(df$Description, levels = df$Description)

        if (max(-log10(df$p.adjust)) < 1.30103) {
            ratio <- (1.30103 - min(-log10(df$p.adjust))) / 20
        } else if (min(-log10(df$p.adjust)) > 1.30103) {
            ratio <- (max(-log10(df$p.adjust)) - 1.30103) / 20
        } else {
            ratio <- (max(-log10(df$p.adjust)) - min(-log10(df$p.adjust))) / 20
        }

        dot <- ggplot(df, aes(x = -log10(p.adjust), y = Description)) +
            geom_point(aes(size = Ratio*2), color = "black", fill = col, shape = 21) +
            scale_size(range = c(4, 14)) +
            geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "darkgrey", size = 1) +
            theme(
                plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
                panel.grid = element_line(color = "#EBEBEB"),
                panel.background = element_rect(fill = "transparent", colour = "Black"),
                legend.background = element_rect(colour = "white", fill = "transparent"),
                axis.title = element_text(size = 36),
                axis.text = element_text(size = 32),
                legend.title = element_text(size = 28),
                legend.text = element_text(size = 24)
            ) +
            scale_y_discrete(labels = label_wrap(100)) +
            coord_fixed(ratio = ratio) +
            labs(
            x = "-log10(p.adjust)",
            y = "",
            size = "Ratio"
            )
        
        pdf_path <- file.path(dir_out, paste0(name, "_dotplot.pdf"))
        cairo_pdf(pdf_path, width = 40, height = 15)
        print(dot)
        dev.off()

        num <- num + 1
    }
}

plot_ora(dir_out, go_df, kegg_df, reactome_df, wiki_df)