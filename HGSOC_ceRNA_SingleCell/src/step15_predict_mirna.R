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
hubgenes_table <- readRDS(args[1])  # Hub genes from hdWGCNA
dir_out <- args[2]  # Output directory
mir_lnc <- read.csv(args[3], header = TRUE)
mir_gene1 <- read.csv(args[4], header = TRUE)
mir_gene2 <- read.csv(args[5], header = TRUE)
mir_gene3 <- read.csv(args[6], header = TRUE)

predict_mirna <- function(
    hubgenes_table, dir_out, mir_lnc, mir_gene1, mir_gene2, mir_gene3) {
    print("Function predict_mirna")

    set.seed(231116)
    
    mir_lnc <- mir_lnc %>%
        select(c("mir_id", "symbol"))

    mir_gene <- rbind(mir_gene1, mir_gene2, mir_gene3) %>%
        select(c("mir_id", "symbol")) %>%
        unique()
        
    cor_df <- data.frame(
        lncRNA = character(), 
        pcg = character(), 
        module = character()
    )

    for (mod in unique(hubgenes_table$module)){
        module_df <- subset(hubgenes_table, module == mod)

        if ("lncRNA" %in% unique(module_df$GeneType) & "pcg" %in% unique(module_df$GeneType)) {
            print(paste0("Processing module ", mod))
            print(module_df)

            lncRNA <- module_df$gene_name[module_df$GeneType == "lncRNA"]
            pcg <- module_df$gene_name[module_df$GeneType == "pcg"]

            for(i in 1:length(lncRNA)) {
                for (j in 1:length(pcg)) {
                    df <- data.frame(lncRNA = lncRNA[i], pcg = pcg[j], module = mod)
                    cor_df <- rbind(cor_df, df)
                }
            }
        }
    }

    csv_path1 <- file.path(dir_out, "possible_pairs.csv")
    write.csv(cor_df, csv_path1, quote = FALSE, row.names = FALSE)  

    lncRNA_miRNA_pcg <- data.frame(
        lncRNA = character(), 
        miRNA = character(), 
        pcg = character(), 
        module = character()
    )

    for (i in 1:dim(cor_df)[1]) {
        pair <- cor_df[i, 1:2]
        mod <- cor_df[i, 3]
        
        lncRNA_miRNA <- mir_lnc$mir_id[mir_lnc$symbol == pair[1]$lncRNA]
        gene_miRNA <- mir_gene$mir_id[mir_gene$symbol == pair[2]$pcg]

        
        if (length(intersect(lncRNA_miRNA, gene_miRNA)) != 0) {
            interaction <- data.frame(
                lncRNA = pair[1]$lncRNA, 
                miRNA = intersect(lncRNA_miRNA, gene_miRNA), 
                pcg = pair[2]$pcg,
                module = mod
            )

        lncRNA_miRNA_pcg <- rbind(lncRNA_miRNA_pcg, interaction)
        }
    }

    csv_path2 <- file.path(dir_out, "interaction_prediction.csv")
    write.csv(lncRNA_miRNA_pcg, csv_path2, quote = FALSE, row.names = FALSE)

    miRNA_per_pair <- lncRNA_miRNA_pcg %>%
        group_by(lncRNA, pcg) %>%
        summarise(miRNA_count = n(), .groups = "drop")

    csv_path3 <- file.path(dir_out, "miRNA_per_pair.csv")
    write.csv(miRNA_per_pair, csv_path3, quote = FALSE, row.names = FALSE)

    return(list(lncRNA_miRNA_pcg, cor_df, miRNA_per_pair))
}

result <- predict_mirna(
    hubgenes_table, dir_out, mir_lnc, mir_gene1, mir_gene2, mir_gene3)
pred_cancer <- result[[1]]
possible_pairs <- result[[2]]
miRNA_per_pair <- result[[3]]

print(paste0("Unique mRNA-lncRNA pairs having shared miRNA targets: ", dim(unique(select(pred_cancer, lncRNA, pcg)))[1]))
print(paste0("Mean miRNA counts per pair: ", mean(miRNA_per_pair$miRNA_count)))

saveRDS(pred_cancer, file.path(dir_out, "pred_cancer.rds"))