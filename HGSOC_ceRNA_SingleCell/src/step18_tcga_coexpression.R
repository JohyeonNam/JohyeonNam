#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(dplyr)
library(tibble)
library(data.table)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
dir_out <- args[1]  # Output directory
dir_tcga <- args[2] # Directory containing TCGA counts
pcg_count <- readRDS(args[3])
lnc_count <- readRDS(args[4])
miRNA_count <- readRDS(args[5])
pred_cancer <- readRDS(args[6])

tcga_coexpression <- function(dir_out, pcg_count, lnc_count, miRNA_count, pred_df) {
    print("Function tcga_coexpression")

    val_df <- data.frame(lncRNA = character(), 
        miRNA = character(), 
        pcg = character(), 
        cor_lnc_mir = numeric(),
        p_lnc_mir = character(),
        cor_mir_pcg = numeric(),
        p_mir_pcg = character(),
        cor_lnc_pcg = numeric(),
        p_lnc_pcg = character(),
        module = character()
    )
    
    for (i in 1:dim(pred_df)[1]) {
        print(i)
        
        lncRNA <- pred_df[i, 1]
        miRNA <- gsub("(-3p|-5p)$", "", pred_df[i, 2])
        miRNA_full <- pred_df[i, 2]
        pcg <- pred_df[i, 3]
        mod <- pred_df[i, 4]

        print(lncRNA)
        print(miRNA)
        print(pcg)
        print(mod)

        if (lncRNA %in% lnc_count$GeneName & pcg %in% pcg_count$GeneName & miRNA %in% miRNA_count$GeneName) {
            lnc_df <- lnc_count %>%
                filter(GeneName == lncRNA) %>%
                column_to_rownames(var = "GeneName") %>%
                t() %>%
                as.data.frame() %>%
                setnames(lncRNA, "lncRNA") %>%
                rownames_to_column(var = "SampleID")

            pcg_df <- pcg_count %>%
                filter(GeneName == pcg) %>%
                column_to_rownames(var = "GeneName") %>%
                t() %>%
                as.data.frame() %>%
                setnames(pcg, "pcg") %>%
                rownames_to_column(var = "SampleID")

            miRNA_df <- miRNA_count %>%
                filter(GeneName == miRNA) %>%
                column_to_rownames(var = "GeneName") %>%
                t() %>%
                as.data.frame() %>%
                setnames(miRNA, "miRNA") %>%
                rownames_to_column(var = "SampleID")
            
            merged_df <- merge(x = lnc_df, y = miRNA_df, by = "SampleID") %>%
                merge(y = pcg_df, by = "SampleID")

            test1 <- cor.test(merged_df$lncRNA, merged_df$miRNA, method = "spearman")
            test2 <- cor.test(merged_df$miRNA, merged_df$pcg, method = "spearman")
            test3 <- cor.test(merged_df$lncRNA, merged_df$pcg, method = "spearman")

            summary_df <- data.frame(
                lncRNA = lncRNA, 
                miRNA = miRNA_full, 
                pcg = pcg,
                cor_lnc_mir = test1$estimate,
                p_lnc_mir = test1$p.value,
                cor_mir_pcg = test2$estimate,
                p_mir_pcg = test2$p.value,
                cor_lnc_pcg = test3$estimate,
                p_lnc_pcg = test3$p.value,
                module = mod
            )

            if (sd(merged_df$lncRNA) != 0 & sd(merged_df$pcg) != 0 & sd(merged_df$miRNA) != 0) {
                val_df <- rbind(val_df, summary_df)
            }
        }
    }
    val_df$fdr_lnc_mir <- p.adjust(val_df$p_lnc_mir, method = "BH")
    val_df$fdr_lnc_pcg <- p.adjust(val_df$p_lnc_pcg, method = "BH")
    val_df$fdr_mir_pcg <- p.adjust(val_df$p_mir_pcg, method = "BH")

    return(val_df)
}

filter_cerna <- function(dir_out, val_df, pcut, fdrcut, coefcut) {
    print("Function filter_cerna")

    fdr_sig <- val_df %>%
        filter(fdr_lnc_pcg < fdrcut) %>%
        filter(fdr_mir_pcg < fdrcut) %>%
        filter(fdr_lnc_mir < fdrcut)

    fdr_coef_sig <- fdr_sig %>%
        filter(cor_lnc_pcg > coefcut) %>%
        filter(abs(cor_mir_pcg) > coefcut) %>%
        filter(abs(cor_lnc_mir) > coefcut) %>%
        filter(abs(cor_lnc_pcg) > coefcut) %>%
        filter(cor_lnc_mir < 0 & cor_mir_pcg < 0)

    csv_path1 <- file.path(dir_out, "tcga_fdr_sig.csv")
    write.csv(fdr_sig, csv_path1, quote = FALSE, row.names = FALSE)
    csv_path2 <- file.path(dir_out, "tcga_fdr_coef_sig.csv")
    write.csv(fdr_coef_sig, csv_path2, quote = FALSE, row.names = FALSE)
    
    return(list(fdr_sig, fdr_coef_sig))
}

val_cancer <- tcga_coexpression(dir_out, pcg_count, lnc_count, miRNA_count, pred_cancer)
result <- filter_cerna(dir_out, val_cancer, 0.05, 0.05, 0.3)
fdr_sig <- result[[1]]
fdr_coef_sig <- result[[2]]

saveRDS(val_cancer, file.path(dir_out, "val_cancer.rds"))
saveRDS(fdr_sig, file.path(dir_out, "fdr_sig.rds"))
saveRDS(fdr_coef_sig, file.path(dir_out, "fdr_coef_sig.rds"))