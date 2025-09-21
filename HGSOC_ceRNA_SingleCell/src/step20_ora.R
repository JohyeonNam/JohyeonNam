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
allgenes_table <- read.csv(args[2])
wgcna_obj_cancer <- readRDS(args[3])

clusterprofiler_ora <- function(dir_out, allgenes_table, wgcna_obj_cancer) {
    print("clusterprofiler_ora")

    gene <- mapIds(org.Hs.eg.db, keys = allgenes_table$gene_name[allgenes_table$module == "Cancer9"], keytype = "SYMBOL", column = "ENTREZID")
    bg <- mapIds(org.Hs.eg.db, keys = rownames(wgcna_obj_cancer), keytype = "SYMBOL", column = "ENTREZID")

    kegg_ora <- enrichKEGG(
        gene = gene,
        organism = "hsa",
        pvalueCutoff = 0.05
    ) %>%
    setReadable("org.Hs.eg.db", "ENTREZID")

    go_ora <- enrichGO(
        gene = gene,
        universe = bg,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pvalueCutoff = 0.05,
    ) %>%
    setReadable("org.Hs.eg.db", "ENTREZID")

    reactome_ora <- enrichPathway(
        gene = gene,
        pvalueCutoff = 0.05,
        readable = TRUE
    ) %>%
    setReadable("org.Hs.eg.db", "ENTREZID")

    wiki_ora <- enrichWP(
        gene = gene,
        organism = "Homo sapiens"
    ) %>%
    setReadable("org.Hs.eg.db", "ENTREZID")

    make_df <- function(df) {
        out_df <- df %>%
        #filter(p.adjust < 0.05) %>%
        filter(pvalue < 0.05) %>%
        arrange(p.adjust) %>%
        #head(n = 20L) %>%
        mutate(Ratio = Count / as.numeric(sub("/.*", "", BgRatio))) %>%
        mutate(p.adjust = p.adjust) %>%
        mutate(Description = factor(Description, levels = Description[order(p.adjust)]))
        return(out_df)
    }

    go_df <- make_df(go_ora@result)
    kegg_df <- make_df(kegg_ora@result)   
    reactome_df <- make_df(reactome_ora@result)
    wiki_df <- make_df(wiki_ora@result)

    table_path1 <- file.path(dir_out, "go_ora.txt")
    write.table(go_df, table_path1, quote = FALSE, row.names = FALSE, sep = "\t")
    table_path2 <- file.path(dir_out, "kegg_ora.txt")
    write.table(kegg_df, table_path2, quote = FALSE, row.names = FALSE, sep = "\t")
    table_path3 <- file.path(dir_out, "reactome_ora.txt")
    write.table(reactome_df, table_path3, quote = FALSE, row.names = FALSE, sep = "\t")
    table_path4 <- file.path(dir_out, "wiki_ora.txt")
    write.table(wiki_df, table_path4, quote = FALSE, row.names = FALSE, sep = "\t")

    return(list(
        head(go_df, n = 20),
        head(kegg_df, n = 20),
        head(reactome_df, n = 20),
        head(wiki_df, n = 20)
    ))
}

result <- clusterprofiler_ora(dir_out, allgenes_table, wgcna_obj_cancer)
go_df <- result[[1]]
kegg_df <- result[[2]]
reactome_df <- result[[3]]
wiki_df <- result[[4]]

saveRDS(go_df, file.path(dir_out, "go_df.rds"))
saveRDS(kegg_df, file.path(dir_out, "kegg_df.rds"))
saveRDS(reactome_df, file.path(dir_out, "reactome_df.rds"))
saveRDS(wiki_df, file.path(dir_out, "wiki_df.rds"))