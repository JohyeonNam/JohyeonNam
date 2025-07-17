#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: Running edgeR

library(devtools)
library(dplyr)
library(edgeR)
library(genefilter)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
dir_in <- args[1]  # featureCounts directory
file_geneinfo <- args[2]  # custom_geneinfo file
file_sampleinfo <- args[3]  # sampleinfo file
dir_out <- args[4]  # output directory

geneinfo    <- read.csv(file_geneinfo, header = TRUE)
sampleinfo  <- read.csv(file_sampleinfo, header = TRUE)

coldata <- data.frame(
  sample = sampleinfo$SampleID,
  condition = sampleinfo$cond
) %>% column_to_rownames(var = "sample")

extract_counts <- function(sample, dir) {
  path <- paste0(dir, sample, ".txt")
  df <- read.table(
    path, skip = 1, header = TRUE,
    col.names = c("Geneid", "Chr", "Start", "End", "Strand", "Length", sample),
    quote = ""
  )
  df %>% select(Geneid, all_of(sample))
}

cts <- extract_counts(rownames(coldata)[1], dir_in)
for (sample in rownames(coldata)[-1]) {
  fcnt <- extract_counts(sample, dir_in)
  cts <- merge(cts, fcnt, by = "Geneid")
}
cts <- column_to_rownames(cts, var = "Geneid")

stopifnot(all(rownames(coldata) %in% colnames(cts)))
stopifnot(all(rownames(coldata) == colnames(cts)))

dge <- DGEList(counts = cts, group = coldata$condition)
design <- model.matrix(~dge$samples$group)
dge$samples$group <- relevel(dge$samples$group, ref = "ctrl")
rownames(design) <- rownames(dge$samples)

dge_disp <- estimateGLMCommonDisp(dge, design)
dge_tagw <- estimateGLMTagwiseDisp(dge_disp, design)
fit <- glmFit(dge_tagw, design, dispersion = dge_tagw$tagwise.dispersion)
lrt <- glmLRT(fit)

norm <- cpm(dge_tagw, normalized.lib.sizes = TRUE)
norm <- norm[rowSums(norm == 0) == 0, ]

res <- topTags(lrt, n = Inf)$table
res$GeneID <- rownames(res)
res <- merge(res, norm, by.x = "GeneID", by.y = "row.names") %>%
  arrange(PValue) %>%
  relocate(logFC, .after = FDR)

output <- merge(geneinfo, res, by = "GeneID", all.y = TRUE)

groups <- unique(coldata$condition)
if (length(groups) != 2) stop("Condition must have exactly two groups.")

ctrl_group <- "ctrl"
case_group <- setdiff(groups, ctrl_group)
if (length(case_group) != 1) stop("Unable to determine case group.")

output_filename <- paste0(case_group, "_vs_", ctrl_group, "_edgeR_results.csv")

write.csv(output, file = file.path(dir_out, output_filename), row.names = FALSE, quote = FALSE)