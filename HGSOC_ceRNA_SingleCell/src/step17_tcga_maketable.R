#!/usr/bin/env Rscript

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: 

set.seed(231116)

library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
dir_out <- args[1]  # Output directory
dir_tcga <- args[2] # Directory containing TCGA counts
sample_sheet <- read.table(args[3], sep = "\t", header = TRUE)
clinical_sheet <- fread(args[4])
idnametype <- read.csv(args[5], header = TRUE) # Custom GTF file

tcga_maketable <- function(dir_tcga, dir_out, sample_sheet, clinical_sheet, idnametype) {
    print("Function tcga_maketable")

    clinical_sheet <- clinical_sheet %>%
        distinct(case_submitter_id, .keep_all = T) %>%
        mutate('case_submitter_id' = gsub('-', '_', .$case_submitter_id)) %>%
        mutate(OS = ifelse(.$vital_status == "Dead", .$days_to_death, .$days_to_last_follow_up)) %>%
        mutate(vital_status = ifelse(.$vital_status == 'Dead', 1, 0))

    idnametype <- idnametype %>%
        filter(!grepl("_.*$", GeneID))

    sample_sheet <- sample_sheet %>%
        filter(Sample.Type == "Primary Tumor") %>%
        mutate(path = paste0(dir_tcga, File.ID, "/", File.Name))

    miRNA_info <- sample_sheet %>%
        subset(Data.Type == "miRNA Expression Quantification") %>%
        mutate(command = paste0("cp ", path,  " ", dir_out, "miRNA/", gsub("-", "_", Sample.ID), ".txt")) %>%
        select(-File.ID, -File.Name, -Data.Category, -Data.Type, -Project.ID, -Case.ID, -Sample.Type, -path) %>%
        mutate(input = paste0(gsub("-", "_", Sample.ID), ".txt"))
    
    RNA_info <- sample_sheet %>%
        subset(Data.Type == "Gene Expression Quantification") %>%
        mutate(command = paste0("cp ", path,  " ", dir_out, "RNA/", gsub("-", "_", Sample.ID), ".txt")) %>%
        select(-File.ID, -File.Name, -Data.Category, -Data.Type, -Project.ID, -Case.ID, -Sample.Type, -path) %>%
        mutate(input = paste0(gsub("-", "_", Sample.ID), ".txt"))
    
    samples_to_include <- intersect(miRNA_info$input, RNA_info$input)

    print(paste0("Including samples with miRNA and RNA counts. Sample N = ", as.character(length(samples_to_include))))

    miRNA_info <- miRNA_info %>%
        filter(input %in% samples_to_include)
    
    RNA_info <- RNA_info %>%
        filter(input %in% samples_to_include)
    
    print(dim(miRNA_info))
    print(dim(RNA_info))
    system(paste0("mkdir ", dir_out, "miRNA"))
    system(paste0("mkdir ", dir_out, "RNA"))

    for (command in miRNA_info$command) {
        system(command)
    }
    for (command in RNA_info$command) {
        system(command)
    }

    count <- 1
    for (sample in RNA_info$input) {
        table <- read.table(paste0(dir_out, 'RNA/', sample), sep = '\t', col.names = c('gene_id', 'gene_name', 'gene_type', 'unstranded', 'stranded_first', 'stranded_second', 'tpm_unstranded',  'fpkm_unstranded', 'fpkm_uq_unstranded'), skip = 6)

        pcg <- table %>%
            filter(gene_type == "protein_coding") %>%
            select(gene_id, tpm_unstranded) %>%
            filter(!grepl("_.*$", gene_id)) %>%
            setNames(c("gene_id", gsub(".txt", "", sample))) %>%
            mutate(gene_id = gsub("\\..*", "", gene_id)) 

        lncRNA <- table %>%
            filter(gene_type == "lncRNA") %>%
            select(gene_id, tpm_unstranded) %>%
            filter(!grepl("_.*$", gene_id)) %>%
            setNames(c("gene_id", gsub(".txt", "", sample))) %>%
            mutate(gene_id = gsub("\\..*", "", gene_id))         

        print(head(pcg))
        print(head(lncRNA))

        if (count == 1) {
            pcg_count <- pcg
            lnc_count <- lncRNA
        }
        
        if (count != 1) {
            print(all(pcg_count$gene_id == pcg$gene_id))
            pcg_count <- cbind(pcg_count, select(pcg, gsub(".txt", "", sample)))
            lnc_count <- cbind(lnc_count, select(lncRNA, gsub(".txt", "", sample)))
        }
                
        print(paste0(as.character(count), " Processing ", gsub(".txt", "", sample)))
        count <- count + 1
    }

    print(head(pcg_count))
    print(head(lnc_count))

    pcg_count <- merge(pcg_count, idnametype[, c("GeneID", "GeneName")], by.x = "gene_id", by.y = "GeneID", all.x = T) %>%
        select(-gene_id) %>%
        select(GeneName, everything())
    colnames(pcg_count) <- substr(colnames(pcg_count), 1, 12)
        
    lnc_count <- merge(lnc_count, idnametype[, c("GeneID", "GeneName")], by.x = "gene_id", by.y = "GeneID", all.x = T) %>%
        select(-gene_id) %>%
        select(GeneName, everything())
    colnames(lnc_count) <- substr(colnames(lnc_count), 1, 12)

    count <- 1
    for (sample in miRNA_info$input) {
        table <- read.table(paste0(dir_out, "miRNA/", sample), sep = "\t", header = T)
        miRNA <- table %>%
            select(miRNA_ID, reads_per_million_miRNA_mapped) %>%
            setNames(c("GeneName", gsub(".txt", "", sample))) %>%
            mutate(GeneName = gsub("(-3p|-5p)$", "", GeneName))
        
        if (count == 1) {
            miRNA_count <- miRNA
        }
                
        if (count != 1) {
            print(all(miRNA_count$GeneName == miRNA$GeneName))
            miRNA_count <- cbind(miRNA_count, select(miRNA, gsub(".txt", "", sample)))
        }
                
        colnames(miRNA_count) <- substr(colnames(miRNA_count), 1, 12)
        print(paste0(as.character(count), " Processing ", gsub(".txt", "", sample)))
        count <- count + 1
    }

    write.csv(pcg_count, paste0(dir_out, "pcg_count.csv"), quote = F, row.names = F)
    write.csv(lnc_count, paste0(dir_out, "lnc_count.csv"), quote = F, row.names = F)
    write.csv(miRNA_count, paste0(dir_out, "miRNA_count.csv"), quote = F, row.names = F)
    return(list(pcg_count, lnc_count, miRNA_count))
}

tcga_maketable_result <- tcga_maketable(
    dir_tcga, dir_out, sample_sheet, clinical_sheet, idnametype)
pcg_count <- tcga_maketable_result[[1]]
lnc_count <- tcga_maketable_result[[2]]
miRNA_count <- tcga_maketable_result[[3]]

saveRDS(pcg_count, file.path(dir_out, "pcg_count.rds"))
saveRDS(lnc_count, file.path(dir_out, "lnc_count.rds"))
saveRDS(miRNA_count, file.path(dir_out, "miRNA_count.rds"))