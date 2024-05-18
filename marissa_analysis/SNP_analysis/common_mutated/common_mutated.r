#***********************************
# Author: Marissa Esteban
# Date: July 24 2023
# Title: SNP analysis pt4 -- looking for the commonly mutated genes among all patients using the files in */top_mutated_filtered_snps
#                                                                                                     and */top_mutated_all_snps
#
#***********************************

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)

all_snps <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/merged_snp_annotations.rds")

# common genes with SNPs for ALL SNPs (including introns) that pass MQ > 40.0
s2_all <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_all_snps/S2_top_mutated_all.rds")
s6_all <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_all_snps/S6_top_mutated_all.rds")
s7_all <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_all_snps/S7_top_mutated_all.rds")
s9_all <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_all_snps/S9_top_mutated_all.rds")
s21_all <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_all_snps/S21_top_mutated_all.rds")
s32_all <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_all_snps/S32_top_mutated_all.rds")
s45_all <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_all_snps/S45_top_mutated_all.rds")

common_gene_ids_all <- Reduce(intersect, list(s2_all$gene_ID, s6_all$gene_ID, s7_all$gene_ID, s9_all$gene_ID, s21_all$gene_ID, s32_all$gene_ID, s45_all$gene_ID))
common_gene_names_all <- Reduce(intersect, list(s2_all$gene_name, s6_all$gene_name, s7_all$gene_name, s9_all$gene_name, s21_all$gene_name, s32_all$gene_name, s45_all$gene_name))
common_gene_names_all <- common_gene_names_all[common_gene_names_all != ""]


# common genes with SNPs for SNPs only in encoding regions that pass MQ > 40.0
s2_filtered <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_encoding_snps/S2_top_mutated_encoding.rds")
s6_filtered <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_encoding_snps/S6_top_mutated_encoding.rds")
s7_filtered <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_encoding_snps/S7_top_mutated_encoding.rds")
s9_filtered <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_encoding_snps/S9_top_mutated_encoding.rds")
s21_filtered <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_encoding_snps/S21_top_mutated_encoding.rds")
s32_filtered <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_encoding_snps/S32_top_mutated_encoding.rds")
s45_filtered <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_encoding_snps/S45_top_mutated_encoding.rds")

common_gene_ids_encoding <- Reduce(intersect, list(s2_filtered$gene_ID, s6_filtered$gene_ID, s7_filtered$gene_ID, s9_filtered$gene_ID, s21_filtered$gene_ID, s32_filtered$gene_ID, s45_filtered$gene_ID))
common_gene_names_encoding <- Reduce(intersect, list(s2_filtered$gene_name, s6_filtered$gene_name, s7_filtered$gene_name, s9_filtered$gene_name, s21_filtered$gene_name, s32_filtered$gene_name, s45_filtered$gene_name))




# ---------------------------------- Making a df with every variant associated with common mutated genes

# making a dataframe with all of the mutations pointing to the common genes between all patients
c_header <- c(colnames(all_snps))
all_common_genes <- data.frame(matrix(nrow = 0, ncol = length(c_header)))
colnames(all_common_genes) <- c_header

all_common_genes <- all_snps %>%
  filter(str_detect(gene_name, paste(common_gene_names_all, collapse = "|"))) %>%
  mutate(common_gene_match = str_extract(gene_name, paste(common_gene_names_all, collapse = "|", sep = "|"))) %>%
  select(common_gene_match, everything())


output_path <- "/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/common_mutated/"

saveRDS(all_common_genes, paste0(output_path, "variants_common_genes_all"))


# ---------------------------------- Now separating by classification (inv vs non inv)

#non_inv <- c("S2")
#inv <- c("S6", "S7", "S9", "S32", "S45") # NOT INCLUDING S21 BECAUSE THIS IS UNK

common_gene_names_inv <- Reduce(intersect, list(s6_all$gene_name, s7_all$gene_name, s9_all$gene_name, s32_all$gene_name, s45_all$gene_name))
common_gene_names_inv <- common_gene_names_inv[common_gene_names_inv != ""]

final_common_values <- setdiff(common_gene_names_inv, s2_all$gene_name) # mutated genes that are unique to invasive group

other_gene_names <- c(s6_all$gene_name, s7_all$gene_name, s9_all$gene_name, s32_all$gene_name, s45_all$gene_name)
unique_to_ninv <- setdiff(s2_all$gene_name, other_gene_names) # Find unique values in s2_all$gene_name that are not present in any of the other data frames

common_genes_inv_df <- all_snps %>%
  filter(grepl(paste(final_common_values, collapse = "|"), gene_name))

common_genes_ninv_df <- all_snps %>%
  filter(grepl(paste(unique_to_ninv, collapse = "|"), gene_name))

output_path <- "/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/common_mutated/"

saveRDS(common_genes_inv_df, paste0(output_path, "variants_common_genes_inv"))
saveRDS(common_genes_ninv_df, paste0(output_path, "variants_common_genes_non_inv"))




# ------------------------ Creating a df of a count matrix of each gene common gene hit by x variants per each sample

gene_variant_counts <- data.frame(Patient = character(0))

# this list is based from brief look through literature
common_cancer_related <- c("KCNJ18", "IRF4", "CRLF2", "MIR663AHG", "CNTNAP2", "ADARB2",
                           "RBMS3", "DUSP22", "KIF26B", "PTPRN2", "NRXN3", "CSMD1", "CTNNA3",
                           "FHIT")

# Iterate through each common gene name
for (common_gene in common_gene_names_all) {
  # Count how many rows match the current gene for each patient
  gene_counts <- all_snps %>%
    filter(str_detect(gene_name, common_gene)) %>%
    group_by(Patient, Chromosome) %>%
    summarise(Variant_Count = n()) %>%
    ungroup()
  
  # Add a new column with the current gene name
  gene_counts$Gene = common_gene
  
  # Merge the gene counts for the current gene with the main results data frame
  gene_variant_counts <-rbind(gene_variant_counts, gene_counts)
}

gene_variant_counts$Chromosome <- sub("chr", "", gene_variant_counts$Chromosome)

# # Pivot the data to create the desired format
# gene_variant_counts_pivot <- gene_variant_counts %>%
#   pivot_wider(names_from = Gene, values_from = Variant_Count, values_fill = 0)
# 
# # Set row names to the Patient column
# gene_variant_counts_df <- as.data.frame(gene_variant_counts_pivot)
# rownames(gene_variant_counts_df) <- gene_variant_counts_df$Patient
# gene_variant_counts_df <- gene_variant_counts_df[, -1]

# # accounting for U6 being %in% RNU6-417P
# colnames(gene_variant_counts_df) <- gsub("-", ".", colnames(gene_variant_counts_df))
# gene_variant_counts_df$U6 <- (gene_variant_counts_df$U6 - gene_variant_counts_df$RNU6.417P)



# ------------------------ Plotting mutation of genes per patient


# gene_variant_counts <- gene_variant_counts %>%
#   filter(Gene != "U6")



gene_variant_counts$Gene <- factor(gene_variant_counts$Gene, levels = unique(gene_variant_counts$Gene[order(gene_variant_counts$Chromosome)]))

ggplot(gene_variant_counts, aes(x = factor(Patient, level = c("S2", "S7", "S9", "S32", "S6", "S45", "S21")), y= Gene, fill = Variant_Count)) +
        geom_tile() + 
        scale_fill_gradient(name = "Variant Count",
                      low = "#FFFFFF",
                      high = "#012345") +
        labs(y = "Common Genes", x = "Patients")


ggsave(paste0(output_path, "variant_counts_per_patient.jpg"))

# ------------------------ Creating a df of a count matrix of each gene common gene hit by x variants per each sample (FOR CANCER RELATED GENES)

gene_variant_counts <- data.frame(Patient = character(0))

# this list is based from brief look through literature
common_cancer_related <- c("KCNJ18", "IRF4", "CRLF2", "MIR663AHG", "CNTNAP2", "ADARB2",
                           "RBMS3", "DUSP22", "KIF26B", "PTPRN2", "NRXN3", "CSMD1", "CTNNA3",
                           "FHIT")

# Iterate through each common gene name
for (common_gene in common_cancer_related) {
  # Count how many rows match the current gene for each patient
  gene_counts <- all_snps %>%
    filter(str_detect(gene_name, common_gene)) %>%
    group_by(Patient, Chromosome) %>%
    summarise(Variant_Count = n()) %>%
    ungroup()
  
  # Add a new column with the current gene name
  gene_counts$Gene = common_gene
  
  # Merge the gene counts for the current gene with the main results data frame
  gene_variant_counts <-rbind(gene_variant_counts, gene_counts)
}

gene_variant_counts$Chromosome <- sub("chr", "", gene_variant_counts$Chromosome)

# # Pivot the data to create the desired format
# gene_variant_counts_pivot <- gene_variant_counts %>%
#   pivot_wider(names_from = Gene, values_from = Variant_Count, values_fill = 0)
# 
# # Set row names to the Patient column
# gene_variant_counts_df <- as.data.frame(gene_variant_counts_pivot)
# rownames(gene_variant_counts_df) <- gene_variant_counts_df$Patient
# gene_variant_counts_df <- gene_variant_counts_df[, -1]

# # accounting for U6 being %in% RNU6-417P
# colnames(gene_variant_counts_df) <- gsub("-", ".", colnames(gene_variant_counts_df))
# gene_variant_counts_df$U6 <- (gene_variant_counts_df$U6 - gene_variant_counts_df$RNU6.417P)



# ------------------------ Plotting mutation of genes per patient


# gene_variant_counts <- gene_variant_counts %>%
#   filter(Gene != "U6")



gene_variant_counts$Gene <- factor(gene_variant_counts$Gene, levels = unique(gene_variant_counts$Gene[order(gene_variant_counts$Chromosome)]))

ggplot(gene_variant_counts, aes(x = factor(Patient, level = c("S2", "S7", "S9", "S32", "S6", "S45", "S21")), y= Gene, fill = Variant_Count)) +
  geom_tile() + 
  scale_fill_gradient(name = "Variant Count",
                      low = "#FFFFFF",
                      high = "#012345") +
  labs(y = "Common Genes", x = "Patients")


ggsave(paste0(output_path, "cancer_variant_counts_per_patient.jpg"))



