#***********************************
# Author: Marissa Esteban
# Date: July 24 2023
# Title: SNP analysis pt2 -- finding top mutated genes for each patient, this is looking at ALL the found SNPs (not filtered by encoding regions)
#
#***********************************

library(biomaRt)
library(dplyr)
library(tidyr)


# Connect to the Ensembl database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") # May take a moment
ensembl_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"), mart = ensembl)

# Makes some kind of data frame as dictionary?
ensembl_converter <- ensembl_genes %>%
    data.frame(name = .$external_gene_name, start = .$start_position, end = .$end_position, row.names = .$ensembl_gene_id)


# Function to get gene length based on gene ID or gene name
get_gene_length <- function(gene_id_or_name) {
  
  if (!(gene_id_or_name %in% rownames(ensembl_converter))) {
    print(gene_id_or_name)
    return(1)
  }
  
  gene_length <- ensembl_converter[gene_id_or_name, "end"] - ensembl_converter[gene_id_or_name,"start"] + 1
  return(gene_length)
}

variants_per_gene <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/merged_snp_variants_per_gene.rds")


# defining df to find the top mutated genes for each patient
mutated_header <- c(colnames(variants_per_gene))

# finding top mutated genes for each patient based on the master list of of all SNPs
patients <- c("S2", "S6", "S7", "S9", "S21", "S32", "S45")
for(patient in patients){

    # filtering by patient (patient_filtered)
    patient_filtered <- variants_per_gene[variants_per_gene$Patient == patient, ]

    # # splitting the variants column into a list
    # patient_filtered$variants <- strsplit(patient_filtered$variants, ";")

    # counting number of variants per gene
    patient_filtered$num_variants <- sapply(patient_filtered$variants, length)
   
    # getting the gene lengths for each gene
    # if gene length doesn't exist, get_gene_length() returns 1 (will filter these genes out)
    patient_filtered$gene_length <- sapply(patient_filtered$gene_ID, get_gene_length)

    patient_filtered <- patient_filtered %>% filter(gene_length != 1)

    # getting the number of mutations per kb
    patient_filtered$mut_per_kb <- patient_filtered$num_variants / (patient_filtered$gene_length / 1000)

    # sorting the df based on the number of variants in decending order
    patient_filtered <- patient_filtered[order(-patient_filtered$mut_per_kb), ]

    patient_filtered$gene_name <- ensembl_converter[patient_filtered$gene_ID, "name"]

    # reordering the columns
    patient_filtered <- patient_filtered[, c("Patient", "gene_ID", "gene_name", "variants", "num_variants", "mut_per_kb")]

    # saving the rds into a new file
    output_path <- "/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_all_snps/"
    saveRDS(patient_filtered, paste0(output_path, patient, "_top_mutated_all.rds"))

}