#***********************************
# Author: Marissa Esteban
# Date: July 24 2023
# Title: SNP analysis pt3 -- finding top mutated genes for each patient, using the SNPs filtered by MQ and non intron variants
#
# First filtering the merged snp annotations file by MQ and coding regions
# Then making a new df consisting of the genes mapped to the variants (that passed the filter)
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

# filtering the master list of ALL SNPs
merged_snp_annotations <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/merged_snp_annotations.rds")
filtered_mq <- subset(merged_snp_annotations, MQ >= 40.0) # MQ >= 40.0 (mapping quality)
no_introns <- subset(filtered_mq, effect != "intron_variant") # filtering so there are no intron variants





# making a variants_per_gene list using annotations filtered by MQ and coding regions
variants_per_gene_filtered <- no_introns %>%
                        dplyr::select(Chromosome, POS, REF, ALT, gene_ID, Patient)  %>%
                        filter(gene_ID != "") %>%
                        separate_longer_delim(c(`gene_ID`), "-") %>%
                        unique.data.frame() %>%
                        mutate(variant_id = paste(`Chromosome`, `POS`, `REF`, `ALT`, sep = "_")) %>%
                        group_by(gene_ID, Patient) %>%
                        summarise(variants = list(variant_id))


variants_per_gene_filtered$gene_name <- ensembl_converter[variants_per_gene_filtered$gene_ID, "name"]




# defining df to find the top mutated genes for each patient
mutated_header <- c(colnames(variants_per_gene_filtered))

# finding top mutated genes for each patient based on the master list of of all SNPs
patients <- c("S2", "S6", "S7", "S9", "S21", "S32", "S45")
for(patient in patients){

    # filtering by patient (patient_filtered)
    patient_filtered <- variants_per_gene_filtered[variants_per_gene_filtered$Patient == patient, ]

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
    output_path <- "/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/top_mutated_encoding_snps/"
    saveRDS(patient_filtered, paste0(output_path, patient, "_top_mutated_encoding.rds"))

}