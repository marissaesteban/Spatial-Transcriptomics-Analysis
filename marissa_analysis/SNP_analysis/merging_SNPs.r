#***********************************
# Author: Marissa Esteban
# Date: July 20 2023
# Title: SNP analysis pt1 -- merging output for strelka SNPs
#
# Merging all of the strelka_somatic_snvs.* into their own master merged files that include ALL patients.

# SNPs defined in these files have only passed the 'PASS' filter. NOT FILTERED BY MQ > 40.0 YET
# Three files saved as
# merged_snp_annotations.rds
# merged_snp_variants_per_gene.rds
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

input_wd <- "/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/sarek/parsed_vcfs4/"
patients <- c("S2", "S6", "S7", "S9", "S21", "S32", "S45")

# this is getting the headers for the combined and variatns per gene files to make the 'merged_*.rds' files
c <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/sarek/parsed_vcfs4/S2/strelka_somatic_snvs.combined.rds")
g <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/sarek/parsed_vcfs4/S6/strelka_somatic_snvs.gene_lists.rds")
v <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/sarek/parsed_vcfs4/S2/strelka_somatic_snvs.variants_per_gene.rds")

c_header <- c(colnames(c), "Patient")
g_header <- c(colnames(g), "Patient")
v_header <- c(colnames(v), "Patient")

merged_snp_annotations <- data.frame(matrix(nrow = 0, ncol = length(c_header)))
colnames(merged_snp_annotations) <- c_header
merged_snp_genes_per_variant <- data.frame(matrix(nrow = 0, ncol = length(g_header)))
colnames(merged_snp_genes_per_variant) <- g_header
merged_snp_variants_per_gene <- data.frame(matrix(nrow = 0, ncol = length(v_header)))
colnames(merged_snp_variants_per_gene) <- v_header

# iterate through every patient directory from the haga parser output
for(patient in patients){

    print(patient)

    path <- paste0(input_wd, patient, "/")
    file_list <- list.files(path, pattern = "strelka_somatic_snvs")

    for (file in file_list){

        # merging in the combined
        if (grepl("combined", file)){
            print("in combined")
            file_path <- paste0(path, "/", file)
            snp_annotations <- readRDS(file_path)
            snp_annotations$Patient <- patient

            common_cols <- intersect(c_header, colnames(snp_annotations))
            non_common_cols <- setdiff(c_header, colnames(snp_annotations))
            snp_annotations <- snp_annotations[, common_cols]

            if (length(non_common_cols) != 0 ) {

                for ( col_name in non_common_cols){
                    snp_annotations[[col_name]] <- NA
                }
            }

            merged_snp_annotations <- rbind(merged_snp_annotations, snp_annotations)
        }
    }
}

output_path <- "/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/"


filtered_mq <- subset(merged_snp_annotations, MQ >= 40.0) # MQ >= 40.0 (mapping quality)
filtered_mq <- filtered_mq %>% distinct(Chromosome, POS, REF, ALT, gene_name, variant_id, Patient, .keep_all = TRUE)
saveRDS(filtered_mq, paste0(output_path, "merged_snp_annotations.rds"))

#saveRDS(merged_snp_genes_per_variant, paste0(output_path, "merged_snp_genes_per_variant.rds"))


# redoing the variants per gene df using the snp annotations with an MQ >= 40
merged_variants_per_gene <- filtered_mq %>%
                        dplyr::select(Chromosome, POS, REF, ALT, gene_ID, Patient)  %>%
                        filter(gene_ID != "") %>%
                        separate_longer_delim(c(`gene_ID`), "-") %>%
                        unique.data.frame() %>%
                        mutate(variant_id = paste(`Chromosome`, `POS`, `REF`, `ALT`, sep = "_")) %>%
                        group_by(gene_ID, Patient) %>%
                        summarise(variants = list(variant_id))


merged_variants_per_gene$gene_name <- ensembl_converter[merged_variants_per_gene$gene_ID, "name"]


saveRDS(merged_variants_per_gene, paste0(output_path, "merged_snp_variants_per_gene.rds"))