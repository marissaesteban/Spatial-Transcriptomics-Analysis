#***********************************
# Author: Marissa Esteban
# Date: July 24 2023
# Title: SNP analysis pt5 -- Visualizatino of SNPs on the genome using CoMutPlotter 
# http://tardis.cgu.edu.tw/comutplotter/CoMutPlotter/

# ** CoMutPlotter is not working trying something different
#
# output file: /disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/CoMutPlotter/all_snps_tsv.tsv
#***********************************

library(stringr)

# input file needs to be a tsv so turning merged_snp_ann into a tsv

all_snps <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/merged_snp_annotations.rds")

header <- c("sample", "chromosome", "chromosome_start", "chromosome_end", "reference_genome_allele", "mutated_to_allele")
snps_tsv <- data.frame(matrix(nrow = nrow(all_snps), ncol = length(header)))
colnames(snps_tsv) <- header

snps_tsv$sample <- str_trim(all_snps$Patient)
snps_tsv$chromosome <- str_trim(sub("chr", "", all_snps$Chromosome))
snps_tsv$chromosome_start <- str_trim(as.character(all_snps$POS))  # Convert POS to character before using str_trim()
snps_tsv$chromosome_end <- str_trim(as.character(all_snps$POS))    # Convert POS to character before using str_trim()
snps_tsv$reference_genome_allele <- str_trim(all_snps$REF)
snps_tsv$mutated_to_allele <- str_trim(all_snps$ALT)

output_path <- "/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/CoMutPlotter/all_snps_tsv.tsv"
write.table(snps_tsv, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
