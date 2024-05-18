#***********************************
# Author: Marissa Esteban
# Date: July 24 2023
# Title: SNP analysis pt8 --  Analysis of mutations within the ERBB2 gene
#
#***********************************

library(dplyr)
library(tidyr)

all_snps <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/merged_snp_annotations.rds")
only_erbb2 <- all_snps %>%
  filter(gene_name == "HMGCS2")