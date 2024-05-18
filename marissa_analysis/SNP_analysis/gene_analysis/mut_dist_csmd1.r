#***********************************
# Author: Marissa Esteban
# Date: July 24 2023
# Title: SNP analysis pt7 --  Analysis of mutations within the CSMD1 gene
#
#***********************************

library(dplyr)
library(tidyr)

all_snps <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/merged_snp_annotations.rds")
only_csmd1 <- all_snps %>%
  filter(gene_name == "CSMD1")
#only_csmd1 = subset(only_csmd1, select = c(Chromosome, POS, REF, ALT, gene_name, variant_id, Patient))

# remove duplicate rows, because some variants were annotated twice...
only_csmd1 <- only_csmd1 %>% distinct(Chromosome, POS, REF, ALT, gene_name, variant_id, Patient, .keep_all = TRUE)



