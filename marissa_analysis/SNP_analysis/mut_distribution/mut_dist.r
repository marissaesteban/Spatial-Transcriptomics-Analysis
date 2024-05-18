#***********************************
# Author: Marissa Esteban
# Date: July 24 2023
# Title: SNP analysis pt6 --  Visualization of SNPs on the genome, and how many SNPs there are per patient
#
#***********************************


library(ggplot2)
library(dplyr)

all_snps <- readRDS("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/merged_snp_annotations.rds")
all_snps$Chromosome <- sub("chr", "", all_snps$Chromosome)

mb <- 1000000

order = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
          "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

# plotting the number of variants there are per mb (this is for ALL patients)
ggplot(all_snps, aes(x = POS)) +
  geom_histogram(binwidth = mb, fill = "steelblue", color = "black") +
  labs(title = "Variant Distribution",
       x = "Position (by mb)",
       y = "Variant Count") +
  facet_grid(~ factor(Chromosome, levels = order), scales = "free_x", space = "free") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
ggsave("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/mut_distribution/var_dist.png")



# --------------------------------------------------------------------------------------------------------------------
# plotting the number of mutations for each patient
# variants_count_df <- as.data.frame(variants_count)
# colnames(variants_count_df) <- c("Patient", "Variant_Count")

variants_count_df <- all_snps %>% count(Patient)

# Create the bar plot with single color bars
ggplot(variants_count_df, aes(x = Patient, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") + 
  labs(title = "Number of Variants per Patient",
       x = "Patient ID",
       y = "Variant Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotates x-axis labels for better readability

ggsave("/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/marissa_analysis/SNP_analysis/mut_distribution/patient_dist.png")