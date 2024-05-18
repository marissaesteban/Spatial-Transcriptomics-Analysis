#***********************************
# Merging the gene lists for each patient into one csv, and merging the gene lists for ALL patients into one master csv
#
#***********************************
library(tibble)

patients <- c("S2", "S6", "S7", "S9", "S21", "S32", "S45")
input_dir <- "/disk2/user/cda/SpatialTranscriptomics/Whole-Genome-Sequencing_Analysis/wholegenomesequencing_breastcancer_summer23/sarek/parsed_vcfs/"

for(patient in patients){

    # setup for the merged dataframe
    header <- c("Patient", "Chromosome", "POS", "REF", "ALT", "names", "ids", "variant_type")
    merged_df <- data.frame(matrix(ncol = length(header), nrow = 0))
    colnames(merged_df) <- header

    patient_dir <- paste0(input_dir, patient, "/")
    file_list <- list.files(patient_dir)

    # iterating through patient directory and looking for the ".gene_lists.csv" files
    for (file in file_list){
        if (grepl("gene_lists", file)){
            file_path <- paste0(patient_dir, "/", file)
            gene_list_df <- read.delim(file_path)
            gene_list_df$Patient <- patient

            merged_df <- rbind(merged_df, gene_list_df)
        }
    }

    # reordering the columns so $patient is first
    merged_df <- merged_df[, c("Patient", "Chromosome", "POS", "REF", "ALT", "names", "ids", "variant_type")]

    outpath <- paste0(patient_dir, "/merged_gene_lists.csv")
    write.table(merged_df, outpath, sep = "\t", row.names = FALSE, quote = FALSE)
}




# making a master merged gene list!!!!!!!
# setup for the GIANT merged dataframe
    header <- c("Patient", "Chromosome", "POS", "REF", "ALT", "names", "ids", "variant_type")
    master_df <- data.frame(matrix(ncol = length(header), nrow = 0))
    colnames(merged_df) <- header

for(patient in patients){

    patient_dir <- paste0(input_dir, patient, "/")
    file_list <- list.files(patient_dir)

    # iterating through patient directory and looking for the "merged_gene_lists.csv" file
    for (file in file_list){
        if (grepl("merged_gene_lists", file)){
            file_path <- paste0(patient_dir, "/", file)
            gene_list_df <- read.delim(file_path)

            master_df <- rbind(master_df, gene_list_df)
        }
    }
}

outpath <- paste0(input_dir, "/master_merged_gene_lists.csv")
write.table(master_df, outpath, sep = "\t", row.names = FALSE, quote = FALSE)