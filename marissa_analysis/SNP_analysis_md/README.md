# WGS Variant Analysis -- SNPs

The structure of the repository is as shown below:

Path to Strelka VCF files: **\*\*\*\*\*\*\*\*\*\***

| Folder                 | Description                                                                                                                                    |
|------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| Processing SNPs        | Taking the Strelka_SNV VCFs (from Sarek) and putting information into .rds files for each patient, and also merged versions with ALL patients. |
| SNP Distribution       | Visualization of the SNPs on the genome                                                                                                        |
| Commonly Mutated Genes | Analyzation of the commonly SNP hit genes                                                                                                      |
| Gene Analysis          | Mutated genes from WGS --> ST data                                                                                                             |
|                        | Diff Expressed genes from ST data --> WGS data                                                                                                 |


### Workflow

**Variant Calling and VCF Parsing** -- The fastq files for all 7 patients were pushed through Sarek. The outputted VCF files are here\*\***.** We extracted all of the useful information from the VCFs with Haga Parser and put them into rds files, the output files are here\*\*\*\*\*\*\*.

**SNP Analysis** -- The tool Strelka called the SNP variants. I merged the rds files of SNP annotations for each patient into one df, and also merged the df of "genes: variants" from each patient. Then looked at SNP distribution, commonly mutated genes, and specific gene analysis.
