# leona knusel, inspired by Chiara Auwerx

### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

### Load Main Data ##############################

# Sample QC; This corresponds to the "ukb_sqc_v2.txt" file available from the UKBB portal
sample_qc <- as.data.frame(fread(snakemake@input[["file_sample_qc"]], header = F, select = c(3:68)))
colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")

print(paste0("Number of individuals in sample_qc file = ", nrow(sample_qc)))

# Sample eids; This corresponds to a list of all sample eids (with sex information), retrieved from a ".fam" file available from the UKBB portal
sample_eid <- as.data.frame(fread(snakemake@input[["file_sample_eid"]], header = F, select = c(1,5), col.names = c("eid", "sex")))

print(paste0("Number of individuals in sample_eid file = ", nrow(sample_eid)))

### Merge Main Data #############################

df <- cbind(sample_eid, sample_qc)
print(paste0("Start: ", nrow(df), " individuals"))

### Sample Filtering ############################

### STEP 1: Exclude related samples (pca_calculation = 1)
df <- df[which(df$pca_calculation == 1), ]
print(paste0("STEP 1: Exclude related samples: ", nrow(df), " individuals"))


### STEP 2: Exclude non-white, non-British ancestry samples (white_british = 1)
df <- df[which(df$white_british == 1), ]  
print(paste0("STEP 2: Exclude non-white, non-British ancestry samples: ", nrow(df), " individuals"))


### STEP 3: Exclude retracted samples
retracted <- as.data.frame(fread(snakemake@input[["file_retracted_samples"]], header = F, col.names = "eid"))
df <- df[!df$eid %in% retracted$eid, ]
print(paste0("STEP 3: Exclude retracted samples: ", nrow(df), " individuals"))


### STEP 4: Exclude samples with non-matching submitted vs. inferred sex
df <- df[which(df$submitted_gender == df$inferred_gender), ]
print(paste0("STEP 4: Exclude sex mismatches: ", nrow(df), " individuals"))


write_delim(df[, "eid", drop = F], snakemake@output[["active_samples"]])
