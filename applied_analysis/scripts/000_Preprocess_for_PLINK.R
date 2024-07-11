# script no. 1 
# to read in files, filter for SNPs with p>10^-4
# to obtain RSIDs, chr, pos, ref and alt
# and save everything PLINK-compatible

library(stringr)
library(dplyr)
library(readr)

data <- read.table(snakemake@input[["expo"]], header = T)

print(paste0("Number of variants in SumStats = ", nrow(data)))

variants_df <- read.table(snakemake@params[["variant_df"]], header = T)

print(paste0("Number of variants in the variants file = ", nrow(variants_df)))

data <- merge(data, variants_df[,c("variant", "chr", "pos", "ref", "alt", "rsid", "varid")], 
              by = "variant", all = FALSE)

print(paste0("Number of variants after merging with variants file = ", nrow(data)))

data <- data[data$pval < 10^-4, ]

print(paste0("Number of variants after removing SNPs with p > 10^-4 = ", nrow(data)))

reference_df <- read.delim(snakemake@params[["file_snp_reference"]])
reference_df$unicol <- paste0(reference_df$X.CHROM, ":", reference_df$POS)

data$unicol <- paste0(data$chr, ":", data$pos)

data <- data[data$unicol %in% reference_df$unicol, ]

print(paste0("Number of variants after filtering for SNPs available in reference file = ",
             nrow(data)))

data <- rename(data, effect_allele = alt)
data <- rename(data, other_allele = ref)
data <- rename(data, SNP = rsid)
data <- rename(data, P = pval)

data <- data[!is.na(data$P), ]

print(paste0("Number of variants after removing NAs from P column = ", nrow(data)))

write_delim(data, snakemake@output[["raw"]])

