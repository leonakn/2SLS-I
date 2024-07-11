# leona knusel
# heavily inspired by Marie Sadler
# select pruned variants from .bgen files

library(rbgen)
library(readr)
library(dplyr)
library(data.table)
library(stringr)


# map people to order
df_sample <- read.table(snakemake@params[["file_sample_order"]], header = T)
df_sample <- df_sample[2:nrow(df_sample),] #get rid of 2. header


print(paste0("Number of samples in file_sample_order: ", length(df_sample$ID_1)))

# add column with sample name from bgen files
df_sample["order"] <- as.character(seq.int(nrow(df_sample)))
df_sample["bgen_ID"] <- paste0("(anonymous_sample_", df_sample$order, ")")

# participants of interest: white british

active_participants <- read.csv(snakemake@input[["file_wb_participants"]])
active_participants <- as.character(active_participants$eid)

print(paste0("Number of samples in file_wb_participants = ", length(active_participants)))

samples <- as.vector(df_sample[as.character(df_sample$ID_1) %in% active_participants, "bgen_ID"]) #BGen ID of white british people
print(paste0("Number of samples with white british ancestry available in file_sample_order: ", 
            length(samples)))

#read SNP rsids

snps <- read.table(snakemake@input[["file_snps"]], header = T)
print(paste0("Number of SNPS = ", nrow(snps)))

chromosomes <- unique(snps$CHR) %>% sort()

files_bgen <- snakemake@input[["files_bgen"]]

gen <- sapply(chromosomes, function(x) bgen.load(files_bgen[x], 
                                                 rsids = paste(snps[snps$CHR == x, "SNP"]),
                                                 samples = samples))

lsGeno <- NULL

lsGeno$variants <- do.call(rbind, gen["variants", ])
lsGeno$data <- as.data.frame(matrix(data = numeric(), nrow = length(gen[,1]$samples)),
                             row.names = df_sample[df_sample$ID_1 %in% active_participants, "ID_1"])

# fill lsGeno with values

chromosomes <- dim(gen)[2]

for (i in 1:chromosomes){
  rsids_chr <- gen[,i]$variants$rsid
  for(j in 1:length(rsids_chr)){
    lsGeno$data[,paste0(rsids_chr[j])] <- gen[,i]$data[paste0(rsids_chr[j]),,2] + 
      gen[,i]$data[paste0(rsids_chr[j]),,3]*2
  }
}

lsGeno$data$eid <- df_sample[df_sample$ID_1 %in% active_participants, "ID_1"]

print(paste0("number of snps in lsGeno$variants = ", length(unique(lsGeno$variants$rsid))))
print(paste0("number of snps in lsGeno$data = ", ncol(lsGeno$data)-1))
print(paste0("number of individuals in lsGeno$data = ", nrow(lsGeno$data)))

print("Genotype list created")

save(lsGeno, file = snakemake@output[["genetics_list"]])

print(paste0("Geno Data saved!"))