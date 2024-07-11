# lknusel
# Obtain GRS from SNPs and SNP effects.
# Compare with raw phenotype and save plot.

# snp effects
library(ggplot2)
library(data.table)
snp_effects_df <- fread(snakemake@input[["file_snp_effects"]], 
                        data.table = F)
snps_clumped <- fread(snakemake@input[["clumped_snps"]],
                      data.table = F, select = "SNP")

# use only independent genome wide significant variants:
snp_effects_df <- merge(snps_clumped, snp_effects_df, by = "SNP",
                        all.x = T, all.y = F)

sign_SNPs <- snp_effects_df$SNP

print(paste0("Number of genome wide significant SNPs = ", length(sign_SNPs)))

# Genetics data:
load(snakemake@input[["file_genotypes"]])

geno_df <- lsGeno$data
rm(lsGeno)

# only keep eid and variants, which were kept in 
# make sure we're having the right variants in the right order
available_snps <- colnames(geno_df)[!colnames(geno_df) %in% "eid"]
snp_order_index <- match(available_snps, snp_effects_df$SNP)
betas <- snp_effects_df[snp_order_index, c("beta", "SNP")]

eid_index <- which(colnames(geno_df) == "eid")

geno_df <- geno_df[,c(betas$SNP, "eid")]

if(all(betas$SNP == colnames(geno_df[,1:length(betas$SNP)]))){
    betas <- betas$beta
}

GRS <- as.matrix(geno_df[,1:length(available_snps)]) %*% betas |> scale()

df_gen_score <- data.frame(eid = geno_df$eid)

exposure_ID <- snakemake@params[["exposure_IDs"]]
colname_pheno_GRS <- paste0("GRS_", exposure_ID)

df_gen_score[,colname_pheno_GRS] <- GRS


#### Now, obtain effect of GRS on X
# correct for the standard stuff, age, age2, sex, medication
load(snakemake@input[["file_phenotypes"]])
load(snakemake@input[["file_covariates"]])
load(snakemake@input[["file_medication"]])


df_gen_score <- df_gen_score[df_gen_score$eid %in% pheno_df$eid, ]

df_gen_score <- df_gen_score[order(df_gen_score$eid), ]
pheno_df <- pheno_df[order(pheno_df$eid), ]

# Correct for medication if appropriate

X_medication <- data.frame("eid" = medication_df$eid)

# Insulin:
if(exposure_ID %in% snakemake@params[["correct_for_insulin"]]){
	X_medication <- merge(X_medication, medication_df[,c("eid", "Insulin")], by = "eid")
}

# Blood pressure medication
if(exposure_ID %in% snakemake@params[["correct_for_bp_medication"]]){
	X_medication <- merge(X_medication, medication_df[,c("eid", "bp_regul")], by = "eid")
}

# Cholesterol lowering medication
if(exposure_ID %in% snakemake@params[["correct_for_cholesterol_lowering"]]){
	X_medication <- merge(X_medication, medication_df[,c("eid", "chol_inhibit")], by = "eid")
}

# Hormone replacement therapy
if(exposure_ID %in% snakemake@params[["correct_for_hormone_replacement"]]){
	X_medication <- merge(X_medication, medication_df[,c("eid", "hormone_repl")], by = "eid")
}

# Oral contraceptive
if(exposure_ID %in% snakemake@params[["correct_for_oral_contraceptive"]]){
	X_medication <- merge(X_medication, medication_df[,c("eid", "oral_contraceptive")], by = "eid")
}

# obtain covariates
# keep: age, age2, sex
# remove sex if its a female_only_trait
if(sum(grepl(exposure_ID, snakemake@params[["female_only_IDs"]])) > 0){
    X <- X[,colnames(X) %in% c("eid", "age_exact", "age_exact2")]
} else {
     X <- X[,colnames(X) %in% c("eid", "sex", "age_exact", "age_exact2")]
}

X <- merge(X, X_medication, by = "eid")

X <- X[X$eid %in% pheno_df$eid,]

pheno_df <- pheno_df[pheno_df$eid %in% X$eid, ]


if(all(X$eid == pheno_df$eid)){
    X <- X[,colnames(X) != "eid"]
} else {
    stop("Houston, we've had a problem!! Eids in covariate matrix do not match with those in pheno df!")
}

X <- as.matrix(X)

df_gen_score <- df_gen_score[df_gen_score$eid %in% pheno_df$eid, ]

if(all(pheno_df$eid == df_gen_score$eid)){
    Y <- paste0(snakemake@params[["exposure_IDs"]], "_IRNT_c")
    # for full GRS (including medication correction)
    colname_pheno_GRS <- paste0("GRS_", snakemake@params[["exposure_IDs"]])
    GRS_on_pheno_model <- lm(pheno_df[,Y] ~ df_gen_score[,colname_pheno_GRS] + X) |> summary()
    save(GRS_on_pheno_model, 
         file = snakemake@output[["file_GRS_on_pheno_model"]])

    # for full GRS (no medication correction)
    minimal_GRS_on_pheno_model <- lm(pheno_df[,Y] ~ df_gen_score[,colname_pheno_GRS]) |> summary()
    save(minimal_GRS_on_pheno_model, file = snakemake@output[["file_minimal_GRS_on_pheno_model"]])

    save(df_gen_score, file = snakemake@output[["file_GRS"]])
    print("Obtained GRS, fitted GRS on X, saved model, plotted association, saved GRS!")

} else {
    warning("eids do not match, please check what went wrong there..")
}



