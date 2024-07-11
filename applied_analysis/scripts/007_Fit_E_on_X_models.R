# lknusel
# obtain estimate for effect of E on X
# corrected for medication if indicated
# corrected for sex
# if E != age, corrected for age, age^2, age*sex

library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)


# needed:
# covariate files (to correct for PCs and sex)
load(snakemake@input[["file_covariates"]])

# medication file, to correct for medication effect
load(snakemake@input[["file_medication"]])

# exposure phenotye (to use as outcome for once)
load(snakemake@input[["file_exposure_phenotypes"]])

# interaction phenotype (to use as exposure)
load(snakemake@input[["file_interaction_phenotype"]])

# Correct for medication if appropriate

exposure_ID <- snakemake@params[["exposure_IDs"]]

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
X <- merge(X, X_medication, by = "eid")

eids_both_index <- which(X[,"eid"] %in% pheno_df$eid)

X <- X[eids_both_index,]

eid_col <- grep("eid", colnames(X))

# rounding is necessary, because
# the annotation differed (scientific vs. standard) and
# as a consequence one value differed minimally.

eids_covariate_matrix <- X[,eid_col] |> as.numeric() |> round()
eids_pheno_df <- pheno_df$eid |> as.numeric() |> round()

if(!all(eids_covariate_matrix == eids_pheno_df)){
    pheno_df <- pheno_df[pheno_df$eid %in% X$eid, ]
    eids_pheno_df <- pheno_df$eid |> as.numeric() |> round()
    eid_order <- match(pheno_df$eid, X[,eid_col])
    X <- X[eid_order, ]
} 

eids_covariate_matrix <- X[,eid_col] |> as.numeric() |> round()
eids_pheno_df <- pheno_df$eid |> as.numeric() |> round()

if(!all(eids_covariate_matrix == eids_pheno_df)){
	stop("HOUSTON, WE'VE HAD A PROBLEM!! Eids in covariate matrix and pheno df do not match!")
} else {
	X <- X[,-eid_col]
}

# obtain Effect of Sex and age on Exposure Phenotypes (for correction)
# remove age column as its duplicated
X <- X[,colnames(X) != "age"]

#remove PC columns
PC_cols <- grep("^PC", colnames(X))
X <- X[, -PC_cols]

# if we're interested in the age effect, remove all other age related parameters!
if(snakemake@params[["interaction_ID"]] == "21022"){
    X <- X[, colnames(X) != "age_exact"]
    X <- X[,colnames(X) != "age_exact_x_sex"]
    X <- X[,colnames(X) != "age_exact2"]
}

X <- as.matrix(X)

# let's check if participant IDs also agree with the interaction df:
interaction_df <- interaction_df[interaction_df$eid %in% eids_pheno_df, ]
eids_interaction_df <- interaction_df$eid

if(!all(eids_interaction_df == eids_pheno_df)){
    new_order <- match(eids_pheno_df, eids_interaction_df)
	interaction_df <- interaction_df[new_order,]
	eids_interaction_df <- interaction_df$eid |> as.numeric() |> round()
}

if(all(eids_interaction_df == eids_pheno_df)){
    interaction_column <- paste0(snakemake@params[["interaction_ID"]], "_z_c")
    interaction_df$E <- interaction_df[,interaction_column]
    z_pheno_column <- paste0(snakemake@params[["exposure_IDs"]], "_z")
    model_E_on_x <- lm(unlist(pheno_df[, z_pheno_column]) ~ interaction_df$E + X) |> summary()
    save(model_E_on_x, file = snakemake@output[["model_E_on_X"]])
    print("Model for E on X effect is fitted and saved!")
} else {
    print("Eids do not align between pheno_df, interactio_df, and X. Please double check!")
}



