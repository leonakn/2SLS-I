# lknusel
# Run MR to obtain causal main- and quadratic effects

outcome_column <- paste0(snakemake@params[["outcome_column"]])

#load(snakemake@input[["file_phenotypes"]]) #pheno_df
load(snakemake@input[["file_GRS"]]) #df_gen_score
load(snakemake@input[["file_outcome_phenotype"]]) # outcome phenotype

trait_id <- snakemake@params[["exposure_IDs"]]

# df <- merge(pheno_df, df_gen_score, by = "eid")
# print(paste0("pheno_df and df_gen_score merged! N individuals = ",
#              nrow(df)))

df <- merge(df_gen_score, outcome_df, by = "eid")

print(paste0("df and outcome_df are merged! N individuals = ",
             nrow(df)))

GRS_column <- snakemake@params[["exposure_GRS_column"]]
GRS2_column <- paste0(GRS_column, "_2")

df[,GRS2_column] <- df[,GRS_column]^2

print(paste0("Class of df is ", class(df)))
print(paste0("GRS_column is called ", GRS_column))
print(paste0("Class of df[,GRS_column] is ", class(df[,GRS_column])))

# Correct for medication
load(snakemake@input[["file_medication"]])

X <- data.frame("eid" = df$eid)

# Insulin:
if(trait_id %in% snakemake@params[["correct_for_insulin"]]){
	X <- merge(X, medication_df[,c("eid", "Insulin")], by = "eid")
}

# Blood pressure medication
if(trait_id %in% snakemake@params[["correct_for_bp_medication"]]){
	X <- merge(X, medication_df[,c("eid", "bp_regul")], by = "eid")
}

# Cholesterol lowering medication
if(trait_id %in% snakemake@params[["correct_for_cholesterol_lowering"]]){
	X <- merge(X, medication_df[,c("eid", "chol_inhibit")], by = "eid")
}

# Hormone replacement therapy
if(trait_id %in% snakemake@params[["correct_for_hormone_replacement"]]){
	X <- merge(X, medication_df[,c("eid", "hormone_repl")], by = "eid")
}

# Oral contraceptive
if(trait_id %in% snakemake@params[["correct_for_oral_contraceptive"]]){
	X <- merge(X, medication_df[,c("eid", "oral_contraceptive")], by = "eid")
}

if(!all(X$eid == df$eid)){
	stop("HOUSTON, WE'VE HAD A PROBLEM!! Eids in covariate matrix and pheno df do not match!")
} else {
    eid_column <- grep("eid", colnames(X))
	X <- X[,-eid_column]
}

X <- as.matrix(X)

if(! trait_id %in% c(snakemake@params[["correct_for_insulin"]], snakemake@params[["correct_for_bp_medication"]], 
                     snakemake@params[["correct_for_cholesterol_lowering"]], snakemake@params[["correct_for_hormone_replacement"]],
                     snakemake@params[["correct_for_oral_contraceptive"]])){
        # Y ~ GRSX + Age + Age2 + Sex
        model_0 <- lm(df[,outcome_column] ~ df[,GRS_column] + 
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex) 

        print("Model 0 is fitted")

        # include quadratic effect Y ~ GRSX + GRSX2 + Age + Age2 + Sex
        model_1 <- lm(df[,outcome_column] ~ df[,GRS_column] + 
                                            df[,GRS2_column] +
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex) 

        print("Model 1 is fitted")

} else {
    # Y ~ GRSX + Age + Age2 + Sex
        model_0 <- lm(df[,outcome_column] ~ df[,GRS_column] + 
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex + X) 

        print("Model 0 is fitted")

        # include quadratic effect Y ~ GRSX + GRSX2 + Age + Age2 + Sex
        model_1 <- lm(df[,outcome_column] ~ df[,GRS_column] + 
                                            df[,GRS2_column] +
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex + X)

        print("Model 1 is fitted")

}


print("Complete model is fitted")

all_models <- list(model_0, model_1)
save(all_models, file = snakemake@output[["mr_models"]])














