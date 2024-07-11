# Run MR and obtain the necessary correction measures to 
# correct the interactions

outcome_column <- paste0(snakemake@params[["outcome_column"]])

load(snakemake@input[["file_GRS"]]) #df_gen_score
load(snakemake@input[["file_outcome_phenotype"]]) # outcome phenotype
load(snakemake@input[["file_interaction_phenotype"]])

trait_id <- snakemake@params[["exposure_IDs"]]

df <- merge(df_gen_score, outcome_df, by = "eid")

print(paste0("df_gen_score and outcome_df are merged! N individuals = ",
             nrow(df)))

df <- merge(df, interaction_df, by = "eid")

print(paste0("df and interaction_df are merged! N individuals = ",
             nrow(df)))

GRS_column <- snakemake@params[["exposure_GRS_column"]]
GRS2_column <- paste0(GRS_column, "_2")

df[,GRS2_column] <- df[,GRS_column]^2

interaction_column <- paste0(snakemake@params[["interaction_ID"]], "_z_c")
interaction2_column <- paste0(interaction_column, "_2")
df[,interaction2_column] <- df[,interaction_column]^2


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
        
    if(snakemake@params[["interaction_ID"]] != "21022"){
            # include interaction Y ~ GRSX + GRSX * Age + Age + Age2 + Sex
        model_1 <- lm(df[,outcome_column] ~ df[,GRS_column] * 
                                            df[, interaction_column] + 
                                            df[,interaction2_column] +
                                            df$age_exact +
                                            df$age_exact2 + 
                                            df$sex) 

        print("Model 1 is fitted")

        # include both Y ~ GRSX + GRSX2 + GRSX * Age + Age + Age2 + Sex
        model_complete <- lm(df[,outcome_column] ~ df[,GRS_column] * 
                                                   df[,interaction_column] + 
                                                   df[,interaction2_column] +
                                                   df[,GRS2_column] +
                                                   df$age_exact +
                                                   df$age_exact2 + 
                                                   df$sex)
        print("Full model is fitted")

    } else {
        model_1 <- lm(df[,outcome_column] ~ df[,GRS_column] * 
                                            df[, interaction_column] + 
                                            df[,interaction2_column] +
                                            df$sex)

        print("Model 1 is fitted")

        # include both Y ~ GRSX + GRSX2 + GRSX * Age + Age + Age2 + Sex
        model_complete <- lm(df[,outcome_column] ~ df[,GRS_column] * 
                                                   df[,interaction_column] + 
                                                   df[,interaction2_column] +
                                                   df[,GRS2_column] +
                                                   df$sex) 
        print("Full model is fitted")
        }

        
} else {
    if(snakemake@params[["interaction_ID"]] != "21022"){
        # include interaction Y ~ GRSX + GRSX * Age + Age + Age2 + Sex
        model_1 <- lm(df[,outcome_column] ~ df[,GRS_column] * 
                                            df[,interaction_column] + 
                                            df[,interaction2_column] +
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex + X)

        print("Model 1 is fitted")

        # include both Y ~ GRSX + GRSX2 + GRSX * Age + Age + Age2 + Sex
        model_complete <- lm(df[,outcome_column] ~ df[,GRS_column] * 
                                                   df[,interaction_column] + 
                                                   df[,interaction2_column] +
                                                   df[,GRS2_column] +
                                                   df$age_exact + 
                                                   df$age_exact2 + 
                                                   df$sex + X)
        print("Full model is fitted")
    } else {
        model_1 <- lm(df[,outcome_column] ~ df[,GRS_column] * 
                                            df[,interaction_column] + 
                                            df[,interaction2_column] +
                                            df$sex + X)

        print("Model 1 is fitted")

        # include both Y ~ GRSX + GRSX2 + GRSX * Age + Age + Age2 + Sex
        model_complete <- lm(df[,outcome_column] ~ df[,GRS_column] * 
                                                   df[,interaction_column] +
                                                   df[,interaction2_column] + 
                                                   df[,GRS2_column] +
                                                   df$sex + X)
        print("Full model is fitted")
    }
}

all_models <- list(model_1, model_complete)
save(all_models, file = snakemake@output[["mr_models"]])














