# Fit MR models including the interaction parameter
# For extended GRS

library(data.table)
library(R.utils)

exposure_ID <- snakemake@params[["exposure_IDs"]]
interaction_column <- paste0(snakemake@params[["interaction_ID"]], "_z_c")
outcome_column <- paste0(snakemake@params[["outcome_column"]])


PRS_overview <- fread(snakemake@input[["file_PRS_overview"]], header = T, data.table = F)
PRS_overview <- PRS_overview[PRS_overview$X_ID %in% exposure_ID, ]
PRS_df <- fread(snakemake@input[["file_PRS"]], header = T, data.table = F)

PRS_column <- paste0(PRS_overview$PGS_ID, "_SUM")
# PRS_df <- PRS_df[,colnames(PRS_df) %in% c("IID", PRS_overview[PRS_overview$Field_ID == exposure_ID, "phenotye"])]
PRS_df <- PRS_df[,colnames(PRS_df) %in% c("IID", PRS_column)]

#PRS_column <- PRS_overview[PRS_overview$Field_ID == exposure_ID, "phenotye"]
PRS_2_column <- "PRS_2_column"
PRS_df$PRS_column <- PRS_df[,PRS_column]
PRS_df$PRS_2_column <- PRS_df$PRS_column^2

load(snakemake@input[["file_outcome_phenotype"]]) # outcome phenotype
load(snakemake@input[["file_interaction_phenotype"]])
interaction2_column <- paste0(interaction_column, "_2")
interaction_df[,interaction2_column] <- interaction_df[,interaction_column]^2

# merge all the dfs to one:
df <- merge(PRS_df, outcome_df, by.x = "IID", by.y = "eid")

print(paste0("PRS_df and outcome_df are merged! N individuals = ",
             nrow(df)))

df <- merge(interaction_df, df, by.x = "eid", by.y = "IID")

print(paste0("df and interaction_df are merged! N individuals = ",
             nrow(df)))

# Correct for medication
load(snakemake@input[["file_medication"]])

X <- data.frame("eid" = medication_df$eid)

# Insulin:
if(exposure_ID %in% snakemake@params[["correct_for_insulin"]]){
	X <- merge(X, medication_df[,c("eid", "Insulin")], by = "eid")
}

# Blood pressure medication
if(exposure_ID %in% snakemake@params[["correct_for_bp_medication"]]){
	X <- merge(X, medication_df[,c("eid", "bp_regul")], by = "eid")
}

# Cholesterol lowering medication
if(exposure_ID %in% snakemake@params[["correct_for_cholesterol_lowering"]]){
	X <- merge(X, medication_df[,c("eid", "chol_inhibit")], by = "eid")
}

# Hormone replacement therapy
if(exposure_ID %in% snakemake@params[["correct_for_hormone_replacement"]]){
	X <- merge(X, medication_df[,c("eid", "hormone_repl")], by = "eid")
}

# Oral contraceptive
if(exposure_ID %in% snakemake@params[["correct_for_oral_contraceptive"]]){
	X <- merge(X, medication_df[,c("eid", "oral_contraceptive")], by = "eid")
}


if(!all(df$eid == X$eid)){
    eid_order <- match(df$eid, X$eid)
    X <- X[eid_order, ]
    if(!all(df$eid == X$eid)){
        stop("Something's wrong with the eids, do double check!")
    }
}

X <- X[,colnames(X) != "eid"]

X <- as.matrix(X)

if(! exposure_ID %in% c(snakemake@params[["correct_for_insulin"]], snakemake@params[["correct_for_bp_medication"]], 
                     snakemake@params[["correct_for_cholesterol_lowering"]], snakemake@params[["correct_for_hormone_replacement"]],
                     snakemake@params[["correct_for_oral_contraceptive"]])){
        
    if(snakemake@params[["interaction_ID"]] != "21022"){
            # include interaction Y ~ PRSX + PRSX * Age + Age + Age2 + Sex
        model_1 <- lm(df[,outcome_column] ~ df[,PRS_column] * 
                                            df[, interaction_column] + 
                                            df[, interaction2_column] +
                                            df$age_exact +
                                            df$age_exact2 + 
                                            df$sex) |> summary()

        print("Model 1 is fitted")

        # include both Y ~ PRSX + PRSX2 + PRSX * Age + Age + Age2 + Sex
        model_complete <- lm(df[,outcome_column] ~ df[,PRS_column] * 
                                                   df[,interaction_column] + 
                                                   df[, interaction2_column] +
                                                   df[,PRS_2_column] +
                                                   df$age_exact +
                                                   df$age_exact2 + 
                                                   df$sex) |> summary()
        print("Full model is fitted")

    } else {
        model_1 <- lm(df[,outcome_column] ~ df[,PRS_column] * 
                                            df[, interaction_column] + 
                                            df[, interaction2_column] +
                                            df$sex) |> summary()

        print("Model 1 is fitted")

        # include both Y ~ PRSX + PRSX2 + PRSX * Age + Age + Age2 + Sex
        model_complete <- lm(df[,outcome_column] ~ df[,PRS_column] * 
                                                   df[,interaction_column] + 
                                                   df[, interaction2_column] +
                                                   df[,PRS_2_column] +
                                                   df$sex) |> summary()
        print("Full model is fitted")
        }

        
} else {
    if(snakemake@params[["interaction_ID"]] != "21022"){
        # include interaction Y ~ PRSX + PRSX * Age + Age + Age2 + Sex
        model_1 <- lm(df[,outcome_column] ~ df[,PRS_column] * 
                                            df[,interaction_column] + 
                                            df[, interaction2_column] +
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex + X) |> summary()

        print("Model 1 is fitted")

        # include both Y ~ PRSX + PRSX2 + PRSX * Age + Age + Age2 + Sex
        model_complete <- lm(df[,outcome_column] ~ df[,PRS_column] * 
                                                   df[,interaction_column] + 
                                                   df[, interaction2_column] +
                                                   df[,PRS_2_column] +
                                                   df$age_exact + 
                                                   df$age_exact2 + 
                                                   df$sex + X) |> summary()
        print("Full model is fitted")
    } else {
        model_1 <- lm(df[,outcome_column] ~ df[,PRS_column] * 
                                            df[,interaction_column] + 
                                            df$sex + X) |> summary()

        print("Model 1 is fitted")

        # include both Y ~ PRSX + PRSX2 + PRSX * Age + Age + Age2 + Sex
        model_complete <- lm(df[,outcome_column] ~ df[,PRS_column] * 
                                                   df[,interaction_column] + 
                                                   df[, interaction2_column] +
                                                   df[,PRS_2_column] + 
                                                   df$sex + X) |> summary()
        print("Full model is fitted")
    }
}

all_models <- list(model_1, model_complete)
save(all_models, file = snakemake@output[["mr_models"]])
