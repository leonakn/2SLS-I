# run MR with the extended GRS

library(data.table)

# load file with PRS
PRS_overview <- fread(snakemake@input[["file_PRS_overview"]], header = T, data.table = F)

exposure_ID <- snakemake@params[["exposure_IDs"]]

PRS_overview <- PRS_overview[PRS_overview$X_ID %in% exposure_ID, ]

PRS_df <- fread(snakemake@input[["file_PRS"]], header = T, data.table = F)

#note typo in "phenotype"!
PRS_column <- paste0(PRS_overview$PGS_ID, "_SUM")
#PRS_df <- PRS_df[,colnames(PRS_df) %in% c("IID", PRS_overview[PRS_overview$Field_ID == exposure_ID, "phenotye"])]
PRS_df <- PRS_df[,colnames(PRS_df) %in% c("IID", PRS_column)]

#PRS_column <- PRS_overview[PRS_overview$Field_ID == exposure_ID, "phenotye"]
PRS_2_column <- "PRS_2_column"
PRS_df$PRS_column <- PRS_df[,PRS_column]
PRS_df$PRS_2_column <- PRS_df$PRS_column^2

load(snakemake@input[["file_outcome_phenotype"]])

df <- merge(PRS_df, outcome_df, by.x = "IID", by.y = "eid")

print(paste0("PRS_df and outcome_df are merged! N individuals = ",
             nrow(df)))

outcome_column <- paste0(snakemake@params[["outcome_column"]])

# load medication file
load(snakemake@input[["file_medication"]])

if(!all(medication_df$eid == df$eid)){
    medication_df <- medication_df[medication_df$eid %in% df$eid, ]
    df <- df[df$eid %in% medication_df$eid, ]

    eids_sorting <- match(df$eid, medication_df$eid)
    medication_df <- medication_df[eids_sorting, ]
}

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

if(!all(df$IID == X$eid)){
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
        # Y ~ GRSX + Age + Age2 + Sex
        model_0 <- lm(df[,outcome_column] ~ df[,PRS_column] + 
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex) |> summary()

        print("Model 0 is fitted")

        # include quadratic effect Y ~ GRSX + GRSX2 + Age + Age2 + Sex
        model_1 <- lm(df[,outcome_column] ~ df[,PRS_column] + 
                                            df[,PRS_2_column] +
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex) |> summary()

        print("Model 1 is fitted")

    } else {
        # Y ~ GRSX + Age + Age2 + Sex
        model_0 <- lm(df[,outcome_column] ~ df[,PRS_column] + 
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex + X) |> summary()

        print("Model 0 is fitted")

        # include quadratic effect Y ~ GRSX + GRSX2 + Age + Age2 + Sex
        model_1 <- lm(df[,outcome_column] ~ df[,PRS_column] + 
                                            df[,PRS_2_column] +
                                            df$age_exact + 
                                            df$age_exact2 + 
                                            df$sex + X) |> summary()

        print("Model 1 is fitted")

}

all_models <- list(model_0, model_1)
save(all_models, file = snakemake@output[["mr_models_PRS"]])


