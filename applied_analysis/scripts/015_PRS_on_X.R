# obtain the models to investigate how much variance the PRS explains on X
# plot PRS vs X

library(data.table)
library(R.utils)
library(ggplot2)

# load file with PRS
PRS_overview <- fread(snakemake@input[["file_PRS_overview"]], 
                      header = T, data.table = F)

exposure_ID <- snakemake@params[["exposure_IDs"]]

PRS_overview <- PRS_overview[PRS_overview$X_ID %in% exposure_ID, ]

PRS_df <- fread(snakemake@input[["file_PRS"]], 
                header = T, data.table = F)

# keep only the two columns of interest:
# IID (= participant ID) and relevant PRS
PRS_column <- paste0(PRS_overview$PGS_ID, "_SUM")
PRS_df <- PRS_df[,colnames(PRS_df) %in% c("IID", PRS_column)]
PRS_df$PRS_column <- PRS_df[,PRS_column]

load(snakemake@input[["file_exposure_phenotype"]])

df <- merge(pheno_df, PRS_df, by.x = "eid", by.y = "IID")

print(paste0("PRS_df and pheno_df are merged! N individuals = ",
             nrow(df)))

load(snakemake@input[["file_covariates"]])
load(snakemake@input[["file_medication"]])

##### Correct for medication if appropriate
X_medication <- data.frame("eid" = df$eid)

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

##### obtain covariates and make sure eids are aligned!
# keep: age, age2, sex
# remove sex if its a female_only_trait
if(sum(grepl(exposure_ID, snakemake@params[["female_only_IDs"]])) > 0){
    X <- X[,colnames(X) %in% c("eid", "age_exact", "age_exact2")]
} else {
     X <- X[,colnames(X) %in% c("eid", "sex", "age_exact", "age_exact2")]
}

X <- merge(X, X_medication, by = "eid")

X <- X[X$eid %in% df$eid,]

df <- df[df$eid %in% X$eid, ]

if(!all(df$eid == X$eid)){
    eid_order <- match(df$eid, X$eid)
    X <- X[eid_order, ]
    if(!all(df$eid == X$eid)){
        stop("Something's wrong with the eids, do double check!")
    }
}

X <- X[,colnames(X) != "eid"]

X <- as.matrix(X)

Y <- paste0(snakemake@params[["exposure_IDs"]], "_IRNT_c")
PRS_on_pheno_model <- lm(df[,Y] ~ df[,PRS_column] + X) |> summary()
save(PRS_on_pheno_model, file = snakemake@output[["file_PRS_on_pheno_model"]])

minimal_PRS_on_pheno_model <- lm(df[,Y] ~ df[,PRS_column]) |> summary()
save(minimal_PRS_on_pheno_model, file = snakemake@output[["file_minimal_PRS_on_pheno_model"]])
