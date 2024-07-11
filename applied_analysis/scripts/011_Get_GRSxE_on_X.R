#lknusel
# X ~ GRS*E to see if there is a level 1 interaction

load(snakemake@input[["file_exposure_phenotypes"]])

load(snakemake@input[["file_GRS"]])

load(snakemake@input[["file_interaction_phenotype"]])

df <- merge(pheno_df, df_gen_score, by = "eid")
df <- merge(df, interaction_df, by = "eid")

Y <- paste0(snakemake@params[["exposure_IDs"]], "_IRNT_c")
E <- paste0(snakemake@params[["interaction_ID"]], "_z_c")
GRS_column <- grep("GRS", colnames(df), value = T)

GRSxE_pheno_model <- lm(df[,Y] ~ df[,GRS_column] * df[,E]) |> summary()

print("Model is fitted!")

print(GRSxE_pheno_model)

save(GRSxE_pheno_model, file = snakemake@output[["file_GRSxE_on_pheno_model"]])
