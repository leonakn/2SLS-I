# lknusel

# obtain the ratio between the observed GRS -> Y effects and the observed GRS -> X effects to
# have an understanding of the causal effect size of X -> Y
library(stringr)
library(dplyr)

path_to_results <- snakemake@params[["path_to_results"]]


df_ids <- read.csv(snakemake@input[["df_ids"]], header = T)
outcome_overview <- read.csv(snakemake@input[["outcome_overview"]],
                             header = T)
interaction_overview <- read.csv(snakemake@input[["interaction_overview"]],
                             header = T)


###########################################################
# GET CAUSAL MAIN EFFECTS AND CAUSAL QUADRATIC EFFECTS
##########################################################
df_causal_main_effects <- data.frame(GRS_Y_Estimate = NA,
                                     GRS_Y_pval = NA,
                                     GRS_Y_SE = NA,
                                     GRS_Y_tval = NA,
                                     GRS_X_Estimate = NA,
                                     X_ID = NA,
                                     Y_ID = NA,
                                     degrees_of_freedom = NA
                                     )

df_causal_quadratic <- data.frame(GRS2_Y_Estimate = NA,
                                          GRS2_Y_pval = NA,
                                          GRS2_Y_SE = NA,
                                          GRS2_Y_tval = NA,
                                          GRS_X_Estimate = NA,
                                          X_ID = NA,
                                          Y_ID = NA,
                                          degrees_of_freedom = NA
                                          )


# columns to add:
# X_Y_causal_Estimate
# X_Y_causal_SE
# X_Y_causal_tval --> should be the same as the GRS_X_tval
# X_Y_causal_pval --> should be the same as the GRS_X_pval

files_MR_models <- snakemake@input[["MR_1_files"]]

for(i in 1:length(files_MR_models)){

    ##### FIRST: MAIN EFFECT!
    load(files_MR_models[i])
    model1 <- all_models[[1]] |> summary()
    model2 <- all_models[[2]] |> summary()
    coefs_MR_1 <- model1 |> coef()

    # get variable IDS
    
    df_causal_main_effects[i, "X_ID"] <- 
        str_remove_all(files_MR_models[i], 
                        paste0(path_to_results, "|012_MR_I_Expo_|_Outc_[0-9]+.Rdata"))
    
    df_causal_main_effects[i, "Y_ID"] <- 
        str_remove_all(files_MR_models[i], 
                        paste0(path_to_results, "|012_MR_I_Expo_[0-9]+_Outc_|.Rdata"))

    
    # get GRS on Y effects
    df_causal_main_effects[i, "GRS_Y_Estimate"] <- 
        coefs_MR_1["df[, GRS_column]", "Estimate"]

    df_causal_main_effects[i, "GRS_Y_SE"] <- 
        coefs_MR_1["df[, GRS_column]", "Std. Error"] 

    df_causal_main_effects[i, "GRS_Y_tval"] <- 
        coefs_MR_1["df[, GRS_column]", "t value"] 

    df_causal_main_effects[i, "GRS_Y_pval"] <- 
        coefs_MR_1["df[, GRS_column]", "Pr(>|t|)"]
    
    # get GRS on X effects
    file_GRS_X_model <- paste0(path_to_results, "010_Model_GRS_on_pheno_", 
                                  df_causal_main_effects[i, "X_ID"],
                                  ".Rdata")
    load(file_GRS_X_model)
    coefs_GRS_X <- coef(GRS_on_pheno_model)
    df_causal_main_effects[i, "GRS_X_Estimate"] <- 
        coefs_GRS_X["df_gen_score[, colname_pheno_GRS]", "Estimate"]
    
    df_causal_main_effects[i, "degrees_of_freedom"] <- 
       model1$df[1] + model1$df[2]

    ##### SECOND: QUADRATIC EFFECT
    coefs_quadratic <- model2 |> coef()
    # get variable IDS
    
    df_causal_quadratic[i, "X_ID"] <- 
        str_remove_all(files_MR_models[i], 
                        paste0(path_to_results, "|012_MR_I_Expo_|_Outc_[0-9]+.Rdata"))
    
    df_causal_quadratic[i, "Y_ID"] <- 
        str_remove_all(files_MR_models[i], 
                        paste0(path_to_results, "|012_MR_I_Expo_[0-9]+_Outc_|.Rdata"))

    # get GRS on Y effects
    df_causal_quadratic[i, "GRS2_Y_Estimate"] <- 
        coefs_quadratic["df[, GRS2_column]", "Estimate"]

    df_causal_quadratic[i, "GRS2_Y_SE"] <- 
        coefs_quadratic["df[, GRS2_column]", "Std. Error"] 

    df_causal_quadratic[i, "GRS2_Y_tval"] <- 
        coefs_quadratic["df[, GRS2_column]", "t value"] 

    df_causal_quadratic[i, "GRS2_Y_pval"] <- 
        coefs_quadratic["df[, GRS2_column]", "Pr(>|t|)"]
    
    df_causal_quadratic[i, "GRS_X_Estimate"] <- 
        coefs_GRS_X["df_gen_score[, colname_pheno_GRS]", "Estimate"]
    
    df_causal_quadratic[i, "degrees_of_freedom"] <- 
       model2$df[1] + model2$df[2]

}


##### FINALIZE df_causal_main_effects
df_causal_main_effects$X_Y_causal_Estimate <- 
    df_causal_main_effects$GRS_Y_Estimate /
    df_causal_main_effects$GRS_X_Estimate

df_causal_main_effects$X_Y_causal_SE <- 
   df_causal_main_effects$GRS_Y_SE /
   df_causal_main_effects$GRS_X_Estimate

df_causal_main_effects$X_Y_causal_tval <-
   df_causal_main_effects$X_Y_causal_Estimate/
   df_causal_main_effects$X_Y_causal_SE

df_causal_main_effects$X_Y_causal_pval <- 
   2 * pt(-abs(df_causal_main_effects$X_Y_causal_tval), 
          df_causal_main_effects$degrees_of_freedom)

df_causal_main_effects$X_Y_causal_LCI <-
   df_causal_main_effects$X_Y_causal_Estimate - 
   (df_causal_main_effects$X_Y_causal_SE * qnorm(0.975))

df_causal_main_effects$X_Y_causal_UCI <-
   df_causal_main_effects$X_Y_causal_Estimate + 
   (df_causal_main_effects$X_Y_causal_SE * qnorm(0.975))

df_causal_main_effects <- merge(df_causal_main_effects,
                                df_ids[,c("X_ID", "X_trait")],
                                by = "X_ID")

# df_causal_main_effects <- dplyr::rename(df_causal_main_effects,
#                                         X_trait = Abbreviation)


df_causal_main_effects <- merge(df_causal_main_effects,
                                outcome_overview, 
                                by.x = "Y_ID",
                                by.y = "Y_ID")

df_causal_main_effects <- dplyr::rename(df_causal_main_effects,
                                        Y_trait = "Outcome")


save(df_causal_main_effects, 
     file = snakemake@output[["file_causal_main_effects"]])

##### FINALIZE df_causal_quadratic

df_causal_quadratic$X2_Y_causal_Estimate <- 
    df_causal_quadratic$GRS2_Y_Estimate /
    (df_causal_quadratic$GRS_X_Estimate)^2

df_causal_quadratic$X2_Y_causal_SE <-
    df_causal_quadratic$GRS2_Y_SE /
    (df_causal_quadratic$GRS_X_Estimate)^2

df_causal_quadratic$X2_Y_causal_tval <-
    df_causal_quadratic$X2_Y_causal_Estimate/
    df_causal_quadratic$X2_Y_causal_SE

df_causal_quadratic$X2_Y_causal_pval <-
    2 * pt(-abs(df_causal_quadratic$X2_Y_causal_tval),
                df_causal_quadratic$degrees_of_freedom)

# df_causal_quadratic <- merge(df_causal_quadratic,
#                                      df_ids[,c("Field_ID", "Abbreviation")],
#                                      by.x = "X_ID", by.y = "Field_ID")
df_causal_quadratic <- merge(df_causal_quadratic,
                                     df_ids[,c("X_ID", "X_trait")],
                                     by = "X_ID")
# df_causal_quadratic <- dplyr::rename(df_causal_quadratic,
#                                              X_trait = Abbreviation)

df_causal_quadratic <- merge(df_causal_quadratic,
                                     outcome_overview, 
                                     by.x = "Y_ID",
                                     by.y = "Y_ID")

df_causal_quadratic <- dplyr::rename(df_causal_quadratic,
                                             Y_trait = "Outcome")

save(df_causal_quadratic, 
     file = snakemake@output[["file_causal_quadratic_effects"]])


df_causal_interaction <- data.frame(GRSxE_Y_Estimate = NA,
                                    GRSxE_Y_pval = NA,
                                    GRSxE_Y_SE = NA,
                                    GRSxE_Y_tval = NA,
                                    GRS_X_Estimate = NA,
                                    X_ID = NA,
                                    Y_ID = NA,
                                    E_ID = NA,
                                    degrees_of_freedom = NA)


files_MR_interaction_models <- snakemake@input[["MR_2_files"]]

for(i in 1:length(files_MR_interaction_models)){
    load(files_MR_interaction_models[i])
    model2 <- all_models[[2]]
    coefs_MR_2 <- model2 |> coef()

    df_causal_interaction[i, "X_ID"] <- 
        str_remove_all(files_MR_interaction_models[i],
        paste0(path_to_results, "|013_MR_II_Expo_|_Outc_[0-9]+_interact_[0-9]+.Rdata"))
    
    df_causal_interaction[i, "Y_ID"] <- 
        str_remove_all(files_MR_interaction_models[i], 
        paste0(path_to_results, "|013_MR_II_Expo_[0-9]+_Outc_|_interact_[0-9]+.Rdata"))
    
    df_causal_interaction[i, "E_ID"] <- 
        str_remove_all(files_MR_interaction_models[i], 
        paste0(path_to_results, "|013_MR_II_Expo_[0-9]+_Outc_[0-9]+_interact_|.Rdata"))
    
    df_causal_interaction[i, "GRSxE_Y_Estimate"] <- 
        coefs_MR_2["df[, GRS_column]:df[, interaction_column]", "Estimate"]
    
    df_causal_interaction[i, "GRSxE_Y_SE"] <- 
        coefs_MR_2["df[, GRS_column]:df[, interaction_column]", "Std. Error"]
    
    df_causal_interaction[i, "GRSxE_Y_tval"] <- 
        coefs_MR_2["df[, GRS_column]:df[, interaction_column]", "t value"]
    
    df_causal_interaction[i, "GRSxE_Y_pval"] <- 
        coefs_MR_2["df[, GRS_column]:df[, interaction_column]", "Pr(>|t|)"]
    
    df_causal_interaction[i, "degrees_of_freedom"] <- model2$df[1] + 
                                                      model2$df[2]
    
    # get GRS on X effect:
    file_GRS_X_model <- paste0(path_to_results, "010_Model_GRS_on_pheno_", 
                                  df_causal_interaction[i, "X_ID"],
                                  ".Rdata")
    load(file_GRS_X_model)
    coefs_GRS_X <- coef(GRS_on_pheno_model)
    df_causal_interaction[i, "GRS_X_Estimate"] <- 
        coefs_GRS_X["df_gen_score[, colname_pheno_GRS]", "Estimate"]

}

df_causal_interaction$XxE_causal_Estimate <- 
    df_causal_interaction$GRSxE_Y_Estimate/
    df_causal_interaction$GRS_X_Estimate

df_causal_interaction$XxE_causal_SE <- 
    df_causal_interaction$GRSxE_Y_SE/
    df_causal_interaction$GRS_X_Estimate

df_causal_interaction$XxE_causal_tval <- 
    df_causal_interaction$XxE_causal_Estimate/
    df_causal_interaction$XxE_causal_SE

df_causal_interaction$XxE_causal_pval <- 
    2 * pt(-abs(df_causal_interaction$XxE_causal_tval),
                df_causal_interaction$degrees_of_freedom)


### Add Trait names in addition to trait Field IDs:
df_causal_interaction <- merge(df_causal_interaction,
                               df_ids[,c("X_ID", "X_trait")],
                               by = "X_ID")
# df_causal_interaction <- dplyr::rename(df_causal_interaction,
#                                        X_trait = Abbreviation)

df_causal_interaction <- merge(df_causal_interaction,
                               outcome_overview, 
                               by.x = "Y_ID",
                               by.y = "Y_ID")

df_causal_interaction <- dplyr::rename(df_causal_interaction,
                                             Y_trait = "Outcome")


df_causal_interaction <- merge(df_causal_interaction, 
                               interaction_overview, 
                               by.x = "E_ID",
                               by.y = "E_ID")
df_causal_interaction <- dplyr::rename(df_causal_interaction,
                                       E_trait = Interaction)

# Also, obtain corrected interaction:

df_causal_interaction$corrected_XxE_estimate <- NA
df_causal_interaction$corrected_XxE_SE <- NA
df_causal_interaction$corrected_XxE_tval <- NA
df_causal_interaction$corrected_XxE_pval <- NA 


for(i in 1:nrow(df_causal_interaction)){
    X_ID <- df_causal_interaction[i, "X_ID"]
    Y_ID <- df_causal_interaction[i, "Y_ID"]
    E_ID <- df_causal_interaction[i, "E_ID"]

    raw_GRSxE_estimate <- df_causal_interaction[i, "GRSxE_Y_Estimate"]
    raw_GRSxE_se <- df_causal_interaction[i, "GRSxE_Y_SE"]
    raw_GRSxE_n <- df_causal_interaction[i, "degrees_of_freedom"]
    raw_GRSxE_var <- (raw_GRSxE_se * sqrt(raw_GRSxE_n))^2

    # GET GRS ON X EFFECT
    load(paste0(path_to_results, "010_Model_GRS_on_pheno_", X_ID, ".Rdata"))

    GRS_x_estimate <- 
        GRS_on_pheno_model$coefficients["df_gen_score[, colname_pheno_GRS]", "Estimate"]
    GRS_x_se <- 
        GRS_on_pheno_model$coefficients["df_gen_score[, colname_pheno_GRS]", "Std. Error"]
    GRS_x_n <- GRS_on_pheno_model$df[1] + GRS_on_pheno_model$df[2]
    GRS_x_variance <- (GRS_x_se * sqrt(GRS_x_n))^2

    # GET X2 ON Y EFFECT
    
    load(paste0(path_to_results, "013_MR_II_Expo_", X_ID, 
                "_Outc_", Y_ID, "_interact_", E_ID, ".Rdata"))
    
    model2 <- all_models[[2]] |> summary()

    GRS2_y_estimate <- model2$coefficients["df[, GRS2_column]","Estimate"]
    GRS2_y_se <- model2$coefficients["df[, GRS2_column]","Std. Error"]
    GRS2_y_n <- raw_GRSxE_n
    GRS2_y_variance <- (GRS2_y_se * sqrt(GRS2_y_n))^2

    # GET E ON X EFFECT
    load(paste0(path_to_results, "007_E_", E_ID, "_X_", X_ID, "_model.Rdata"))
    E_x_estimate <- model_E_on_x$coefficients["interaction_df$E", "Estimate"]
    E_x_se <- model_E_on_x$coefficients["interaction_df$E", "Std. Error"]
    E_x_n <- model_E_on_x$df[1] + model_E_on_x$df[2]
    E_x_variance <- (E_x_se * sqrt(E_x_n))^2

    # correct GRS*E -> Y estimate
    corrected_interaction_estimate <- raw_GRSxE_estimate -
                                      ((GRS2_y_estimate/GRS_x_estimate) *
                                      2 * E_x_estimate)
    # get causal X*E -> Y effect
    corrected_interaction_estimate <- corrected_interaction_estimate/
                                      GRS_x_estimate
    
    # correct GRS*E -> variance
    corrected_interaction_variance <- raw_GRSxE_var +
                                      (GRS2_y_estimate^2/GRS_x_estimate^2) *
                                      ((GRS2_y_variance/GRS2_y_estimate^2) +
                                      (GRS_x_variance/GRS_x_estimate^2)) *
                                       E_x_estimate^2 * 2^2
   
    corrected_interaction_se <- sqrt(corrected_interaction_variance)/
                                sqrt(raw_GRSxE_n)
    corrected_interaction_se <- corrected_interaction_se/GRS_x_estimate
    corrected_interaction_tval <- corrected_interaction_estimate/
                                  corrected_interaction_se

    degrees_of_freedom <- raw_GRSxE_n - 2
    corrected_interaction_pval <- 2 * pt(-abs(corrected_interaction_tval), 
                                              degrees_of_freedom)
    
    df_causal_interaction[i, "corrected_XxE_estimate"] <- 
        corrected_interaction_estimate
    df_causal_interaction[i, "corrected_XxE_SE"] <- 
        corrected_interaction_se
    df_causal_interaction[i, "corrected_XxE_tval"] <- 
        corrected_interaction_tval
    df_causal_interaction[i, "corrected_XxE_pval"] <- 
        corrected_interaction_pval

}

save(df_causal_interaction, file = snakemake@output[["file_causal_interaction_effects"]])




   
      
