# obtain summary files for main effects, quadratic main effects, 
# and interaction effects (both raw and corrected)
library(stringr)
library(dplyr)
library(data.table)

path_to_results <- snakemake@params[["path_to_results"]]

df_ids <- fread(snakemake@input[["df_ids"]], header = T, data.table = F)
outcome_overview <- fread(snakemake@input[["outcome_overview"]],
                             header = T, data.table = F)
interaction_overview <- fread(snakemake@input[["interaction_overview"]],
                             header = T, data.table = F)
df_causal_main_effects <- data.frame(PRS_Y_Estimate = NA,
                                     PRS_Y_pval = NA,
                                     PRS_Y_SE = NA,
                                     PRS_Y_tval = NA,
                                     PRS_X_Estimate = NA,
                                     X_ID = NA,
                                     Y_ID = NA,
                                     degrees_of_freedom = NA
                                     )

df_causal_quadratic <- data.frame(PRS2_Y_Estimate = NA,
                                          PRS2_Y_pval = NA,
                                          PRS2_Y_SE = NA,
                                          PRS2_Y_tval = NA,
                                          PRS_X_Estimate = NA,
                                          X_ID = NA,
                                          Y_ID = NA,
                                          degrees_of_freedom = NA
                                          )


files_MR_models <- snakemake@input[["MR_1_files"]]

for(i in 1:length(files_MR_models)){

    ##### FIRST: MAIN EFFECT!
    load(files_MR_models[i])
    model1 <- all_models[[1]] 
    coefs_MR_1 <-  model1 |> coef()

    # get variable IDS
    
    df_causal_main_effects[i, "X_ID"] <- 
        str_remove_all(files_MR_models[i], 
                        paste0(path_to_results, "|016_MR_I_PRS_Expo_|_Outc_[0-9]+.Rdata"))
    
    df_causal_main_effects[i, "Y_ID"] <- 
        str_remove_all(files_MR_models[i], 
                        paste0(path_to_results, "|016_MR_I_PRS_Expo_[0-9]+_Outc_|.Rdata"))

    
    # get PRS on Y effects
    df_causal_main_effects[i, "PRS_Y_Estimate"] <- 
        coefs_MR_1["df[, PRS_column]", "Estimate"]

    df_causal_main_effects[i, "PRS_Y_SE"] <- 
        coefs_MR_1["df[, PRS_column]", "Std. Error"] 

    df_causal_main_effects[i, "PRS_Y_tval"] <- 
        coefs_MR_1["df[, PRS_column]", "t value"] 

    df_causal_main_effects[i, "PRS_Y_pval"] <- 
        coefs_MR_1["df[, PRS_column]", "Pr(>|t|)"]
    
    # get PRS on X effects
    file_PRS_X_model <- paste0(path_to_results, "015_Model_PRS_on_pheno_", 
                                  df_causal_main_effects[i, "X_ID"],
                                  ".Rdata")
    load(file_PRS_X_model)
    coefs_PRS_X <- coef(PRS_on_pheno_model)
    df_causal_main_effects[i, "PRS_X_Estimate"] <- 
        coefs_PRS_X["df[, PRS_column]", "Estimate"]
    
    df_causal_main_effects[i, "degrees_of_freedom"] <- 
       model1$df[1] + model1$df[2]

    ##### SECOND: QUADRATIC EFFECT
    model2 <- all_models[[2]]
    coefs_quadratic <- model2 |> coef()
    # get variable IDS
 
        df_causal_quadratic[i, "X_ID"] <- 
            str_remove_all(files_MR_models[i], 
                           paste0(path_to_results, "|016_MR_I_PRS_Expo_|_Outc_[0-9]+.Rdata"))
    
        df_causal_quadratic[i, "Y_ID"] <- 
            str_remove_all(files_MR_models[i], 
                           paste0(path_to_results, "|016_MR_I_PRS_Expo_[0-9]+_Outc_|.Rdata"))

    
    # get PRS on Y effects
    df_causal_quadratic[i, "PRS2_Y_Estimate"] <- 
        coefs_quadratic["df[, PRS_2_column]", "Estimate"]

    df_causal_quadratic[i, "PRS2_Y_SE"] <- 
        coefs_quadratic["df[, PRS_2_column]", "Std. Error"] 

    df_causal_quadratic[i, "PRS2_Y_tval"] <- 
        coefs_quadratic["df[, PRS_2_column]", "t value"] 

    df_causal_quadratic[i, "PRS2_Y_pval"] <- 
        coefs_quadratic["df[, PRS_2_column]", "Pr(>|t|)"]
    
    df_causal_quadratic[i, "PRS_X_Estimate"] <- 
        coefs_PRS_X["df[, PRS_column]", "Estimate"]
    
    df_causal_quadratic[i, "degrees_of_freedom"] <- 
       model2$df[1] + model2$df[2]

}


##### FINALIZE df_causal_main_effects
df_causal_main_effects$X_Y_causal_Estimate <- 
    df_causal_main_effects$PRS_Y_Estimate /
    df_causal_main_effects$PRS_X_Estimate

df_causal_main_effects$X_Y_causal_SE <- 
   df_causal_main_effects$PRS_Y_SE /
   df_causal_main_effects$PRS_X_Estimate

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
                                by = "X_ID", all.x = T)

df_causal_main_effects <- merge(df_causal_main_effects,
                                outcome_overview, 
                                by.x = "Y_ID",
                                by.y = "Y_ID")


save(df_causal_main_effects, 
     file = snakemake@output[["file_PRS_main_effects"]])

print("df main effects saved!")

##### FINALIZE df_causal_quadratic

df_causal_quadratic$X2_Y_causal_Estimate <- 
    df_causal_quadratic$PRS2_Y_Estimate /
    (df_causal_quadratic$PRS_X_Estimate)^2

df_causal_quadratic$X2_Y_causal_SE <-
    df_causal_quadratic$PRS2_Y_SE /
    (df_causal_quadratic$PRS_X_Estimate)^2

df_causal_quadratic$X2_Y_causal_tval <-
    df_causal_quadratic$X2_Y_causal_Estimate/
    df_causal_quadratic$X2_Y_causal_SE

df_causal_quadratic$X2_Y_causal_pval <-
    2 * pt(-abs(df_causal_quadratic$X2_Y_causal_tval),
                df_causal_quadratic$degrees_of_freedom)

df_causal_quadratic <- merge(df_causal_quadratic,
                             df_ids[,c("X_ID", "X_trait")],
                             by = "X_ID",
                             all.x = T)

df_causal_quadratic <- merge(df_causal_quadratic,
                             outcome_overview, 
                             by.x = "Y_ID",
                             by.y = "Y_ID")


save(df_causal_quadratic, 
     file = snakemake@output[["file_PRS_quadratic_effects"]])

print("df quadratic effects saved!")

df_causal_interaction <- data.frame(PRSxE_Y_Estimate = NA,
                                    PRSxE_Y_pval = NA,
                                    PRSxE_Y_SE = NA,
                                    PRSxE_Y_tval = NA,
                                    PRS_X_Estimate = NA,
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
                       paste0(path_to_results, "|017_MR_II_PRS_Expo_|_Outc_[0-9]+_interact_[0-9]+.Rdata"))
    
    df_causal_interaction[i, "Y_ID"] <- 
        str_remove_all(files_MR_interaction_models[i], 
                           paste0(path_to_results, "|017_MR_II_PRS_Expo_[0-9]+_Outc_|_interact_[0-9]+.Rdata"))
    
    df_causal_interaction[i, "E_ID"] <- 
        str_remove_all(files_MR_interaction_models[i], 
                           paste0(path_to_results, "|017_MR_II_PRS_Expo_[0-9]+_Outc_[0-9]+_interact_|.Rdata"))
    
    df_causal_interaction[i, "PRSxE_Y_Estimate"] <- 
        coefs_MR_2["df[, PRS_column]:df[, interaction_column]", "Estimate"]
    
    df_causal_interaction[i, "PRSxE_Y_SE"] <- 
        coefs_MR_2["df[, PRS_column]:df[, interaction_column]", "Std. Error"]
    
    df_causal_interaction[i, "PRSxE_Y_tval"] <- 
        coefs_MR_2["df[, PRS_column]:df[, interaction_column]", "t value"]
    
    df_causal_interaction[i, "PRSxE_Y_pval"] <- 
        coefs_MR_2["df[, PRS_column]:df[, interaction_column]", "Pr(>|t|)"]
    
    df_causal_interaction[i, "degrees_of_freedom"] <- model2$df[1] + 
                                                      model2$df[2]
    
    # get PRS on X effect:
    file_PRS_X_model <- paste0(path_to_results, "015_Model_PRS_on_pheno_", 
                               df_causal_interaction[i, "X_ID"],
                               ".Rdata")
    load(file_PRS_X_model)
    coefs_PRS_X <- coef(PRS_on_pheno_model)
    df_causal_interaction[i, "PRS_X_Estimate"] <- 
        coefs_PRS_X["df[, PRS_column]", "Estimate"]

}

df_causal_interaction$XxE_causal_Estimate <- 
    df_causal_interaction$PRSxE_Y_Estimate/
    df_causal_interaction$PRS_X_Estimate

df_causal_interaction$XxE_causal_SE <- 
    df_causal_interaction$PRSxE_Y_SE/
    df_causal_interaction$PRS_X_Estimate

df_causal_interaction$XxE_causal_tval <- 
    df_causal_interaction$XxE_causal_Estimate/
    df_causal_interaction$XxE_causal_SE

df_causal_interaction$XxE_causal_pval <- 
    2 * pt(-abs(df_causal_interaction$XxE_causal_tval),
                df_causal_interaction$degrees_of_freedom)


df_causal_interaction <- merge(df_causal_interaction,
                               df_ids[,c("X_ID", "X_trait")],
                               by = "X_ID",
                               all.x = T)

df_causal_interaction <- merge(df_causal_interaction,
                               outcome_overview, 
                               by.x = "Y_ID",
                               by.y = "Y_ID")


df_causal_interaction <- merge(df_causal_interaction, 
                               interaction_overview, 
                               by.x = "E_ID",
                               by.y = "E_ID")


# Also, obtain corrected interaction:

df_causal_interaction$corrected_XxE_estimate <- NA
df_causal_interaction$corrected_XxE_SE <- NA
df_causal_interaction$corrected_XxE_tval <- NA
df_causal_interaction$corrected_XxE_pval <- NA 


for(i in 1:nrow(df_causal_interaction)){
    X_ID <- df_causal_interaction[i, "X_ID"]
    Y_ID <- df_causal_interaction[i, "Y_ID"]
    E_ID <- df_causal_interaction[i, "E_ID"]

    raw_PRSxE_estimate <- df_causal_interaction[i, "PRSxE_Y_Estimate"]
    raw_PRSxE_se <- df_causal_interaction[i, "PRSxE_Y_SE"]
    raw_PRSxE_n <- df_causal_interaction[i, "degrees_of_freedom"]
    raw_PRSxE_var <- (raw_PRSxE_se * sqrt(raw_PRSxE_n))^2

    # GET PRS ON X EFFECT
    load(paste0(path_to_results, "015_Model_PRS_on_pheno_", 
                X_ID, ".Rdata"))

    PRS_x_estimate <- 
        PRS_on_pheno_model$coefficients["df[, PRS_column]", "Estimate"]
    PRS_x_se <- 
        PRS_on_pheno_model$coefficients["df[, PRS_column]", "Std. Error"]
    PRS_x_n <- PRS_on_pheno_model$df[1] + PRS_on_pheno_model$df[2]
    PRS_x_variance <- (PRS_x_se * sqrt(PRS_x_n))^2

    
    load(paste0(path_to_results, "017_MR_II_PRS_Expo_", X_ID, 
                "_Outc_", Y_ID, "_interact_", E_ID, ".Rdata"))
    

    PRS2_y_estimate <- model2$coefficients["df[, PRS_2_column]","Estimate"]
    PRS2_y_se <- model2$coefficients["df[, PRS_2_column]","Std. Error"]
    PRS2_y_n <- raw_PRSxE_n
    PRS2_y_variance <- (PRS2_y_se * sqrt(PRS2_y_n))^2

    # GET E ON X EFFECT
    load(paste0(path_to_results, "007_E_", E_ID, "_X_", X_ID, "_model.Rdata"))
    E_x_estimate <- model_E_on_x$coefficients["interaction_df$E", "Estimate"]
    E_x_se <- model_E_on_x$coefficients["interaction_df$E", "Std. Error"]
    E_x_n <- model_E_on_x$df[1] + model_E_on_x$df[2]
    E_x_variance <- (E_x_se * sqrt(E_x_n))^2

    # perform correction
    corrected_interaction_estimate <- raw_PRSxE_estimate -
                                      ((PRS2_y_estimate/PRS_x_estimate) *
                                      2 * E_x_estimate)
    
    #obtain causal effect
    corrected_interaction_estimate <- corrected_interaction_estimate/PRS_x_estimate

    # perform correction
    corrected_interaction_variance <- raw_PRSxE_var +
                                      (PRS2_y_estimate^2/PRS_x_estimate^2) *
                                      ((PRS2_y_variance/PRS2_y_estimate^2) +
                                      (PRS_x_variance/PRS_x_estimate^2)) *
                                       E_x_estimate^2 * 2^2
    
   
    # obtain SE, get ratio.
    corrected_interaction_se <- sqrt(corrected_interaction_variance)/
                                sqrt(raw_PRSxE_n)
    corrected_interaction_se <- corrected_interaction_se/PRS_x_estimate
    corrected_interaction_tval <- corrected_interaction_estimate/
                                  corrected_interaction_se

    degrees_of_freedom <- raw_PRSxE_n - 2
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

save(df_causal_interaction, file = snakemake@output[["file_PRS_interaction_effects"]])

