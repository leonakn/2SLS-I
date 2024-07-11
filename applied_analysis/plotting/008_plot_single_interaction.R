# get two plots for each interaction with tier > 0
# first, only split by environment,
# second, split by environment, presenting raw predicted effect and corrected predicted effect

path_to_results <- "" # where your preprocessed phenotypes and fitted MR models are
path_to_data <- "" # where your data frames to connect trait IDs with trait labels is
path_out <- "" # where you want your plots to be saved

# if(!dir.exists(path_out)){
#   dir.create(path_out)
# }

# libraries
library(ggplot2)
library(data.table)
library(MASS)
library(stringr)

# translate Y_ID to Y_trait
outcome_overview <- fread(paste0(path_to_data, "outcome_overview.csv"),
                                 data.table = F)

# translate X_ID to X_trait
df_IDs <- fread(paste0(path_to_data,
                       "Pheno_search_PRS_avail_v2.csv"), 
                data.table = F, 
                select = c("X_ID", "X_trait"))

# Translate E_ID to E_trait
interaction_overview <- fread(paste0(path_to_data, "interaction_overview.csv"),
                                     data.table = F)

# relevant interactions
load(paste0(path_to_results, "019_df_interactions.Rdata")) #df_interaction
load(paste0(path_to_results, "014_causal_main_effects.Rdata")) # df_causal_main_effects
load(paste0(path_to_results, "014_causal_corrected_interaction.Rdata")) # df_causal_interaction
load(paste0(path_to_results, "014_causal_quadratic_effects.Rdata")) #df_causal_quadratic

# only plot relevant interactions:
df_interaction <- df_interaction[df_interaction$tier > 0, ]
df_interaction <- df_interaction[df_interaction$X_ID != df_interaction$Y_ID, ]

obs_to_plot <- df_interaction

print(paste0("number of obs to plot = ", nrow(obs_to_plot)))

for(i in 1:nrow(obs_to_plot)){

    X_ID <- obs_to_plot[i, "X_ID"]
    E_ID <- obs_to_plot[i, "E_ID"]
    Y_ID <- obs_to_plot[i, "Y_ID"]

    X_trait <- df_IDs[df_IDs$X_ID == X_ID, "X_trait"]
    E_trait <- interaction_overview[interaction_overview$E_ID == E_ID, "E_trait"]
    Y_trait <- outcome_overview[outcome_overview$Y_ID == Y_ID, "Y_trait"]


    load(paste0(path_to_results, "004_Exposure_phenotype_", X_ID, ".Rdata")) #pheno_df
    load(paste0(path_to_results, "006_interaction_phenotype_", E_ID, ".Rdata")) #interaction_df
    
    tier <- obs_to_plot[i, "tier"]

    raw_cor <- obs_to_plot[i, "raw_cor_PRS"]

    # Get range GRS
    x_range <- seq(from = min(pheno_df[,paste0(X_ID, "_IRNT_c")], na.rm = T), 
                   to = max(pheno_df[,paste0(X_ID, "_IRNT_c")], na.rm = T),
                   length.out = 1000)
    
    # Get range environment
    e_levels <- quantile(interaction_df[,paste0(E_ID, "_z_c")], 
                         probs = c(0.1, 0.5, 0.9), na.rm = T)

    #load(paste0(path_to_results, "003_Covariates.Rdata"))

    pred_data <- expand.grid(X = x_range, environment = e_levels)

    pred_data$X2 <- pred_data$X^2
    pred_data$environment2 <- pred_data$environment^2

    pred_data$ExX <- pred_data$X * pred_data$environment

    ####################### get model coefficients:
    load(paste0(path_to_results, "013_MR_II_Expo_", X_ID,
                "_Outc_", Y_ID, "_interact_", E_ID, ".Rdata"))
    
    coefs <- model |> summary() |> coefficients()
    # BETAS 
    b_intercept <- coefs["(Intercept)", "Estimate"]
    b_X <- df_causal_main_effects[df_causal_main_effects$X_ID == X_ID &
                                  df_causal_main_effects$Y_ID == Y_ID,
                                  "X_Y_causal_Estimate"]
    b_X2 <- df_causal_quadratic[df_causal_quadratic$X_ID == X_ID &
                                df_causal_quadratic$Y_ID == Y_ID,
                                "X2_Y_causal_Estimate"]
    b_E <- coefs["df[, interaction_column]", "Estimate"]
    b_E2 <- coefs["df[, interaction2_column]", "Estimate"]
    b_X_E_raw <- df_causal_interaction[df_causal_interaction$X_ID == X_ID &
                                       df_causal_interaction$Y_ID == Y_ID &
                                       df_causal_interaction$E_ID == E_ID,
                                       "XxE_causal_Estimate"]
    
    b_X_E_corrected <- df_causal_interaction[df_causal_interaction$X_ID == X_ID &
                                             df_causal_interaction$Y_ID == Y_ID &
                                             df_causal_interaction$E_ID == E_ID,
                                             "corrected_XxE_estimate"]

    pred_data$pred_Y_raw <- with(pred_data, b_intercept + b_X * X + b_X2 * X2 +
                                         b_E * environment + b_E2 * environment2 + 
                                         b_X_E_raw * ExX)
    
    pred_data$pred_Y_corr <- with(pred_data, b_intercept + b_X * X + b_X2 * X2 +
                                         b_E * environment + b_E2 * environment2 + 
                                         b_X_E_corrected * ExX)
    
    # lets aim for the high art of confidence intervals

    # load the model that was fitted
    load(paste0(path_to_results, "013_MR_II_Expo_", X_ID, 
                "_Outc_", Y_ID, "_interact_", E_ID, ".Rdata"))

    model <- all_models[[2]]

    cov_mat <- vcov(model)

    cov_mat <- cov_mat[c("(Intercept)", "GRS_column",
                         "interaction_column",
                         "interaction2_column",
                         "GRS2_column",
                         "GRS_column:interaction_column"),
                       c("(Intercept)", "GRS_column",
                         "interaction_column",
                         "interaction2_column",
                         "GRS2_column",
                         "GRS_column:interaction_column")]
      
    colnames(cov_mat) <- c("intercept", "X", "E", "E2", "X2", "XxE")
    rownames(cov_mat) <- c("intercept", "X", "E", "E2", "X2", "XxE")

    b_vector_raw <- c(b_intercept, b_X, b_E, b_E2, b_X2, b_X_E_raw)
    b_vector_corrected <- c(b_intercept, b_X, b_E, b_E2, b_X2, b_X_E_corrected)

    beta_sim_raw <- mvrnorm(n = 10000, mu = b_vector_raw,
                            Sigma = cov_mat)
    beta_sim_corrected <- mvrnorm(n = 10000, mu = b_vector_corrected,
                                  Sigma = cov_mat)

    df_raw <- data.frame(r1 = rep(NA, 3000))
    df_corrected <- data.frame(r1 = rep(NA, 3000))

    # takes long, but CIs may not look smooth with significantly fewer reps
    for(j in 1:10000){ 
        y_raw <- with(pred_data, beta_sim_raw[j, "intercept"] + 
                      beta_sim_raw[j, "X"] * X + 
                      beta_sim_raw[j, "X2"] * X2 +
                      beta_sim_raw[j, "E"] * environment + 
                      beta_sim_raw[j, "E2"] * environment2 + 
                      beta_sim_raw[j, "XxE"] * ExX)

        df_raw <- cbind(df_raw, y_raw)
        
        y_corrected <-  with(pred_data, beta_sim_corrected[j, "intercept"] + 
                             beta_sim_corrected[j, "X"] * X + 
                             beta_sim_corrected[j, "X2"] * X2 +
                             beta_sim_corrected[j, "E"] * environment + 
                             beta_sim_corrected[j, "E2"] * environment2 + 
                             beta_sim_corrected[j, "XxE"] * ExX)
        
        df_corrected <- cbind(df_corrected, y_corrected)
    }

    df_raw <- df_raw[,-1]
    df_corrected <- df_corrected[,-1]


    for(j in 1:nrow(pred_data)){

      predictions_raw <- as.numeric(df_raw[j, ])

      pred_data[j, "ymin_raw"] <- sort(predictions_raw)[250] #LCI (raw)
      pred_data[j, "ymax_raw"] <- sort(predictions_raw)[9750] # UCI (raw)

      predictions_corrected <- as.numeric(df_corrected[j, ])

      pred_data[j, "ymin_cor"] <- sort(predictions_corrected)[250] #LCI (corrected)
      pred_data[j, "ymax_cor"] <- sort(predictions_corrected)[9750] #UCI (corrected)
    }

    pred_data$environment_factor <- factor(pred_data$environment, labels = c("p10", "p50", "p90"))
    
    plot_interaction_raw <- ggplot(data = pred_data[pred_data$environment_factor != "p50", ]) +
        geom_ribbon(alpha = 0.5,
                    aes(x = X, y = pred_Y_raw,
                        ymin = ymin_raw, 
                        ymax = ymax_raw,
                        group = environment_factor, 
                        fill = environment_factor)) +
        geom_line(aes(x = X, y = pred_Y_raw, 
                  color = environment_factor), lwd = 1.2) +
        labs(color = "% environment", x = X_trait, y = paste0("Predicted ", Y_trait), fill = NULL) +
        scale_color_manual(values = c("p10" = "#a6611a", "p90" = "#018571"),
                          name = paste0("% ", E_trait)) +
        scale_fill_manual(values = c("p10" = "#dfc27d", "p90" = "#80cdc1"),
                          name = paste0("% ", E_trait)) +
        theme_minimal(base_size = 14) +
        theme(panel.grid.minor = element_blank(),
              axis.text = element_text(face = "bold"),
              axis.ticks = element_blank())

    ggsave(filename = paste0(path_out, "008_X_", X_ID, "_E_", E_ID, "_Y_", Y_ID, "_raw.png"),
           plot_interaction_raw, width = 6, height = 5)

    df_pred_long <- data.frame(Y = c(pred_data$pred_Y_raw, 
                                     pred_data$pred_Y_corr),
                               environment = rep(pred_data$environment_factor, 
                                                 times = 2),
                                ymin = c(pred_data$ymin_raw,
                                          pred_data$ymin_corrected),
                                ymax = c(pred_data$ymax_raw,
                                          pred_data$ymax_corrected),
                                X = rep(pred_data$X, times = 2),
                                raw_cor = rep(c("raw", "cor"),
                                              each = length(pred_data$pred_Y_raw)))
    
    df_pred_long$group <- paste0(df_pred_long$environment,
                                    df_pred_long$raw_cor)
    df_pred_long$group <- as.factor(df_pred_long$group)
    df_pred_long <- df_pred_long[df_pred_long$environment != "p50", ]
    
    plot_interaction_raw_cor <- 
        ggplot(data = df_pred_long) +
        geom_ribbon(alpha = 0.5,
                    aes(x = X, y = Y,
                        ymin = ymin, 
                        ymax = ymax,
                        group = group, 
                        fill = group)) +
        geom_line(aes(x = X, y = Y, 
                  color = group), lwd = 1.2) +
        labs(color = "% environment\nraw/corrected", x = X_trait, 
             y = paste0("Predicted", Y_trait), fill = NULL) +
        scale_color_manual(values = c("p10raw" = "#8c510a", 
                                      "p10cor" = "#bf812d",
                                      "p90raw" = "#01665e",
                                      "p90cor" = "#35978f"),
                          name = paste0("% ", E_trait),
                          labels = c("p10raw" = "p10 raw", 
                                     "p10cor" = "p10 corrected",
                                     "p90raw" = "p90 raw",
                                     "p90cor" = "p90 corrected")) +
        scale_fill_manual(values = c("p10raw" = "#dfc27d", 
                                      "p10cor" = "#f6e8c3",
                                      "p90raw" = "#c7eae5",
                                      "p90cor" = "#c7eae5"),
                          name = paste0("% ", E_trait),
                          labels = c("p10raw" = "p10 raw", 
                                     "p10cor" = "p10 corrected",
                                     "p90raw" = "p90 raw",
                                     "p90cor" = "p90 corrected")) +
        theme_minimal(base_size = 14) +
        theme(panel.grid.minor = element_blank(),
              axis.text = element_text(face = "bold"),
              axis.ticks = element_blank())

    ggsave(filename = paste0(path_out, "008_X_", X_ID, 
                             "_E_", E_ID, "_Y_", Y_ID, "_raw_cor.png"),
           plot_interaction_raw_cor,
           width = 6, height = 5)

    print(paste0("done with round ", i))

}   
