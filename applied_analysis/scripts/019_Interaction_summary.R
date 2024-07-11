# lknusel

################################################################################
# QUADRATIC EFFECTS: when to consider the corrected interaction?
################################################################################
file_quadratic_GRS <- snakemake@input[["file_quadratic_GRS"]]

load(file_quadratic_GRS)

df_causal_quadratic <- df_causal_quadratic[df_causal_quadratic$X_ID != 
             df_causal_quadratic$Y_ID, ] 
quadratic_index <-  which(df_causal_quadratic$X2_Y_causal_pval < 0.05) 


check_corrected_interaction <- 
  df_causal_quadratic[quadratic_index, c("X_trait", "Y_trait")]

if(length(quadratic_index) > 0){
  check_corrected_interaction$unicol <- 
     paste0(check_corrected_interaction$X_trait,
            "_",
            check_corrected_interaction$Y_trait)
}

################################################################################
# INTRACITONS: GET RELEVANT ESTIMATES
################################################################################

GRS_interaction_file <- snakemake@input[["file_interaction_GRS"]]

load(GRS_interaction_file)

df_GRS <- df_causal_interaction

df_GRS <- df_GRS[,c("X_trait", "Y_trait", "E_trait", 
                    "X_ID", "Y_ID", "E_ID",
                    "XxE_causal_Estimate", "XxE_causal_SE", 
                    "XxE_causal_tval", "XxE_causal_pval",
                    "corrected_XxE_estimate", "corrected_XxE_SE", 
                    "corrected_XxE_tval", "corrected_XxE_pval",
                    "degrees_of_freedom")]

df_GRS$unicol <- paste0(df_GRS$X_trait, "_", df_GRS$Y_trait)

# choose either raw or corrected interaction, depending on quadratic effects:
df_GRS$estimate <- NA
df_GRS$se <- NA
df_GRS$pval <- NA
df_GRS$raw_cor <- NA

for(i in 1:nrow(df_GRS)){
  df_GRS$estimate[i] <- ifelse(df_GRS$unicol[i] %in% check_corrected_interaction$unicol,
                               yes = df_GRS$corrected_XxE_estimate[i],
                               no = df_GRS$XxE_causal_Estimate[i])
  df_GRS$se[i] <- ifelse(df_GRS$unicol[i] %in% check_corrected_interaction$unicol,
                         yes = df_GRS$corrected_XxE_SE[i],
                         no = df_GRS$XxE_causal_SE[i])
  df_GRS$pval[i] <- ifelse(df_GRS$unicol[i] %in% check_corrected_interaction$unicol,
                           yes = df_GRS$corrected_XxE_pval[i],
                           no = df_GRS$XxE_causal_pval[i])
  df_GRS$tval[i] <- ifelse(df_GRS$unicol[i] %in% check_corrected_interaction$unicol,
                           yes = df_GRS$corrected_XxE_tval[i],
                           no = df_GRS$XxE_causal_pval[i])
  df_GRS$raw_cor[i] <- ifelse(df_GRS$unicol[i] %in% check_corrected_interaction$unicol,
                              yes = "cor",
                              no = "raw")
}

df_GRS <- df_GRS[,c("X_trait", "Y_trait", "E_trait",
                    "X_ID", "Y_ID", "E_ID",
                    "degrees_of_freedom", "unicol", 
                    "estimate", "se", "tval", "pval", "raw_cor")] # only keep relevant observations

df_GRS$PRS_GRS <- "GRS"

PRS_interaction_file <- snakemake@input[["file_interaction_PRS"]]

load(PRS_interaction_file)

df_PRS <- df_causal_interaction

df_PRS <- df_PRS[,c("X_trait", "Y_trait", "E_trait", 
                    "X_ID", "Y_ID", "E_ID",
                    "XxE_causal_Estimate", "XxE_causal_SE", 
                    "XxE_causal_tval", "XxE_causal_pval",
                    "corrected_XxE_estimate", "corrected_XxE_SE", 
                    "corrected_XxE_tval", "corrected_XxE_pval",
                    "degrees_of_freedom")]

# create unique identifier column
# for(i in 1:nrow(check_corrected_interaction)){
#   check_corrected_interaction$unicol[i] <- 
#     paste0(check_corrected_interaction$X_trait[i], "_", 
#            check_corrected_interaction$Y_trait[i])
# }

for(i in 1:nrow(df_PRS)){
  df_PRS$unicol[i] <- paste0(df_PRS$X_trait[i], "_", df_PRS$Y_trait[i])
}

# choose either raw or corrected interaction, depending on quadratic effects:
df_PRS$estimate <- NA
df_PRS$se <- NA
df_PRS$pval <- NA
df_PRS$raw_cor <- NA

for(i in 1:nrow(df_PRS)){
  df_PRS$estimate[i] <- ifelse(df_PRS$unicol[i] %in% check_corrected_interaction$unicol,
                               yes = df_PRS$corrected_XxE_estimate[i],
                               no = df_PRS$XxE_causal_Estimate[i])
  df_PRS$se[i] <- ifelse(df_PRS$unicol[i] %in% check_corrected_interaction$unicol,
                         yes = df_PRS$corrected_XxE_SE[i],
                         no = df_PRS$XxE_causal_SE[i])
  df_PRS$pval[i] <- ifelse(df_PRS$unicol[i] %in% check_corrected_interaction$unicol,
                           yes = df_PRS$corrected_XxE_pval[i],
                           no = df_PRS$XxE_causal_pval[i])
  df_PRS$tval[i] <- ifelse(df_PRS$unicol[i] %in% check_corrected_interaction$unicol,
                           yes = df_PRS$corrected_XxE_tval[i],
                           no = df_PRS$XxE_causal_pval[i])
  df_PRS$raw_cor[i] <- ifelse(df_PRS$unicol[i] %in% check_corrected_interaction$unicol,
                              yes = "cor",
                              no = "raw")
}

df_PRS <- df_PRS[,c("X_trait", "Y_trait", "E_trait", 
                    "X_ID", "Y_ID", "E_ID",
                    "degrees_of_freedom", "unicol", 
                    "estimate", "se", "tval", "pval", "raw_cor")] # only keep relevant observations


df_PRS$PRS_GRS <- "PRS"

adjust_names <- grep(pattern = "_trait|_ID", names(df_PRS), invert = T)
names(df_PRS)[adjust_names] <- paste0(names(df_PRS)[adjust_names], "_PRS")
adjust_names <- grep(pattern = "_trait|_ID", names(df_GRS), invert = T)
names(df_GRS)[adjust_names] <- paste0(names(df_GRS)[adjust_names], "_GRS")


df_interaction <- merge(df_PRS, df_GRS, by = c("X_trait", "Y_trait", "E_trait", 
                                               "X_ID", "Y_ID", "E_ID"), 
                        all = T)


# Prepare tier computation:
df_interaction$both_sig <- NA
df_interaction$only_GRS_sig <- NA
df_interaction$only_PRS_sig <- NA
df_interaction$effect_agrees <- NA



for(i in 1:nrow(df_interaction)){
  
  X_trait <- df_interaction$X_trait[i]
  Y_trait <- df_interaction$Y_trait[i]
  
  df_interaction$both_sig[i] <- ifelse(df_interaction$pval_GRS[i] < 10^-3 & df_interaction$pval_PRS[i] < 10^-3,
                                       yes = 1,
                                       no = 0)
  df_interaction$only_GRS_sig[i] <- ifelse(df_interaction$pval_GRS[i] < 10^-3 & df_interaction$pval_PRS[i] > 10^-3,
                                           yes = 1,
                                           no = 0)
  df_interaction$only_PRS_sig[i] <- ifelse(df_interaction$pval_GRS[i] > 10^-3 & df_interaction$pval_PRS[i] < 10^-3,
                                           yes = 1,
                                           no = 0)
  std_err <- sqrt(df_interaction$se_GRS[i]^2 + df_interaction$se_PRS[i]^2)
  tval <- (df_interaction$estimate_GRS[i] - df_interaction$estimate_PRS[i])/std_err
  degfree <- mean(df_interaction$degrees_of_freedom_GRS[i], df_interaction$degrees_of_freedom_PRS[i])
  pval <- 2 * pt(-abs(tval), degfree)
  
  df_interaction$effect_agrees[i] <- ifelse(pval > 0.05  & 
                                                 (sign(df_interaction[i, "estimate_GRS"]) == sign(df_interaction[i, "estimate_PRS"])),
                                            yes = 1, 
                                            no = 0)

  df_interaction$effect_direction_agrees[i] <- ifelse(sign(df_interaction[i, "estimate_GRS"]) == sign(df_interaction[i, "estimate_PRS"]),
                                                      yes = 1, 
                                                      no = 0)
  
}

df_interaction$any_sig <- df_interaction$only_GRS_sig + df_interaction$only_PRS_sig

df_interaction$tier <- 3 * (df_interaction$both_sig * df_interaction$effect_direction_agrees) + 
  (2 * ((df_interaction$any_sig) * df_interaction$effect_agrees)) +
  df_interaction$only_GRS_sig

save(df_interaction, file = snakemake@output[["df_interaction"]])