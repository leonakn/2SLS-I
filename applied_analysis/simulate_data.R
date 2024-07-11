# simulate data for example workflow

set.seed(1505)
setwd("/Users/lknusel/Documents/05_strati_polyMR/08_public_github/2SLS-I/applied_analysis/data/")
# assume X_ID = 37, Y_ID = 11, and E_ID = 1204
n_sample <- 100000
b_G_X <- 0.2
b_GxE_X <- 0.01
b_E_X <- 0.2
b_X_Y <- 0.2
b_X2_Y <- 0.05 #preetty strong!
b_XxE_Y <- 0.02

E <- rnorm(n = n_sample) # the environmental variable

G <- rnorm(n = n_sample) #the GRS

# define X (=the exposure)
#random noise in X, variance fixed so that var(X) will be ~ 1
eps_X <- rnorm(n = n_sample) * c(sqrt(1-var(b_G_X * G +
                                              b_GxE_X * G * E +
                                              b_E_X * E)))
X <- b_G_X * G +
  b_GxE_X * G * E +
  b_E_X * E +
  eps_X

pheno_df <- data.frame(eid = c(1:n_sample),
                       '37-0.0' = X,
                       '37_z' = X,
                       '37_IRNT' = X,
                       "37_IRNT_c" = X,
                       check.names = F)
save(x = pheno_df, file = "004_Exposure_phenotype_37.Rdata")

X2 <- X^2

# define Y (=the outcome)
eps_Y <- rnorm(n = n_sample) * c(sqrt(1-var(b_X_Y * X +
                                              b_X2_Y * X2 +
                                              b_XxE_Y * X * E)))
Y <- b_X_Y * X +
  b_X2_Y * X2 +
  b_XxE_Y * E * X +
  eps_Y

sex <- sample(x = c(0, 1),
              size = n_sample,
              prob = c(0.5, 0.5),
              replace = T)

age <- rnorm(n_sample)

outcome_df <- data.frame(eid = c(1:n_sample),
                         '11-0.0' = Y,
                         sex = sex,
                         age_exact = age,
                         age_exact2 = age^2,
                         age_exact_x_sex = age*sex,
                         '11_z' = Y,
                         '11_z_c' = Y,
                         '11_IRNT' = Y,
                         check.names = F)

save(x = outcome_df, file = "005_outcome_phenotype_11.Rdata")

interaction_df <- data.frame(eid = c(1:n_sample),
                             '1204-0.0' = E,
                             '1204_z' = E,
                             '1204_z_c' = E,
                             '1204_c' = E,
                             check.names = F)

save(x = interaction_df, file = "006_interaction_phenotype_1204.Rdata")

df_gen_score <- data.frame(eid = c(1:n_sample),
                           GRS_37 = G)

save(x = df_gen_score, file = "010_individual_GRS_37.Rdata")

# get GRS on Exposure model (for causal effects)

Y <- "37_IRNT_c"
GRS_on_pheno_model <- lm(pheno_df[,Y] ~ df_gen_score[,"GRS_37"]) |> summary()
save(GRS_on_pheno_model, file = "010_Model_GRS_on_pheno_37.Rdata")


# brilliantly, the covariates df is called "X"
X <- data.frame(eid = c(1:n_sample),
                sex = sex,
                age = age,
                age_exact = age,
                age_exact2 = age^2,
                age_exact_x_sex = age*sex)

PCs <- matrix(data = rnorm(40*n_sample), nrow = n_sample, ncol = 40)
colnames(PCs) <- paste0("PC", c(1:40))
PCs <- as.data.frame(PCs)
X <- cbind(X, PCs)

save(x = X, file = "003_Covariates.Rdata")

medication_df <- data.frame(eid = c(1:n_sample),
                            chol_inhibit = sample(x = c(0, 1),
                                                  size = n_sample,
                                                  prob = c(0.8, 0.2),
                                                  replace = T),
                            bp_regul = sample(x = c(0, 1),
                                              size = n_sample,
                                              prob = c(0.8, 0.2),
                                              replace = T),
                            Insulin = sample(x = c(0, 1),
                                             size = n_sample,
                                             prob = c(0.95, 0.05),
                                             replace = T),
                            hormone_repl = sample(x = c(0, 1),
                                                  size = n_sample,
                                                  prob = c(0.97, 0.03),
                                                  replace = T),
                            oral_contraceptive = sample(x = c(0, 1),
                                                        size = n_sample,
                                                        prob = c(0.9, 0.1),
                                                        replace = T))

save(x = medication_df, file = "003_Medication.Rdata")

df_PGS <- data.frame(IID = c(1:n_sample),
                     PGS_37_SUM = G)
library(readr)
library(R.utils)
write_delim(x = df_PGS, file = "../data/aggregated_scores.txt", delim = "\t")
gzip("../data/aggregated_scores.txt")
