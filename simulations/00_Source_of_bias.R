# lknusel
# simulation script that accesses Functions.R to simulate and process data as indicated

library(lubridate)
library(stringr)
library(ggplot2)
library(tidyr)
library(grid)
library(magrittr)
library(gridExtra)
library(tictoc)
library(readr)
require(plyr)
library(dplyr)
library(reshape2)
library(gridGraphics)
library(MetBrewer)


Date <- today() %>% str_remove_all(pattern = "^20|-")

path_sim_sett <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/lknusel/",
                        "results/age_stratified_polymr/BMIAnalysis/Simulations/sim_sett/")

pathRFile <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/lknusel/script/",
                    "age_stratified_polymr/00_simulation_paper_figures/00_Source_of_bias.R")

path_out <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/lknusel/results/",
                   "age_stratified_polymr/BMIAnalysis/Simulations/",
                   "00_Source_of_bias/", 
                   Date, "/")


if(!dir.exists(path_out)){
  dir.create(path_out, recursive = T)
}

# Create Simulation Settings
sim_sett <- data.frame(nSample = rep(10^4, times = 27),
                       a0 = rep(c(0, 0.1, 0.3), each = 9), # effects of E (envir) on X (expo)
                       a1 = rep(0, 27), # effects of E on Y (outcome)
                       a2 = rep(0.2, 27), # effect of X on Y
                       a3 = rep(c(0, 0.1, 0.3), each = 3, times = 3), # strength of E*X on Y
                       a4 = rep(sqrt(0.1), 27), # amount of variance GRS explains in X
                       a5 = rep(0, 27), # strength of E*G on X
                       b0 = rep(0, 27), # effects of E^2 on X
                       b1 = rep(0, 27), # effects of E^2 on Y
                       b2 = rep(c(0, 0.05, 0.15), times = 9), # effects of X^2 on Y
                       qX = rep(0.3, 27), # strength of confounding on X
                       qY = rep(0.5, 27), # strength of confounding on Y
                       qY2 = rep(0.1, 27), # strength of non-linear confounding on Y
                       Szenario = c(1:27)
                       )


write.csv(sim_sett, paste0(path_sim_sett, "00_", Date, "_settings_source_of_bias.csv"))


n_snp <- 324
n_repet <- 500

#### Simulate data, fit models
source("/data/FAC/FBM/DBC/zkutalik/default_sensitive/lknusel/script/age_stratified_polymr/00_simulation_paper_figures/Functions.R")


# Create data frame for simulation settings

df <- list()
# same parameteres for each row, except from "Seed" (random number)

for (i in 1:nrow(sim_sett)) {
  df[[i]] <- data.frame(
    a0 = rep(sim_sett[i, "a0"], n_repet),
    a1 = sim_sett[i, "a1"],
    a2 = sim_sett[i, "a2"],
    a3 = sim_sett[i, "a3"],
    a4 = sim_sett[i, "a4"],
    a5 = sim_sett[i, "a5"],
    b0 = sim_sett[i, "b0"],
    b1 = sim_sett[i, "b1"],
    b2 = sim_sett[i, "b2"],
    qX = sim_sett[i, "qX"],
    qY = sim_sett[i, "qY"],
    qY2 = sim_sett[i, "qY2"],
    Seed = seq(from = 1, to = n_repet, by = 1),
    LEVEL = sim_sett[i, "Szenario"],
    n_sample = sim_sett[i, "nSample"]
  )
}


# Check whether parameter setting allows for sufficient variance
# in the residuals and adjust settings if it does not.
# commented out because already run

for (i in 1:length(df)) {
  df[[i]] <- obtain_maximum_parameters(df[[i]], nObs = 5, 
                                       path_out = path_out, n_repet = n_repet)
  print(paste("Analysed szenario ", unique(df[[i]]$LEVEL)))
}

save(df, file = paste0(path_out, "DFSim.Rdata"))

load(file = paste0(path_out, "DFSim.Rdata"))

# Create Model Formulae
#ModelFormulae <- data.frame("Model_no" = c(1:4), "Formula" = NA)
ModelFormulae <- data.frame("Model_no" = c(1, 2), "Formula" = NA)
 
ModelFormulae[1, "Formula"] <- as.character("Y ~ GRS + E + GRS2 + E:GRS")
ModelFormulae[2, "Formula"] <- as.character("Y ~ GRS + E + GRS2 + E2 + E:GRS")


# Run The Magic Function

for (i in 1:nrow(sim_sett)){
  ModellingXY(df = df[[i]], n_snp = n_snp, n_repet = n_repet, 
              ModelFormulae = ModelFormulae, path_out = path_out)
}

files <- list.files(path_out, pattern = "Model_2_Data.csv", full.names = T)
all(file.exists(files))


data_m2 <- lapply(files, read.csv)
df_m2 <- do.call(rbind, data_m2)
df_m2$deviation_true_a3 <- df_m2$ExGRS_Estimate - df_m2$a3
median_deviation_a3 <- aggregate(deviation_true_a3 ~ LEVEL, df_m2, median)
df_m2 <- merge(df_m2, median_deviation_a3, by = "LEVEL", all = T)
df_m2$deviation_corrected <- df_m2$CorrectedEstimate_ExGRS - df_m2$a3
median_deviation_corrected <- aggregate(deviation_corrected ~ LEVEL, df_m2, median)
df_m2 <- merge(df_m2, median_deviation_corrected, by = "LEVEL", all = T)

df <- df_m2[ ,c("a0", "b2", "a3", "ExGRS_Estimate", "CorrectedEstimate_ExGRS", "deviation_true_a3.y", 
                "deviation_corrected.y", "LEVEL")]

df1 <- df_m2[,c("a0", "b2", "a3", "ExGRS_Estimate", "deviation_true_a3.y", "LEVEL")]
df1$raw_vs_corrected <- "raw"
to_rename_Estim <- which(colnames(df1) == "ExGRS_Estimate")
to_rename_deviation <- which(colnames(df1) == "deviation_true_a3.y")
names(df1)[c(to_rename_Estim, to_rename_deviation)] <- c("Estimate", "Deviation")

df2 <- df_m2[,c("a0", "b2", "a3", "CorrectedEstimate_ExGRS", "deviation_corrected.y", "LEVEL")]
df2$raw_vs_corrected <- "corrected"
to_rename_Estim <- which(colnames(df2) == "CorrectedEstimate_ExGRS")
to_rename_deviation <- which(colnames(df2) == "deviation_corrected.y")
names(df2)[c(to_rename_Estim, to_rename_deviation)] <- c("Estimate", "Deviation")

df <- rbind(df1, df2)

df$a0 <- factor(df$a0,
                levels = c(0.3, 0.1, 0),
                labels = paste0("beta[E %->% X]", ": ", c(0.3, 0.1, 0)))
df$b2 <- factor(df$b2,
                labels = paste0("beta[X^2 %->% Y]", ": ", levels(factor(df$b2))))
df$a3 <- factor(df$a3)
df$raw_vs_corrected <- factor(df$raw_vs_corrected, 
                              levels = c("raw", "corrected"))
df$group <- paste0(df$raw_vs_corrected, "_", df$a3) 
df$group <- factor(df$group, 
                   levels = c("raw_0", "raw_0.1", "raw_0.3",
                              "corrected_0", "corrected_0.1", 
                              "corrected_0.3"))

bias_plot <- ggplot(data = df, 
                       mapping = aes(x = a3, y = Estimate, group = group)) +
  geom_violin(aes(fill = raw_vs_corrected, group = group),
              color = "lightgrey",
              position = position_dodge(0.9))+
  geom_boxplot(mapping = aes(group = group),
               color = "lightgrey",
               fill = "transparent",
               position = position_dodge(0.9),
               width = 0.3) +
  scale_fill_viridis_d(begin = 0.3, end = 0.8) +
  labs(fill = "measurement") +
  geom_segment(mapping = aes(x = as.numeric(a3) - 0.55,
                             xend = as.numeric(a3) + 0.55,
                             y = as.numeric(as.character(a3)),
                             yend = as.numeric(as.character(a3))), 
               color = "darkgrey", lwd = 0.3)+

  ylab(label = expression(hat(beta)[E~"*"~GRS[X] %->% Y]))+
  xlab(label = expression(beta[E~"*"~GRS[X] %->% Y]))+
  theme_light(base_size = 12) +
  theme(panel.grid.major.x = element_blank())+
  facet_grid(a0 ~ b2,
             labeller = label_parsed) +
  theme(strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black"),
        strip.text.y = element_text(angle = 0)) 

ggsave(filename = paste0(path_out, "00_Source_of_bias_correction_plot.png"),
       plot = bias_plot, width = 10, height = 6)
