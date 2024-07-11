# lknusel

# perform power analysis for 2SLS-I

# note: now a5 = GRS*E


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
#library(ggiraphExtra)
#require(ggiraph)

Date <- today() %>% str_remove_all(pattern = "^20|-")


path_sim_sett <- paste0("where_you_want_the_simulation_settings/")

pathRFile <- paste0("where_is_this_file_saved/file.R")

path_out <- paste0("/where_you_want_the_results/", 
                   Date, "/")

print(path_out)

if(!dir.exists(path_out)){
  dir.create(path_out, recursive = T)
}

# Create Simulation Settings
sim_sett <- data.frame(nSample = rep(c(10^4, 5 * 10^4, 10^5), each = 18), # vary sample size
                       a0 = rep(0.2, 54), #keep effect of E on X robust at 0.2 (could also be varied)
                       a1 = rep(0, 54),
                       a2 = rep(0.2, 54), 
                       a3 = rep(c(0, 0.05, 0.1), each = 6, times = 3), # vary strength of the true interaction
                       a4 = rep(c(sqrt(0.05), sqrt(0.1), sqrt(0.2)), each = 2, times = 9), # vary amount of variance GRS explains in X
                       a5 = rep(0, 54),
                       b0 = rep(0, 54),
                       b1 = rep(0, 54),
                       b2 = rep(c(0, 0.1), times = 27), # vary strength of quadratic effect (i.e. bias)
                       qX = rep(0.3, 54),
                       qY = rep(0.5, 54),
                       qY2 = rep(0.1, 54),
                       Szenario = c(1:54)
                       )


write.csv(sim_sett, paste0(path_sim_sett, "02_", Date, "_settings_poweranalysis.csv"))

n_snp <- 324
n_repet <- 500

# load in functions for simulation
source("scripts/Functions.R")

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
ModelFormulae <- data.frame("Model_no" = c(1), "Formula" = NA)

ModelFormulae[1, "Formula"] <- as.character("Y ~ GRS + E + GRS2 + E2 + E:GRS")


# Run The Magic Function
for (i in 1:nrow(sim_sett)){
  ModellingXY(df = df[[i]], n_snp = n_snp, n_repet = n_repet, 
              ModelFormulae = ModelFormulae, path_out = path_out)
}

df_power_analysis <- data.frame(Level = rep(NA, times = nrow(sim_sett)*2), 
                                TrueInteract = NA, 
                                X2 = NA, nSample = NA, 
                                Heritability = NA, 
                                PercentSigP = NA,
                                Parameter = "ExGRS",
                                Correction = rep(c("Raw", "'2SLS-I'"), times = nrow(sim_sett)))


files <- list.files(path_in, pattern = "Model_1_Data.csv", full.names = T)

index <- seq(from = 1, to = nrow(sim_sett)*2, by = 2)

for (i in 1:nrow(sim_sett)){
  
  data <- read.csv(files[i], header = T)
  
  df_power_analysis$Level[c(index[i], index[i]+1)] <- unique(data$LEVEL)
  df_power_analysis$TrueInteract[c(index[i], index[i]+1)] <- unique(data$a3)
  df_power_analysis$X2[c(index[i], index[i]+1)] <- unique(data$b2)
  df_power_analysis$nSample[c(index[i], index[i]+1)] <- unique(data$n_sample)
  df_power_analysis$Heritability[c(index[i], index[i]+1)] <- unique(data$a4)^2
  df_power_analysis$PercentSigP[index[i]] <-  sum(data$ExGRS_pval < 0.05)/5
  df_power_analysis$PercentSigP[index[i]+1] <- sum(data$CorrectedPval_ExGRS < 0.05)/5
}

df_power_analysis$error_rate <- ifelse(df_power_analysis$TrueInteract == 0,
                                       yes = df_power_analysis$PercentSigP/100,
                                       no = -(100 - df_power_analysis$PercentSigP)/100)


df_power_analysis$fpr_power <- ifelse(df_power_analysis$TrueInteract == 0,
                                      yes = df_power_analysis$PercentSigP/100,
                                      no = -(df_power_analysis$PercentSigP/100))

##### FALSE POSITIVE RATE
# False positive rate of raw interaction in absence of bias:
summary(df_power_analysis[df_power_analysis$TrueInteract == 0 & 
                            df_power_analysis$X2 == 0 &
                            df_power_analysis$Correction == "Raw", "PercentSigP"])
sd(df_power_analysis[df_power_analysis$TrueInteract == 0 & 
                       df_power_analysis$X2 == 0 &
                       df_power_analysis$Correction == "Raw", "PercentSigP"])

# False positive rate of 2SLS-I in absence of bias:
summary(df_power_analysis[df_power_analysis$TrueInteract == 0 & 
                            df_power_analysis$X2 == 0 &
                            df_power_analysis$Correction != "Raw", "PercentSigP"])
sd(df_power_analysis[df_power_analysis$TrueInteract == 0 & 
                       df_power_analysis$X2 == 0 &
                       df_power_analysis$Correction != "Raw", "PercentSigP"])

# False positive rate of 2SLS-I in presence of bias:
summary(df_power_analysis[df_power_analysis$TrueInteract == 0 & 
                            df_power_analysis$X2 == 0.1 &
                            df_power_analysis$Correction != "Raw", "PercentSigP"])
sd(df_power_analysis[df_power_analysis$TrueInteract == 0 & 
                       df_power_analysis$X2 == 0.1 &
                       df_power_analysis$Correction != "Raw", "PercentSigP"])
# False positive rate of raw interaction in presence of bias
summary(df_power_analysis[df_power_analysis$TrueInteract == 0 & 
                            df_power_analysis$X2 == 0.1 &
                            df_power_analysis$Correction == "Raw", "PercentSigP"])
sd(df_power_analysis[df_power_analysis$TrueInteract == 0 & 
                       df_power_analysis$X2 == 0.1 &
                       df_power_analysis$Correction == "Raw", "PercentSigP"])

##### TRUE POSITIVE RATE IF INTERACTION = 0.05
# TPR no bias raw
summary(df_power_analysis[df_power_analysis$TrueInteract == 0.05 &
                            df_power_analysis$X2 == 0 &
                            df_power_analysis$Correction == "Raw", "PercentSigP"])
sd(df_power_analysis[df_power_analysis$TrueInteract == 0.05 &
                       df_power_analysis$X2 == 0 &
                       df_power_analysis$Correction == "Raw", "PercentSigP"])
# TPR no bias 2SLS-I
summary(df_power_analysis[df_power_analysis$TrueInteract == 0.05 &
                            df_power_analysis$X2 == 0 &
                            df_power_analysis$Correction != "Raw", "PercentSigP"])
sd(df_power_analysis[df_power_analysis$TrueInteract == 0.05 &
                       df_power_analysis$X2 == 0 &
                       df_power_analysis$Correction != "Raw", "PercentSigP"])

##### TRUE POSITIVE RATE IF INTERACTION = 0.1
# TPR no bias raw
summary(df_power_analysis[df_power_analysis$TrueInteract == 0.1 &
                            df_power_analysis$X2 == 0 &
                            df_power_analysis$Correction == "Raw", "PercentSigP"])
sd(df_power_analysis[df_power_analysis$TrueInteract == 0.1 &
                       df_power_analysis$X2 == 0 &
                       df_power_analysis$Correction == "Raw", "PercentSigP"])
# TPR no bias 2SLS-I
summary(df_power_analysis[df_power_analysis$TrueInteract == 0.1 &
                            df_power_analysis$X2 == 0 &
                            df_power_analysis$Correction != "Raw", "PercentSigP"])
sd(df_power_analysis[df_power_analysis$TrueInteract == 0.1 &
                       df_power_analysis$X2 == 0 &
                       df_power_analysis$Correction != "Raw", "PercentSigP"])


df_power_analysis$Correction <- factor(df_power_analysis$Correction, 
                                       levels = c("Raw", "'2SLS-I'"))
df_power_analysis$X2 <- factor(df_power_analysis$X2, 
                               labels = paste0("beta [X^2 %->% Y]", ": ", 
                                               levels(factor(df_power_analysis$X2))))
df_power_analysis$Heritability <- factor(df_power_analysis$Heritability)
df_power_analysis$nSample <- factor(df_power_analysis$nSample)
df_power_analysis$TrueInteract <- factor(df_power_analysis$TrueInteract,
                                         labels = paste0('beta ["E*X" %->% Y]',
                                                         ':', 
                                                         levels(factor(df_power_analysis$TrueInteract))))


FilenamePoweranalysis <- paste0(path_out, "Poweranalysis.csv")
write.csv(df_power_analysis, file = FilenamePoweranalysis)

# colors obtained from:
# https://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3
# 231203: plot power instead of false negative rate for readability

plot_power_analysisExGRS <- ggplot(data = df_power_analysis, 
                                      mapping = aes(x = Heritability, 
                                                    y = nSample,
                                                    fill = error_rate,
                                                    label = abs(round(fpr_power, 2))))+
  geom_tile()+
  #scale_fill_viridis_c(begin = 0.2, end = 0.8)+
  # scale_fill_gradient2(low = "#67a9cf", high = "#ef8a62", 
  #                      mid = "#f7f7f7", midpoint = 0)+
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365", 
                       mid = "#f5f5f5", midpoint = 0)+
  geom_text(color = "#8F361F")+
  theme_minimal(base_size = 12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(TrueInteract ~ X2 + Correction, 
             labeller = label_parsed)+
  ylab(label = "Sample Size") +
  labs(fill = "error rate") +
  theme(strip.text.y = element_text(angle = 0)) 


ggsave(filename = paste0(path_out, "PlotPowerExGRS.png"), plot = plot_power_analysisExGRS_v3, 
       width = 10, height = 7)

file.copy(from = pathRFile,
          to = path_out, overwrite = T)


