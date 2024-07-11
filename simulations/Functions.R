# Selection of the functions needed to perform simulations

# basically just rnorm for E and GRS,
# considering the heritability for the GRS as sd
# make sure to reload existing data instead of creating new data each time
GRSFun <- function(n_sample, n_snp, n_repet, path_out, h_X, LEVEL = "X") {
  
  if(file.exists(paste0(path_out, "GRS_n", n_sample, "nRep", n_repet, "h2_", round(h_X, 2), "Level_", LEVEL, ".Rdata"))){
    
    load(paste0(path_out, "GRS_n", n_sample, "nRep", n_repet, "h2_", round(h_X, 2), "Level_", LEVEL, ".Rdata"))
    
    print("Genetic data exists already")
    
  } else {
    
    GRS <- list()
    
    for (i in 1:n_repet){
      GRS[[i]] <- rnorm(n_sample, sd = h_X)
    }
    
    save(GRS, file = paste0(path_out, "GRS_n", n_sample, "nRep", n_repet, "h2_", round(h_X, 2), "Level_", LEVEL, ".Rdata"))
  }
  
  if(file.exists(paste0(path_out, "E_n", n_sample, "nRep", n_repet, "Level_", LEVEL, ".Rdata"))){
    
    load(paste0(path_out, "E_n", n_sample, "nRep", n_repet, "Level_", LEVEL, ".Rdata"))
    
    print("Environment variable exists already")
    
  } else {

    E <- list()
    
    for(i in 1:n_repet){
      E[[i]] <- rnorm(n_sample)
    }
    
    save(E, file = paste0(path_out, "E_n", n_sample, "nRep", n_repet, "Level_", LEVEL, ".Rdata"))
  }
}

# function that allows to obtain the maximal height of the parameters of intersest
# will provide a df with the parameter levels.
# if those allow for sufficient variance, they will be left as is.
# otherwise, all parameters except the heritability will be adjusted accordingly.

obtain_maximum_parameters <- function(df, nObs, path_out, n_repet, LEVEL = "obtain_max_par") {
  
  n_sample <- unique(df$n_sample)
  h_X <- unique(df$a4)
  GRSFun(n_sample, n_snp, n_repet, path_out, h_X, LEVEL = LEVEL)
  
  load(file = paste0(path_out, "GRS_n", n_sample, "nRep", n_repet, "h2_", round(h_X, 2), "Level_", LEVEL, ".Rdata"))
  load(file = paste0(path_out, "E_n", n_sample, "nRep", n_repet, "Level_", LEVEL, ".Rdata"))
  
  i <- sample(x = c(1:n_repet), size = 1)
  
  GRS <- GRS[[i]]
  E <- E[[i]]
  
  ExceededVarianceDf <- data.frame(
    VarParX = rep(NA, nObs), VarParY = rep(NA, nObs),
    ReductionRateX = rep(NA, nObs), ReductionRateY = rep(NA, nObs),
    sdX1 = rep(NA, nObs), sdY1 = rep(NA, nObs)
  )
  
  for (i in 1:nObs) {
    dftemp <- df[i, ]
    a0 <- dftemp$a0
    a0.1 <- dftemp$a0
    a1 <- dftemp$a1
    a1.1 <- dftemp$a1
    a2 <- dftemp$a2
    a3 <- dftemp$a3
    a5 <- dftemp$a5
    #a4 <- dftemp$a4
    b0 <- dftemp$b0
    b1 <- dftemp$b1
    b2 <- dftemp$b2
    qX <- dftemp$qX
    qY <- dftemp$qY
    qY2 <- dftemp$qY2
    SeedToSet <- dftemp$Seed
    n_snp <- 324
    set.seed(SeedToSet)
    U <- rnorm(n_sample)
    
    # if var(allpars(X)) > 0.68, multiply all parameters by a constant to find
    # max par possible...
    
    #if (var(a0 * E + b0 * E^2 + a4 * GRS + U * qX) > 0.68) {
    if (var(a0 * E + b0 * E^2 + a5 * GRS * E + GRS + U * qX) > 0.68) {
      ExceededVarianceDf[i, "VarParX"] <- 1
      #  while (var(a0 * E + b0 * E^2 + a4 * GRS + U * qX) >= 0.68) {
      while (var(a0 * E + b0 * E^2 + a5 * GRS * E + GRS + U * qX) >= 0.68) {
        a0 <- a0 * 0.99
        b0 <- b0 * 0.99
        qX <- qX * 0.99
        a5 <- a5 * 0.99
      }
      ExceededVarianceDf[i, "ReductionRateX"] <- a0 / a0.1
    } else {
      ExceededVarianceDf[i, "VarParX"] <- 0
      ExceededVarianceDf[i, "ReductionRateX"] <- 1
    }
    #eps0 <- rnorm(n_sample) * c(sqrt(1 - var(a0 * E + b0 * E^2 + a4 * GRS + U * qX)))
    eps0 <- rnorm(n_sample) * c(sqrt(1 - var(a0 * E + b0 * E^2 + a5 * GRS * E + GRS + U * qX)))
    #X <- a0 * E + b0 * E^2 + a4 * GRS + U * qX + eps0
    X <- a0 * E + b0 * E^2 + a5 * GRS * E + GRS + U * qX + eps0
    
    # little quality check
    if (sd(X) > 1.03 | sd(X) < 0.97) {
      ExceededVarianceDf[i, "sdX1"] <- 0
    } else {
      ExceededVarianceDf[i, "sdX1"] <- 1
    }
    a1 <- a1 * ExceededVarianceDf[i, "ReductionRateX"]
    a2 <- a2 * ExceededVarianceDf[i, "ReductionRateX"]
    a3 <- a3 * ExceededVarianceDf[i, "ReductionRateX"]
    #a4 <- a4 * ExceededVarianceDf[i, "ReductionRateX"]
    a5 <- a5 * ExceededVarianceDf[i, "ReductionRateX"]
    b1 <- b1 * ExceededVarianceDf[i, "ReductionRateX"]
    b2 <- b2 * ExceededVarianceDf[i, "ReductionRateX"]
    qY <- qY * ExceededVarianceDf[i, "ReductionRateX"]
    qY2 <- qY2 * ExceededVarianceDf[i, "ReductionRateX"]
    
    if (var(a1 * E + a2 * X + a3 * X * E + b1 * E^2 + b2 * X^2 + U * qY + U^2 * qY2) > 0.68) {
      ExceededVarianceDf[i, "VarParY"] <- 1
      # keep reducing parameters until they reach a level where var(allpars(Y)) <= 0.68
      while (var(a1 * E + a2 * X + a3 * X * E + b1 * E^2 + b2 * X^2 + U * qY + U^2 * qY2) >= 0.68) {
        a1 <- a1 * 0.99
        a2 <- a2 * 0.99
        a3 <- a3 * 0.99
        a7 <- a7 * 0.99
        b1 <- b1 * 0.99
        b2 <- b2 * 0.99
        qY <- qY * 0.99
        qY2 <- qY2 * 0.99
      }
      ExceededVarianceDf[i, "ReductionRateY"] <- a1 / a1.1
    } else {
      ExceededVarianceDf[i, "VarParY"] <- 0
      ExceededVarianceDf[i, "ReductionRateY"] <- 1
    }
    eps1 <- rnorm(n_sample) * c(sqrt((1 - var(a1 * E + a2 * X + a3 * X * E + b1 * E^2 + b2 * X^2 + U * qY + U^2 * qY2))))
    Y <- a1 * E + b1 * E^2 + a2 * X + a3 * X * E + b2 * X^2 + U * qY + U^2 * qY2 + eps1 # in all scenarios where b1 and b2 = 0 only linear effect of E on X
    # check if the resulting sd(Y) is ~ 1
    if (sd(Y) > 1.03 | sd(Y) < 0.97) {
      ExceededVarianceDf[i, "sdY1"] <- 0
    } else {
      ExceededVarianceDf[i, "sdY1"] <- 1
    }
  }
  if (sum(ExceededVarianceDf[, "VarParX"]) > 0) {
    df[, c("a0", "a1", "a2", "a3", "a5", "b0", "b1", "b2", "qX", "qY", "qY2")] <- 
      df[, c("a0", "a1", "a2", "a3", "a5", "b0", "b1", "b2", "qX", "qY", "qY2")] * 
      min(ExceededVarianceDf[, "ReductionRateX"])
    print(paste0("df ", i, " has been adjusted"))
  }
  if (sum(ExceededVarianceDf[, "VarParY"]) > 0) {
    df[, c("a0", "a1", "a2", "a3", "a5", "b0", "b1", "b2", "qX", "qY", "qY2")] <- 
      df[, c("a0", "a1", "a2", "a3", "a5", "b0", "b1", "b2", "qX", "qY", "qY2")] * 
      min(ExceededVarianceDf[, "ReductionRateY"])
    print(paste0("df ", i, " has been adjusted"))
  }
  if (sum(ExceededVarianceDf[, "sdX1"]) != nrow(ExceededVarianceDf)) {
    print("oh-ooh, correction for param X was not sufficent")
  }
  if (sum(ExceededVarianceDf[, "sdY1"]) != nrow(ExceededVarianceDf)) {
    print("oh-ooh, correction for param Y was not sufficent")
  }
  if (sum(ExceededVarianceDf[, "VarParX"]) > 0 | sum(ExceededVarianceDf[, "VarParY"]) > 0) {
    print(paste0(
      "Parameters were adjusted by ", min(ExceededVarianceDf[, "ReductionRateX"]), "; and ", min(ExceededVarianceDf[, "ReductionRateX"]),
      ". Leading to the following new parameters: a0 = ", unique(df$a0),
      "; a0 = ", unique(df$a0),
      "; a1 = ", unique(df$a1),
      "; a2 = ", unique(df$a2),
      "; a3 = ", unique(df$a3),
      "; a4 = ", unique(df$a4),
      "; a5 = ", unique(df$a5),
      "; b0 = ", unique(df$b0),
      "; b1 = ", unique(df$b1),
      "; b2 = ", unique(df$b2),
      "; qX = ", unique(df$qX),
      "; qY = ", unique(df$qY),
      "; qY2 = ", unique(df$qY2)
    ))
  } else {
    print("Parameters remain as planned")
  }
  return(df)
}

# Function that generates all raw data points
# GRS, E, E2, E2, X, X2, Y, and U


generate_raw_data <- function(a0, a1, a2, a3, a4, a5, b0, b1, b2, qX, qY, qY2, SeedToSet, n_sample, n_snp, LEVEL, path_out) {
  
  #If data for GRS and E does not yet exist, create data first
  
  h_X <- a4
  
  if(!file.exists(paste0(path_out, "GRS_n", n_sample, "nRep", n_repet, "h2_", round(h_X, 2), "Level_", LEVEL, ".Rdata"))){
    #    GRSFun(n_sample = n_sample, n_snp = n_snp, n_repet = n_repet, path_out = path_out)}
    GRSFun(n_sample = n_sample, n_snp = n_snp, n_repet = n_repet, path_out = path_out, h_X = a4, LEVEL = LEVEL)}
  
  load(file = paste0(path_out, "GRS_n", n_sample, "nRep", n_repet, "h2_", round(h_X, 2), "Level_", LEVEL, ".Rdata"))
  load(file = paste0(path_out, "E_n", n_sample, "nRep", n_repet, "Level_", LEVEL, ".Rdata"))
  
  set.seed(SeedToSet)
  U <- rnorm(n_sample)
  
  # QC1:
  #  if (var(a0 * E + b0 * E^2 + a4 * c(GRS[[SeedToSet]]) + U * qX) > 0.7) {
  if (var(a0 * c(E[[SeedToSet]]) + b0 * c(E[[SeedToSet]])^2 + c(GRS[[SeedToSet]]) + a5 * c(GRS[[SeedToSet]]) * c(E[[SeedToSet]]) + U * qX) > 0.7) {
    stop(paste0(
      "WARNING: for a0 = ", a0, " b0 = ", b0, ", and qX = ", qX,
      "the variance of all parameter equals ", 
      #      var(a0 * E + b0 * E^2 + a4 * c(GRS[[SeedToSet]]) + U * qX), 
      var(a0 * c(E[[SeedToSet]]) + b0 * c(E[[SeedToSet]])^2 + c(GRS[[SeedToSet]]) + a5 * c(GRS[[SeedToSet]]) * c(E[[SeedToSet]]) + U * qX),
      "! Please adjust Parameters."
    ))
  }
  # define more parameters
  #  eps0 <- rnorm(n_sample) * c(sqrt(1 - var(a0 * E + b0 * E^2 + a4 * c(GRS[[SeedToSet]]) + U * qX))) 
  eps0 <- rnorm(n_sample) * c(sqrt(1 - var(a0 * c(E[[SeedToSet]]) + b0 * c(E[[SeedToSet]])^2 + c(GRS[[SeedToSet]])+ a5 * c(GRS[[SeedToSet]]) * c(E[[SeedToSet]]) + U * qX))) 
  #  X <- a0 * E + b0 * E^2 + a4 * c(GRS[[SeedToSet]]) + U * qX + eps0
  X <- a0 * c(E[[SeedToSet]]) + b0 * c(E[[SeedToSet]])^2 + c(GRS[[SeedToSet]])+ a5 * c(GRS[[SeedToSet]]) * c(E[[SeedToSet]]) + U * qX + eps0
  
  # QC2:
  if (sd(X) > 1.02 | sd(X) < 0.98) {
    print(paste0("The sd of X for LEVEL ", LEVEL, " is ", sd(X)))
  }
  
  # define more parameters
  var_allpars_Y <- var(a1 * c(E[[SeedToSet]]) + a2 * X + a3 * X * c(E[[SeedToSet]]) + b1 * c(E[[SeedToSet]])^2 + b2 * X^2 + U * qY + U^2 * qY2)
  
  # QC3:
  if (var_allpars_Y > 0.7) {
    stop(paste0(
      "WARNING: for a1 = ", a1, ", a2 = ", a2, ", b1 = ", b1, ", b2 = ", b2, ", and qY = ", qY,
      "the variance of all parameter equals ", var_allpars_Y, "! Please adjust Parameters."
    ))
  }
  
  # define more parameters
  eps1 <- rnorm(n_sample) * c(sqrt(1 - var_allpars_Y))
  Y <- a1 * c(E[[SeedToSet]]) + b1 * c(E[[SeedToSet]])^2 + a2 * X + b2 * X^2 + a3 * X * c(E[[SeedToSet]]) + U * qY + U^2 * qY2 + eps1
  
  # QC4:
  if (sd(Y) > 1.02 | sd(Y) < 0.98) {
    print(paste0("The sd of Y for LEVEL ", LEVEL, " is = ", sd(Y)))
  }
  
  # Put this all together to a new data frame
  dftemp <- data.frame(X = X, Y = Y, E = c(E[[SeedToSet]]), 
                       GRS = c(GRS[[SeedToSet]]), U = U, eps0 = eps0, eps1 = eps1)
  
  dftemp$GRS2 <- dftemp$GRS^2
  dftemp$X2 <- dftemp$X^2
  dftemp$E2 <- dftemp$E^2
  dftemp$E3 <- dftemp$E^3
  dftemp$E4 <- dftemp$E^4
  dftemp$resXGRS <- residuals(lm(dftemp$X ~ dftemp$GRS))
  dftemp$resXGRS2 <- dftemp$resXGRS^2
  dftemp$X_corrected_for_E <- residuals(lm(dftemp$X ~ dftemp$E))
  dftemp$Y_corrected_for_E <- residuals(lm(dftemp$Y ~ dftemp$E))
  
  return(dftemp)
}




runModelFun <- function(df, ModelFormula) {
  model_summary <- eval(bquote(lm(formula = .(ModelFormula), data = df))) %>% summary()
  return(model_summary)
}

#### Get Overview DF of the parameters of interest

createOverviewFun <- function(Models, dfModel) {
  
  dfOut <- dfModel
  
  Coeff_Names_readable <- gsub(pattern = "\\(|\\)", rownames(Models[[1]]$coefficients), replacement = "") %>%
    gsub(pattern = ":", ., replacement = "x")
  
  if("GRSxE" %in% Coeff_Names_readable){
    to_rename <- which(Coeff_Names_readable == "GRSxE")
    Coeff_Names_readable[to_rename] <- "ExGRS"
  }
  
  dfOut[, paste0(Coeff_Names_readable, "_pval")] <- NA
  dfOut[, paste0(Coeff_Names_readable, "_Estimate")] <- NA
  dfOut[, paste0(Coeff_Names_readable, "_SE")] <- NA
  dfOut[, paste0(Coeff_Names_readable, "_tval")] <- NA
  dfOut[, "ModelPval"] <- NA
  
  
  for (ind in 1:nrow(dfModel)) {
    Model <- Models[[ind]]
    Coeffs <- rownames(Model$coefficients)
    ModelPval <- pf(Model$fstatistic[1],
                    Model$fstatistic[2],
                    Model$fstatistic[3],
                    lower.tail = F
    )
    
    dfOut[ind, "ModelPval"] <- ModelPval
    
    for (i in 1:length(Coeff_Names_readable)) {
      dfOut[ind, paste0(Coeff_Names_readable[[i]], "_pval")] <- Model$coefficients[Coeffs[i], "Pr(>|t|)"]
      dfOut[ind, paste0(Coeff_Names_readable[[i]], "_Estimate")] <- Model$coefficients[Coeffs[i], "Estimate"]
      dfOut[ind, paste0(Coeff_Names_readable[[i]], "_SE")] <- Model$coefficients[Coeffs[i], "Std. Error"]
      dfOut[ind, paste0(Coeff_Names_readable[[i]], "_tval")] <- Model$coefficients[Coeffs[i], "t value"]
    }
  }
  return(dfOut)
}

#### Create pretty qqplots

qqplotFunction <- function(expect, observ, col) {
  observ <- sort(observ, decreasing = T)
  dftemp <- data.frame(expect, observ, col)
  plottemp <- ggplot(data = dftemp) +
    geom_point(mapping = aes(x = expect, y = observ), colour = unique(dftemp$col), size = 3) +
    geom_abline(intercept = 0, slope = 1) +
    ylim(0, max(dftemp$observ)) +
    xlim(range(dftemp$expect)) +
    theme_grey(base_size = 22) +
    theme(axis.title = element_blank(), legend.position = "none")
  return(plottemp)
}

# the wrapper function

ModellingXY <- function(df, n_snp, n_repet, ModelFormulae, path_out) {
  
  n_sample <- unique(df$n_sample)
  
  AllParametersX <- c("reg_X_E", "reg_X_E2", "reg_X_GRS")
  
  LEVEL <- unique(df$LEVEL)
  # check: has individual level data of according or smaller size already been created?
  # if yes: load data and downsize to the selected size (n_repet and n_sample)
  tic("Load or create individual level data")
  if(file.exists(paste0(path_out, "individualSimulations_nS", n_sample, "_nR", n_repet, "_LEVEL_", LEVEL, ".Rdata"))) {
    load(paste0(path_out, "individualSimulations_nS", n_sample, "_nR", n_repet, "_LEVEL_", LEVEL, ".Rdata"), verbose = T)
    
    print("Existing Data was loaded!")
    
  } else {
    
    rawSimData <- lapply(c(1:nrow(df)), function(x) {
      generate_raw_data(
        a0 = df$a0[x], a1 = df$a1[x], a2 = df$a2[x], a3 = df$a3[x], 
        a4 = df$a4[x], a5 = df$a5[x], 
        b0 = df$b0[x], b1 = df$b1[x], b2 = df$b2[x],
        qX = df$qX[x], qY = df$qY[x], qY2 = df$qY2[x],
        SeedToSet = df$Seed[x], n_sample = n_sample,
        n_snp = n_snp, LEVEL = df$LEVEL[x], path_out = path_out 
      )
    })
    save(rawSimData, file = paste0(path_out, "individualSimulations_nS", n_sample, "_nR", n_repet, "_LEVEL_", LEVEL, ".Rdata"))
    print("Data was created and saved")
    
  }
  
  toc()
  
  if(unique(df$b0 == 0)){ #when E2 is not in X
    ModelsX <- lapply(c(1:length(rawSimData)), function(x) {
      runModelFun(df = rawSimData[[x]], ModelFormula = "X ~ E + GRS")
    })
  } else {
    ModelsX <- lapply(c(1:length(rawSimData)), function(x) {
      runModelFun(df = rawSimData[[x]], ModelFormula = "X ~ E + E2 + GRS")
    })
  }
  
  dfOverviewX <- createOverviewFun(Models = ModelsX, dfModel = df)
  
  ind_cols_pvalsX <- grep("_pval", colnames(dfOverviewX), value = FALSE) #get indices of cols with pvals
  cols_pvalsX <- grep("_pval", colnames(dfOverviewX), value = TRUE) #get names of cols with pvals
  names(dfOverviewX)[ind_cols_pvalsX] <- paste0("reg_X_", cols_pvalsX)
  cols_pvalsX <- names(dfOverviewX)[ind_cols_pvalsX]
  cols_pvalsX <- cols_pvalsX[!cols_pvalsX %in% c("reg_X_Intercept_pval")]
  
  
  ind_cols_estimatesX <- grep("_Estimate", colnames(dfOverviewX), value = FALSE)
  cols_estimatesX <- grep("_Estimate", colnames(dfOverviewX), value = TRUE)
  names(dfOverviewX)[ind_cols_estimatesX] <- paste0("reg_X_", cols_estimatesX)
  cols_estimatesX <- names(dfOverviewX)[ind_cols_estimatesX]
  cols_estimatesX <- cols_estimatesX[!cols_estimatesX %in% c("reg_X_Intercept_Estimate")]
  
  numberofestimatesX <- length(cols_estimatesX)
  
  ParamX <- str_replace_all(cols_pvalsX, "_pval", "")
  ParamValX <- c("a0", "b0", "", "qX")
  
  IndexcolsX <- rep(NA, times = length(ParamX))
  
  # mainly for accurate plot colors: get index of the analysed parameters within all potential parameters
  for (i in 1:length(ParamX)) {
    IndexcolsX[i] <- grep(paste0("^", ParamX[i], "$"), AllParametersX)
  }
  
  ParamValX <- ParamValX[IndexcolsX]
  
  logpval <- lapply(dfOverviewX[, cols_pvalsX], function(x) log10(x) %>% multiply_by(-1))
  cols_pval_log10X <- paste0(names(logpval), "_log10")
  dfOverviewX[, cols_pval_log10X] <- logpval
  

  ## FOR Y
  AllParameters <- c("E", "GRS", "X", "E2", "GRS2", "X2", "E3", "E4", 
                     "ExGRS", "ExX", "GRSxE2", "ExGRS2", "E2xGRS2")
  

  for(index1 in 1:nrow(ModelFormulae)){
    
    ModelFormula <- ModelFormulae[index1, "Formula"]
    
    tic("Fit model for according data")
    Models <- lapply(c(1:length(rawSimData)), function(x) {
      runModelFun(df = rawSimData[[x]], ModelFormula = ModelFormula)
    })
    toc()
    
    dfOverview <- createOverviewFun(Models = Models, dfModel = dfOverviewX)
    
    cols_pvals <- grep("_pval", colnames(dfOverview), value = TRUE) %>% grep(pattern = "reg_X_", x = ., invert = T, value = T)
    cols_pvals <- cols_pvals[!cols_pvals %in% c("resXGRS_pval", "resXGRS2_pval", "Intercept_pval")]
    Param <- str_replace_all(cols_pvals, "_pval", "")
    
    logpval <- lapply(dfOverview[, cols_pvals], function(x) log10(x) %>% multiply_by(-1))
    cols_pval_log10 <- paste0(names(logpval), "_log10")
    dfOverview[, cols_pval_log10] <- logpval
    
    cols_estimates <- grep("_Estimate", colnames(dfOverview), value = TRUE) %>% grep(pattern = "reg_X_", x = ., invert = T, value = T)
    cols_estimates <- cols_estimates[!cols_estimates %in% c("resXGRS_Estimate", "resXGRS2_Estimate", "Intercept_Estimate")]
    numberofestimates <- length(cols_estimates)
    
    Indexcols <- rep(NA, times = length(Param))
    
    for (i in 1:length(Param)) {
      Indexcols[i] <- grep(paste0("^", Param[i], "$"), AllParameters)
    }
    
    # Get estimates for the correction terms, correction terms, and corrected parameters
    
    if(grepl("E:GRS", ModelFormula) & grepl("GRS2", ModelFormula)){ 
      #some models are observational, others do not include any interaction or are missing the quadratic estimate
      
      print("Correction is being performed!")
      newcols <- paste0(rep(c("a0", "b0", "a4", "b2"), each = 3), c("_Estimate", "_SE", "_var"))
      dfOverview[,newcols] <- NA
      
      
      for (index2 in 1:n_repet){
        
        coefs <- lm(X ~ E + E2, data = rawSimData[[index2]]) %>% summary() %>% coef()
        dfOverview$a0_Estimate[index2] <- coefs["E", "Estimate"]
        dfOverview$a0_SE[index2] <- coefs["E", "Std. Error"]
        dfOverview$a0_var[index2] <- (coefs["E", "Std. Error"]*sqrt(n_sample))^2
        
        dfOverview$b0_Estimate[index2] <- coefs["E2", "Estimate"]
        dfOverview$b0_SE[index2] <- coefs["E2", "Std. Error"]
        dfOverview$b0_var[index2] <- (coefs["E", "Std. Error"]*sqrt(n_sample))^2
        
        coefs <- lm(X ~ GRS, data = rawSimData[[index2]]) %>% summary() %>% coef()
        dfOverview$a4_Estimate[index2] <- coefs["GRS", "Estimate"]
        dfOverview$a4_SE[index2] <- coefs["GRS", "Std. Error"]
        dfOverview$a4_var[index2] <- (coefs["GRS", "Std. Error"]*sqrt(n_sample))^2
        
      }
      
      
      #### Correction of parameters
      
      # obtain variances of ExGRS and GRS2 parameters from SE
      if("ExGRS" %in% AllParameters){
        dfOverview$ExGRS_var <- (dfOverview$ExGRS_SE*sqrt(n_sample))^2
      } else if("GRSxE" %in% AllParameters){
        dfOverview$ExGRS_var <- (dfOverview$GRSxE_SE*sqrt(n_sample))^2
      }
      
      dfOverview$GRS2_var <- (dfOverview$GRS2_SE*sqrt(n_sample))^2
      
      #dfOverview$E2xGRS_var <- (dfOverview$GRSxE2_SE * sqrt(n_sample))^2
      
      # ExGRS correction
      
      dfOverview$CorrectionVar_ExGRS <- (dfOverview$GRS2_Estimate^2/dfOverview$a4_Estimate^2) *
        ((dfOverview$GRS2_var/dfOverview$GRS2_Estimate^2)+
           (dfOverview$a4_var/dfOverview$a4_Estimate^2)) *
        dfOverview$a0_Estimate^2 *2^2
    
      
      if("ExGRS" %in% AllParameters){
        ExGRS_col <- "ExGRS_Estimate"
      } else if("GRSxE" %in% AllParameters){
        ExGRS_col <- "GRSxE_Estimate"
      }
      
      dfOverview[,"CorrectedEstimate_ExGRS"] <- dfOverview[, ExGRS_col] - 
        (2 * dfOverview[,"a0_Estimate"] *
           (dfOverview[,"GRS2_Estimate"]/dfOverview[,"a4_Estimate"]))
      
      dfOverview[,"CorrectedVar_ExGRS"] <- (dfOverview$ExGRS_var + 
                                              dfOverview$CorrectionVar_ExGRS)
      
      dfOverview[,"CorrectedSE_ExGRS"] <- sqrt(dfOverview$CorrectedVar_ExGRS)/sqrt(n_sample)
      mu = 0
      degfree = n_sample - 2
      
      dfOverview[,"CorrectedTstat_ExGRS"] <- (dfOverview$CorrectedEstimate_ExGRS - mu)/
        dfOverview$CorrectedSE_ExGRS
      dfOverview[,"CorrectedPval_ExGRS"] <- 2 * pt(-abs(dfOverview$CorrectedTstat_ExGRS), degfree)
      
      dfOverview$CorrectedPval_ExGRS_log10 <- log10(dfOverview$CorrectedPval_ExGRS) %>% 
        multiply_by(-1)
      
      
      if("ExGRS" %in% AllParameters){
        cols_estimates_old <- "ExGRS_Estimate"
      } else if("GRSxE" %in% AllParameters){
        cols_estimates_old <- "GRSxE_Estimate"
      }
      
      cols_estimates_new <- "CorrectedEstimate_ExGRS"
      
      df$expectedpvals <- -log10(seq(from = 1 / nrow(dfOverview), to = 1, length.out = nrow(dfOverview)))
      
      if("ExGRS" %in% AllParameters){
        cols_pval_log10_old <- c("ExGRS_pval_log10")
      } else if("GRSxE" %in% AllParameters){
        cols_pval_log10_old <- c("GRSxE_pval_log10")
      }
      
    }
    
    CSVFileName <- paste0(path_out, "Level_", LEVEL, "nS", n_sample, "_Model_", ModelFormulae[index1, "Model_no"], "_Data.csv")
    write_csv(dfOverview, file = CSVFileName) # on Jura: path, on private laptop: file
  }
}

