# script to obtain interaction between E and SNP on X
# lknusel

#load packages
library(dplyr)
library(stringr)
library(ggplot2)
library(magrittr)

# define_function

get_snp_estiamtes <- function(df = NA, covariates = NA, 
                              Y = NA, interaction_parameter = NA){

    coefficients <- as.data.frame(matrix(data = NA, ncol = 8, nrow = 10))
    colnames(coefficients) <- c("SNP", "Parameter", "Estimate", "SE", "t-val", 
			 	                "p-val", "nSample", "PhenoFormat")

    df$E <- df[,interaction_parameter]
    df$E2 <- df$E^2
    df$Y <- df[,Y]
	
	snps <- grep("^rs", colnames(df), value = T)

	Out <- lapply(1:length(snps), function(i){

        # with a little help from your friends:
        # https://stackoverflow.com/questions/71737719/is-there-a-way-to-loop-through-column-names-not-numbers-in-r-for-linear-models
        if(!is.na(covariates)){
                model_formula <- formula(paste("Y ~ ", snps[i], " * E + ", snps[i], " * E2 + ", covariates))
        } else {
                model_formula <- formula(paste("Y ~ ", snps[i], " * E + ", snps[i], " * E2"))
        }
        
        model <- lm(model_formula, data = df) |> summary()

        parameters <- grep("PC|Intercept", rownames(model$coefficients), invert = T, value = T)
        relevant_model_coefs <- model$coefficients[parameters,]
		
        if(any(class(relevant_model_coefs) == "matrix")){
            number_of_parameters <- nrow(relevant_model_coefs)
            coefficients[1:number_of_parameters, 3:6] <- relevant_model_coefs[1:number_of_parameters,]
        } else if(class(relevant_model_coefs) == "numeric"){
            number_of_parameters <- 1
            coefficients[1, 3:6] <- relevant_model_coefs
        } else {
            warning(paste0("For snp ", snps[i], " there was an issue with the obtained model coefficients."))
        }
		
		coefficients[c(1:number_of_parameters), "SNP"] <- snps[i]
		coefficients[c(1:number_of_parameters),"Parameter"] <- parameters
		coefficients$PhenoFormat <- Y
    	coefficients$nSample <- nrow(df)
			
		coefficients <- coefficients[!is.na(coefficients$SNP), ]
        coefficients
	})
	Out <- do.call(rbind, Out)
	return(Out)
}

qqplotFunction <- function(expect, observ) {
  observ <- sort(observ, decreasing = T)
  dftemp <- data.frame(expect, observ)
  plottemp <- ggplot(data = dftemp) +
    geom_point(mapping = aes(x = expect, y = observ), size = 3) +
    geom_abline(intercept = 0, slope = 1) +
    ylim(0, max(dftemp$observ)) +
    xlim(range(dftemp$expect)) +
    theme_grey(base_size = 22) +
    theme(axis.title = element_blank(), legend.position = "none")
  return(plottemp)
}

interaction_ID <- snakemake@params[["interaction_ID"]]
exposure_IDs <- snakemake@params[["exposure_IDs"]]

load(snakemake@input[["file_exposure_phenotypes"]])

print(paste0("file_exposure_phenotype loaded. Number of individuals = ", 
             nrow(pheno_df)))

load(snakemake@input[["file_genotypes"]])

print(paste0("file_genotype loaded. Number of individuals = ", 
             nrow(lsGeno$data)))

load(snakemake@input[["file_covariates"]])

print(paste0("file_covariates loaded. Number of individuals = ", 
             nrow(X)))

load(snakemake@input[["file_interaction_phenotype"]])

print(paste0("file interaction pheno loaded. Number of individuals = ", 
             nrow(interaction_df)))

geno_df <- lsGeno$data

df_model <- merge(pheno_df, geno_df, by = "eid") |> as.data.frame()
rm(lsGeno, pheno_df)
print(paste0("Geno and pheno df are loaded and merged! Number of individuals = ", nrow(df_model)))

df_model <- merge(df_model, interaction_df, by = "eid")
rm(interaction_df)
print(paste0("df_model and interaction_df are merged! Number of individuals = ", 
             nrow(df_model)))

duplicated_column <- grep("^age$", colnames(X))
X <- X[,-duplicated_column]

# note: for some unknown reason, the number of participants in the updated
# pheno file is *smaller* than the number of wb participants, and thus the 
# covariate file. So here, remove all participants witch are in covariate 
# file but not in pheno file...

eids_both_index <- which(X[,"eid"] %in% df_model$eid)

X <- X[eids_both_index,]

eid_col <- grep("eid", colnames(X))

# rounding is necessary, because for some wierd reason,
# the annotation differed (scientific vs. standard) and
# as a consequence one value differed minimally.

eids_covariate_matrix <- X[,eid_col] |> as.character() |> as.numeric() |> round()
eids_df <- df_model$eid |> as.character() |> as.numeric() |> round()

if(!all(eids_covariate_matrix == eids_df)){
	stop("HOUSTON, WE'VE HAD A PROBLEM!! Eids in covariate matrix and df do not match!")
} else {
    if(interaction_ID != 21022){ 
        #make sure not to have age twice in your model...
        Covariates <- X[,c(paste0("PC", c(1:40)), "age_exact", "age_exact2", "sex")]
    } else {
        Covariates <- X[,c(paste0("PC", c(1:40)), "sex")]
    }
}

df_model <- cbind(df_model, Covariates)

Y <- paste0(exposure_IDs, "_IRNT_c")
print(paste0("Exposure (X) column is = ", Y))

E <- paste0(interaction_ID, "_z_c")
print(paste0("Interaction (E) column is = ", E))

# define covariates to include in the model...
# Always: 40 PCs & sex
# Only if interaction parameter != age: age and age^2
# age and age2 will be included automatically as part of the
# interaction, so avoid duplication!

PCs <- paste0("PC", c(1:40)) |> paste0(collapse = " + ")
if(interaction_ID != 21022){
    model_covariates <- paste0("sex + age_exact + ", 
                               "age_exact2 + ",
                               PCs)
} else {
    model_covariates <- paste0("sex + ", 
                               PCs)
}

#### NOW CHECK FOR INTERACTION OF SNP WITH E ON PHENOTYPE!!                           
snp_E_interactions <- get_snp_estiamtes(df = df_model, Y = Y, 
                                        covariates = model_covariates,
                                        interaction_parameter = E)

save(snp_E_interactions, file = snakemake@output[["snp_E_interactions"]])

print("SNP*E interactions models fitted and saved")

# obtain only pvalues from the interaction of SNP with age
pvals_SNP_E <- snp_E_interactions[grep("rs[0-9]+:E$", 
                                       snp_E_interactions$Parameter), 
                                       "p-val"] |> 
        log10() |> 
        multiply_by((-1))


pvals_expect <- seq(from = 1 / length(pvals_SNP_E), to = 1, 
                    length.out = length(pvals_SNP_E)) |> 
    log10() |> 
    multiply_by((-1))

qq_E_x_SNP <- qqplotFunction(pvals_expect, pvals_SNP_E)

ggsave(snakemake@output[["qq_plot_SNP_x_E_effects"]], plot = qq_E_x_SNP)

print("QQ plot for E * SNP effect plotted and saved!")

# same for age2
pvals_SNP_E2 <- snp_E_interactions[grep("rs[0-9]+:E2$", 
                                        snp_E_interactions$Parameter), 
                                        "p-val"] |> 
        log10() |> 
        multiply_by((-1))

pvals_expect <- seq(from = 1 / length(pvals_SNP_E2), to = 1, 
                    length.out = length(pvals_SNP_E2)) |> 
    log10() |> 
    multiply_by((-1))

qq_E2_x_SNP <- qqplotFunction(pvals_expect, pvals_SNP_E2)

ggsave(snakemake@output[["qq_plot_SNP_x_E2_effects"]], plot = qq_E2_x_SNP)

print("QQ plot for E2 * SNP effect plotted and saved!")
