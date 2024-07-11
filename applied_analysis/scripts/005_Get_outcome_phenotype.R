# lknusel
# obtain outcome phenotypes


library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)

# step 1: get trait ids.
trait_id <- snakemake@params[["outcome_ID"]]

# step 2: define functions
find_file_function <- function(ID){
	path_to_search <- snakemake@params[["path_ukbb_files"]]
	library(stringr)
	# List all UKBB pheno files
	All_files <- list.files(path = path_to_serach, pattern = 'ukb.*[.]csv', full.names = T)
	# Get colnames for each file
	get_colnames <- function(file){
	  colnames(data.table::fread(file, nrow = 1))
	}
	All_colnames <- sapply(All_files, function(x) get_colnames(x))
	# Look if the variable is in the colnames of any of the files
	find_column <- function(list_colnames, ID){
	  grep(paste0('^', ID, '-0.0'), list_colnames)
	}
	File = sapply(All_colnames, function(x) find_column(x, ID))
	File = File[sapply(File, length)>0]
    if(length(File) > 0){
	# By default, use the most recent file version
	file_number <- str_remove_all(names(File), 
								  paste0(path_to_search, "/ukb|.csv")) |> 
		as.numeric()
	max_file_number_index <- which(file_number == max(file_number))
	my_File = names(File)[max_file_number_index]
	# By default, use 1st visit (baseline characteristics)
	my_Column = File[[max_file_number_index]][1]
	return(c(my_File, my_Column))
     } else {
	stop("Variable not found")
	}
 }

correct_outcome_function <- function(Variable, CovariateMatrix){

   VariableAvailable <- !is.na(Variable)
   CorrectedMeasures <- rep(NA, length(Variable))
   CorrectedMeasures[VariableAvailable] <- residuals(lm(Variable ~ CovariateMatrix, na.action = na.omit))
   CorrectedMeasures <- scale(CorrectedMeasures)

   return(CorrectedMeasures)
}

correct_pheno_get_IRNT <- function(Variable, CovariateMatrix){

   VariableAvailable <- !is.na(Variable)
   CorrectedMeasures <- rep(NA, length(Variable))
   CorrectedMeasures[VariableAvailable] <- residuals(lm(Variable ~ CovariateMatrix, na.action = na.omit))
   CorrectedMeasures <- qnorm((rank(CorrectedMeasures, na.last = "keep") - 0.5)/sum(!is.na(CorrectedMeasures)))

   return(CorrectedMeasures)
}


# step 3: obtain phenotype data
pheno_file <- find_file_function(trait_id)

outcome_df <- fread(pheno_file[1], select = c(1, as.numeric(pheno_file[2])), data.table = F)

print(paste0("Number of participants at start = ", nrow(outcome_df)))

# step 3: keep only relevant participants

participants <- read.csv(snakemake@params[["file_wb_participants"]])

outcome_df <- outcome_df[outcome_df$eid %in% participants$eid,]

rm(participants)

print(paste0("Filtering for active participants, check. Number of active wb participants = ", 
             nrow(outcome_df)))

# step 4: load covariate matrix, merge relevant demographic data with pheno df

load(snakemake@input[["file_covariates"]])

demographics_columns <- grep("[A|a]ge|[S|s]ex|[E|e]id", colnames(X))

demographics_df <- as.data.frame(X[,demographics_columns])

#remove the raw age column to only correct for the more accurate age effect.
duplicated_column <- grep("^age$", colnames(demographics_df))
demographics_df <- demographics_df[,-duplicated_column]

outcome_df <- merge(outcome_df, demographics_df, by = "eid")

print(paste0("outcome_df and demographics_df are merged, n individuals = ", 
             nrow(outcome_df)))

# Step 5: check if eids from X and pheno df match and get rid of eid in X

eids_both_index <- which(X[,"eid"] %in% outcome_df$eid)

X <- X[eids_both_index,]

eid_col <- grep("eid", colnames(X))

# rounding is necessary, because for some wierd reason,
# the annotation differed (scientific vs. standard) and
# as a consequence one value differed minimally.

eids_covariate_matrix <- X[,eid_col] |> as.numeric() |> round()
eids_outcome_df <- outcome_df$eid |> as.numeric() |> round()

if(!all(eids_covariate_matrix == eids_outcome_df)){
	X <- X[X[,eid_col] %in% outcome_df$eid, ]
	outcome_df <- outcome_df[outcome_df$eid %in% X[,eid_col], ]
	eid_order <- match(outcome_df$eid, X[,eid_col])
	X <- X[eid_order, ]

    if(!all(X[,eid_col] == outcome_df$eid)){
      stop("HOUSTON, WE'VE HAD A PROBLEM!! Eids in covariate matrix and pheno df do not match!")
	}
}

X <- X[,-eid_col]

# Step 6: standardize data, obtain residuals after including covariates in the model

raw_pheno_column <- grep(paste0("^", trait_id, "-"), 
	                     colnames(outcome_df), value = T)

print(raw_pheno_column)
print(head(outcome_df))

if(length(raw_pheno_column) != 1){
	stop("not one single phenotype column selected!")
}

outcome_df[,paste0(trait_id, "_z")] <- scale(outcome_df[,raw_pheno_column])

outcome_df[,paste0(trait_id, "_z_c")] <- 
	    correct_outcome_function(Variable = unlist(outcome_df[,raw_pheno_column]),
			                         CovariateMatrix = X)

outcome_df[,paste0(trait_id, "_IRNT_c")] <- 
        correct_pheno_get_IRNT(Variable = unlist(outcome_df[,raw_pheno_column]),
		                       CovariateMatrix = X)

save(outcome_df, file = snakemake@output[["file_outcome_phenotype"]])
print("pheno data corrected and saved!")
