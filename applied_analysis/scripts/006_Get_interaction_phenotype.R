# lknusel
# get environment phenotype

library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)

# step 1: get trait ids.
trait_id <- snakemake@params[["interaction_ID"]]

# step 2: define functions
find_file_function <- function(ID){
	path_to_search <- snakemake@params[["path_ukbb_files"]]
	library(stringr)
	# List all UKBB pheno files
	All_files <- list.files(path = path_to_search, pattern = 'ukb.*[.]csv', full.names = T)
	# Get colnames for each file
	get_colnames <- function(file){
	  colnames(data.table::fread(file, nrow = 1))
	}
	All_colnames <- sapply(All_files, function(x) get_colnames(x))
	# Look if the variable is in the colnames of any of the files
	find_column <- function(list_colnames, ID){
	  grep(paste0('^', ID, '-*'), list_colnames)
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

correct_interaction_function <- function(Variable, CovariateMatrix, zstd = TRUE){

   VariableAvailable <- !is.na(Variable)
   CorrectedMeasures <- rep(NA, length(Variable))
   CorrectedMeasures[VariableAvailable] <- residuals(lm(Variable ~ CovariateMatrix, na.action = na.omit))
   if(zstd){
	CorrectedMeasures <- scale(CorrectedMeasures)
   }
   return(CorrectedMeasures)
}


# step 3: obtain phenotype data
pheno_file <- find_file_function(trait_id)

interaction_df <- fread(pheno_file[1], select = c(1, as.numeric(pheno_file[2])), data.table = F)

print(paste0("Number of participants at start = ", nrow(interaction_df)))

# step 3: keep only relevant participants

participants <- read.csv(snakemake@params[["file_wb_participants"]])

interaction_df <- interaction_df[interaction_df$eid %in% participants$eid,]

rm(participants)

print(paste0("Filtering for active participants, check. Number of active wb participants = ", 
             nrow(interaction_df)))

load(snakemake@input[["file_covariates"]])


# Step 5: check if eids from X and pheno df match and get rid of eid in X

eids_both_index <- which(X[,"eid"] %in% interaction_df$eid)

X <- X[eids_both_index,]

eid_col <- grep("eid", colnames(X))

# rounding is necessary, because for some reason,
# the annotation differed (scientific vs. standard) and
# as a consequence one value differed minimally.

eids_covariate_matrix <- X[,eid_col] |> as.numeric() |> round()
eids_interaction_df <- interaction_df$eid |> as.numeric() |> round()

if(!all(eids_covariate_matrix == eids_interaction_df)){
	new_order <- match(eids_covariate_matrix, eids_interaction_df)
	interaction_df <- interaction_df[new_order,]
	eids_interaction_df <- interaction_df$eid |> as.numeric() |> round()
	}

if(!all(eids_covariate_matrix == eids_interaction_df)){
	stop("HOUSTON, WE'VE HAD A PROBLEM!! Eids in covariate matrix and pheno df do not match!")
} 

X <- X[,-eid_col]
PC_cols <- grep("^PC", colnames(X))
X <- X[,-PC_cols]

#make sure you don't end up correcting age for itself...
if(trait_id == "21022"){
	age_related_columns <- grep("[A|a]ge", colnames(X))
	X <- X[,-age_related_columns]
} else {
	duplicated_columns <- grep("^age$", colnames(X))
	X <- X[,-duplicated_columns]
}

# Step 6: standardize data, obtain residuals after regressing out covariates

raw_pheno_column <- grep(paste0("^", trait_id, "-"), 
	                     colnames(interaction_df), value = T)

if(length(raw_pheno_column) != 1){
	stop("not one single phenotype column selected!")
}

# do some preprocessing:
if(trait_id == "40049" | trait_id == "40048" | trait_id == "40047"){
	# note: accelerometer data requires some quality control
	# 90015 is 0 for everyone with total wear time < 72 h and who did not manage to
	# cover every hour of the day at least once
	correction_file <- find_file_function("90015")
	correction_df <- fread(correction_file[1], select = c(1, as.numeric(correction_file[2])))
	eids_to_set_to_NA <- correction_df[correction_df$'90015-0.0' == 0, "eid"] #this works
	print(paste0("Initial number of not NA individuals of interaction phenotype = ",  
				 sum(!is.na(interaction_df[, raw_pheno_column]))))
	interaction_df[interaction_df$eid %in% eids_to_set_to_NA$eid, raw_pheno_column] <- NA
	print(paste0("Number of not NA individuals of interaction phenotype after quality control = ",  
				 sum(!is.na(interaction_df[, raw_pheno_column]))))
}

if(trait_id == "1070"){
	# set -1 (do not know) and -3 (prefer not to answer) to NA
	min_1_eids <- interaction_df[interaction_df[, raw_pheno_column] == -1, "eid"]
	min_3_eids <- interaction_df[interaction_df[, raw_pheno_column] == -3, "eid"]
	eids_to_set_to_NA <- c(min_1_eids, min_3_eids)
	interaction_df[interaction_df$eid %in% eids_to_set_to_NA, raw_pheno_column] <- NA

	# set -10 (less than one hour) to 0.3
	min_10_eids <- interaction_df[interaction_df[, raw_pheno_column] == -10, "eid"]
	interaction_df[interaction_df$eid %in% min_10_eids, raw_pheno_column] <- 0.3

}

if(trait_id == "1239"){ #binarize for simplicity
	# set -3 to NA (prefer not to answer)
	eids_to_set_to_NA <- interaction_df[interaction_df[, raw_pheno_column] == -3, "eid"]
	interaction_df[interaction_df$eid %in% eids_to_set_to_NA, raw_pheno_column] <- NA
	# set 2 to 1 (only occasionally)
	eids_to_set_to_1 <- interaction_df[interaction_df[, raw_pheno_column] == 2, "eid"]
	interaction_df[interaction_df$eid %in% eids_to_set_to_1, raw_pheno_column] <- 1
}

if(trait_id == "845"){
	# replace NAs of educated people with 19 years, 
	# replace never went to school with 0, do not know with NA 
	# and prefer not to answer also with NA
	pheno_file_2 <- find_file_function(6138)
    interaction_df_2 <- fread(pheno_file_2[1], select = c(1, as.numeric(pheno_file_2[2])), data.table = FALSE)
    print(paste0("Number of participants at start = ", nrow(interaction_df_2)))

    # make sure you're working on the same individuals
	# in the same order
	interaction_df_2 <- interaction_df_2[interaction_df_2$eid %in% interaction_df$eid, ]
	sort_ids <- match(interaction_df$eid, interaction_df_2$eid)
	interaction_df_2 <- interaction_df_2[sort_ids, ]

	if(!all(interaction_df$eid == interaction_df_2$eid)){
		stop("Phenotype education: something went wrong with eids in df1 and df2")
	}
    
	# replace NAs in 845 due to "University degree" by 19 years
	interaction_df[, raw_pheno_column] <- ifelse(interaction_df_2$'6138-0.0' == "1" , 
	                                       19, interaction_df[,raw_pheno_column])
	# replace never went to school with 0
    interaction_df[, raw_pheno_column] <- ifelse(interaction_df[,raw_pheno_column] == "-2", 
	                                       0, interaction_df[,raw_pheno_column])
	# replace "do not know" with NA
	interaction_df[, raw_pheno_column] <- ifelse(interaction_df[,raw_pheno_column] == "-1", 
	                                       NA, interaction_df[,raw_pheno_column])
	# replace prefer not to answer with NA
	interaction_df[, raw_pheno_column] <- ifelse(interaction_df[,raw_pheno_column] == "-3", 
	                                       NA, interaction_df[,raw_pheno_column])

}

interaction_df[,paste0(trait_id, "_z")] <- scale(interaction_df[, raw_pheno_column])
interaction_df[,paste0(trait_id, "_z_c")] <- 
	    correct_interaction_function(Variable = unlist(interaction_df[, raw_pheno_column]),
			                         CovariateMatrix = X, zstd = TRUE)
interaction_df[,paste0(trait_id, "_c")] <- 
	    correct_interaction_function(Variable = unlist(interaction_df[, raw_pheno_column]),
			                         CovariateMatrix = X, zstd = FALSE)

save(interaction_df, file = snakemake@output[["file_interaction_phenotype"]])
print("pheno data corrected and saved!")
