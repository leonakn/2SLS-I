# lknusel
# obtain exposure phenotype data

library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)


# step 1: get trait ids.
trait_id <- snakemake@params[["exposure_IDs"]]

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
	  grep(paste0('^', ID, '-0*'), list_colnames)
	}
	File = sapply(All_colnames, function(x) find_column(x, ID))
	File = File[sapply(File, length)>0]
    if(length(File) > 0){
		# By default, use the most recent file version
		file_number <- str_remove_all(names(File), paste0(path_to_search, "/ukb|.csv")) |> as.numeric()
		max_file_number_index <- which(file_number == max(file_number))
		my_File = names(File)[max_file_number_index]
		# By default, use 1st visit (baseline characteristics)
		my_Column = File[[max_file_number_index]][1]

		if(grepl(paste0("^", ID, "-0*"), All_colnames[[my_File]][my_Column])){
			return(c(my_File, my_Column))
		} else {
			print("Something went wrong, try again!")
		}
	} else {
	stop("Variable not found")
	}
 }

correct_pheno_get_IRNT <- function(Variable, CovariateMatrix){

   VariableAvailable <- !is.na(Variable)
   CorrectedMeasures <- rep(NA, length(Variable))
   CorrectedMeasures[VariableAvailable] <- residuals(lm(Variable ~ CovariateMatrix, na.action = na.omit))
   CorrectedMeasures <- qnorm((rank(CorrectedMeasures, na.last = "keep") - 0.5)/sum(!is.na(CorrectedMeasures)))

   return(CorrectedMeasures)
}


# step 3: obtain phenotype data

# if the trait is categoric, make sure you only keep the UKBB ID and not the
# according categoric value, otherwise won't find the file.

if(grepl("_", trait_id)){
	UKBB_trait_ID <- str_split(trait_id, pattern = "_")[[1]][1]
	pheno_file <- find_file_function(UKBB_trait_ID)
} else {
	pheno_file <- find_file_function(trait_id)
}


pheno_df <- fread(pheno_file[1], select = c(1, as.numeric(pheno_file[2])), data.table = FALSE)
# 502616 participants

print(paste0("Number of participants at start = ", nrow(pheno_df)))

# step 3: keep only relevant participants

participants <- read.csv(snakemake@params[["file_wb_participants"]])

# down to 337392
pheno_df <- pheno_df[pheno_df$eid %in% participants$eid,]

rm(participants)

print(paste0("Filtering for active participants, check. Number of active wb participants = ", 
             nrow(pheno_df)))

# step 4: load covariate matrix, merge relevant demographic data with pheno df

load(snakemake@input[["file_covariates"]])

#remove the raw age column to only correct for the more accurate age effect.
duplicated_column <- grep("^age$", colnames(X))
X <- X[,-duplicated_column]


# Step 5: check if eids from X and pheno df match and get rid of eid in X

eids_both_index <- which(X[,"eid"] %in% pheno_df$eid)

X <- X[eids_both_index,]

eid_col <- grep("eid", colnames(X))

# rounding is necessary, because the annotation 
# differed (scientific vs. standard) and
# as a consequence one value differed minimally.

eids_covariate_matrix <- X[,eid_col] |> as.numeric() |> round()
eids_pheno_df <- pheno_df$eid |> as.numeric() |> round()

if(!all(eids_covariate_matrix == eids_pheno_df)){
	X <- X[X[,eid_col] %in% pheno_df$eid, ]
	pheno_df <- pheno_df[pheno_df$eid %in% X[,eid_col], ]
	eid_order <- match(pheno_df$eid, X[,eid_col])
	X <- X[eid_order, ]

    if(!all(X[,eid_col] == pheno_df$eid)){
      stop("HOUSTON, WE'VE HAD A PROBLEM!! Eids in covariate matrix and pheno df do not match!")
	}
}

X <- X[,-eid_col]


# Step 6: standardize data, obtain residuals after including covariates in the model
if(grepl("_", trait_id)){
	UKBB_trait_ID <- str_split(trait_id, pattern = "_")[[1]][1]
	raw_pheno_column <- grep(paste0("^", UKBB_trait_ID, "-"), 
	                     colnames(pheno_df), value = T)
	new_raw_pheno_column <- paste0(trait_id, "-0.0")
	pheno_df[, new_raw_pheno_column] <- pheno_df[,raw_pheno_column]
}

raw_pheno_column <- grep(paste0("^", trait_id, "-"), 
	                     colnames(pheno_df), value = T)

print(paste0("The raw pheno column is = ", raw_pheno_column))
print("Head of the raw_pheno_column is:")

head(pheno_df[,raw_pheno_column])

if(length(raw_pheno_column) != 1){
	stop("not one single phenotype column selected!")
}

# removed all categoric phenotypes, but this would work for 
# categoric data, too.
if(grepl("_", trait_id)){

	print("Working on a categoric phenotype!")
	print(paste0("Number of NAs initially in phenotype = ", sum(is.na(pheno_df[, raw_pheno_column]))))
	
	remove_missing_values_index <- which(pheno_df[,raw_pheno_column] < 0)
	#remove_missing_values_index <- remove_missing_values_index[,1]

	print(paste0("Class of remove_missing_values_index is = ", 
				 class(remove_missing_values_index)))
	
	pheno_df[remove_missing_values_index, raw_pheno_column] <- NA

	print(paste0("Number of NAs after removing prefer not to answer = ", sum(is.na(pheno_df[, raw_pheno_column]))))

	effect_category <- str_split(trait_id, pattern = "_")[[1]][2] |> as.numeric() # get the number of the effect category

	effect_category_index <- pheno_df[, raw_pheno_column] == effect_category
	effect_category_index <- effect_category_index[,1]

	print(paste0("Class of effect_category_index is = ",
	             class(effect_category_index), ". And the length is = ",
				 length(effect_category_index)))

	pheno_df[, raw_pheno_column] <- ifelse(effect_category_index, 
	                                         yes = 1, no = 0)
	
	pheno_df[,paste0(trait_id, "_z")] <- scale(pheno_df[, raw_pheno_column])
	pheno_df[,paste0(trait_id, "_IRNT")] <- qnorm((rank(pheno_df[, raw_pheno_column])-0.5)/
	                                            sum(!is.na(pheno_df[, raw_pheno_column])))
	pheno_df[,paste0(trait_id, "_IRNT_c")] <- 
	    correct_pheno_get_IRNT(Variable = unlist(pheno_df[, raw_pheno_column]),
			                     CovariateMatrix = X)

} else if(trait_id %in% snakemake@params[["remove_negative_code"]]){

	print("Working on a phenotype with negative answer categories!")
	print(paste0("Number of NAs initially in phenotype = ", sum(is.na(pheno_df[, raw_pheno_column]))))

    # in some phenos, there is a code -10 indicating less than 1
	# replaced this by 0.3 (guess)
	set_min_10_to_0.3 <- which(pheno_df[, raw_pheno_column] == -10)
	pheno_df[set_min_10_to_0.3, raw_pheno_column] <- 0.3
	
	# replace the "do not know"/"prefer not to answer" options with NAs
	remove_missing_values_index <- which(pheno_df[, raw_pheno_column] < 0)
	pheno_df[remove_missing_values_index, raw_pheno_column] <- NA

	print(paste0("Number of NAs after removing prefer not to answer = ", sum(is.na(pheno_df[,raw_pheno_column]))))

	pheno_df[,paste0(trait_id, "_z")] <- scale(pheno_df[,raw_pheno_column])
	pheno_df[,paste0(trait_id, "_IRNT")] <- qnorm((rank(pheno_df[,raw_pheno_column])-0.5)/
	                                            sum(!is.na(pheno_df[,raw_pheno_column])))
	pheno_df[,paste0(trait_id, "_IRNT_c")] <- 
	    correct_pheno_get_IRNT(Variable = unlist(pheno_df[,raw_pheno_column]),
			                     CovariateMatrix = X)

} else if(trait_id == "845"){ 

	pheno_file_2 <- find_file_function(6138)
    pheno_df_2 <- fread(pheno_file_2[1], select = c(1, as.numeric(pheno_file_2[2])), data.table = FALSE)
    print(paste0("Number of participants at start = ", nrow(pheno_df_2)))

    # make sure you're working on the same individuals
	# in the same order
	pheno_df_2 <- pheno_df_2[pheno_df_2$eid %in% pheno_df$eid, ]
	sort_ids <- match(pheno_df$eid, pheno_df_2$eid)
	pheno_df_2 <- pheno_df_2[sort_ids, ]

	if(!all(pheno_df$eid == pheno_df_2$eid)){
		stop("Phenotype education: something went wrong with eids in df1 and df2")
	}
    
	# replace NAs in 845 due to "University degree" by 19 years
	pheno_df[, raw_pheno_column] <- ifelse(pheno_df_2$'6138-0.0' == "1" , 
	                                       19, pheno_df[,raw_pheno_column])
	# replace never went to school with 0
    pheno_df[, raw_pheno_column] <- ifelse(pheno_df[,raw_pheno_column] == "-2", 
	                                       0, pheno_df[,raw_pheno_column])
	# replace "do not know" with NA
	pheno_df[, raw_pheno_column] <- ifelse(pheno_df[,raw_pheno_column] == "-1", 
	                                       NA, pheno_df[,raw_pheno_column])
	# replace prefer not to answer with NA
	pheno_df[, raw_pheno_column] <- ifelse(pheno_df[,raw_pheno_column] == "-3", 
	                                       NA, pheno_df[,raw_pheno_column])

	pheno_df[,paste0(trait_id, "_z")] <- scale(pheno_df[,raw_pheno_column])
	pheno_df[,paste0(trait_id, "_IRNT")] <- qnorm((rank(pheno_df[,raw_pheno_column])-0.5)/
	                                            sum(!is.na(pheno_df[,raw_pheno_column])))
	pheno_df[,paste0(trait_id, "_IRNT_c")] <- 
	    correct_pheno_get_IRNT(Variable = unlist(pheno_df[,raw_pheno_column]),
			                     CovariateMatrix = X)

} else {

	pheno_df[,paste0(trait_id, "_z")] <- scale(pheno_df[,raw_pheno_column])
	pheno_df[,paste0(trait_id, "_IRNT")] <- qnorm((rank(pheno_df[,raw_pheno_column])-0.5)/
	                                            sum(!is.na(pheno_df[,raw_pheno_column])))
	pheno_df[,paste0(trait_id, "_IRNT_c")] <- 
	    correct_pheno_get_IRNT(Variable = unlist(pheno_df[,raw_pheno_column]),
			                     CovariateMatrix = X)
}

print("Columns in pheno_df are:")
print(colnames(pheno_df))

# get rid of demographic columns, as they are in the outcome phenotype
cols_to_keep <- c(grep(trait_id, colnames(pheno_df), value = T), "eid")

pheno_df <- pheno_df[,cols_to_keep]

print(paste0("After exclusion of demographic info, the following columns remain:"))
colnames(pheno_df)

save(pheno_df, file = snakemake@output[["file_phenotypes"]])
print("pheno data corrected and saved!")
