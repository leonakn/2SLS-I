# leona knusel
# get relevant covariates

# load packages
library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)


# step 2: obtain demographics data
header <- colnames(fread(snakemake@input[["file_demographics"]], nrows = 0, header = T))
# "/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb21067.csv"

#select only time point 1 measurments
column_age <- grep("eid|21003-0.0|^31-0.0", header)

age_df <- as.data.frame(fread(snakemake@input[["file_demographics"]], select = (column_age), header = T))

#select parameteres to get a more precise measure of Age:
columns_age_detailed <- grep("eid|^52-0.0|^34-0.0|^53-0.0", header)
age_df_2 <- fread(snakemake@input[["file_demographics"]], 
                  select = (columns_age_detailed), header = T) |> as.data.frame()

months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

age_df_2$month_of_birth <- months[age_df_2$"52-0.0"]
age_df_2$day_of_birth <- 15

age_df_2 <- mutate(age_df_2, DateOfBirth = paste(age_df_2$"34-0.0", 
                   age_df_2$month_of_birth, age_df_2$day_of_birth, sep = "-"))

age_df_2 <- mutate(age_df_2, age_exact = as.numeric(difftime(as.Date(age_df_2$"53-0.0"), 
                   as.Date(age_df_2$DateOfBirth), units = "weeks"))/52.25)

age_df_2 <- age_df_2[,c("eid", "age_exact")]

age_df <- merge(age_df, age_df_2, by = "eid")
age_df <- rename(age_df, sex = '31-0.0', age = '21003-0.0')

print("STEP 2: Age and sex selected")



# Step 3: obtain quality control data and create covariates matrix
qc_data <- data.table(fread(snakemake@params[["file_qc"]], header = F, select = c(3:68)))

colnames(qc_data) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", 
                       "dna_concentration", "submitted_gender", 
                        "inferred_gender", "X_intensity", "Y_intensity", 
						"submitted_plate", "submitted_well", 
						"missing_rate", "heterozygosity", "heterozygosity_pc_corrected", 
						"heterozygosity_missing_outlier", "PSCA", "in_kinship", 
						"excluded_kinship_inference", "excess_relatives", "white_british", 
						"pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", 
						"phasing_X", "phasing_Y")

eids <- data.table(fread(snakemake@params[["file_eids"]], 
                   header = F, select = c(1), 
				   col.names = c("eid")))

qc_data <- cbind(qc_data, eids)

print("STEP 3: Quality control file loaded and merged with eids")

# step 4: filter for active participants
participants <- read.csv(snakemake@params[["file_wb_participants"]])
# correct number in participants

qc_data <- qc_data[qc_data$eid %in% participants$eid,]
# correct number in qc data

print(paste0("Covariates file created, number of individuals = ", 
             nrow(qc_data)))

age_df$age_exact2 <- age_df$age_exact^2
age_df$age_exact_x_sex <- age_df$age_exact * age_df$sex

X <- data.matrix(merge(age_df, 
                       qc_data[,c("eid", paste0("PC", seq(1:40)))], 
					   by = "eid"))

print("Covariate matrix X is created!")
print(paste0("Number of individuals in covariate matrix = ", nrow(X)))
print(colnames(X))

save(X, file = snakemake@output[["file_covariates"]])

print("Covariance matrix X is saved!")


## Obtain medication parameters


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


file_medication_fem <- find_file_function("6153-0.") # medication file for females
file_medication_male <- find_file_function("6177-0.") # same file, different columns

medication_df <- fread(file_medication_fem[1], 
	       	  		   select = c(1, as.numeric(file_medication_fem[2]), 
				   	   as.numeric(file_medication_male[2])))

# check out if NA in 0.0 by definition mans NA in 0.1-0.x
NAin6153.0.x <- medication_df[is.na(medication_df$'6153-0.0'), "eid"] #these should all be from males

NAin6177.0.x <- medication_df[is.na(medication_df$'6177-0.0'), "eid"] # these should all be from females


female_NA_index <- which(medication_df$'6153-0.0' %in% c("-3", "-1")) # 
#-1 = do not know, -3 = prefer not to answer
#note: if participant chose any of these two, multiple answers were not possible
male_NA_index <- which(medication_df$'6177-0.0' %in% c("-3", "-1"))
NA_index <- c(female_NA_index, male_NA_index) 

#Thus, these are the "true" NAs (i.e. people did not (want to) provide information)

#Create meaningful columns out of this information
medication_df$chol_inhibit <- 0
chol_inhibit_male <- which(medication_df$'6177-0.0' == 1)
chol_inhibit_fem <- which(medication_df$'6153-0.0' == 1)
index_chol_inhibit <- c(chol_inhibit_male, chol_inhibit_fem)
medication_df$chol_inhibit[index_chol_inhibit] <- 1
medication_df$chol_inhibit[NA_index] <- NA

medication_df$bp_regul <- 0
bp_regul_male <- which(medication_df$'6177-0.0' == 2)
bp_regul_female <- which(medication_df$'6153-0.0' == 2)
index_bp_regul <- c(bp_regul_male, bp_regul_female)
medication_df$bp_regul[index_bp_regul] <- 1
medication_df$bp_regul[NA_index] <- NA

medication_df$Insulin <- 0
insulin_male <- which(medication_df$'6177-0.0' == 3)
insulin_fem <- which(medication_df$'6153-0.0' == 3)
index_insulin <- c(insulin_male, insulin_fem)
medication_df$Insulin[index_insulin] <- 1
medication_df$Insulin[NA_index] <- NA

# Females only (will have to decide if this needs to be included based on its effect in the female only sample). 
# as of now, men are treated as 0, not NA

medication_df$hormone_repl <- 0
hormone_repl_fem <- which(medication_df$'6153-0.0' == 4)
medication_df$hormone_repl[hormone_repl_fem] <- 1
medication_df$hormone_repl[female_NA_index] <- NA

medication_df$oral_contraceptive <- 0
oral_contra_fem <- which(medication_df$'6153-0.0' == 5)
medication_df$oral_contraceptive[oral_contra_fem] <- 1
medication_df$oral_contraceptive[female_NA_index] <- NA

rawcols <- c(grep("6177", colnames(medication_df), value = T), grep("6153", colnames(medication_df), value = T))
medication_df[, rawcols] <- list(NULL)

print(paste0("Initial number of individuals = ", nrow(medication_df)))

medication_df <- medication_df[medication_df$eid %in% X[,"eid"]]

print(paste0("Number of individuals after filtering as in covariates = ", nrow(medication_df)))

save(medication_df, file = snakemake@output[["file_medication"]])

print("Medication data created and saved!")



