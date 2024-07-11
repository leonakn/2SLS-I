# 2SLS-I

Welcome to the GitHub page of 2SLS-I, a method to investigate how environmental variables modify causal effects.

!Please note, 2SLS-I has not yet been peer reviewed and this repository is under development!

# Usage
This Github contains the workflow that was used to investigate how environmental variables may modify causal effects. The code presented here has been tested on Linux combined with R 4.2.1 and Snakemake 7.25.3.
The workflow is organized as a snakemake pipeline. This allows for a reproducible analysis and parallel computations - which is especially useful since the analysis has been conducted on range of exposure- environment- and outcome traits. However, individual scripts can also be run without the snakemake workflow manager by replacing snakemake variables by hard-coded input and output paths.

# Input data
The application of 2SLS-I was performed on data of the UK Biobank.
Summary statistics were obtained from Neale et al. (2017), and extended Genetic risk scores (here referred to as "PRS"), were obtained using weights and methods from Privé et al. (2022).

# Simulations
Scripts to replicate the simulation analysis are provided in the folder "R/simulations".

# Application
Please note that the full pipeline is provided, even though it is likely that some of the preprocessing steps for the phenotypes will require adequate attention depending on the structure of the cluster that is used.
To account for this, simulated data is provided to run the second half of the analysis, assuming the data has been prepared accordingly. Thus, all elements that would require access to the UKBB to run have been commented out. The working example takes less than 10 minutes to run.

