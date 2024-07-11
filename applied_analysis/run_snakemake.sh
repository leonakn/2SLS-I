#!/bin/bash

#SBATCH --job-name='Applied_Interaction_MR'
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition XXX
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=2GB

snakemake --keep-going -j 1000 --cluster-config cluster_config.yaml --cluster "sbatch -p urblauna --job-name='Applied Interaction MR' --time={cluster.time} --nodes=1 --cpus-per-task=1 --mem-per-cpu={cluster.mem} --output={cluster.log}"

