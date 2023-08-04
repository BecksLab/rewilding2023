#!/bin/bash
#
#SBATCH --job-name=rewild
# Set number of iteration
# Amount of RAM requested per job
#SBATCH --mem=64G
# Nb of threads requested per job (smp = shared memory)
#SBATCH --cpus-per-task=15
#SBATCH --ntasks=1

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/users/${USER}/rewilding2023
DIR_SCRIPT=scripts

# Load modules
#module load apps/julia/1.8.5/binary

# Put 1.8.5
JULIA="/users/bi1ahd/.juliaup/bin/julialauncher +release"

cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/00_generate_foodweb.jl
cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/01_predator_present.jl
