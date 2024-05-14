#!/bin/bash
#
#SBATCH --job-name=rewild
# Set number of iteration
# Amount of RAM requested per job
#SBATCH --mem=78G
# Nb of threads requested per job (smp = shared memory)
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
#SBATCH --time=96:00:00

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/users/${USER}/rewilding2023
DIR_SCRIPT=script

# Load modules
#module load apps/julia/1.8.5/binary

# Put 1.8.5
JULIA="/users/bi1ahd/.juliaup/bin/julialauncher +release"

#cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/00_generate_foodweb.jl
cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/01_connectance_richness.jl
#cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/02_Z_predator.jl
#cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/03_propagule_size.jl
#cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/04_generality.jl
