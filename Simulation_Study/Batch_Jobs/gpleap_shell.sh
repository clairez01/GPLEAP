#!/bin/bash

#SBATCH --job-name=gpleaprun
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH	--array=1-510
#SBATCH --output=/work/users/c/l/clairez1/Paper2/GPLEAP/Log/slurmLogFiles%a.out
#SBATCH --error=/work/users/c/l/clairez1/Paper2/GPLEAP/Error/%a.err
#SBATCH --constraint=rhel8

## add R module
##module add gcc/11.2.0
##module add r/4.3.1
module add gcc/12.2.0
module add r/4.4.0

R CMD BATCH --no-restore /nas/longleaf/home/clairez1/Paper2/Programs/gpleap_sims.R /work/users/c/l/clairez1/Paper2/GPLEAP/sim_$SLURM_ARRAY_TASK_ID.Rout


