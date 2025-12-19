#!/bin/bash

# I.Tarter 10/21/2025

#SBATCH --time 05:00:00                             #Time for the job to run
#SBATCH --job-name=Tyler-Escape
#SBATCH --nodes=1				#To allocate all the requested nodes on the same compute node ( each node has 128 cores)
#SBATCH --ntasks=24                                   #Number of cores needed for the job
#SBATCH --partition=normal                      #Name of the GPU queue
#SBATCH -e slurm-%j.err         # Error file for this job.
#SBATCH -o slurm-%j.out         # Output file for this job.
#SBATCH --mem=256gb                              # RAM
#SBATCH --account=coa_jdbr310_uksr               #Name of account to run under
#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user ilta223@uky.edu              #Where to send email

#Module needed for this Python job
module load ccs/Miniforge3
source activate

conda activate /project/jdbr310_uksr/ilta223/my_conda/Tyler-escape

cd /scratch/ilta223/Tyler-ZIKV-Escape

snakemake --cores 24 -s Snakefile.txt
