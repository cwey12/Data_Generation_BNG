#!/bin/bash
#SBATCH -J Hello_World # Job Name
#SBATCH -o out.txt # Output log name
#SBATCH -n 16 # Number of processors 
#SBATCH -N 1  # Number of nodes
#SBATCH -t 24:00:00 # Run time (hh:mm:ss)
#SBATCH -A elread_lab  # Account to "charge" credits to

##### Constants #####
JOBID=$SLURM_JOBID

# Load anaconda module if not installed on home directory
module load anaconda/2020.07

#Define variables to be used in the script
num_cells=1000

# Command to submit job here
python Run_Simulation.py $num_cells 

# Sucessful finishing note
echo This will write to the log when the job is finished
