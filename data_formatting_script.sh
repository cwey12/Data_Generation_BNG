#!/bin/bash
#SBATCH -J Hello_World # Job Name
#SBATCH -o out.txt # Output log name
#SBATCH -n 16 # Number of processors 
#SBATCH -N 1  # Number of nodes
#SBATCH -t 24:00:00 # Run time (hh:mm:ss)
#SBATCH -A cwey  # Account to "charge" credits to

##### Constants #####
JOBID=$SLURM_JOBID

# Load anaconda module if not installed on home directory
module load anaconda/2020.07

#Define variable arguments for combining data script
lower_time_limit=0
upper_time_limit=6
num_timepoints=10
gdat_dir='/data/homezvol2/cwey/Working_Dir/Data/22-10-2020_20:09:10/gdat'

# Command to submit job here
python Combine_Data_for_Sincerities.py $lower_time_limit $upper_time_limit $num_timepoints $gdat_dir 

# Sucessful finishing note
echo This will write to the log when the job is finished
