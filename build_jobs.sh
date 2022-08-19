#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=1:00:00
#SBATCH --account=rrg-najmanov
source /home/thomasd/ENV/bin/activate

python3 /home/thomasd/projects/rrg-najmanov/thomasd/New_binding_software/FastAID_Py/build_jobs_all_target.py False
