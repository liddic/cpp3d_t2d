#!/bin/bash -e

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=1-0
#SBATCH --mem=8000M
#SBATCH --cpus-per-task=2
python ft2d-CHN_2_clean_fastp.py
