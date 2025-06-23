#!/bin/bash -e
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=16000
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
fastqc -o /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1b_qc -t 4 /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw/SRR413618_pass_2.fastq
