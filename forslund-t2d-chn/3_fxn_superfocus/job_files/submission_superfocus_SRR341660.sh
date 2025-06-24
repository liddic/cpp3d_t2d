#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4-0
#SBATCH --mem=64000M
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
export TMPDIR=/cluster/jobs/lidd0026/temp/SRR341660
superfocus -q /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_2_clean_fastp/SRR341660_R1.good.fastq -dir /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_3_fxn_superfocus/superfocus_out_SRR341660 -a diamond -db DB_100 --alternate_directory /home/lidd0026/miniconda3/lib/python3.8/site-packages/superfocus_app
