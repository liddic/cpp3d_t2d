#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=16000
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
fastp --length_required 50 --average_qual 25 --n_base_limit 1 --dedup \
    --trim_front1 4 --trim_tail1 2 --trim_front2 4 --trim_tail2 2 \
    --cut_mean_quality 30 --cut_window_size 5 --cut_front --cut_tail \
    --out1 /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_2_clean_fastp/SRR413580_R1.good.fastq --unpaired1 /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_2_clean_fastp/SRR413580_R1.single.fastq \
    --out2 /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_2_clean_fastp/SRR413580_R2.good.fastq --unpaired2 /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_2_clean_fastp/SRR413580_R2.single.fastq \
    --in1 /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw/SRR413580_pass_1.fastq --in2 /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw/SRR413580_pass_2.fastq
    --adapter_fasta /cluster/jobs/lidd0026/forslund-t2d-chn/IlluminaAdapters.fa \
    --html  /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_2_clean_fastp/QC_fastq_qual_files/SRR413580_fastp.html --json /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_2_clean_fastp/QC_fastq_qual_files/SRR413580_fastp.json \
