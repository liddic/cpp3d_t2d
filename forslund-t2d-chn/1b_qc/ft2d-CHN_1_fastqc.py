# # # # # # # # # # # # # # # #
# # FastQC (without MultiQC) - for Forslund et al 2015 CHN selected metagenomics
# # # # # # # # # # # # # # # #
#$ sinfo -N -o "%N %C %e %m"   # get cluster usage info
#$ cd <check run folder!!>

import os
import sys
import time
import pandas as pd


print(sys.path)

# directory for raw fastq files
READDIR = '/cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw'

# directory for QC outputs
OUTDIR = '/cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1b_qc'

# set working directory
workDir = OUTDIR
os.chdir(workDir)
print("Current Working Directory: ",os.getcwd())

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


job_directory = "%s/job_files" % os.getcwd()
mkdir_p(job_directory)

## read in manifest of sampleID and corresponding raw fastq R1/R2 files
table = pd.read_csv('/cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw/ft2d-CHN-samples.csv', sep=',')

# extract info from pandas data table
r1_files = table["R1_filenames"]
r2_files = table["R2_filenames"]



# iterate through creating jobs for fastp QC of each sample # len(r1_files) = 285
n1 = len(r1_files)
n2 = len(r2_files)

for i in range(n1):
    #i=0
    job_file = os.path.join(job_directory, "slurm_fastqc_r1_%s.sh" % i)
    
    infile1 = os.path.join(READDIR,"%s" % r1_files[i] )
        
    #with open(job_file) as fh:
    fh = open(job_file, "w")
    fh.writelines("#!/bin/bash -e\n")
    fh.writelines("#SBATCH --ntasks=1\n")
    fh.writelines("#SBATCH --time=1-0\n")
    fh.writelines("#SBATCH --mem=16000\n")
    fh.writelines("#SBATCH --cpus-per-task=4\n")
    fh.writelines("#SBATCH --nodes=1\n")
    
    fh.writelines("fastqc -o %s -t 4 %s\n" % (OUTDIR,infile1))
    
    fh.close()
    time.sleep(2)
    os.system("sbatch %s" % job_file)


for i in range(n2):
    #i=0
    job_file = os.path.join(job_directory, "slurm_fastqc_r2_%s.sh" % i)
    
    infile2 = os.path.join(READDIR,"%s" % r2_files[i] )
        
    #with open(job_file) as fh:
    fh = open(job_file, "w")
    fh.writelines("#!/bin/bash -e\n")
    fh.writelines("#SBATCH --ntasks=1\n")
    fh.writelines("#SBATCH --time=1-0\n")
    fh.writelines("#SBATCH --mem=16000\n")
    fh.writelines("#SBATCH --cpus-per-task=4\n")
    fh.writelines("#SBATCH --nodes=1\n")
    
    fh.writelines("fastqc -o %s -t 4 %s\n" % (OUTDIR,infile2))
    
    fh.close()
    time.sleep(2)
    os.system("sbatch %s" % job_file)

## END
