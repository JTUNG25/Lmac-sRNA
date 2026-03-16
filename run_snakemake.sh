#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=shortstack_test
#SBATCH --output=sm.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8
cd /QRISdata/Q9140/lmac/lmac_srna

snakemake -s shortstack.smk --profile profiles/bunya/ --use-singularity