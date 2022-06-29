#!/bin/bash -l
#SBATCH -A snic2022-5-34
#SBATCH -p node -n 5
#SBATCH -J juvireID
#SBATCH -t 200:00:00
#SBATCH --mail-user xuejing.hu.7095@student.uu.se
#SBATCH --mail-type=ALL
conda activate snakemake
snakemake -s juvireID.py --cores 100  --use-conda --resources mem_gb=256 --rerun-incomplete
	
