#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000
#SBATCH --time=02:00:00
#SBATCH --job-name=bamdownsample
#SBATCH --mail-user=stavrosmakrodi
#SBATCH --mail-type=ALL

ml use /opt/insy/modules/DBL/modulefiles/
ml load samtools
ml load bedtools/2.27.0

python sonia_downsampling.py
