#!/bin/sh
#SBATCH --job-name=kallisto_new.sh
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH -t 6:00:00
#SBATCH -o /home/eyscott/scratch/FAIZ_2/FAIZ/run_fixer.sh.out
#SBATCH -e /home/eyscott/scratch/FAIZ_2/FAIZ/kalboom.sh.err
#SBATCH --mem 100000

module load StdEnv/2020
module load  gcc/9.3.0
module load intel/2020.1.217

for f in $(ls *R1.fastq | awk -F'[-_.]' '{print $1 "_" $2}');
do /home/eyscott/scratch/FAIZ_2/FAIZ/kallisto-0.48.0/build/src/kallisto bus \
--index /home/eyscott/scratch/FAIZ_2/FAIZ/Mouse/transcriptome.idx -o /home/eyscott/scratch/FAIZ_2/FAIZ/${f}_new --technology=BDWTA --threads=16 --unstranded ${f}_R1.fastq ${f}_R2.fastq \
-g /home/eyscott/scratch/FAIZ_2/FAIZ/Mouse/Mus_musculus.GRCm38.96.gtf; done
