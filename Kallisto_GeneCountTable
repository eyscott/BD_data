#!/bin/sh
#SBATCH --job-name=bustool_correct.sh
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH -t 1:00:00
#SBATCH -o /home/eyscott/scratch/FAIZ_2/FAIZ/Faiz_olig1_new/bus_roll.sh.out
#SBATCH -e /home/eyscott/scratch/FAIZ_2/FAIZ/Faiz_olig1_new/bus_crash.sh.err
#SBATCH --mem 100000

module load StdEnv/2020
module load  gcc/9.3.0
module load intel/2020.1.217

/home/eyscott/scratch/FAIZ_2/FAIZ/bustools/build/src/bustools sort -o sorted_out.bus output.bus
/home/eyscott/scratch/FAIZ_2/FAIZ/bustools/build/src/bustools correct -o /home/eyscott/scratch/FAIZ_2/FAIZ/corr_out.bus \
--whitelist /home/eyscott/scratch/FAIZ_2/FAIZ/BD_WTA_whitelist_final.txt sorted_out.bus
/home/eyscott/scratch/FAIZ_2/FAIZ/bustools/build/src/bustools count --genecounts -g /home/eyscott/scratch/FAIZ_2/FAIZ/Mouse/transcripts_to_genes.txt \
-t transcripts.txt -e matrix.ec -o /home/eyscott/scratch/FAIZ_2/FAIZ/out_counts corr_out.bus
