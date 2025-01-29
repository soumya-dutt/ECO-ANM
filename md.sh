#!/usr/bin/bash
#SBATCH -N 1
#SBATCH -p general
#SBATCH -q grp_asinghar
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH -t 7-00:00:00
#SBATCH --mem=90000

source ./setup.sh
python3 ANM_highthroughput_pipeline.py > eco_anm.log

