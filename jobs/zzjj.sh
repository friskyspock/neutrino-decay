#!/bin/bash
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=25
#SBATCH --mem=10000M
#SBATCH --time=01-18:00
#SBATCH --requeue
#SBATCH --mail-user=kalp.shah@research.iiit.ac.in
#SBATCH --mail-type=ALL
python ~/neutrinoDecay/codes/multi_run/zzjj.py >> ~/neutrinoDecay/logs/zzjj.txt
