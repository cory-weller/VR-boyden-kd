#!/usr/bin/env bash
#SBATCH --time 10:00:00
#SBATCH --mem 120G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6

module load python/3

cd /data/CARD_ARDIS/users/wellerca/VR-boyden-kd/APA/APA-Scan

python ./APA-Scan.py