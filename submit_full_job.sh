#!/usr/bin/env bash
#SBATCH --job-name=cdglap
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4096M
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --account=rrg-jeon-ac

echo $OMP_NUM_THREADS

time ./microjets cdglap_input
