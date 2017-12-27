#!/bin/bash
#SBATCH -J trainA
#SBATCH -p normal_q
#SBATCH -N 1
#SBATCH -t 20:10:00
#SBATCH --mem=200G
#SBATCH --gres=gpu:pascal:2

echo "Allocated GPU with ID $CUDA_VISIBLE_DEVICES"
# echo "Activate virtual environment: "

# source /work/newriver/gustavo1/deepLearning/bitPredEnv/bin/activate
# export PYTHONNOUSERSITE=True
source activate gustavo1
# module load gcc cuda theano 

cd /work/newriver/gustavo1/deepLearning/deepGene/deepGene/

python DLModel1.py
