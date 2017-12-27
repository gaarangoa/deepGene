#!/bin/bash
#SBATCH -J trainA
#SBATCH -p normal_q
#SBATCH -N 1
#SBATCH -t 20:10:00
#SBATCH --mem=200G
#SBATCH --gres=gpu:pascal:1

echo "Allocated GPU with ID $CUDA_VISIBLE_DEVICES"
# echo "Activate virtual environment: "
export PYTHONNOUSERSITE=True
source /work/newriver/gustavo1/deepLearning/bitPredEnv/bin/activate

cd /work/newriver/gustavo1/deepLearning/deepGene/deepGene/
module load cuda gcc/5.4.0 theano

python DLModel1.py