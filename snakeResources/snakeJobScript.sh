#!/bin/bash -x
#$ -N TEST
#$ -cwd
#$ -pe smp 1
#$ -l mem=8G,time=8::
#$ -V
echo "Using snakeJobScript.sh to submit job"
echo "module load conda"
module load conda
echo "source activate atac"
source activate atac