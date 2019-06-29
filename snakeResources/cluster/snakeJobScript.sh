#!/bin/bash -x
#$ -N TEST
#$ -cwd
#$ -pe smp 1
#$ -l mem=8G,time=8::
#$ -V
#
echo "Using snakeJobScript.sh to submit job"