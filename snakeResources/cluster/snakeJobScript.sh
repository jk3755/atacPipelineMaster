#!/bin/bash
#$ -N TEST
#$ -cwd
#$ -pe smp 1
#$ -l mem=8G,time=8::

module load conda
source activate atac