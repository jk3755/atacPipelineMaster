#!/bin/bash
#$ -N TEST
#$ -cwd
#$ -pe smp 1
#$ -l mem=8G,time=8::

module load conda
source activate atac
snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 100 rawFP_lncap_ex01 --latency-wait=30 --cluster-config snakeResources/cluster/qsubConfig.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0"