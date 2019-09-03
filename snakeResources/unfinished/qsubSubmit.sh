#!/bin/bash -x
#$ -N TEST
#$ -cwd
#$ -pe smp 1
#$ -l mem=4G,time=15::
#$ -V
echo "Running qsubSubmit.sh script"
#
#echo "module load conda"
#module load conda
#
#echo "source activate atac"
#source activate atac
#
snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 1000 rawFP_lncap --cluster-config snakeResources/cluster/qsubConfig.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V" --use-conda --restart-times 12 --latency-wait 30 --verbose
