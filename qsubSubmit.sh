#!/bin/bash -x
#$ -N TEST
#$ -cwd
#$ -pe smp 1
#$ -l mem=4G,time=15::
echo "Running qsubSubmit.sh script"
#
echo "module load conda"
module load conda
#
echo "source activate atac"
source activate atac
#
echo "snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 100 rawFP_lncap_ex01 --cluster-config snakeResources/cluster/qsubConfig.json --cluster 'qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0'-V" --latency-wait 60
snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 100 rawFP_lncap_ex01 --cluster-config snakeResources/cluster/qsubConfig.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V" --latency-wait 60