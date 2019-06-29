#!/bin/bash -x
#$ -N TEST
#$ -cwd
#$ -pe smp 1
#$ -l mem=8G,time=8::
echo "Running qsubSubmit.sh script"
#
echo "module load conda"
module load conda
#
echo "source activate atac"
source activate atac
#
echo "snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 1 rawFP_lncap_ex01 --cluster-config snakeResources/cluster/qsubConfig.json --cluster 'qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0'-V"
snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 1 rawFP_lncap_ex01 --cluster-config snakeResources/cluster/qsubConfig.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0 -V"