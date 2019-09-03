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
#echo "snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 100 rawFP --cluster-config qsubConfig.json --cluster 'qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0'-V" --latency-wait 60
snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 100 rawFP_lncap --cluster-config qsubConfig.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V" --use-conda --restart-times 5 --latency-wait 60 --verbose