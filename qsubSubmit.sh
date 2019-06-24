#!/bin/bash
#$ -N TEST
#$ -cwd
#$ -pe smp 1
#$ -l mem=8G,time=8::

snakemake -j 100 full_analysis_test --latency-wait=30 --cluster-config snakeResources/cluster/qsubConfig.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0"