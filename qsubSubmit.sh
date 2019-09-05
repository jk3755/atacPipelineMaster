#!/bin/bash -x
#$ -N TEST
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac
#$ -pe smp 25
#$ -l mem=50G,time=96::
echo "Running qsubSubmit.sh script"
#
echo "module load conda"
module load conda
#
echo "source activate atac"
source activate atac
#
#echo "snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 999 rawFP --cluster-config qsubConfig.json --cluster 'qsub -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0'-V" --latency-wait 60
# always unlock first
echo "unlocking"
snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 999 rawFP_lncap --cluster-config qsubConfig.json --cluster "qsub -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V" --use-conda --conda-prefix /ifs/scratch/c2b2/ac_lab/jk3755/atac/conda --restart-times 12 --latency-wait 60 --rerun-incomplete --unlock
echo "unlocked"
snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 999 rawFP_lncap --cluster-config qsubConfig.json --cluster "qsub -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V" --use-conda --conda-prefix /ifs/scratch/c2b2/ac_lab/jk3755/atac/conda --restart-times 12 --latency-wait 60 --rerun-incomplete