#!/bin/bash -x
#$ -N snakemakeATACseqFP
#$ -j y
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac
#$ -pe smp 4
#$ -l mem=40G,time=168::
##
## NOTES ##
## -x
## -N
## -wd
## -pe
## -l
## -j
##
## Load the conda environment for running snakemake on the master job
echo "Running qsubSubmit.sh script"
echo "Loading conda module"
module load conda
echo "Activating atac conda env"
source activate atac
echo "Conda env activated"
##
## Set up the desired variables for running the job
echo "Setting up variables"
SNAKEFILE="../../snakefileATACseqWorkflow.snakefile"
SNAKEJOB="full_lncap"
CORES="2000"
JOBRESTARTS="5"
LATENCYWAIT="5"
CLUSTCONFIG="qsubConfig.json"
WORKDIR="/ifs/scratch/c2b2/ac_lab/jk3755/atac"
CONDADIR="/ifs/scratch/c2b2/ac_lab/jk3755/atac/conda"
CLUSTSUBMIT="qsub -wd {$WORKDIR} -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V"
##
## Echo variable settings
echo "Echoing current variables"
echo "Snakefile: $SNAKEFILE"
echo "Snakejob: $SNAKEJOB"
echo "Cores: $CORES"
echo "Job restarts: $JOBRESTARTS"
echo "Filesystem latency wait: $LATENCYWAIT"
echo "Cluster config: $CLUSTCONFIG"
echo "Cluster submit: $CLUSTSUBMIT"
echo "Working directory: $WORKDIR"
echo "Conda directory: $CONDADIR"
##
## Always unlock the working directory first and use the --rerun-incomplete flag, in case a previous job crashed
echo "Unlocking snakemake directory"
snakemake \
--snakefile $SNAKEFILE \
--cores $CORES \
$SNAKEJOB \
--cluster-config $CLUSTCONFIG \
--cluster "qsub -j y -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete \
--keep-going \
--unlock
echo "Snakemake directory unlocked"
##
## Run the job
echo "Spooling the snakemake job"
snakemake \
--snakefile $SNAKEFILE \
--cores $CORES \
$SNAKEJOB \
--cluster-config $CLUSTCONFIG \
--cluster "qsub -j y -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete \
--keep-going