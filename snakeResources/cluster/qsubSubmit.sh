#!/bin/bash -x
#$ -N snakemakeATACseq
#$ -j y
#$ -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac
#$ -pe smp 1
#$ -l mem=20G,time=72::
# Load the conda environment for running snakemake on the master job
echo "Running qsubSubmit.sh script"
echo "Loading conda module"
module load conda
echo "Activating atac conda env"
source activate atac
echo "Conda env activated"
# Set up the desired variables for running the job
echo "Setting up variables"
TARGETRULE="all_lncap"
SNAKEFILE="/ifs/scratch/c2b2/ac_lab/jk3755/atac/snakefileATACseqWorkflow.snakefile"
WORKDIR="/ifs/scratch/c2b2/ac_lab/jk3755/atac"
CONDADIR="conda"
CLUSTCONFIG="/ifs/scratch/c2b2/ac_lab/jk3755/atac/snakeResources/cluster/qsubConfig.json"
CORES="999"
JOBRESTARTS="10"
LATENCYWAIT="60"
# Echo variable settings
echo "Current cluster settings for snakemake run:"
echo "Target rule: $TARGETRULE"
echo "Snakefile: $SNAKEFILE"
echo "Working directory: $WORKDIR"
echo "Conda environment directory: $CONDADIR"
echo "Cluster configuration file: $CLUSTCONFIG"
echo "Provided cluster cores: $CORES"
echo "Job restart attempts: $JOBRESTARTS"
echo "Filesystem latency wait: $LATENCYWAIT"
# Always unlock the working directory first
echo "Unlocking snakemake directory"
snakemake \
--snakefile $SNAKEFILE \
--cores $CORES \
$TARGETRULE \
--cluster-config $CLUSTCONFIG \
--cluster "qsub -terse -j y -o pipeline_log.txt -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete \
--unlock
echo "Snakemake directory unlocked"
# Run the pipeline
echo "Spooling the snakemake pipeline"
snakemake \
--snakefile $SNAKEFILE \
--cores $CORES \
$TARGETRULE \
--cluster-config $CLUSTCONFIG \
--cluster "qsub -terse -j y -o pipeline_log.txt -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -wd /ifs/scratch/c2b2/ac_lab/jk3755/atac -V" \
--use-conda \
--conda-prefix $CONDADIR \
--restart-times $JOBRESTARTS \
--latency-wait $LATENCYWAIT \
--rerun-incomplete
echo "Finished processing!"
