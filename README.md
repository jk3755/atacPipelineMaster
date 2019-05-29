atacPipelineMaster

## Snakemake execution guide
# A dry run of the pipeline can be run with:
# snakemake -np h508go
# On one of the virtualization servers, run the pipeline with the following to allocate 20 threads and 90 gb max memory (to avoid crashing the process)
# snakemake -j 20 h508go --resources hg38align=1 fp_by_chr=5 raw_fp_graph=2 parse_fp=2 make_bigwig=1
#
## To initiate the pipeline, you must create a directory called "preprocessing" with a subdirectory called "1gz" containing the fastq.gz files
## Rename the files as described below before spooling the pipeline
# 
## Raw file info
# H508-1_S3_L001_R1_001.fastq.gz - Sample 1
# H508-2_S2_L001_R1_001.fastq.gz - Sample 2
# H508-3_S1_L001_R1_001.fastq.gz - Sample 3
## Before running pipeline, if you have three replicates, rename these to:
# H508-WT-01_REP1of3_L1_R1.fastq.gz
# H508-WT-01_REP2of3_L1_R1.fastq.gz
# H508-WT-01_REP3of3_L1_R1.fastq.gz
## If you only have two reps, rename files to:
# H508-WT-01_REP1of2_L1_R1.fastq.gz
# H508-WT-01_REP2of2_L1_R1.fastq.gz
## If you only have one replicate, rename files to:
# H508-WT-01_REP1of1_L1_R1.fastq.gz
#
# Note that additional rule definitions (large group aggregators for panTF, scanPWM) are defined and imported into this main
# Script from the files located in snakeModules directory

## Note that even though this will be sped up by making 20 redundant copies of the bam file,
## There is still a chance two processes will access the same file the way it is currently written
## This will happen if two processes are launched with the same hard coded bam file
## Note sure how to fix this, its probably fine for now
## This code needs some work. Something is tripping it up if I try to run all TFs at once (gets stuck),
## And I have also not been able to enforce strict group ordering in the execution
## For now, I can run each group sequencially by using the shell command:
## for i in {1..62}; do snakemake -j 20 run_group$i; done

## Potential observation when writing/testing this block of code:
## If I put all the TF targets into 62 target rule groups of 20 each,
## And then attempt to run the pipeline by pulling an aggregator tool
## That collects all 62 groups at once, it doesn't crash but stalls 
## and does not run. This may be because the pipeline is pulling target
## TFs from all 62 groups at once, so the entire cohort is available
## to start new processes as soon as one finishes. What this means is,
## FP targets that have very little computational requirements will finish
## Quickly and then that thread will move on to a new target - until it reaches
## One that has a heavy memory/computational load. All threads will do this until
## Eventually all 20 processes are stuck on targets that have serious comp. requirements
## And the pipeline will stall.
## If, alternatively, you run the pipeline so that each group must finish completely before
## the next one starts, this will not be a problem, as all the processes will sync up at
## Each step and wait for the heavier ones to finish.

## Footprinting:
## Run with a terminal command like: for i in {1..62}; do snakemake --config group=$i -j 20 run_pantf_ls1034wt01; done
