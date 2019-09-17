##################################################################
#### Important notes when running on SGE cluster #################
##################################################################
For this to work properly, you need to give snakemake some way to check the job status
If not, it will submit the original number of jobs specified in -j, but then get stuck
and fail to submit new jobs, as it seems to be waiting for the original jobs to end


## NOTES ##
## -x
## -N
## -wd
## -pe
## -l
## -j
##


##################################################################
#### Description of files and directories ########################
##################################################################

cluster - 
config -
envs -
modules -
scripts -


##################################################################
#### Description of pipeline execution on a local environment ####
##################################################################

# snakemake
# --snakefile snakefileATACseqWorkflow.snakefile
# -j 20 
# [rule]
# --resources mem_mb=50000
# --use-conda
# --restart-times 3

# --snakefile: specify the file where the snakemake workflow is contained
# j: specifies the number of threads the run will use
# restart-times: sets the number of times snakemake will attempt to restart a failed job
# --rerun-incomplete: can be used to regenerate incomple files
# --unlock: unlocks the snakemake working directory, such as after a power loss or kill signal
#
## Resource definitions
# mem_mb: specifies the total memory limit of the pipeline. This relies on mem_mb being user defined in individual rules
#
## Cluster configuration
# When running snakemake workflows on a cluster, initiate an interactivate session from which you will run the pipeline
# Note, limits seem to be weird on interactive, cant seem to initiate more than 6 hours?
# qlogin -l mem=4G,time=8:0:0
# make sure to cd to the appropriate directory
# for example, /ifs/scratch/c2b2/ac_lab/jk3755/atac/
# Activate the appropriate conda environment
# module load conda
# source activate atac
# you can then use a command similar to this to execute the workflow from the interactive session
# snakemake --snakefile snakefileATACseqWorkflow.snakefile -j 1000 rawFP_lncap --cluster-config qsubConfig.json
# --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time={cluster.runtime}:0:0 -V" --use-conda --restart-times 5 --latency-wait 60
#
## Options that may be useful when running
# --list-conda-envs - lists conda envs and their locations
# --cleanup-conda - cleans up unused conda envs
# --create-envs-only - allows you to run the pipeline but only create the conda envs necessary (first time runs)
# --verbose - set verbose output
# --debug-dag - print the dag workflow
# --printshellcmds - print shell commands that will be executed
# --reason - print the reason for execution of each rule

# Spool the pipeline with the following parameters:
# snakemake -j 20 [rule] --resources hg38align=1 rawFPanalysisLarge=1 purgeduplicates=10 mem_mb=95000 --restart-times=3
#
# Parameters:
# j: specifies the number of threads the run will use
# hg38align:
# rawFPanalysisLarge:
# purgeDuplicates:
# mem_mb: specifies the global memory limit the snakemake run can use
# restart-times: sets the number of times snakemake will attempt to restart a failed job

##################################################################
#### Description pipeline execution on an SGE cluster ############
##################################################################






##################################################################
####  Description of rules and options/parameters ################
##################################################################

## Sample correlation
# This analysis relies on the deepTools package: https://deeptools.readthedocs.io/en/develop/index.html
#
# parameters:
# -b input bam files
# -o output file name
# -bs set the bin size used for comparison, default is 10000 bp
# -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
# -p set the number of computing processors to use
# -v verbose mode

rule CORRELATION_spearman_8samples:
    # parameters:
    # -b input bam files
    # -o output file name
    # -bs set the bin size used for comparison, default is 10000 bp
    # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
    # -p set the number of computing processors to use
    # -v verbose mode

Build dirs
## -p ignore error if existing
## -v verbose

## 
rule STEP1_gunzip:
    # -k keep original files
    # -c write to standard output

    ## Fastp fastq QC Filtering
rule STEP2_fastp_filtering:
    # -i, -I specify paired end input
    # -o, -O specifies paired end output
    # -w specifies the number of threads to use (16 max)
    # -h specifies output for the html QC report
    # -j specifies the json output QC report

    rule STEP3_mycoalign:
    # -q fastq input file format
    # -p num threads to use
    # -X1000 align to a maximum of 2000 bp frag length
    # -1 is read 1 input fastq file
    # -2 is read 2 input fastq file
    # -S output file path
    # 2> bowtie2 outputs alignment metrics to STDERR, 2> will allow redirect to a text file

rule STEP4_hg38align:
    # use 'snakemake --resources hg38align=1' to limit the number of parallel instances of this rule
    # -q fastq input file format
    # -p num threads to use
    # -X1000 align to a maximum of 2000 bp frag length
    # -1 is read 1 input fastq file
    # -2 is read 2 input fastq file
    # -S output file path
    # 2> bowtie2 outputs alignment metrics to STDERR, 2> will allow redirect to a text file

 rule STEP5_coordsort_sam:
    # -o output file path
    # -O output file format

rule STEP6_blacklistfilter_bamconversion:
    # -b output in bam format
    # -h include header in output file
    # -o specify output file path
    # -L only output alignments that overlap with the provided BED file
    # -U write the alignments NOT selected by other parameters to the specified file
    # -@ specify number of threads

rule STEP7_chrM_contamination:
    # remove mitochondrial reads
    # params:
    # -b input file is in bam format
    # -h keep the sam header. important downstream
    # -o output filepath for reads NOT matching to blacklist region
    # -L path to the blacklist BED file
    # -U output filepath for reads matching blacklist region
    # -@ number of threads to use

rule STEP8_addrgandcsbam:
    # refer to https://software.broadinstitute.org/gatk/documentation/article.php?id=6472 for information on read group tags
    # note - proper specification of RG tags is critical for downstream analysis and unique sample identification when submitting for publication
    # specification of the lane allows optical duplicates to be detected(?)
    # required by GATK standards
    # also important for identifying batch effects/technical artifacts(?)
    # see: https://software.broadinstitute.org/gatk/documentation/article.php?id=6472
    # Required @RG parameter specifications:
    # RGID (read group ID) - this must be a globally unique string. for illumina data, use flowcell + lane
    # RGLB (read group library) - This is used by MarkDuplicates to collect reads from the same library on different lanes, so it must be common to all files from the same library
    # RGPL (read group platform) - ILLUMINA
    # RGPU (read group platform unit) - The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
    # RGSM (read group sample name) - the name of the sample sequenced in this file. should be consistent across different files from different lanes
    # params:
    # I specifies the input file
    # O specifies the output file

rule STEP9_cleansam:
    # soft-clips bases aligned past the end of the ref sequence
    # soft-clipping retains the bases in the SEQ string, but they are not displayed or used in downstream data analysis
    # sets MAPQ score to 0 for unmapped reads
    # I specifies the input file
    # O specifies the output file

rule STEP10_mergelanes:
    # Merge files for individual lanes
    # I specifies input files for each lane
    # O specifies the output files
    # SORT_ORDER/ASSUME_SORTED specify the type of sorting in the input files
    # MERGE_SEQUENCE_DICTIONARIES will combine the sequence dictionaries from the individual files
    # a sequence dictionary contains information about sequence name, length, genome assembly ID, etc
    # USE_THREADING allows multithreadded operation for the bam compression

rule STEP10b_clean_intermediate_data:
    # -rm removes files
    # -f forces the removal

rule STEP11_purgeduplicates:
	# -Xmx50g specifies a 50 gb memory limit per process
    # I specifies the input file
    # O specifies the output file
    # M specifies the duplication metrics output file
    # REMOVE_DUPLICATES enables removal of duplicate reads from the output file
    # ASSUME_SORTED indicates the input file is already sorted

rule STEP12_mapqfilter:
    # Remove multimapping reads
    # for an explanation of how bowtie2 calculates mapq scores:
    # http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
    # for bowtie2, mapq higher than 2 is a uniquely mapped read
    # params:
    # -h include the header in the output
    # -q only include reads with mapping quality X or higher
    # -b output as a bam file

rule STEP13_build_index:
    # creates a bai index for the bam files
    # this is required for many downstream operations
    # the bai index allows other processes to access specific reads in the bam file without having to read through the entire bam contents to find them (its like a table of contents)
    # I specifies the input bam file
    # O specifies the output index file

rule STEP14_makebigwig_bamcov:
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)

rule STEP15_MACS2_peaks_global_normilization:
    # notes:
    # because we are going to use the TCGA data downstream likely as a reference point,
    # we will need to call the peaks in the exact same way as they did in this paper:
    # http://science.sciencemag.org/content/sci/suppl/2018/10/24/362.6413.eaav1898.DC1/aav1898_Corces_SM.pdf
    # which is "macs2 callpeak --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
    # params:
    # -t input bam file (treatment)
    # -n base name for output files
    # --outdir output directory
    # --shift find all tags in the bam, and shift them by 75 bp
    # --extsize extend all shifted tags by 150 bp (should be roughly equal to avg frag size in lib)
    # --nomodel do not use the macs2 function to determine shifting model
    # --call-summits call the peak summits, detect subpeaks within a peaks
    # --nolambda do not use local bias correction, use background nolambda
    # --keep-dup all keep all duplicate reads (bam should be purged of PCR duplicates at this point)
    # -p set the p-value cutoff for peak calling

rule STEP16_MACS2_peaks_local_normalization:
    # call peaks with MACS2 local normalization (+/- 1000 bp) enabled
    # params:
    # -t input bam file (treatment)
    # -n base name for output files
    # --outdir output directory
    # --shift find all tags in the bam, and shift them by 75 bp
    # --extsize extend all shifted tags by 150 bp (should be roughly equal to avg frag size in lib)
    # --nomodel do not use the macs2 function to determine shifting model
    # --call-summits call the peak summits, detect subpeaks within a peaks
    # --keep-dup all keep all duplicate reads (bam should be purged of PCR duplicates at this point)
    # -p set the p-value cutoff for peak calling

rule STEP17a_percent_peak_genome_coverage_globalnorm:
    # returns a fraction value of the basepairs of the genome covered by the merged peak file. multiple by 100 for percentages
    # parameters:
    # --echo output will be at least a three-column bed file
    # --bases-uniq the number of distinct bases from ref covered by overlap bed file
    # --delim change output delimeter from '|' to <delim>, e.g. '\t'

rule STEP17b_percent_peak_genome_coverage_localnorm:
    # returns a fraction value of the basepairs of the genome covered by the merged peak file. multiple by 100 for percentages
    # parameters:
    # --echo output will be at least a three-column bed file
    # --bases-uniq the number of distinct bases from ref covered by overlap bed file
    # --delim change output delimeter from '|' to <delim>, e.g. '\t'

rule PEAKS_differential_peak_calling_2samples:
    # notes:
    # because we are going to use the TCGA data downstream likely as a reference point,
    # we will need to call the peaks in the exact same way as they did in this paper:
    # http://science.sciencemag.org/content/sci/suppl/2018/10/24/362.6413.eaav1898.DC1/aav1898_Corces_SM.pdf
    # which is "macs2 callpeak --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
    # params:
    # -t input bam file (treatment)
    # -n base name for output files
    # --outdir output directory
    # --shift find all tags in the bam, and shift them by 75 bp
    # --extsize extend all shifted tags by 150 bp (should be roughly equal to avg frag size in lib)
    # --nomodel do not use the macs2 function to determine shifting model
    # --call-summits call the peak summits, detect subpeaks within a peaks
    # --nolambda do not use local bias correction, use background nolambda
    # --keep-dup all keep all duplicate reads (bam should be purged of PCR duplicates at this point)
    # -p set the p-value cutoff for peak calling