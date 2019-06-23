########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################
## Spool the pipeline with the following parameters:
# snakemake -j 20 [rule] --resources hg38align=1 rawFPanalysisSectored=1 purgeduplicates=10 mem_mb=95000 --restart-times=3
#
## Parameters:
# j: specifies the number of threads the run will use
# restart-times: sets the number of times snakemake will attempt to restart a failed job
# --rerun-incomplete: can be used to regenerate incomple files
# --unlock: unlocks the snakemake working directory, such as after a power loss or kill signal
#
## Resource definitions
# mem_mb: specifies the total memory limit of the pipeline. This relies on mem_mb being user defined in individual rules
# bowtie2align: jobs that use bowtie2 for aligning reads
# purgeDuplicates: rules that use picard to remove duplicate reads
# rawFPanalysis: rules that analyze raw footprinting data
# mergeRawFPSectors: rule that merges to sectored data for the larger footprint analysis
#
#
########################################################################################################################################
#### IMPORT MODULES AND CONFIG #########################################################################################################
########################################################################################################################################
configfile: "snakeResources/config/config.yaml"
#
include: "snakeResources/modules/generateSites.snakefile"
include: "snakeResources/modules/spoolPreprocessing.snakefile"
include: "snakeResources/modules/spoolFootprinting.snakefile"
include: "snakeResources/modules/spoolSampleCorrelation.snakefile"
include: "snakeResources/modules/spoolFullAnalysis.snakefile"

########################################################################################################################################
#### FULL ANALYSIS AGGREGATOR ##########################################################################################################
########################################################################################################################################
# This rule determines what is run in the full analysis spooling option
rule AGGREGATOR_full_analysis:
    input:
        "{path}operations/preprocessing/{sample}-REP{repnum}.preprocessing.complete",
        "{path}operations/footprints/{sample}-REP{repnum}.footprinting_raw_analysis.complete",
        "{path}operations/footprints/{sample}-REP{repnum}.sectored_footprinting_analysis_raw.complete"
    output:
        "{path}operations/modules/{sample}-REP{repnum}.full_analysis.finished"
    shell:
        "touch {output}"

########################################################################################################################################
#### CREATE LOCAL PWM SCAN DATABASE ####################################################################################################
########################################################################################################################################
# Run this rule to generate all needed data for scanning the genome for matches to PWMs
# Will generate data for all annotated genes in motifDB, for all unique motifs
rule run_PWMscan:
    input:
        "snakeResources/sites/operations/PWMscan.allgroups.done"

########################################################################################################################################
#### PREPROCESSING AGGREGATOR ##########################################################################################################
########################################################################################################################################
rule AGGREGATOR_preprocessing:
    input:
        #"snakeResources/sites/operations/PWMscan.allgroups.done",
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai",
        "{path}preprocessing/11bigwig/{sample}-REP{repnum}.bw",
        "{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak",
        "{path}peaks/localnorm/{sample}-REP{repnum}_localnorm_peaks.narrowPeak",
        "{path}metrics/{sample}-REP{repnum}.peak.globalnorm.genomecov.txt",
        "{path}metrics/{sample}-REP{repnum}.peak.localnorm.genomecov.txt",
        "{path}metrics/{sample}-REP{repnum}.fragsizes.svg",
        "{path}operations/preprocessing/{sample}-REP{repnum}.globalpeak.annotations.done",
        "{path}operations/preprocessing/{sample}-REP{repnum}.localpeak.annotations.done",
        "{path}metrics/{sample}-REP{repnum}.totalreads.Rdata",
        "{path}operations/saturation/{sample}-REP{repnum}.saturation_analysis.done"
    output:
        "{path}operations/preprocessing/{sample}-REP{repnum}.preprocessing.complete"
    shell:
        "touch {output}"

########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################
# Build the directory structure
rule PREP_builddirstructure:
    # params: -p ignore error if existing, make parent dirs, -v verbose
    output:
        "{path}operations/preprocessing/dirtree.built"
    shell:
        """
        ####################################################################################################################################################################
        mkdir -p -v {wildcards.path}benchmark
        mkdir -p -v {wildcards.path}benchmark/preprocessing
        mkdir -p -v {wildcards.path}benchmark/correlation
        mkdir -p -v {wildcards.path}benchmark/saturation
        #
        mkdir -p -v {wildcards.path}benchmark/footprints
        mkdir -p -v {wildcards.path}benchmark/footprints/raw {wildcards.path}benchmark/footprints/parsed
        mkdir -p -v {wildcards.path}benchmark/footprints/processed {wildcards.path}benchmark/footprints/merge
        ####################################################################################################################################################################
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}operations/modules
        mkdir -p -v {wildcards.path}operations/preprocessing
        #
        mkdir -p -v {wildcards.path}operations/footprints
        mkdir -p -v {wildcards.path}operations/footprints/raw {wildcards.path}operations/footprints/parsed {wildcards.path}operations/footprints/processed
        #
        mkdir -p -v {wildcards.path}operations/footprints/temp
        mkdir -p -v {wildcards.path}operations/footprints/merged
        #
        mkdir -p -v {wildcards.path}operations/saturation
        mkdir -p -v {wildcards.path}operations/saturation/footprints {wildcards.path}operations/saturation/footprints/raw
        ####################################################################################################################################################################
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}preprocessing/2fastq
        mkdir -p -v {wildcards.path}preprocessing/3goodfastq
        mkdir -p -v {wildcards.path}preprocessing/4mycoalign
        mkdir -p -v {wildcards.path}preprocessing/5hg38align
        #
        mkdir -p -v {wildcards.path}preprocessing/6rawbamm
        mkdir -p -v {wildcards.path}preprocessing/6rawbam/mitochondrial {wildcards.path}preprocessing/6rawbam/blacklist {wildcards.path}preprocessing/6rawbam/nonblacklist
        #
        mkdir -p -v {wildcards.path}preprocessing/7rgsort
        mkdir -p -v {wildcards.path}preprocessing/8merged
        mkdir -p -v {wildcards.path}preprocessing/9dedup
        mkdir -p -v {wildcards.path}preprocessing/10unique
        mkdir -p -v {wildcards.path}preprocessing/11bigwig
        ####################################################################################################################################################################
        mkdir -p -v {wildcards.path}preprocessing/saturation
        mkdir -p -v {wildcards.path}preprocessing/saturation/complexity
        mkdir -p -v {wildcards.path}preprocessing/saturation/peaks 
        mkdir -p -v {wildcards.path}preprocessing/saturation/downsampled
        #
        mkdir -p -v {wildcards.path}preprocessing/saturation/downsampled
        mkdir -p -v {wildcards.path}preprocessing/saturation/downsampled/raw {wildcards.path}preprocessing/saturation/downsampled/cs {wildcards.path}preprocessing/saturation/downsampled/md
        #
        mkdir -p -v {wildcards.path}preprocessing/saturation/footprints
        mkdir -p -v {wildcards.path}preprocessing/saturation/footprints/raw {wildcards.path}preprocessing/saturation/footprints/parsed {wildcards.path}preprocessing/saturation/footprints/processed
        ####################################################################################################################################################################
        mkdir -p -v {wildcards.path}footprints
        #
        mkdir -p -v {wildcards.path}footprints/data 
        mkdir -p -v {wildcards.path}footprints/data/raw {wildcards.path}footprints/data/parsed 
        mkdir -p -v {wildcards.path}footprints/data/processed {wildcards.path}footprints/data/aggregated
        #
        mkdir -p -v {wildcards.path}footprints/data/temp
        #
        mkdir -p -v {wildcards.path}footprints/graphs
        mkdir -p -v {wildcards.path}footprints/graphs/insprob {wildcards.path}footprints/graphs/heatmaps
        ####################################################################################################################################################################
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}peaks/localnorm {wildcards.path}peaks/globalnorm
        ####################################################################################################################################################################
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}metrics/saturation
        ####################################################################################################################################################################
        mkdir -p -v {wildcards.path}correlation
        ####################################################################################################################################################################
        touch {output}
        """

# Gunzip the fastq files
rule STEP1_gunzip:
    # -k keep original files
    # -c write to standard output
    input:
        a="{path}preprocessing/1gz/{sample}-REP{repnum}_L{lane}_R{read}.fastq.gz",
        b="{path}operations/preprocessing/dirtree.built"
    output:
        c="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R{read}.fastq"
    shell:
        "gunzip -k -c {input.a} > {output.c}"

# Fastq QC Filtering
rule STEP2_afterqc_fastqfiltering:
    # -1 specifies read 1 fastq file
    # -2 specifies read 2 fastq file
    # -g specifies the output directory for the good fastq files
    # -b specifies the output directory for the bad fastq files
    # -f -1 autodetects number of bases to trim at front
    # -t -1 autodetects number of bases to trim at tail
    # -s is the shortest trimmed read length allowed past QC filter
    input:
        a="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R1.fastq",
        b="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R2.fastq"
    output:
        c="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        d="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.{lane}.fastqfilter.benchmark.txt'
    shell:
        "after.py -1 {input.a} -2 {input.b} -g {wildcards.path}preprocessing/3goodfastq -b {wildcards.path}preprocessing/3goodfastq -f -1 -t -1 -s 15"
    
# Check for mycoplasma contamination
rule STEP3_mycoalign:
    # -q fastq input file format
    # -p num threads to use
    # -X1000 align to a maximum of 2000 bp frag length
    # -1 is read 1 input fastq file
    # -2 is read 2 input fastq file
    # -S output file path
    # 2> bowtie2 outputs alignment metrics to STDERR, 2> will allow redirect to a text file
    input:
        a="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq"
    output:
        "{path}preprocessing/4mycoalign/{sample}-REP{repnum}_L{lane}.myco.sam"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.{lane}.mycoalign.benchmark.txt'
    threads:
        20
    resources:
        bowtie2align=1,
        mem_mb=5000
    shell:
        "bowtie2 -q -p 20 -X2000 -x genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}_L{wildcards.lane}.myco.alignment.txt"
    
# Align reads to human hg38 build
rule STEP4_hg38align:
    # use 'snakemake --resources hg38align=1' to limit the number of parallel instances of this rule
    # -q fastq input file format
    # -p num threads to use
    # -X1000 align to a maximum of 2000 bp frag length
    # -1 is read 1 input fastq file
    # -2 is read 2 input fastq file
    # -S output file path
    # 2> bowtie2 outputs alignment metrics to STDERR, 2> will allow redirect to a text file
    input:
        a="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq",
        c="{path}preprocessing/4mycoalign/{sample}-REP{repnum}_L{lane}.myco.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.sam"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.{lane}.hg38align.benchmark.txt'
    threads:
        20
    resources:
        bowtie2align=1,
        mem_mb=10000
    shell:
        "bowtie2 -q -p 20 -X2000 -x genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}_L{wildcards.lane}.hg38.alignment.txt"
    
# Coordinate sort the aligned reads. This is required for blacklist filtering
rule STEP5_coordsort_sam:
    # -o output file path
    # -O output file format
    input:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.cs.sam"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.{lane}.coordsort.benchmark.txt'
    shell:
        "samtools sort {input} -o {output} -O sam"
    
# Remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
rule STEP6_blacklistfilter_bamconversion:
    # -b output in bam format
    # -h include header in output file
    # -o specify output file path
    # -L only output alignments that overlap with the provided BED file
    # -U write the alignments NOT selected by other parameters to the specified file
    # -@ specify number of threads
    input:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.cs.sam"
    output:
        a="{path}preprocessing/6rawbam/blacklist/{sample}-REP{repnum}_L{lane}.hg38blacklist.bam",
        b="{path}preprocessing/6rawbam/nonblacklist/{sample}-REP{repnum}_L{lane}.blrm.bam"
    threads:
        20
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.{lane}.bamconvert.benchmark.txt'
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 20 {input}"
    
# Remove reads mapping to mitochondrial DNA
rule STEP7_chrM_contamination:
    # remove mitochondrial reads
    # params:
    # -b input file is in bam format
    # -h keep the sam header. important downstream
    # -o output filepath for reads NOT matching to blacklist region
    # -L path to the blacklist BED file
    # -U output filepath for reads matching blacklist region
    # -@ number of threads to use
    input:
        "{path}preprocessing/6rawbam/nonblacklist/{sample}-REP{repnum}_L{lane}.blrm.bam"
    output:
        a="{path}preprocessing/6rawbam/mitochondrial/{sample}-REP{repnum}_L{lane}.mitochondrial.bam",
        b="{path}preprocessing/6rawbam/{sample}-REP{repnum}_L{lane}.goodbam"
    threads:
        20
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.{lane}.chrMfilter.benchmark.txt'
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 20 {input}"
    
# Add @RG tags to the reads and perform coordinate sorting
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
    input:
        "{path}preprocessing/6rawbam/{sample}-REP{repnum}_L{lane}.goodbam"
    output:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.rg.cs.bam"
    resources:
        mem_mb=50000
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.{lane}.addRGtag.benchmark.txt'
    shell:
        "java -jar snakeResources/programs/picard/picard.jar AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID=H5YHHBGX3.{wildcards.lane} \
        RGLB={wildcards.sample} \
        RGPL=ILLUMINA \
        RGPU=H5YHHBGX3.{wildcards.lane}.{wildcards.sample} \
        RGSM={wildcards.sample}"
    
# Clean the bam file
rule STEP9_cleansam:
    # soft-clips bases aligned past the end of the ref sequence
    # soft-clipping retains the bases in the SEQ string, but they are not displayed or used in downstream data analysis
    # sets MAPQ score to 0 for unmapped reads
    # I specifies the input file
    # O specifies the output file
    input:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.rg.cs.bam"
    output:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.clean.bam"
    resources:
        mem_mb=50000
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.{lane}.cleansam.benchmark.txt'
    shell:
        "java -jar snakeResources/programs/picard/picard.jar CleanSam \
        I={input} \
        O={output}"
    
# Merge reads from different NextSeq lanes
rule STEP10_mergelanes:
    # Merge files for individual lanes
    # I specifies input files for each lane
    # O specifies the output files
    # SORT_ORDER/ASSUME_SORTED specify the type of sorting in the input files
    # MERGE_SEQUENCE_DICTIONARIES will combine the sequence dictionaries from the individual files
    # a sequence dictionary contains information about sequence name, length, genome assembly ID, etc
    # USE_THREADING allows multithreadded operation
    input:
        a="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L1.clean.bam",
        b="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L2.clean.bam",
        c="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L3.clean.bam",
        d="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L4.clean.bam"
    output:
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    resources:
        mem_mb=50000
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.mergelanes.benchmark.txt'
    shell:
        "java -jar snakeResources/programs/picard/picard.jar MergeSamFiles \
        I={input.a} \
        I={input.b} \
        I={input.c} \
        I={input.d} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

# Clean up intermediate data to this point
# Also copy the fastq filtering QC files to the metrics folder
rule STEP10b_clean_intermediate_data:
    # -rm removes files
    # -f forces the removal
    # -a is recursive option for cp
    # /. causes cp to copy all contents of folder including hidden items
    input:
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    output:
        "{path}operations/preprocessing/clean10b.{sample}.{repnum}.done"
    shell:
        """
        rm -f {wildcards.path}preprocessing/2fastq/*REP{wildcards.repnum}*.fastq
        rm -f {wildcards.path}preprocessing/3goodfastq/*REP{wildcards.repnum}*.fq
        rm -f {wildcards.path}preprocessing/4mycoalign/*REP{wildcards.repnum}*.sam
        rm -f {wildcards.path}preprocessing/5hg38align/*REP{wildcards.repnum}*.sam
        rm -f {wildcards.path}preprocessing/6rawbam/*REP{wildcards.repnum}*.goodbam
        rm -f {wildcards.path}preprocessing/6rawbam/blacklist/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/6rawbam/mitochondrial/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/6rawbam/nonblacklist/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/7rgsort/*REP{wildcards.repnum}*.bam
        ##
        cp -a {wildcards.path}preprocessing/QC/. {wildcards.path}metrics/
        touch {output}
        """

# Purge PCR duplicate reads
rule STEP11_purgeduplicates:
	# -Xmx50g specifies a 50 gb memory limit per process
    # I specifies the input file
    # O specifies the output file
    # M specifies the duplication metrics output file
    # REMOVE_DUPLICATES enables removal of duplicate reads from the output file
    # ASSUME_SORTED indicates the input file is already sorted
    input:
        a="{path}operations/preprocessing/clean10b.{sample}.{repnum}.done",
        b="{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    output:
        "{path}preprocessing/9dedup/{sample}-REP{repnum}.dp.bam"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.purgeduplicates.benchmark.txt'
    resources:
        purgeDuplicates=1,
        mem_mb=50000
    shell:
        "java -Xmx50g -jar snakeResources/programs/picard/picard.jar MarkDuplicates \
        I={input.b} \
        O={output} \
        M={wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}.duplication.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"
    
# Filter reads for only uniquely mapping
rule STEP12_mapqfilter:
    # Remove multimapping reads
    # for an explanation of how bowtie2 calculates mapq scores:
    # http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
    # for bowtie2, mapq higher than 2 is a uniquely mapped read
    # params:
    # -h include the header in the output
    # -q only include reads with mapping quality X or higher
    # -b output as a bam file
    input:
        "{path}preprocessing/9dedup/{sample}-REP{repnum}.dp.bam"
    output:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.mapqfilter.benchmark.txt'
    shell:
        "samtools view -h -q 2 -b {input} > {output}"
    
# Clean up intermediate data to this point
rule STEP12b_clean_intermediate_data:
    input:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam"
    output:
        "{path}operations/preprocessing/clean12b.{sample}.{repnum}.done"
    shell:
        """
        rm -f {wildcards.path}preprocessing/9dedup/*REP{wildcards.repnum}*.bam
        touch {output}
        """

# Build the .bai index for the processed bam file
rule STEP13_build_index:
    # creates a bai index for the bam files
    # this is required for many downstream operations
    # the bai index allows other processes to access specific reads in the bam file without having to read through the entire bam contents to find them (its like a table of contents)
    # I specifies the input bam file
    # O specifies the output index file
    input:
        "{path}operations/preprocessing/clean12b.{sample}.{repnum}.done",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam"
    output:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    resources:
        mem_mb=50000
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.buildindex.benchmark.txt'
    shell:
        "java -jar snakeResources/programs/picard/picard.jar BuildBamIndex \
        I={input.b} \
        O={output}"
    
# Make a bigwig file from the bam file
rule STEP14_makebigwig_bamcov:
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}preprocessing/11bigwig/{sample}-REP{repnum}.bw"
    threads:
        20
    resources:
        mem_mb=25000
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.makebigwig.benchmark.txt'
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"
    
# Call peaks using global normalization
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
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.callpeaks.globalnorm.benchmark.txt'
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}_globalnorm --outdir {wildcards.path}peaks/globalnorm --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
    
# Call peaks using local normalization
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
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}peaks/localnorm/{sample}-REP{repnum}_localnorm_peaks.narrowPeak"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.callpeaks.localnorm.benchmark.txt'
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}_localnorm --outdir {wildcards.path}peaks/localnorm --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"
    
# Calculate percent genome coverage from peaks with global normalization
rule STEP17a_percent_peak_genome_coverage_globalnorm:
    # returns a fraction value of the basepairs of the genome covered by the merged peak file. multiple by 100 for percentages
    # parameters:
    # --echo output will be at least a three-column bed file
    # --bases-uniq the number of distinct bases from ref covered by overlap bed file
    # --delim change output delimeter from '|' to <delim>, e.g. '\t'
    input:
        a="{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/{sample}-REP{repnum}.peak.globalnorm.genomecov.txt"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.genomecov.globalnorm.benchmark.txt'
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"
    
# Calculate percent genome coverage from peaks with local normalization
rule STEP17b_percent_peak_genome_coverage_localnorm:
    # returns a fraction value of the basepairs of the genome covered by the merged peak file. multiple by 100 for percentages
    # parameters:
    # --echo output will be at least a three-column bed file
    # --bases-uniq the number of distinct bases from ref covered by overlap bed file
    # --delim change output delimeter from '|' to <delim>, e.g. '\t'
    input:
        a="{path}peaks/localnorm/{sample}-REP{repnum}_localnorm_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/{sample}-REP{repnum}.peak.localnorm.genomecov.txt"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.genomecov.localnorm.benchmark.txt'
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"
    
# Generate the fragment size distribution graph
rule STEP18_fragment_size_distribution:
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}metrics/{sample}-REP{repnum}.fragsizes.svg"
    resources:
        mem_mb=20000
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.fragsizes.benchmark.txt'
    script:
        "snakeResources/scripts/QC/snakeFragSizeDist.R"
    
# Annotate the peaks with global normalization
rule STEP19_annotate_peaks_global:
    input:
        "{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak"
    output:
        "{path}operations/preprocessing/{sample}-REP{repnum}.globalpeak.annotations.done"
    resources:
        mem_mb=40000
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.globalpeak.annotations.benchmark.txt'
    script:
        "snakeResources/scripts/QC/snakeAnnotatePeaks.R"

# Annotate the peaks with local normalization
rule STEP19_annotate_peaks_local:
    input:
        "{path}peaks/localnorm/{sample}-REP{repnum}_localnorm_peaks.narrowPeak"
    output:
        "{path}operations/preprocessing/{sample}-REP{repnum}.localpeak.annotations.done"
    resources:
        mem_mb=40000
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.localpeak.annotations.benchmark.txt'
    script:
        "snakeResources/scripts/QC/snakeAnnotatePeaks.R"

# Count the total number of reads in the sample
rule STEP20_sample_total_reads:
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}metrics/{sample}-REP{repnum}.totalreads.Rdata"
    benchmark:
        '{path}benchmark/preprocessing/{sample}-REP{repnum}.totalreads.benchmark.txt'
    script:
        "snakeResources/scripts/QC/snakeCountSampleReads.R"

# Analyze saturation of library in terms of library complexity, called peaks, and selected TF footprints
rule STEP21_saturation_analysis:
    input:
        "{path}operations/saturation/{sample}-REP{repnum}.downsample.done",
        "{path}operations/preprocessing/saturation/clean1.{sample}.{repnum}.done",
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_globalnorm_numpeaks.txt",
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_localnorm_numpeaks.txt",
        "{path}operations/saturation/footprints/{sample}-REP{repnum}.allgenes.footprint.downsampled.done"
    output:
        "{path}operations/saturation/{sample}-REP{repnum}.saturation_analysis.aggregator"
    shell:
    	"touch {output}"

########################################################################################################################################
#### DOWNSAMPLE RULES ##################################################################################################################
########################################################################################################################################

# When all saturation analysis is done, remove the intermediate data, as it is quite large
rule SATURATION_clean_data_final:
    input:
        "{path}operations/saturation/{sample}-REP{repnum}.saturation_analysis.aggregator"
    output:
        "{path}operations/saturation/{sample}-REP{repnum}.saturation_analysis.done"
    shell:
        """
        rm -f {wildcards.path}preprocessing/saturation/downsampled/md/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/md/*REP{wildcards.repnum}*.bai
        rm -f {wildcards.path}preprocessing/saturation/peaks/*REP{wildcards.repnum}*
        touch {output}
        """

# This rule determines which genes will be analyzed for the footprinting saturation analysis
rule AGGREGATOR_saturation_footprints_genes:
    input:
        "{path}operations/saturation/footprints/{sample}-REP{repnum}.CTCF.footprint.downsampled.done"
    output:
        "{path}operations/saturation/footprints/{sample}-REP{repnum}.allgenes.footprint.downsampled.done"
    shell:
        """
        rm -f {wildcards.path}preprocessing/8merged/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/raw/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/cs/*REP{wildcards.repnum}*.bam
        touch {output}
        """

# Determines the downsampling levels of the libraries
rule AGGREGATOR_saturation_downsample:
    input:
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.9.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.8.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.7.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.6.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.5.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.4.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.3.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.2.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.1.md.bai"      
    output:
        "{path}operations/saturation/{sample}-REP{repnum}.downsample.done"
    shell:
        "touch {output}"

# Downsample the processed but NOT duplicate purged .bam files
rule SATURATION_downsample_bam:
    input:
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    output:
        "{path}preprocessing/saturation/downsampled/raw/{sample}-REP{repnum}.{prob}.bam"
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.downsample.bam.benchmark.txt'
    shell:
        "java -jar snakeResources/programs/picard/picard.jar DownsampleSam \
        I={input} \
        O={output} \
        PROBABILITY=0.{wildcards.prob}"

# Coordinate sort the downsampled .bam files
rule SATURATION_sort_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/raw/{sample}-REP{repnum}.{prob}.bam"
    output:
        "{path}preprocessing/saturation/downsampled/cs/{sample}-REP{repnum}.{prob}.cs.bam"
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.sort.downsampled.benchmark.txt'
    shell:
        "java -jar snakeResources/programs/picard/picard.jar SortSam \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate"

# Purge duplicates from the downsampled .bam files
rule SATURATION_purge_duplicates_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/cs/{sample}-REP{repnum}.{prob}.cs.bam"
    output:
        a="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam",
        b="{path}metrics/saturation/{sample}-REP{repnum}.{prob}.duplication-metrics.txt"
    resources:
        purgeDuplicates=1
        mem_mb=50000
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.purgeduplicates.downsampled.benchmark.txt'
    shell:
        "java -Xmx50g -jar snakeResources/programs/picard/picard.jar MarkDuplicates \
        I={input} \
        O={output.a} \
        M={output.b} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"

# Generate .bai index for each downsampled .bam file
rule SATURATION_index_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam"
    output:
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bai"
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.index.downsampled.benchmark.txt'
    shell:
        "java -jar snakeResources/programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"

# Determine the library complexity of the downsampled libraries and output to metrics
rule SATURATION_analyze_complexity_downsampled:
    input:
        a="{path}metrics/saturation/{sample}-REP{repnum}.9.duplication-metrics.txt",
        b="{path}metrics/saturation/{sample}-REP{repnum}.8.duplication-metrics.txt",
        c="{path}metrics/saturation/{sample}-REP{repnum}.7.duplication-metrics.txt",
        d="{path}metrics/saturation/{sample}-REP{repnum}.6.duplication-metrics.txt",
        e="{path}metrics/saturation/{sample}-REP{repnum}.5.duplication-metrics.txt",
        f="{path}metrics/saturation/{sample}-REP{repnum}.4.duplication-metrics.txt",
        g="{path}metrics/saturation/{sample}-REP{repnum}.3.duplication-metrics.txt",
        h="{path}metrics/saturation/{sample}-REP{repnum}.2.duplication-metrics.txt",
        i="{path}metrics/saturation/{sample}-REP{repnum}.1.duplication-metrics.txt",
        j="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.9.md.bai",
        k="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.8.md.bai",
        l="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.7.md.bai",
        m="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.6.md.bai",
        n="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.5.md.bai",
        o="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.4.md.bai",
        p="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.3.md.bai",
        q="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.2.md.bai",
        r="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.1.md.bai"
    output:
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_library_size.txt"
    shell:
        """
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.a} >> {output}
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.b} >> {output}
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.c} >> {output}
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.d} >> {output}
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.e} >> {output}
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.f} >> {output}
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.g} >> {output}
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.h} >> {output}
        awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.i} >> {output}
        """

# Cleanup the uneeded intermediate files
rule SATURATION_clean_intermediate_data:
    input:
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_library_size.txt"
    output:
        "{path}operations/preprocessing/saturation/clean1.{sample}.{repnum}.done"
    shell:
        """
        rm -f {wildcards.path}preprocessing/8merged/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/raw/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/cs/*REP{wildcards.repnum}*.bam
        touch {output}
        """

# Call peaks with global normalization from downsampled libraries
rule SATURATION_peaks_globalnorm:
    input:
        a="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam",
        b="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bai"
    output:
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.{prob}_globalnorm_peaks.xls"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}.{wildcards.prob}_globalnorm --outdir {wildcards.path}preprocessing/saturation/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# Call peaks with local normalization from downsampled libraries
rule SATURATION_peaks_localnorm:
    input:
        a="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam",
        b="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bai"
    output:
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.{prob}_localnorm_peaks.xls"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}.{wildcards.prob}_localnorm --outdir {wildcards.path}preprocessing/saturation/peaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

# Count the number of peaks with global normalization from downsampled libraries and output to metrics
rule SATURATION_analyze_peak_saturation_globalnorm:
    input:
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.9_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.8_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.7_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.6_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.5_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.4_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.3_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.2_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.1_globalnorm_peaks.xls"
    output:
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_globalnorm_numpeaks.txt"
    shell:
        "wc -l < {input} >> {output}"

# Count the number of peaks with local normalization from downsampled libraries and output to metrics
rule SATURATION_analyze_peak_saturation_localnorm:
    input:
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.9_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.8_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.7_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.6_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.5_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.4_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.3_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.2_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.1_localnorm_peaks.xls"
    output:
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_localnorm_numpeaks.txt"
    shell:
        "wc -l < {input} >> {output}"

# An aggregator for the footprinting saturation analysis
rule AGGREGATOR_saturation_footprints:
    input:
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.9.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.8.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.7.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.6.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.5.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.4.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.3.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.2.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.1.done"       
    output:
        "{path}operations/saturation/footprints/{sample}-REP{repnum}.{gene}.footprint.downsampled.done"
    shell:
        "touch {output}"

# Perform the footprint saturation analysis
rule SATURATION_analyze_raw_footprint_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bai",
        "snakeResources/sites/data/genes/{gene}.bindingSites.Rdata",
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_library_size.txt"
    output:
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.{prob}.done"
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.rawfootprint.downsampled.benchmark.txt'
    script:
        "snakeResources/scripts/saturation/snakeAnalyzeRawFootprintSaturation.R"

########################################################################################################################
#### FOOTPRINTING ######################################################################################################
########################################################################################################################

# This rule initiates the raw footprint analysis for all genes found in the config file
rule AGGREGATOR_footprinting_raw_analysis:
    input:
        "{path}operations/preprocessing/{sample}-REP{repnum}.preprocessing.complete",
        expand("{{path}}operations/footprints/raw/{{sample}}-REP{{repnum}}.{genename}.rawFPanalysis.done", genename=config["geneNames"])
    output:
        "{path}operations/footprints/{sample}-REP{repnum}.footprinting_raw_analysis.complete"
    shell:
        "touch {output}"

# Generate the raw data used for downstream footprint analysis
rule FOOTPRINTING_raw_analysis:
    input:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai",
        "snakeResources/sites/data/genes/{gene}.bindingSites.Rdata"
    output:
        "{path}operations/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.done"
    resources:
        rawFPanalysis=1,
        mem_mb=25000
    benchmark:
        '{path}benchmark/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.benchmark.txt'
    script:
        "snakeResources/scripts/footprints/snakeAnalyzeRawFootprint.R"

# # Parse the sites based on the signals present
# rule FOOTPRINTING_parse_sites:
#     input:
#         "{path}operations/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.bamcopy{bamcopy}.done",
#         "{path}metrics/{sample}-REP{repnum}.totalreads.Rdata",
#         "{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak"
#     output:
#         "{path}footprints/operations/parse/{sample}-REP{repnum}.{gene}.parseFP.bamcopy{bamcopy}.done"
#     resources:
#         parseFootprint=1
#     benchmark:
#         '{path}footprints/benchmark/parse/{sample}-REP{repnum}.{gene}.bamcopy{bamcopy}.parseFP.txt'
#     script:
#         "scripts/panTF/snakeParseAndGenerateFootprintStats.R"

# # Process the sites and perform data analysis
# rule FOOTPRINTING_process_sites:
#     input:
#         "{path}footprints/operations/parse/{sample}-REP{repnum}.{gene}.parseFP.bamcopy{bamcopy}.done"
#     output:
#         "{path}footprints/operations/processed/{sample}-REP{repnum}.{gene}.processFP.bamcopy{bamcopy}.done"
#     resources:
#         processFootprint=1
#     benchmark:
#         '{path}footprints/benchmark/processed/{sample}-REP{repnum}.{gene}.bamcopy{bamcopy}.parseFP.txt'
#     script:
#         "scripts/panTF/snakeProcessFootprint.R"

# # Generate the footprinting figures from the data
# rule FOOTPRINTING_generate_figures:
#     input:
#         "{path}footprints/operations/processed/{sample}-REP{repnum}.{gene}.processFP.bamcopy{bamcopy}.done"
#     output:
#         "{path}footprints/operations/graphs/{sample}-REP{repnum}.{gene}.graphFP.bamcopy{bamcopy}.done"
#     resources:
#         graphFootprint=1
#     script:
#         "scripts/panTF/snakeGenerateFootprintGraphs.R"

# # Aggregate the footprinting data into a single R object
# rule FOOTPRINTING_aggregate_data:
#     input:
#         "{path}footprints/data/processed/"
#     output:
#         "{path}footprints/operations/aggregated/{sample}-REP{repnum}.aggregated.done"
#     script:
#         "scripts/panTF/snakeAggregateProcessedFootprintData.R"

########################################################################################################################
#### SECTORED FOOTPRINTING  ############################################################################################
########################################################################################################################

# This rule initiates the raw footprint analysis for all genes found in the config file that require sectored analysis
rule AGGREGATOR_sectored_footprinting_raw_analysis:
    input:
        "{path}operations/preprocessing/{sample}-REP{repnum}.preprocessing.complete",
        "{path}operations/footprints/{sample}-REP{repnum}.footprinting_raw_analysis.complete"
        #expand("{{path}}operations/footprints/merged/{{sample}}-REP{{repnum}}.{genename}.rawFPsectored.merged", genename=config["geneNamesSectored"])
    output:
        "{path}operations/footprints/{sample}-REP{repnum}.sectored_footprinting_analysis_raw.complete"
    shell:
        "touch {output}"

# Perform the raw fp analysis by sector. Output temporary sectored .Rdata files
rule FOOTPRINTING_raw_analysis_sectored:
    input:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai",
        "snakeResources/sites/data/genes/{gene}.bindingSites.Rdata"
    output:
        "{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector{sector}.done"
    resources:
        rawFPanalysis=1,
        mem_mb=200000
    benchmark:
        '{path}benchmark/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.sector{sector}.benchmark.txt'
    script:
        "snakeResources/scripts/footprints/snakeAnalyzeRawFootprintSectored.R"

# This rule aggregates the spooling for all 20 sectors for the larger footprints for raw analysis
rule AGGREGATOR_raw_analysis_sectored:
	input:
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector1.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector2.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector3.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector4.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector5.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector6.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector7.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector8.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector9.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector10.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector11.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector12.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector13.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector14.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector15.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector16.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector17.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector18.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector19.done",
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFP.sector20.done"
	output:
		"{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFPsectored.done"
	shell:
		"touch {output}"

# Merge the sectored raw footprint analysis .Rdata files
rule FOOTPRINTING_merge_sectored_raw_footprints:
    input:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai",
        "snakeResources/sites/data/genes/{gene}.bindingSites.Rdata",
        "{path}operations/footprints/temp/{sample}-REP{repnum}.{gene}.rawFPsectored.done"
    output:
        "{path}operations/footprints/merged/{sample}-REP{repnum}.{gene}.rawFPsectored.merged"
    resources:
        mergeRawFPSectors=1,
        mem_mb=400000
    benchmark:
        '{path}benchmark/footprints/merge/{sample}-REP{repnum}.{gene}.mergeRawFP.sector{sector}.benchmark.txt'
    script:
        "snakeResources/scripts/footprints/snakeMergeRawFootprintSectored.R"