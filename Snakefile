########################################################################################################################################
################################ GENERAL INFO ##########################################################################################
########################################################################################################################################
## Snakemake execution guide
# A dry run of the pipeline can be run with:
# snakemake -np h508go
# On one of the virtualization servers, run the pipeline with the following to allocate 20 threads and 90 gb max memory (to avoid crashing the process)
# snakemake -j 20 h508go --resources mem_gb=90
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
########################################################################################################################################
#### SPOOL PREPROCESSING ###############################################################################################################
########################################################################################################################################

rule run_h508wt01:
    input:
        "h508/wt01/preprocessing/logs/H508-WT-01.preprocessing.cleaning.done.txt"

rule run_ls1034wt01:
    input:
        "ls1034/wt01/preprocessing/logs/LS1034-WT-01.preprocessing.cleaning.done.txt"

rule run_snu61wt01:
    input:
        "snu61/wt01/preprocessing/logs/SNU61-WT-01.preprocessing.cleaning.done.txt"

rule run_mdst8wt01:
	input:
		"mdst8/wt01/operations/MDST8-WT-01-pipeline.complete.txt"

rule AGGREGATOR_pipeline:
	input:
		"{path}operations/{mergedsample}-correlation.done.txt",
		"{path}operations/{mergedsample}-peaks.done.txt",
		"{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt",
		"{path}operations/{mergedsample}-REP1of2-downsample.done.txt",
		"{path}operations/{mergedsample}-REP2of2-downsample.done.txt",
	output:
		"{path}operations/{mergedsample}-pipeline.complete.txt"
	shell:
		"touch {output}"

########################################################################################################################################
#### SPOOL FOOTPRINTING ################################################################################################################
########################################################################################################################################

#rule run_fp_coadmr_mdst8_wt01:
#    input:
#        "mdst8/wt01/footprints/parsed/MDST8-WT-01.coadmr.parsed.done.txt"

# rule run_fp_coadmr_h508wt02a:
#     input:
#         "h508/wt02a/footprints/parsed/H508A-WT-02.CDX2.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.TCF7.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.HOXA3.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.MNX1.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.POU5F1B.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.OVOL1.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.ESRRA.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.ASCL2.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.HNF4A.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.GMEB2.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.ZSWIM1.parsed.done.txt",
#         "h508/wt02a/footprints/parsed/H508A-WT-02.CBFA2T2.parsed.done.txt"

########################################################################################################################################
#### SPOOL CROSS SAMPLE CORRELATION ####################################################################################################
########################################################################################################################################

# rule run_xsample_corr_h508_snu61_ls1034:
#     input:
#         "xsample_analysis/correlation/H508-wt-01.LS1034-wt-01.SNU61-wt-01.spearman.heatmap.svg"

# rule run_xsample_corr_replicates_h508_snu61_ls1034:
#     input:
#         "xsample_analysis/correlation/H508-wt-01.LS1034-wt-01.SNU61-wt-01.spearman.heatmap.svg"

########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################

rule PREP_builddirstructure:
    # params: -p ignore error if existing, make parent dirs, -v verbose
    output:
        "{path}preprocessing/operations/dirtree.built.done"
    shell:
        """
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}preprocessing/2fastq {wildcards.path}preprocessing/3goodfastq {wildcards.path}preprocessing/4mycoalign {wildcards.path}preprocessing/5hg38align
        mkdir -p -v {wildcards.path}preprocessing/6rawbam {wildcards.path}preprocessing/7rgsort {wildcards.path}preprocessing/8merged {wildcards.path}preprocessing/9dedup
        mkdir -p -v {wildcards.path}preprocessing/10unique {wildcards.path}preprocessing/11repmerged {wildcards.path}preprocessing/12bigwig {wildcards.path}preprocessing/operations
        mkdir -p -v {wildcards.path}preprocessing/6rawbam/mitochondrial {wildcards.path}preprocessing/6rawbam/blacklist {wildcards.path}preprocessing/6rawbam/nonblacklist
        mkdir -p -v {wildcards.path}saturation
        mkdir -p -v {wildcards.path}saturation/complexity {wildcards.path}saturation/footprints {wildcards.path}saturation/peaks {wildcards.path}saturation/downsampled
        mkdir -p -v {wildcards.path}saturation/footprints/data {wildcards.path}saturation/footprints/graphs
        mkdir -p -v {wildcards.path}saturation/footprints/data/merged {wildcards.path}saturation/footprints/data/motifmerge {wildcards.path}saturation/footprints/data/parsed {wildcards.path}saturation/footprints/data/bychr
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}footprints/graphs {wildcards.path}footprints/heatmaps {wildcards.path}footprints/data
        mkdir -p -v {wildcards.path}footprints/data/merged {wildcards.path}footprints/data/motifmerge {wildcards.path}footprints/parsed {wildcards.path}footprints/bychr {wildcards.path}footprints/operations
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}peaks/genrich {wildcards.path}peaks/macs2 {wildcards.path}peaks/macs2/individual {wildcards.path}peaks/macs2/merged
        mkdir -p -v {wildcards.path}correlation
        mkdir -p -v {wildcards.path}metrics
        touch {output}
        """

rule STEP1_gunzip:
	# params:
	# -k keep original files
	# -c write to standard output
	input:
		a="{path}preprocessing/1gz/{sample}-REP{repnum}of{reptot}_L{lane}_R{read}.fastq.gz",
		b="{path}preprocessing/operations/dirtree.built.done"
	output:
		c="{path}preprocessing/2fastq/{sample}-REP{repnum}of{reptot}_L{lane}_R{read}.fastq"
	shell:
		"gunzip -k -c {input.a} > {output.c}"

rule STEP2_afterqc_fastqfiltering:
    # params:
    # -1 specifies read 1 fastq file
    # -2 specifies read 2 fastq file
    # -g specifies the output directory for the good fastq files
    # -b specifies the output directory for the bad fastq files
    # -f -1 autodetects number of bases to trim at front
    # -t -1 autodetects number of bases to trim at tail
    # -s is the shortest trimmed read length allowed past QC filter
    input:
        a="{path}preprocessing/2fastq/{sample}-REP{repnum}of{reptot}_L{lane}_R1.fastq",
        b="{path}preprocessing/2fastq/{sample}-REP{repnum}of{reptot}_L{lane}_R2.fastq"
    output:
        c="{path}preprocessing/3goodfastq/{sample}-REP{repnum}of{reptot}_L{lane}_R1.good.fq",
        d="{path}preprocessing/3goodfastq/{sample}-REP{repnum}of{reptot}_L{lane}_R2.good.fq"
    shell:
        "after.py -1 {input.a} -2 {input.b} -g {wildcards.path}preprocessing/3goodfastq -b {wildcards.path}preprocessing/3goodfastq -f -1 -t -1 -s 15"

rule STEP3_mycoalign:
    # params:
    # -q fastq input file format
    # -p num threads to use
    # -X1000 align to a maximum of 2000 bp frag length
    # -1 is read 1 input fastq file
    # -2 is read 2 input fastq file
    # -S output file path
    # 2> bowtie2 outputs alignment metrics to STDERR, 2> will allow redirect to a text file
    input:
        a="{path}preprocessing/3goodfastq/{sample}-REP{repnum}of{reptot}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}-REP{repnum}of{reptot}_L{lane}_R2.good.fq"
    output:
        "{path}preprocessing/4mycoalign/{sample}-REP{repnum}of{reptot}_L{lane}.myco.sam"
    shell:
        "bowtie2 -q -p 20 -X2000 -x genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}of{wildcards.reptot}_L{wildcards.lane}.myco.alignment.txt"

rule STEP4_hg38align:
    # params:
    # -q fastq input file format
    # -p num threads to use
    # -X1000 align to a maximum of 2000 bp frag length
    # -1 is read 1 input fastq file
    # -2 is read 2 input fastq file
    # -S output file path
    # 2> bowtie2 outputs alignment metrics to STDERR, 2> will allow redirect to a text file
    input:
        a="{path}preprocessing/3goodfastq/{sample}-REP{repnum}of{reptot}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}-REP{repnum}of{reptot}_L{lane}_R2.good.fq",
        c="{path}preprocessing/4mycoalign/{sample}-REP{repnum}of{reptot}_L{lane}.myco.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}of{reptot}_L{lane}.hg38.sam"
    shell:
        "bowtie2 -q -p 20 -X2000 -x genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}of{wildcards.reptot}_L{wildcards.lane}.hg38.alignment.txt"

rule STEP5_coordsort_sam:
    # coordinate sorting the sam files is required for blacklist filtering
    # params:
    # -o output file path
    # -O output file format
    input:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}of{reptot}_L{lane}.hg38.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}of{reptot}_L{lane}.hg38.cs.sam"
    shell:
        "samtools sort {input} -o {output} -O sam"

rule STEP6_blacklistfilter_bamconversion:
    # remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
    input:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}of{reptot}_L{lane}.hg38.cs.sam"
    output:
        a="{path}preprocessing/6rawbam/blacklist/{sample}-REP{repnum}of{reptot}_L{lane}.hg38blacklist.bam",
        b="{path}preprocessing/6rawbam/nonblacklist/{sample}-REP{repnum}of{reptot}_L{lane}.blrm.bam"
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 20 {input}"

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
        "{path}preprocessing/6rawbam/nonblacklist/{sample}-REP{repnum}of{reptot}_L{lane}.blrm.bam"
    output:
        a="{path}preprocessing/6rawbam/mitochondrial/{sample}-REP{repnum}of{reptot}_L{lane}.mitochondrial.bam",
        b="{path}preprocessing/6rawbam/{sample}-REP{repnum}of{reptot}_L{lane}.goodbam"
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 20 {input}"

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
    # I specifies the input file
    # O specifies the output file
    input:
        "{path}preprocessing/6rawbam/{sample}-REP{repnum}of{reptot}_L{lane}.goodbam"
    output:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}of{reptot}_L{lane}.rg.cs.bam"
    shell:
        "java -jar programs/picard/picard.jar AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID=H5YHHBGX3.{wildcards.lane} \
        RGLB={wildcards.sample} \
        RGPL=ILLUMINA \
        RGPU=H5YHHBGX3.{wildcards.lane}.{wildcards.sample} \
        RGSM={wildcards.sample}"

rule STEP9_cleansam:
    # params:
    # soft-clips bases aligned past the end of the ref sequence
    # soft-clipping retains the bases in the SEQ string, but they are not displayed or used in downstream data analysis
    # sets MAPQ score to 0 for unmapped reads
    # I specifies the input file
    # O specifies the output file
    input:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}of{reptot}_L{lane}.rg.cs.bam"
    output:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}of{reptot}_L{lane}.clean.bam"
    shell:
        "java -jar programs/picard/picard.jar CleanSam \
        I={input} \
        O={output}"

rule STEP10_mergelanes:
	# Merge files for individual lanes
	# I specifies input files for each lane
	# O specifies the output files
	# SORT_ORDER/ASSUME_SORTED specify the type of sorting in the input files
	# MERGE_SEQUENCE_DICTIONARIES will combine the sequence dictionaries from the individual files
	# a sequence dictionary contains information about sequence name, length, genome assembly ID, etc
	# USE_THREADING allows multithreadded operation
    input:
        a="{path}preprocessing/7rgsort/{sample}-REP{repnum}of{reptot}_L1.clean.bam",
        b="{path}preprocessing/7rgsort/{sample}-REP{repnum}of{reptot}_L2.clean.bam",
        c="{path}preprocessing/7rgsort/{sample}-REP{repnum}of{reptot}_L3.clean.bam",
        d="{path}preprocessing/7rgsort/{sample}-REP{repnum}of{reptot}_L4.clean.bam"
    output:
    	"{path}preprocessing/8merged/{sample}-REP{repnum}of{reptot}.lanemerge.bam"
    shell:
        "java -jar programs/picard/picard.jar MergeSamFiles \
        I={input.a} \
        I={input.b} \
        I={input.c} \
        I={input.d} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

rule STEP11_purgeduplicates:
    # params:
    # I specifies the input file
    # O specifies the output file
    # M specifies the duplication metrics output file
    # REMOVE_DUPLICATES enables removal of duplicate reads from the output file
    # ASSUME_SORTED indicates the input file is already sorted
    input:
        "{path}preprocessing/8merged/{sample}-REP{repnum}of{reptot}.lanemerge.bam"
    output:
        "{path}preprocessing/9dedup/{sample}-REP{repnum}of{reptot}.dp.bam"
    shell:
        "java -jar programs/picard/picard.jar MarkDuplicates \
        I={input} \
        O={output} \
        M={wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}of{wildcards.reptot}.duplication.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"

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
        "{path}preprocessing/9dedup/{sample}-REP{repnum}of{reptot}.dp.bam"
    output:
        "{path}preprocessing/10unique/{sample}-REP{repnum}of{reptot}.u.bam"
    shell:
        "samtools view -h -q 2 -b {input} > {output}"

rule STEP13_buildindex:
    # creates a bai index for the bam files
    # this is required for many downstream operations
    # the bai index allows other processes to access specific reads in the bam file without having to read through the entire bam contents to find them (its like a table of contents)
    # I specifies the input bam file
    # O specifies the output index file
    input:
        "{path}preprocessing/10unique/{sample}-REP{repnum}of{reptot}.u.bam"
    output:
        "{path}preprocessing/10unique/{sample}-REP{repnum}of{reptot}.u.bai"
    shell:
        "java -jar programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"

rule STEP14_merge_1_replicate:
    # If only one replicate is present, you can just copy the previous bam file to the next directory
    input:
        "{path}preprocessing/10unique/{mergedsample}-REP1of1.u.bam"
    output:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam"
    shell:
        "cp {input} {output}"

rule STEP14_merge_2_replicates:
    # This rule will be called when there are two input replicates
    # Merges the bam files from the infividual replicates
    # I specifies the input files for individual replicates
    # O specifies the merged output file
    # SORT_ORDER/ASSUME_SORTED specify the type of sorting in the input files
	# MERGE_SEQUENCE_DICTIONARIES will combine the sequence dictionaries from the individual files
	# a sequence dictionary contains information about sequence name, length, genome assembly ID, etc
	# USE_THREADING allows multithreadded operation
    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bam"
    output:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam"
    shell:
        "java -jar programs/picard/picard.jar MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

rule STEP14_merge_3_replicates:
    # This rule will be called when there are three input replicates
    # Merges the bam files from the infividual replicates
    # I specifies the input files for individual replicates
    # O specifies the merged output file
    # SORT_ORDER/ASSUME_SORTED specify the type of sorting in the input files
	# MERGE_SEQUENCE_DICTIONARIES will combine the sequence dictionaries from the individual files
	# a sequence dictionary contains information about sequence name, length, genome assembly ID, etc
	# USE_THREADING allows multithreadded operation
    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bam",
        c="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bam"
    output:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam"
    shell:
        "java -jar programs/picard/picard.jar MergeSamFiles \
        I={input.a} \
        I={input.b} \
        I={input.c} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

rule STEP15_index_replicate_merged:
    # creates a bai index for the bam files
    # this is required for many downstream operations
    # the bai index allows other processes to access specific reads in the bam file without having to read through the entire bam contents to find them (its like a table of contents)
    # I specifies the input bam file
    # O specifies the output index file
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam"
    output:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai"
    shell:
        "java -jar programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"

rule STEP16_makebigwig_bamcov_individual:
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}of{reptot}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}of{reptot}.u.bai"
    output:
        "{path}preprocessing/12bigwig/{sample}-REP{repnum}of{reptot}.bw"
    shell:
    	"bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

rule STEP17_makebigwig_bamcov_merged_1replicate:
	# This rule will be used when only one replicate is present
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
    input:
        a="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        c="{path}preprocessing/12bigwig/{mergedsample}-REP1of1.bw"
    output:
        "{path}preprocessing/12bigwig/{mergedsample}-repmerged.bw"
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

rule STEP17_makebigwig_bamcov_merged_2replicates:
	# This rule will be used when two replicates are present
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
    input:
        a="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        c="{path}preprocessing/12bigwig/{mergedsample}-REP1of2.bw",
        d="{path}preprocessing/12bigwig/{mergedsample}-REP2of2.bw"
    output:
        "{path}preprocessing/12bigwig/{mergedsample}-repmerged.bw"
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

rule STEP17_makebigwig_bamcov_merged_3replicates:
	# This rule will be used when three replicates are present
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
    input:
        a="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        c="{path}preprocessing/12bigwig/{mergedsample}-REP1of3.bw",
        d="{path}preprocessing/12bigwig/{mergedsample}-REP2of3.bw",
        e="{path}preprocessing/12bigwig/{mergedsample}-REP2of3.bw"
    output:
        "{path}preprocessing/12bigwig/{mergedsample}-repmerged.bw"
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

rule STEP18_preprocessing_metrics_and_delete_intermediate_files:
	# gather and determine the various preprocessing metrics, record to output text file
	# delete unnecessary intermediate preprocessing files
	# -f option will ignore nonexistent files
	input:
		"{path}preprocessing/operations/{mergedsample}-preprocessing.aggregator.txt"
	output:
		"{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt"
	shell:
		"""
      	rm -f {wildcards.path}preprocessing/2fastq/*.fastq
      	rm -f {wildcards.path}preprocessing/3goodfastq/*.fq
       	rm -f {wildcards.path}preprocessing/4mycoalign/*.sam
       	rm -f {wildcards.path}preprocessing/5hg38align/*.sam
       	rm -f {wildcards.path}preprocessing/6rawbam/*.goodbam
       	rm -f {wildcards.path}preprocessing/6rawbam/blacklist/*.bam
       	rm -f {wildcards.path}preprocessing/6rawbam/mitochondrial/*.bam
       	rm -f {wildcards.path}preprocessing/6rawbam/nonblacklist/*.bam
       	rm -f {wildcards.path}preprocessing/7rgsort/*.bam
        touch {output}
        """

rule AGGREGATOR_preprocessing_steps:
    input:
        "{path}preprocessing/12bigwig/{mergedsample}-repmerged.bw"
    output:
        "{path}preprocessing/operations/{mergedsample}-preprocessing.aggregator.txt"
    shell:
        "touch {output}"

########################################################################################################################################
#### PEAK CALLING RULES ################################################################################################################
########################################################################################################################################

rule STEP19_MACS2_peaks_individual_global_normilization:
    # notes:
    # because we are going to use the TCGA data downstream likely as a reference point,
    # we will need to call the peaks in the exact same way as they did in this paper:
    # http://science.sciencemag.org/content/sci/suppl/2018/10/24/362.6413.eaav1898.DC1/aav1898_Corces_SM.pdf
    # which is "macs2 callpeak --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
    ## params:
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
        a="{path}preprocessing/10unique/{mergedsample}-REP{repnum}of{reptot}.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP{repnum}of{reptot}.u.bai",
        c="{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt"
    output:
        "{path}peaks/macs2/individual/{mergedsample}-REP{repnum}of{reptot}_global_normalization_peaks.xls"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.mergedsample}-REP{wildcards.repnum}of{wildcards.reptot}_global_normalization --outdir {wildcards.path}peaks/macs2/individual --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule STEP20_MACS2_peaks_individual_local_normalization:
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
        a="{path}preprocessing/10unique/{mergedsample}-REP{repnum}of{reptot}.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP{repnum}of{reptot}.u.bai",
        c="{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt"
    output:
        "{path}peaks/macs2/individual/{mergedsample}-REP{repnum}of{reptot}_local_normalization_peaks.xls"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.mergedsample}-REP{wildcards.repnum}of{wildcards.reptot}_local_normalization --outdir {wildcards.path}peaks/macs2/individual --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

rule STEP21_MACS2_peaks_merged_global_normilization_1replicate:
    # see above for notes applicable to MACS2 peak calling
    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP1of1.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP1of1.u.bai",
        c="{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt",
        d="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        e="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        f="{path}peaks/macs2/individual/{mergedsample}-REP1of1_global_normalization_peaks.xls",
        g="{path}peaks/macs2/individual/{mergedsample}-REP1of1_local_normalization_peaks.xls"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.xls"
    shell:
        "macs2 callpeak -t {input.d} -n {wildcards.mergedsample}-merged_global_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule STEP21_MACS2_peaks_merged_global_normilization_2replicates:
    # see above for notes applicable to MACS2 peak calling
    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bai",
        c="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bam",
        d="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bai",
        e="{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt",
        f="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        g="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        h="{path}peaks/macs2/individual/{mergedsample}-REP1of2_global_normalization_peaks.xls",
        i="{path}peaks/macs2/individual/{mergedsample}-REP1of2_local_normalization_peaks.xls",
        j="{path}peaks/macs2/individual/{mergedsample}-REP2of2_global_normalization_peaks.xls",
        k="{path}peaks/macs2/individual/{mergedsample}-REP2of2_local_normalization_peaks.xls"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.xls"
    shell:
        "macs2 callpeak -t {input.f} -n {wildcards.mergedsample}-merged_global_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule STEP21_MACS2_peaks_merged_global_normilization_3replicates:
    # see above for notes applicable to MACS2 peak calling
    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bai",
        c="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bam",
        d="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bai",
        e="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bam",
        f="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bai",
        g="{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt",
        h="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        i="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        j="{path}peaks/macs2/individual/{mergedsample}-REP1of3_global_normalization_peaks.xls",
        k="{path}peaks/macs2/individual/{mergedsample}-REP1of3_local_normalization_peaks.xls",
        l="{path}peaks/macs2/individual/{mergedsample}-REP2of3_global_normalization_peaks.xls",
        m="{path}peaks/macs2/individual/{mergedsample}-REP2of3_local_normalization_peaks.xls",
        n="{path}peaks/macs2/individual/{mergedsample}-REP3of3_global_normalization_peaks.xls",
        o="{path}peaks/macs2/individual/{mergedsample}-REP3of3_local_normalization_peaks.xls"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.xls"
    shell:
        "macs2 callpeak -t {input.h} -n {wildcards.mergedsample}-merged_global_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule STEP22_MACS2_peaks_merged_local_normilization:
    # see above for notes applicable to MACS2 peak calling
    input:
        a="{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.xls",
        b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        c="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_local_normalization_peaks.xls"
    shell:
        "macs2 callpeak -t {input.b} -n {wildcards.mergedsample}-merged_local_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule AGGREGATOR_peaks:
	input:
		"{path}peaks/macs2/merged/{mergedsample}-merged_local_normalization_peaks.xls"
	output:
		"{path}operations/{mergedsample}-peaks.done.txt"
	shell:
		"touch {output}"

########################################################################################################################################
#### SAMPLE CORRELATION ANALYSIS RULES #################################################################################################
########################################################################################################################################

rule STEP23_sample_correlation_spearman_2replicates:
    # parameters:
    # -b input bam files
    # -o output file name
    # -bs set the bin size used for comparison, default is 10000 bp
    # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
    # -p set the number of computing processors to use
    # -v verbose mode
    input:
    	a="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bam",
    	b="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bam",
    	c="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bai",
    	d="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bai"
    output:
        "{path}correlation/{mergedsample}.spearman.corrTest"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} -o {output} -bs 10000 -p 20 -v"

rule STEP23_sample_correlation_spearman_3replicates:
    # parameters:
    # -b input bam files
    # -o output file name
    # -bs set the bin size used for comparison, default is 10000 bp
    # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
    # -p set the number of computing processors to use
    # -v verbose mode
    input:
    	a="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bam",
    	b="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bam",
    	c="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bam",
    	d="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bai",
    	e="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bai",
    	f="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bai"
    output:
        "{path}correlation/{mergedsample}.spearman.corrTest"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} {input.c} -o {output} -bs 10000 -p 20 -v"

rule STEP24_makecorrheatmap:
    input:
        "{path}correlation/{sample}.spearman.corrTest"
    output:
        "{path}correlation/{sample}.spearman.heatmap.svg"
    shell:
        "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"

rule AGGREGATOR_correlation:
	input:
		"{path}correlation/{mergedsample}.spearman.heatmap.svg"
	output:
		"{path}operations/{mergedsample}-correlation.done.txt"
	shell:
		"touch {output}"

########################################################################################################################################
#### SATURATION ANALYSIS RULES #########################################################################################################
########################################################################################################################################

rule STEP25_downsample_bam:
    input:
        "{path}preprocessing/8merged/{mergedsample}-REP{repnum}of{reptot}.lanemerge.bam"
    output:
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.bam"
    shell:
        "java -jar programs/picard/picard.jar DownsampleSam \
        I={input} \
        O={output} \
        PROBABILITY=0.{wildcards.prob}"

rule STEP26_sort_downsampled:
    input:
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.bam"
    output:
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.cs.bam"
    shell:
        "java -jar programs/picard/picard.jar SortSam \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate"

rule STEP27_purge_duplicates_downsampled:
    input:
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.cs.bam"
    output:
        a="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bam",
        b="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.{prob}.duplication-metrics.txt",
    shell:
        "java -Xmx5g -jar programs/picard/picard.jar MarkDuplicates \
        I={input} \
        O={output.a} \
        M={output.b} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"

rule STEP28_index_downsampled:
    input:
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bam"
    output:
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bai"
    shell:
        "java -jar programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"

rule STEP29_analyze_complexity_downsampled:
    input:
        a="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.9.duplication-metrics.txt",
        b="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.8.duplication-metrics.txt",
        c="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.7.duplication-metrics.txt",
        d="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.6.duplication-metrics.txt",
        e="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.5.duplication-metrics.txt",
        f="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.4.duplication-metrics.txt",
        g="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.3.duplication-metrics.txt",
        h="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.2.duplication-metrics.txt",
        i="{path}metrics/{mergedsample}-REP{repnum}of{reptot}.1.duplication-metrics.txt",
        j="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.9.md.bai",
        k="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.8.md.bai",
        l="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.7.md.bai",
        m="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.6.md.bai",
        n="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.5.md.bai",
        o="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.4.md.bai",
        p="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.3.md.bai",
        q="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.2.md.bai",
        r="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.1.md.bai",
    output:
        "{path}metrics/{mergedsample}-REP{repnum}of{reptot}.downsampled_library_size.txt"
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
    
rule STEP30_MACS2_peaks_downsampled:
    input:
    	a="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bam",
        b="{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bai"
    output:
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.{prob}_global_normalization_peaks.xls"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.mergedsample}-REP{wildcards.repnum}of{wildcards.reptot}.{wildcards.prob}_global_normalization --outdir {wildcards.path}saturation/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule STEP31_analyze_peak_saturation_downsampled:
    input:
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.9_global_normalization_peaks.xls",
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.8_global_normalization_peaks.xls",
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.7_global_normalization_peaks.xls",
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.6_global_normalization_peaks.xls",
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.5_global_normalization_peaks.xls",
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.4_global_normalization_peaks.xls",
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.3_global_normalization_peaks.xls",
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.2_global_normalization_peaks.xls",
        "{path}saturation/peaks/{mergedsample}-REP{repnum}of{reptot}.1_global_normalization_peaks.xls"
    output:
        "{path}metrics/{mergedsample}-REP{repnum}of{reptot}.downsampled_numpeaks.txt"
    shell:
        "wc -l < {input} >> {output}"

# FP selections can be added here, for example:
#"{path}operations/{mergedsample}-REP{repnum}of{reptot}.allprob.MNX1.done.graphs.txt"
rule AGGREGATOR_saturation:
	input:
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.9.md.bai",
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.8.md.bai",
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.7.md.bai",
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.6.md.bai",
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.5.md.bai",
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.4.md.bai",
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.3.md.bai",
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.2.md.bai",
		"{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.1.md.bai",
		"{path}metrics/{mergedsample}-REP{repnum}of{reptot}.downsampled_library_size.txt",
		"{path}metrics/{mergedsample}-REP{repnum}of{reptot}.downsampled_numpeaks.txt"		
	output:
		"{path}operations/{mergedsample}-REP{repnum}of{reptot}-downsample.done.txt"
	shell:
		"touch {output}"

# rule STEP32_make_footprint_by_chr_downsampled:
#     input:
#         "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bam",
#         "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bai",
#         "sites/{gene}.sites.Rdata"
#     output:
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.{chr}.done.txt"
#     script:
#         "scripts/snakeMakeFPbyChrDownsampled.R"

# rule STEP33_merge_footprint_by_chr_downsampled:
#     input:
#         "sites/{gene}.sites.Rdata",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr1.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr2.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr3.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr4.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr5.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr6.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr7.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr8.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr9.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr10.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr11.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr12.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr13.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr14.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr15.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr16.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr17.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr18.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr19.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr20.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr21.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr22.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chrY.done.txt",
#         "{path}saturation/footprints/data/bychr/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chrX.done.txt"
#     output:
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.done.merged.txt"
#     script:
#         "scripts/snakeMergeFPbyChrDownsampled.R"

# rule STEP34_make_footprint_graph_downsampled:
#     input:
#         "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bam",
#         "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bai",
#         "sites/{gene}.sites.Rdata",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.done.merged.txt"
#     output:
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.graphs.done.txt"
#     script:
#         "scripts/snakeGenerateMergedFPGraphDownsampled.R"

# rule AGGREGATOR_saturation_footprint_graphs:
#     input:
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.9.{gene}.graphs.done.txt",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.8.{gene}.graphs.done.txt",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.7.{gene}.graphs.done.txt",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.6.{gene}.graphs.done.txt",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.5.{gene}.graphs.done.txt",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.4.{gene}.graphs.done.txt",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.3.{gene}.graphs.done.txt",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.2.{gene}.graphs.done.txt",
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.1.{gene}.graphs.done.txt"
#     output:
#         "{path}operations/{mergedsample}-REP{repnum}of{reptot}.allprob.{gene}.done.graphs.txt"
#     shell:
#         "touch {output}"

########################################################################################################################################
#### FOOTPRINT ANALYSIS RULES ##########################################################################################################
########################################################################################################################################

rule STEP35_make_footprint_signal_by_chr:
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        "sites/{gene}.sites.Rdata"
    output:
        "{path}footprints/operations/{mergedsample}.{gene}.{chr}.done.bychr.txt"
    script:
        "scripts/snakeMakeFPbyChr.R"

rule STEP36_merge_footprint_signal_by_chr:
    input:
        "sites/{gene}.sites.Rdata",
        "{path}footprints/operations/{mergedsample}.{gene}.chr1.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr2.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr3.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr4.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr5.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr6.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr7.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr8.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr9.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr10.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr11.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr12.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr13.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr14.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr15.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr16.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr17.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr18.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr19.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr20.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr21.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chr22.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chrX.done.bychr.txt",
        "{path}footprints/operations/{mergedsample}.{gene}.chrY.done.bychr.txt"
    output:
        "{path}footprints/operations/{mergedsample}.{gene}.merged.done.txt"
    script:
        "scripts/snakeMergeFPbyChr.R"

rule STEP37_make_raw_footprint_graphs:
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        "sites/{gene}.sites.Rdata",
        "{path}footprints/operations/{mergedsample}.{gene}.merged.done.txt"
    output:
        "{path}footprints/operations/{mergedsample}.{gene}.graphs.done.txt"
    script:
        "scripts/snakeGenerateMergedFPGraph.R"

rule STEP38_parse_footprint_signals:
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        "sites/{gene}.sites.Rdata",
        "{path}footprints/operations/{mergedsample}.{gene}.merged.done.txt",
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
    output:
        "{path}footprints/operations/{mergedsample}.{gene}.parsed.done.txt"
    script:
        "scripts/snakeParseFP.R"

rule STEP39_make_parsed_footprint_heatmaps:
    input:
        "{path}footprints/parsed/{mergedsample}.{gene}.motif{motif}.info.Rdata",
    output:
        "{path}footprints/heatmaps/{mergedsample}.{gene}.motif{motif}.heatmap.svg"
    script:
        "scripts/snakeFootprintHeatmaps.R"

rule STEP40_merge_footprint_signal_motifs:
    input:
        "{path}footprints/operations/{mergedsample}.{gene}.parsed.done.txt"
    output:
        "{path}footprints/data/motifmerge/{mergedsample}.{gene}.{nummotif}.mergedmotif.Rdata"
    script:
        "scripts/snakeMergeMotifs.R"

rule AGGREGATOR_COADMR_footprinting:
	input:
		"mdst8/wt01/footprints/operations/MDST8-WT-01.CDX2.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.TCF7.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.HOXA3.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.MNX1.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.POU5F1B.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.OVOL1.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.ESRRA.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.ASCL2.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.HNF4A.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.GMEB2.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.ZSWIM1.parsed.done.txt",
		"mdst8/wt01/footprints/operations/MDST8-WT-01.CBFA2T2.parsed.done.txt"
	output:
		"mdst8/wt01/footprints/operations/MDST8-WT-01.coadmr.parsed.done.txt"
	shell:
		"touch {output}"

rule test:
	input:
		"mdst8/wt01/footprints/operations/MDST8-WT-01.coadmr.parsed.done.txt"