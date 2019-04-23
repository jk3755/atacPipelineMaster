
########################################################################################################################################
################################ GENERAL INFO ##########################################################################################
########################################################################################################################################
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
#
########################################################################################################################################
#### IMPORT MODULES AND CONFIG #########################################################################################################
########################################################################################################################################

include: "snakeModules/panTF.snakefile"

include: "snakeModules/scanPWM.snakefile"

#configfile: "snakeModules/config.yaml"

########################################################################################################################################
#### SPOOL PREPROCESSING ###############################################################################################################
########################################################################################################################################

rule run_h508wt01:
    input:
        "h508/wt01/preprocessing/logs/H508-WT-01.preprocessing.cleaning.done.txt"

rule run_h508wt02a:
    input:
        "h508/wt02a/operations/H508A-WT-02-pipeline.complete.txt"

rule run_ls1034wt01:
    input:
        "ls1034/wt01/operations/LS1034-WT-01-pipeline.complete.txt"

rule run_snu61wt01:
    input:
        "snu61/wt01/operations/SNU61-WT-01-pipeline.complete.txt"

rule run_mdst8wt01:
    input:
        "mdst8/wt01/operations/MDST8-WT-01-pipeline.complete.txt"

rule AGGREGATOR_pipeline:
    input:
        "{path}operations/{mergedsample}-correlation.done.txt",
        "{path}operations/{mergedsample}-peaks.done.txt",
        "{path}footprints/operations/{mergedsample}.footprints.coadmr.done.txt",
        "{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt",
        "{path}operations/{mergedsample}-downsample.final.txt",
        "{path}operations/{mergedsample}.metrics.annotations.done.txt"
    output:
        "{path}operations/{mergedsample}-pipeline.complete.txt"
    shell:
        "touch {output}"

########################################################################################################################################
#### SPOOL INDIVIDUAL OPERATIONS #######################################################################################################
########################################################################################################################################

rule run_PWMscan:
    # Run this rule to generate all needed data for scanning the genome for matches to PWMs
    # Will generate data for all annotated genes in motifDB, for all unique motifs
    input:
        "sites/motifData.Rdata",
        "sites/geneNames.txt",
        "sites/operations/groups/PWMscan.allgroups.done"

########################################################################################################################################
#### SPOOL TF saturation analysis ######################################################################################################
########################################################################################################################################

rule SATURATION_footprint_analysis:
    input:
        "snu61/wt01/operations/footsat/SNU61-WT-01-REP1of3.allprob.MNX1.done.parsed.txt"
    output:
        "snu61/wt01/operations/footsat/SNU61-WT-01-REP1of3.MNX1.footprint.satanalysis.done.txt"
    shell:
        "touch {output}"

########################################################################################################################################
#### SPOOL CROSS SAMPLE CORRELATION ####################################################################################################
########################################################################################################################################

rule run_xsample_corr_h508_snu61_ls1034:
    input:
        "xsample_analysis/correlation/H508-wt-01.LS1034-wt-01.SNU61-wt-01.spearman.heatmap.svg"

rule run_xsample_corr_replicates_h508_snu61_ls1034:
    input:
        "xsample_analysis/correlation/H508-wt-01.LS1034-wt-01.SNU61-wt-01.spearman.heatmap.svg"

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
        #
        mkdir -p -v {wildcards.path}saturation
        mkdir -p -v {wildcards.path}saturation/footprints
        mkdir -p -v {wildcards.path}saturation/footprints/data {wildcards.path}saturation/footprints/graphs {wildcards.path}saturation/footprints/operations {wildcards.path}saturation/footprints/benchmark
        mkdir -p -v {wildcards.path}saturation/complexity
        mkdir -p -v {wildcards.path}saturation/peaks
        mkdir -p -v {wildcards.path}saturation/downsampled
        #
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}footprints/benchmark
        mkdir -p -v {wildcards.path}footprints/benchmark/parse {wildcards.path}footprints/benchmark/raw {wildcards.path}footprints/benchmark/processed
        mkdir -p -v {wildcards.path}footprints/data 
        mkdir -p -v {wildcards.path}footprints/data/parsed {wildcards.path}footprints/data/raw {wildcards.path}footprints/data/processed
        mkdir -p -v {wildcards.path}footprints/graphs
		mkdir -p -v {wildcards.path}footprints/graphs/bf {wildcards.path}footprints/graphs/heatmaps {wildcards.path}footprints/graphs/peaks {wildcards.path}footprints/graphs/processed
		mkdir -p -v {wildcards.path}footprints/operations
        mkdir -p -v {wildcards.path}footprints/operations/groups {wildcards.path}footprints/operations/parse {wildcards.path}footprints/operations/raw {wildcards.path}footprints/operations/processed
        #
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}peaks/genrich
        mkdir -p -v {wildcards.path}peaks/macs2
        mkdir -p -v {wildcards.path}peaks/macs2/individual {wildcards.path}peaks/macs2/merged
        #
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}correlation
        #
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
    resources:
        hg38align=1
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
    resources:
        make_bigwig=1
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
    resources:
        make_bigwig=1
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
    resources:
        make_bigwig=1
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
    resources:
        make_bigwig=1
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
        "{path}peaks/macs2/individual/{mergedsample}-REP{repnum}of{reptot}_global_normalization_peaks.narrowPeak"
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
        "{path}peaks/macs2/individual/{mergedsample}-REP{repnum}of{reptot}_local_normalization_peaks.narrowPeak"
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
        f="{path}peaks/macs2/individual/{mergedsample}-REP1of1_global_normalization_peaks.narrowPeak",
        g="{path}peaks/macs2/individual/{mergedsample}-REP1of1_local_normalization_peaks.narrowPeak"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
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
        h="{path}peaks/macs2/individual/{mergedsample}-REP1of2_global_normalization_peaks.narrowPeak",
        i="{path}peaks/macs2/individual/{mergedsample}-REP1of2_local_normalization_peaks.narrowPeak",
        j="{path}peaks/macs2/individual/{mergedsample}-REP2of2_global_normalization_peaks.narrowPeak",
        k="{path}peaks/macs2/individual/{mergedsample}-REP2of2_local_normalization_peaks.narrowPeak"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
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
        j="{path}peaks/macs2/individual/{mergedsample}-REP1of3_global_normalization_peaks.narrowPeak",
        k="{path}peaks/macs2/individual/{mergedsample}-REP1of3_local_normalization_peaks.narrowPeak",
        l="{path}peaks/macs2/individual/{mergedsample}-REP2of3_global_normalization_peaks.narrowPeak",
        m="{path}peaks/macs2/individual/{mergedsample}-REP2of3_local_normalization_peaks.narrowPeak",
        n="{path}peaks/macs2/individual/{mergedsample}-REP3of3_global_normalization_peaks.narrowPeak",
        o="{path}peaks/macs2/individual/{mergedsample}-REP3of3_local_normalization_peaks.narrowPeak"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
    shell:
        "macs2 callpeak -t {input.h} -n {wildcards.mergedsample}-merged_global_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule STEP22_MACS2_peaks_merged_local_normilization:
    # see above for notes applicable to MACS2 peak calling
    input:
        a="{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak",
        b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        c="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_local_normalization_peaks.narrowPeak"
    shell:
        "macs2 callpeak -t {input.b} -n {wildcards.mergedsample}-merged_local_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule AGGREGATOR_peaks:
    input:
        "{path}peaks/macs2/merged/{mergedsample}-merged_local_normalization_peaks.narrowPeak"
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

rule FINISH_saturation_1rep:
    input:
        "{path}operations/{mergedsample}-REP1of1-downsample.done.txt"
    output:
        "{path}operations/{mergedsample}-downsample.final.txt"
    shell:
        "touch {output}"

rule FINISH_saturation_2rep:
    input:
        "{path}operations/{mergedsample}-REP1of2-downsample.done.txt",
        "{path}operations/{mergedsample}-REP2of2-downsample.done.txt",
    output:
        "{path}operations/{mergedsample}-downsample.final.txt"
    shell:
        "touch {output}"

rule FINISH_saturation_3rep:
    input:
        "{path}operations/{mergedsample}-REP1of3-downsample.done.txt",
        "{path}operations/{mergedsample}-REP2of3-downsample.done.txt",
        "{path}operations/{mergedsample}-REP3of3-downsample.done.txt",
    output:
        "{path}operations/{mergedsample}-downsample.final.txt"
    shell:
        "touch {output}"

rule STEP32_make_footprint_by_chr_downsampled:
    input:
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bam",
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bai",
        "sites/{gene}.sites.Rdata"
    output:
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.{chr}.footprint.bychr.sat.done.txt"
    script:
        "scripts/saturation/snakeMakeFPbyChrDownsampled.R"

rule STEP33_merge_footprint_by_chr_downsampled:
    input:
        "sites/{gene}.sites.Rdata",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr1.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr2.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr3.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr4.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr5.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr6.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr7.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr8.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr9.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr10.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr11.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr12.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr13.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr14.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr15.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr16.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr17.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr18.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr19.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr20.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr21.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chr22.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chrX.footprint.bychr.sat.done.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.chrY.footprint.bychr.sat.done.txt"
    output:
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.done.merged.txt"
    script:
        "scripts/saturation/snakeMergeFPbyChrDownsampled.R"

rule STEP34_parse_footprint_downsampled:
    # note that you are using the non-downsampled peaks here
    input:
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bam",
        "{path}saturation/downsampled/{mergedsample}-REP{repnum}of{reptot}.{prob}.md.bai",
        "sites/{gene}.sites.Rdata",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.done.merged.txt",
        "{path}peaks/macs2/individual/{mergedsample}-REP{repnum}of{reptot}_global_normalization_peaks.narrowPeak",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.done.merged.txt" 
    output:
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.{prob}.{gene}.parsed.finished.txt"
    resources:
        parse_fp=1
    script:
        "scripts/saturation/snakeParseFPDownsampled.R"

rule AGGREGATOR_saturation_footprints:
    input:
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.9.{gene}.parsed.finished.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.8.{gene}.parsed.finished.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.7.{gene}.parsed.finished.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.6.{gene}.parsed.finished.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.5.{gene}.parsed.finished.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.4.{gene}.parsed.finished.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.3.{gene}.parsed.finished.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.2.{gene}.parsed.finished.txt",
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.1.{gene}.parsed.finished.txt"
    output:
        "{path}operations/footsat/{mergedsample}-REP{repnum}of{reptot}.allprob.{gene}.done.parsed.txt"
    shell:
        "touch {output}"

########################################################################################################################################
#### METRICS AND ANNOTATIONS RULES #####################################################################################################
########################################################################################################################################

rule METRICS_percent_peak_genome_coverage:
    # returns a fraction value of the basepairs of the genome covered by the merged peak file. multiple by 100 for percentages
    # parameters:
    # --echo output will be at least a three-column bed file
    # --bases-uniq the number of distinct bases from ref covered by overlap bed file
    # --delim change output delimeter from '|' to <delim>, e.g. '\t'
    input:
        a="{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/{mergedsample}.peak.genome.coverage.txt"
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"

rule METRICS_fragment_size_distributions:
    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP{repnum}of{reptot}.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP{repnum}of{reptot}.u.bai"
    output:
        "{path}metrics/{mergedsample}-REP{repnum}of{reptot}.u.fragsizes.svg"
    script:
        "scripts/snakeFragSizeDist.R"

rule AGGREGATOR_fragsize_1rep:
    input:
        "{path}metrics/{mergedsample}-REP1of1.u.fragsizes.svg"
    output:
        "{path}operations/{mergedsample}.fragsizes.done.txt"
    shell:
        "touch {output}"

rule AGGREGATOR_fragsize_2reps:
    input:
        "{path}metrics/{mergedsample}-REP1of2.u.fragsizes.svg",
        "{path}metrics/{mergedsample}-REP2of2.u.fragsizes.svg"
    output:
        "{path}operations/{mergedsample}.fragsizes.done.txt"
    shell:
        "touch {output}"

rule AGGREGATOR_fragsize_3reps:
    input:
        "{path}metrics/{mergedsample}-REP1of3.u.fragsizes.svg",
        "{path}metrics/{mergedsample}-REP2of3.u.fragsizes.svg",
        "{path}metrics/{mergedsample}-REP3of3.u.fragsizes.svg"
    output:
        "{path}operations/{mergedsample}.fragsizes.done.txt"
    shell:
        "touch {output}"

rule METRICS_annotate_peaks_merged:
    input:
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
    output:
        "{path}operations/{mergedsample}.mergedpeak.annotations.done.txt"
    script:
        "scripts/snakeAnnotateATAC.R"

rule AGGREGATOR_metrics_and_annotations:
    input:
        "{path}metrics/{mergedsample}.peak.genome.coverage.txt",
        "{path}operations/{mergedsample}.fragsizes.done.txt",
        "{path}operations/{mergedsample}.mergedpeak.annotations.done.txt"
    output:
        "{path}operations/{mergedsample}.metrics.annotations.done.txt"
    shell:
        "touch {output}"

########################################################################################################################################
#### XSAMPLE ANALYSIS RULES ############################################################################################################
########################################################################################################################################

rule xsample_footprint_direct_comparison:
    input:
        a="{path1}footprints/parsed/{sample1}.{gene}.motif{motif}.info.Rdata",
        b="{path2}footprints/parsed/{sample2}.{gene}.motif{motif}.info.Rdata"
    output:
        "xsample/footprints/{sample1}.{sample2}.{gene}.{motif}.xfootprint.Rdata"
    script:
        "script/snakeXsampleCompareFootprint.R"

########################################################################################################################################
#### PAN TF FOOTPRINTING ANALYSIS ######################################################################################################
########################################################################################################################################

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

## Note - this section utilizes rules defined in an auxillary snakefile called 'panTF.snakefile'

## Spooling commands ###################################################################################################################
# Run this with a terminal command like: for i in {1..62}; do snakemake --config group=$i -j 20 run_pantf_ls1034wt01; done
rule run_pantf_h508wt01:
    input:
        expand("h508/wt01/footprints/operations/H508-WT-01.parseTF.group{param}.done", param=config["group"])

rule run_pantf_ls1034wt01:
    input:
        expand("ls1034/wt01/footprints/operations/LS1034-WT-01.parseTF.group{param}.done", param=config["group"])


# You can also use this rule to run everything at once
# For some reason, this runs MUCH slower/gets stuck. I don't know why
# Will need to spend some time troubleshooting at some point
# For now, use the other method 
#rule run_pantf_ls1034wt01:
#    input:
#        "ls1034/wt01/footprints/operations/LS1034-WT-01.parseFP.allgroups.done"



## Pipeline rules ###################################################################################################################
rule PANTF_run_group:
	input:
		"h508/wt01/footprints/operations/H508-WT-01.parseTF.group{config.group}.done"

rule PANTF_copy_bam:
    # The TF analysis script runs in 20 simultaneous processes
    # Each process will need to access the bam file individually
    # To significantly speed this analysis up, temporarily make 20 copies of the bam file
    # And assign each individual process a unique file to access
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam"
    output:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.{bamcopy}.bam"
    shell:
        "cp {input} {output}"

rule PANTF_copy_bai:
    # The TF analysis script runs in 20 simultaneous processes
    # Each process will need to access the bam file individually
    # To significantly speed this analysis up, temporarily make 20 copies of the bam file
    # And assign each individual process a unique file to access
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai"
    output:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.{bamcopy}.bai"
    shell:
        "cp {input} {output}"

## The 'raw' footprint analysis involves pulling the reads from the bam files and generating insertion matrices
rule PANTF_raw_footprint_analysis:
    input:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.{bamcopy}.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.{bamcopy}.bai",
        "sites/data/{gene}.bindingSites.Rdata",
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.bamcopy.done"
    output:
        "{path}footprints/operations/{mergedsample}.{gene}.rawFPanalysis.bamcopy{bamcopy}.done"
    resources:
        analyzeRawFP=1
    benchmark:
        '{path}footprints/benchmark/{mergedsample}.{gene}.rawFPanalysis.bamcopy{bamcopy}.txt'
    script:
        "scripts/panTF/snakeAnalyzeRawFootprint.R"

## Parsing the raw footprints involves identifying which genomic loci have a TF bound
rule PANTF_parse_footprint_analysis:
    input:
        "{path}footprints/operations/{mergedsample}.{gene}.rawFPanalysis.bamcopy{bamcopy}.done"
    output:
    	"{path}footprints/operations/{mergedsample}.{gene}.parseFP.bamcopy{bamcopy}.done"
    resources:
        parseFootprint=1
    benchmark:
        '{path}footprints/benchmark/{mergedsample}.{gene}.bamcopy{bamcopy}.parseFP.txt'
    script:
    	"scripts/panTF/snakeParseFootprint.R"

## Remove the extra copies of the bam files once they are no longer needed
rule PANTF_remove_bamcopy:
    input:
        "{path}footprints/operations/{mergedsample}.rawTF.allgroups.done"
    output:
        "{path}footprints/operations/{mergedsample}.rawTF.analysis.done"
    shell:
         """
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bam
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bai
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bamcopy.done
         touch {output}
         """


########################################################################################################################################
#### CREATE LOCAL PWM SCAN DATABASE ####################################################################################################
########################################################################################################################################

## Note - this section utilizes rules defined in an auxillary snakefile

rule generate_motifData:
    output:
        "sites/motifData.Rdata"
    script:
        "scripts/scanPWM/generateMotifData.R"

rule generate_geneNames:
    output:
        "sites/geneNames.txt"
    script:
        "scripts/scanPWM/generateNames.R"

rule scanPWM:
    input:
        "sites/motifData.Rdata"
    output:
        'sites/operations/scans/{gene}.PWMscan.done'
    resources:
        scanPWM=1
    script:
        'scripts/scanPWM/snakeScanPWM.R'
