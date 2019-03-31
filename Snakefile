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
        mkdir -p -v {wildcards.path}saturation
        mkdir -p -v {wildcards.path}saturation/footprints
        mkdir -p -v {wildcards.path}saturation/footprints/data {wildcards.path}saturation/footprints/graphs 
        mkdir -p -v {wildcards.path}saturation/footprints/data/merged {wildcards.path}saturation/footprints/data/motifmerge {wildcards.path}saturation/footprints/data/parsed {wildcards.path}saturation/footprints/data/bychr
        mkdir -p -v {wildcards.path}saturation/complexity
        mkdir -p -v {wildcards.path}saturation/peaks
        mkdir -p -v {wildcards.path}saturation/downsampled
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}footprints/graphs {wildcards.path}footprints/heatmaps {wildcards.path}footprints/data
        mkdir -p -v {wildcards.path}footprints/data/merged {wildcards.path}footprints/data/motifmerge {wildcards.path}footprints/parsed {wildcards.path}footprints/data/bychr {wildcards.path}footprints/operations
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}peaks/genrich {wildcards.path}peaks/macs2 {wildcards.path}peaks/macs2/individual {wildcards.path}peaks/macs2/merged
        mkdir -p -v {wildcards.path}correlation
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}operations {wildcards.path}operations/footsat
        mkdir -p -v {wildcards.path}pantf
        mkdir -p -v {wildcards.path}pantf/parsed
        mkdir -p -v {wildcards.path}pantf/operations
        mkdir -p -v {wildcards.path}pantf/graphs
        mkdir -p -v {wildcards.path}pantf/data {wildcards.path}pantf/data/bychr {wildcards.path}pantf/data/merged {wildcards.path}pantf/data/motifmerge
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

rule snu61wt01_pantf_analysis:
    input:
        "snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.analysis.done"

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

rule AGGREGATOR_copy_bam:
    input:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.1.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.2.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.3.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.4.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.5.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.6.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.7.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.8.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.9.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.10.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.11.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.12.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.13.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.14.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.15.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.16.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.17.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.18.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.19.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.20.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.1.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.2.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.3.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.4.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.5.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.6.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.7.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.8.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.9.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.10.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.11.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.12.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.13.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.14.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.15.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.16.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.17.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.18.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.19.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.20.bai"
    output:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.bamcopy.done"
    shell:
        "touch {output}"

rule PANTF_analyze_footprint:
    input:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.{bamcopy}.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.{bamcopy}.bai",
        "sites/data/{gene}.bindingSites.Rdata",
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.bamcopy.done"
    output:
        "{path}footprints/operations/{mergedsample}.{gene}.rawFPanalysis.bamcopy{bamcopy}.done"
    resources:
        analyzeRawTF=1
    script:
        "scripts/panTF/snakeAnalyzeFootprint.R"

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

# rule PANTF_parse_footprint:
#     input:
#         "{path}footprints/operations/{mergedsample}.{gene}.rawFPanalysis.done"
#     output:

#     script:



# rule test:
#     input:
#         "snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.analysis.done"


rule PANTF_TFgroup_aggregator:
    input:
        '{path}footprints/operations/{mergedsample}.rawTF.group1.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group2.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group3.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group4.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group5.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group6.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group7.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group8.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group9.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group10.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group11.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group12.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group13.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group14.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group15.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group16.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group17.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group18.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group19.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group20.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group21.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group22.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group23.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group24.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group25.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group26.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group27.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group28.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group29.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group30.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group31.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group32.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group33.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group34.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group35.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group36.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group37.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group38.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group39.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group40.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group41.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group42.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group43.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group44.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group45.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group46.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group47.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group48.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group49.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group50.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group51.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group52.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group53.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group54.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group55.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group56.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group57.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group58.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group59.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group60.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group61.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group62.done'
    output:
        "{path}footprints/operations/{mergedsample}.rawTF.allgroups.done"
    shell:
        "touch {output}"

rule rawTF_group1:
    input:
        '{path}footprints/operations/{mergedsample}.TFAP2A.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.NFIL3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.HLF.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.NHLH1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.MAX.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.USF1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.CEBPA.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.EBF1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.CEBPB.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.FOS.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.FOSL1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.FOSL2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.JUN.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.JUNB.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.JUND.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.MAFF.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.MAFK.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.TFAP2C.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.USF2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.SREBF1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group1.done'
    shell:
        'touch {output}'
rule rawTF_group2:
    input:
        '{path}footprints/operations/{mergedsample}.SREBF2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.AHR.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.TFAP4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ARNT.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ATF6.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.BACH1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.BACH2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.CREB1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ATF2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.TCF3.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.XBP1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.ARID5B.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.MYOD1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.NFE2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.MYCN.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.NFE2L1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.TEF.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ATF3.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.BATF.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.TCF12.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group2.done'
    shell:
        'touch {output}'
rule rawTF_group3:
    input:
        '{path}footprints/operations/{mergedsample}.MYC.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.MXI1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.BHLHE40.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ARNTL.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ATF4.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ATF7.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.BATF3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.BHLHA15.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.BHLHE41.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.BHLHE22.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.BHLHE23.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.CEBPD.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.CEBPE.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.CEBPG.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.CLOCK.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.CREB3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.CREB3L1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.DBP.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.FIGLA.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.HES5.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group3.done'
    shell:
        'touch {output}'
rule rawTF_group4:
    input:
        '{path}footprints/operations/{mergedsample}.HES7.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.HEY1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.HEY2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ID4.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.JDP2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.MAFG.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.MESP1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.MGA.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.MLX.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.MLXIPL.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.MNT.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.MSC.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.MYF6.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.NEUROD2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.NEUROG2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.NRL.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.OLIG1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.OLIG2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.OLIG3.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.TCF4.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group4.done'
    shell:
        'touch {output}'
rule rawTF_group5:
    input:
        '{path}footprints/operations/{mergedsample}.TFAP2B.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.TFE3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.TFEB.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.TFEC.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.TFAP2D.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ARID3A.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ARNT2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ATF1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ATF5.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.CREM.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.DDIT3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.EPAS1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.FOSB.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.HAND1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.HES1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.HIF1A.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HMGA1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.HMGA2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.MAFA.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.MAFB.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group5.done'
    shell:
        'touch {output}'
rule rawTF_group6:
    input:
        '{path}footprints/operations/{mergedsample}.MAF.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.MITF.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.MYOG.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.NEUROD1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.NFE2L2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.PTF1A.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.TAL1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.TWIST1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.AIRE.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ALX1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ALX3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.ALX4.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ANDR.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.AP2A.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.AP2B.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.AP2C.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.AP2D.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ARI3A.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ARI5B.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ARX.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group6.done'
    shell:
        'touch {output}'
rule rawTF_group7:
    input:
        '{path}footprints/operations/{mergedsample}.ASCL2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ATF6A.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ATOH1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.BARH1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.BARH2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.BARX1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.BARX2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.BC11A.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.BCL6B.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.BCL6.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.BHA15.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.BHE22.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.BHE23.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.BHE40.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.BHE41.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.BMAL1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.BPTF.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.BRAC.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.BRCA1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.BSH.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group7.done'
    shell:
        'touch {output}'
rule rawTF_group8:
    input:
        '{path}footprints/operations/{mergedsample}.CDC5L.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.CDX1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.CDX2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.CEBPZ.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.CENPB.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.COE1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.COT1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.COT2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.CPEB1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.CR3L1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.CR3L2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.CREB5.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.CRX.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.CTCFL.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.CTCF.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.CUX1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.CUX2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.CXXC1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.DLX1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.DLX2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group8.done'
    shell:
        'touch {output}'
rule rawTF_group9:
    input:
        '{path}footprints/operations/{mergedsample}.DLX3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.DLX4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.DLX5.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.DLX6.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.DMBX1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.DPRX.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.DRGX.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.DUXA.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.E2F1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.E2F2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.E2F3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.E2F4.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.E2F5.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.E2F6.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.E2F7.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.E2F8.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.E4F1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.EGR1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.EGR2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.EGR3.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group9.done'
    shell:
        'touch {output}'
rule rawTF_group10:
    input:
        '{path}footprints/operations/{mergedsample}.EGR4.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.EHF.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ELF1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ELF2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ELF3.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ELF5.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ELK1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ELK3.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ELK4.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.EMX1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.EMX2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.EOMES.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ERF.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.ERG.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.ERR1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.ERR2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.ERR3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ESR1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ESR2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ESX1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group10.done'
    shell:
        'touch {output}'
rule rawTF_group11:
    input:
        '{path}footprints/operations/{mergedsample}.ETS1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ETS2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ETV1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ETV2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ETV3.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ETV4.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ETV5.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ETV6.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ETV7.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.EVI1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.EVX1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.EVX2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.FEV.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.FLI1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.FOXA1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.FOXA2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.FOXA3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.FOXB1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.FOXC1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.FOXC2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group11.done'
    shell:
        'touch {output}'
rule rawTF_group12:
    input:
        '{path}footprints/operations/{mergedsample}.FOXD1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.FOXD2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.FOXD3.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.FOXF1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.FOXF2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.FOXG1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.FOXH1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.FOXI1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.FOXJ2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.FOXJ3.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.FOXK1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.FOXL1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.FOXM1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.FOXO1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.FOXO3.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.FOXO4.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.FOXO6.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.FOXP2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.FOXP3.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.FOXQ1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group12.done'
    shell:
        'touch {output}'
rule rawTF_group13:
    input:
        '{path}footprints/operations/{mergedsample}.FUBP1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.GABP1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.GABPA.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.GATA1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.GATA2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.GATA3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.GATA4.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.GATA5.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.GATA6.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.GBX1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.GBX2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.GCM1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.GCM2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.GCR.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.GFI1B.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.GFI1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.GLI1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.GLI2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.GLI3.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.GLIS1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group13.done'
    shell:
        'touch {output}'
rule rawTF_group14:
    input:
        '{path}footprints/operations/{mergedsample}.GLIS2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.GLIS3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.GMEB2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.GRHL1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.GSC2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.GSC.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.GSX1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.GSX2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.HBP1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.HEN1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.HESX1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.HIC1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.HIC2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.HINFP.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.HLTF.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.HMBX1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HME1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.HME2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.HMX1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.HMX2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group14.done'
    shell:
        'touch {output}'
rule rawTF_group15:
    input:
        '{path}footprints/operations/{mergedsample}.HMX3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.HNF1A.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.HNF1B.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.HNF4A.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.HNF4G.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.HNF6.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.HOMEZ.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.HSF1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.HSF2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.HSF4.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.HSFY1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.HTF4.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.HXA10.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.HXA11.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.HXA13.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.HXA1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HXA2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.HXA5.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.HXA7.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.HXA9.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group15.done'
    shell:
        'touch {output}'
rule rawTF_group16:
    input:
        '{path}footprints/operations/{mergedsample}.HXB13.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.HXB1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.HXB2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.HXB3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.HXB6.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.HXB7.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.HXB8.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.HXC10.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.HXC11.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.HXC12.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.HXC13.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.HXC6.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.HXC8.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.HXD10.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.HXD11.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.HXD12.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HXD13.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.HXD3.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.HXD4.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.HXD8.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group16.done'
    shell:
        'touch {output}'
rule rawTF_group17:
    input:
        '{path}footprints/operations/{mergedsample}.HXD9.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.IKZF1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.INSM1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.IRF1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.IRF2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.IRF3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.IRF4.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.IRF5.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.IRF7.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.IRF8.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.IRF9.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.IRX2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.IRX3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.ISL1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.ISL2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.ISX.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.ITF2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.KAISO.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.KLF12.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.KLF13.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group17.done'
    shell:
        'touch {output}'
rule rawTF_group18:
    input:
        '{path}footprints/operations/{mergedsample}.KLF14.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.KLF15.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.KLF16.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.KLF1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.KLF3.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.KLF4.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.KLF6.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.KLF8.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.LBX2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.LEF1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.LHX2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.LHX3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.LHX4.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.LHX6.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.LHX8.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.LHX9.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.LMX1A.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.LMX1B.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.MAZ.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.MBD2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group18.done'
    shell:
        'touch {output}'
rule rawTF_group19:
    input:
        '{path}footprints/operations/{mergedsample}.MCR.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.MECP2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.MEF2A.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.MEF2B.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.MEF2C.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.MEF2D.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.MEIS1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.MEIS2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.MEIS3.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.MEOX1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.MEOX2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.MGAP.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.MIXL1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.MLXPL.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.MNX1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.MSX1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.MSX2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.MTF1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.MUSC.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.MYBA.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group19.done'
    shell:
        'touch {output}'
rule rawTF_group20:
    input:
        '{path}footprints/operations/{mergedsample}.MYBB.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.MYB.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.MZF1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.NANOG.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.NDF1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.NDF2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.NF2L1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.NF2L2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.NFAC1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.NFAC2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.NFAC3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.NFAC4.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.NFAT5.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.NFIA.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.NFIC.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.NFKB1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.NFKB2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.NFYA.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.NFYB.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.NFYC.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group20.done'
    shell:
        'touch {output}'
rule rawTF_group21:
    input:
        '{path}footprints/operations/{mergedsample}.NGN2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.NKX21.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.NKX22.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.NKX23.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.NKX25.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.NKX28.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.NKX31.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.NKX32.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.NKX61.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.NKX62.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.NOBOX.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.NOTO.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.NR0B1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.NR1D1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.NR1H2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.NR1H4.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.NR1I2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.NR1I3.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.NR2C1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.NR2C2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group21.done'
    shell:
        'touch {output}'
rule rawTF_group22:
    input:
        '{path}footprints/operations/{mergedsample}.NR2E1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.NR2E3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.NR2F6.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.NR4A1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.NR4A2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.NR4A3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.NR5A2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.NR6A1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.NRF1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ONEC2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ONEC3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.OTX1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.OTX2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.OVOL1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.P53.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.P5F1B.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.P63.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.P73.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.PAX1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.PAX2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group22.done'
    shell:
        'touch {output}'
rule rawTF_group23:
    input:
        '{path}footprints/operations/{mergedsample}.PAX3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.PAX4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.PAX5.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.PAX6.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.PAX7.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.PAX8.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.PBX1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.PBX2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.PBX3.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.PDX1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.PEBB.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.PHX2A.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.PHX2B.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.PIT1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.PITX1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.PITX2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.PITX3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.PKNX1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.PKNX2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.PLAG1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group23.done'
    shell:
        'touch {output}'
rule rawTF_group24:
    input:
        '{path}footprints/operations/{mergedsample}.PLAL1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.PO2F1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.PO2F2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.PO2F3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.PO3F1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.PO3F2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.PO3F3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.PO3F4.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.PO4F1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.PO4F2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.PO4F3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.PO5F1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.PO6F1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.PO6F2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.PPARA.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.PPARD.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.PPARG.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.PRD14.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.PRDM1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.PRDM4.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group24.done'
    shell:
        'touch {output}'
rule rawTF_group25:
    input:
        '{path}footprints/operations/{mergedsample}.PRGR.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.PROP1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.PROX1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.PRRX1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.PRRX2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.PURA.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.RARA.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.RARB.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.RARG.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.RAX2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.RELB.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.REL.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.REST.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.RFX1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.RFX2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.RFX3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.RFX4.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.RFX5.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.RHXF1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.RORA.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group25.done'
    shell:
        'touch {output}'
rule rawTF_group26:
    input:
        '{path}footprints/operations/{mergedsample}.RORG.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.RREB1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.RUNX1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.RUNX2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.RUNX3.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.RXRA.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.RXRB.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.RXRG.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.RX.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.SCRT1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.SCRT2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.SHOX2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.SHOX.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.SMAD1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.SMAD2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.SMAD3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.SMAD4.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.SMRC1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.SNAI1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.SNAI2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group26.done'
    shell:
        'touch {output}'
rule rawTF_group27:
    input:
        '{path}footprints/operations/{mergedsample}.SOX10.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.SOX11.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.SOX13.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.SOX15.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.SOX17.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.SOX18.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.SOX1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.SOX21.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.SOX2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.SOX3.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.SOX4.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.SOX5.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.SOX7.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.SOX8.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.SOX9.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.SP1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.SP2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.SP3.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.SP4.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.SPDEF.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group27.done'
    shell:
        'touch {output}'
rule rawTF_group28:
    input:
        '{path}footprints/operations/{mergedsample}.SPI1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.SPIB.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.SPIC.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.SPZ1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.SRBP1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.SRBP2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.SRF.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.SRY.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.STA5A.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.STA5B.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.STAT1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.STAT2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.STAT3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.STAT4.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.STAT6.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.STF1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.SUH.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.TBP.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.TBR1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.TBX15.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group28.done'
    shell:
        'touch {output}'
rule rawTF_group29:
    input:
        '{path}footprints/operations/{mergedsample}.TBX19.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.TBX1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.TBX20.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.TBX21.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.TBX2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.TBX3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.TBX4.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.TBX5.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.TCF7.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.TEAD1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.TEAD3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.TEAD4.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.TF2LX.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.TF65.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.TF7L1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.TF7L2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.TFCP2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.TFDP1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.TFE2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.TGIF1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group29.done'
    shell:
        'touch {output}'
rule rawTF_group30:
    input:
        '{path}footprints/operations/{mergedsample}.TGIF2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.THAP1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.THA.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.THB.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.TLX1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.TWST1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.TYY1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.TYY2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.UBIP1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.UNC4.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.VAX1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.VAX2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.VDR.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.VENTX.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.VSX1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.VSX2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.WT1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.YBOX1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ZBED1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ZBT18.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group30.done'
    shell:
        'touch {output}'
rule rawTF_group31:
    input:
        '{path}footprints/operations/{mergedsample}.ZBT49.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ZBT7A.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ZBT7B.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB4.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB6.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ZEB1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ZEP1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ZEP2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ZFHX3.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ZFX.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ZIC1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.ZIC2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ZIC3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.ZIC4.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.ZKSC1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.ZKSC3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.ZN143.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ZN148.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ZN219.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ZN232.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group31.done'
    shell:
        'touch {output}'
rule rawTF_group32:
    input:
        '{path}footprints/operations/{mergedsample}.ZN282.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ZN333.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ZN350.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ZN384.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ZN410.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ZN423.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ZN524.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ZN589.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ZN639.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ZN652.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ZN713.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.ZN740.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ZN784.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.ZSC16.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.ZSCA4.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.ABCF2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.A1CF.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ACO1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ADARB1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.AFF4.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group32.done'
    shell:
        'touch {output}'
rule rawTF_group33:
    input:
        '{path}footprints/operations/{mergedsample}.AGGF1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.AKR1A1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ANXA1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ANXA11.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.APEX2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ARFGAP1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ASCC1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ASPSCR1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.AVEN.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.BAD.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.GPANK1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.BAX.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.BCL11A.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.BOLL.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.CELF4.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.CELF5.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.CELF6.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.C19orf25.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.C19orf40.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.EXO5.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group33.done'
    shell:
        'touch {output}'
rule rawTF_group34:
    input:
        '{path}footprints/operations/{mergedsample}.LINC00471.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.C9orf156.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.CANX.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.CAT.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.CBFA2T2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.CBFB.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.CBX7.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ZNF830.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.CCDC25.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.CD59.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.CDK2AP1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.AGAP2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.CFL2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.FOXN3.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.CKMT1B.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.CLK1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.CNOT6.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.NELFB.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.CPSF4.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.CSNK2B.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group34.done'
    shell:
        'touch {output}'
rule rawTF_group35:
    input:
        '{path}footprints/operations/{mergedsample}.CSTF2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.CYB5R1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.CYCS.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.DAB2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.DAZAP1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ASAP3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.DDX20.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.DDX4.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.DDX43.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.DDX53.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.DGCR8.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.DHX36.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.DIABLO.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.DIS3.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.DNMT3A.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.DTL.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.DUS3L.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.DUSP22.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.DUSP26.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ECSIT.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group35.done'
    shell:
        'touch {output}'
rule rawTF_group36:
    input:
        '{path}footprints/operations/{mergedsample}.EDN1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.EEF1D.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.EIF5A2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ENO1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ESRRA.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ETFB.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.EWSR1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.EXOSC3.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.METTL21B.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.FAM127B.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.FEZ1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.FEZF2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.FGF19.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.FHL2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.FIP1L1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.SRRM3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.FOXP4.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.GADD45A.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.GIT2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.GLYCTK.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group36.done'
    shell:
        'touch {output}'
rule rawTF_group37:
    input:
        '{path}footprints/operations/{mergedsample}.GOT1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.GPAM.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.GPD1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.GRHPR.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.GTF2B.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.GTF2H3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.GTF3C2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.GTF3C5.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.GTPBP1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.GTPBP6.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.H1FX.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.H2AFY.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.H2AFZ.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.HCFC2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.HCLS1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.HDAC8.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HHAT.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.HHEX.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.UBE2K.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.HIRIP3.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group37.done'
    shell:
        'touch {output}'
rule rawTF_group38:
    input:
        '{path}footprints/operations/{mergedsample}.HIST1H2BN.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.HIST2H2AB.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.HIST2H2BE.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.HLCS.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.HMG20A.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.HNRNPA0.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.HNRNPA1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.HNRNPC.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.HNRNPH3.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.HNRNPLL.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.HOXB13.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.HOXB9.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.HOXD3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.HP1BP3.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.HSPA1L.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.HSPA5.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HTATIP2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ID2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.IL24.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ING3.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group38.done'
    shell:
        'touch {output}'
rule rawTF_group39:
    input:
        '{path}footprints/operations/{mergedsample}.IRF6.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.IVD.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.KDM5A.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.KDM5D.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.KCNIP1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.KIAA0907.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.KIF22.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.LARP1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.LARP4.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.LAS1L.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.CERS4.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.UBXN1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.CBX3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.LRRFIP1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.LSM6.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.LUZP1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.LUZP2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.MAGEA8.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.MAGED4B.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.MAGEF1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group39.done'
    shell:
        'touch {output}'
rule rawTF_group40:
    input:
        '{path}footprints/operations/{mergedsample}.MAGOH.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.MAP4K2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.MAPK1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.MBTPS2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.MCTP2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.MDM2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.GLTPD1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.RBM42.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.WDR83.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.MORN1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.MRPL1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.MRPL2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.MRPS25.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.MSI1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.MSI2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.MSRA.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.MSRB3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.MTHFD1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.MXD4.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.MYEF2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group40.done'
    shell:
        'touch {output}'
rule rawTF_group41:
    input:
        '{path}footprints/operations/{mergedsample}.MYLK.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.NANOS1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.NAP1L1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.NCALD.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.NCBP2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.NFATC3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.NFATC4.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.NFIB.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.NFIX.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.NME1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.NMI.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.NMRAL1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.NNT.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.NOC2L.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.GAR1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.NONO.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.NR2F1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.NUCB1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.NUP107.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.NUP133.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group41.done'
    shell:
        'touch {output}'
rule rawTF_group42:
    input:
        '{path}footprints/operations/{mergedsample}.NXPH3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ODC1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.OTUD4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.P4HB.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.PAXIP1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.PCK2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.PDCD11.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.PDE6H.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.PDLIM5.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.PGAM2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.PHLDA2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.PHOX2A.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.PHTF1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.PICK1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.PIK3C3.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.PIR.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.PKM.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.PKNOX2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.PLAGL1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.PLG.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group42.done'
    shell:
        'touch {output}'
rule rawTF_group43:
    input:
        '{path}footprints/operations/{mergedsample}.POLE3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.POLI.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.POU3F2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.POU4F3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.PPP1R10.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.PPP2R3B.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.PPP5C.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.PQBP1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.PRDX5.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.PRKRIR.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.PRNP.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.PSMA6.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.PSMC2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.PTCD1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.PTPMT1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.PURG.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.R3HDM2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.RAB14.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.RAB18.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.RAB2A.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group43.done'
    shell:
        'touch {output}'
rule rawTF_group44:
    input:
        '{path}footprints/operations/{mergedsample}.RAB7A.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.RAN.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.RAX.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.RBBP5.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.RBBP9.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.RBM17.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.RBM22.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.RBM3.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ESRP1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ESRP2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.RBM7.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.RBM8A.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.RBFOX2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.RBMS1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.RFC2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.RFC3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.RFXANK.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.RIOK2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.MEX3C.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.RNASEH2C.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group44.done'
    shell:
        'touch {output}'
rule rawTF_group45:
    input:
        '{path}footprints/operations/{mergedsample}.RNF138.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.RPL35.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.RPL6.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.RPP25.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.RPS10.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.RPS4X.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.RPS6KA5.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.RUFY3.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.RUVBL1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.SCAND2P.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.PDS5A.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.SCMH1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.SEMA4A.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.SF1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.SF3B1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.SFT2D1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.SLC18A1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.SMAP2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.SMCR7L.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.SMPX.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group45.done'
    shell:
        'touch {output}'
rule rawTF_group46:
    input:
        '{path}footprints/operations/{mergedsample}.SMUG1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.SNAPC4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.SNAPC5.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.SND1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.SNRNP70.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.SNRPB2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.SOCS4.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.SOD1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.SOX14.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.SPAG7.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.SPATS2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.SPR.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.SRBD1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.SRP9.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.SSBP3.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.SSX2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.SSX3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.STAU2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.STUB1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.SUCLG1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group46.done'
    shell:
        'touch {output}'
rule rawTF_group47:
    input:
        '{path}footprints/operations/{mergedsample}.TAF1A.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.TAF9.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.TAGLN2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.TBPL1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.TCEAL2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.TCEAL6.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.TFAM.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.TGIF2LX.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.THAP5.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.THRA.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.MED30.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.TIA1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.TIMELESS.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.TIMM44.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.TIMM8A.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.TMSB4XP8.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.TOB2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.TP73.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.TPI1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.TPPP.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group47.done'
    shell:
        'touch {output}'
rule rawTF_group48:
    input:
        '{path}footprints/operations/{mergedsample}.TRIM21.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.TRIM69.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.TRIP10.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.TRMT1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.TROVE2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.TSC22D4.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.TSN.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.TSNAX.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.TULP1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.U2AF1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.UBB.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.UBE2V1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.UGP2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.UQCRB.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.USP39.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.UTP18.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.VAMP3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.EZR.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.VPS4B.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.NELFA.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group48.done'
    shell:
        'touch {output}'
rule rawTF_group49:
    input:
        '{path}footprints/operations/{mergedsample}.WISP2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.XG.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.XRCC1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.YEATS4.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.YWHAE.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.YWHAZ.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB12.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB25.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB43.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB46.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ZC3H7A.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.ZCCHC14.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ZCCHC17.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.ZDHHC15.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.ZDHHC5.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.ZFP3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.ZHX3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ZMAT2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ZMAT4.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ZNF124.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group49.done'
    shell:
        'touch {output}'
rule rawTF_group50:
    input:
        '{path}footprints/operations/{mergedsample}.ZNF131.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ZNF160.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ZKSCAN8.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ZSCAN9.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ZNF205.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ZNF207.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB18.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ZNF250.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ZNF26.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ZNF3.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ZNF304.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.RNF114.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ZSCAN31.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.ZNF326.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.ZNF385A.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.ZNF503.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.ZNF510.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ZNF655.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ZNF671.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ZNF695.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group50.done'
    shell:
        'touch {output}'
rule rawTF_group51:
    input:
        '{path}footprints/operations/{mergedsample}.ZNF706.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ZNF71.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ZNF720.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ZNF76.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.ZNF766.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.ZRSR2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.ZSWIM1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.Myf.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.Pax6.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.RORA_1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.RORA_2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.YY1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.TP53.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.RELA.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.ZNF354C.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.MIZF.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.AP1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.DUX4.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.FOXP1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.POU2F2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group51.done'
    shell:
        'touch {output}'
rule rawTF_group52:
    input:
        '{path}footprints/operations/{mergedsample}.TCF7L2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.TP63.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB33.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ZNF263.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.AR.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.KLF5.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.T.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.EN1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ZNF143.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.NR3C1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ESRRB.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.HOXA5.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.DMRT3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.LBX1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.POU6F1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.BARHL2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.ELF4.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.EN2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.HOXA13.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.HOXC11.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group52.done'
    shell:
        'touch {output}'
rule rawTF_group53:
    input:
        '{path}footprints/operations/{mergedsample}.ONECUT1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.POU4F2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB7B.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB7C.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.RHOXF1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.UNCX.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.NR3C2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.SP8.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.YY2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB7A.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ZNF410.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.ZNF740.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ONECUT2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.ONECUT3.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.MYBL1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.MYBL2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.PAX9.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.PKNOX1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.POU1F1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.POU2F1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group53.done'
    shell:
        'touch {output}'
rule rawTF_group54:
    input:
        '{path}footprints/operations/{mergedsample}.POU3F1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.POU3F3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.POU3F4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.POU4F1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.POU5F1B.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.POU6F2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.HOXD12.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.BSX.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.HMBOX1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.HOXA10.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.HOXA2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.HOXB2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.HOXB3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.HOXC10.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.HOXC12.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.HOXC13.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HOXD11.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.HOXD13.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.NFATC2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ASCL1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group54.done'
    shell:
        'touch {output}'
rule rawTF_group55:
    input:
        '{path}footprints/operations/{mergedsample}.FOXK2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.GRHL2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.KLF9.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.NR2F2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.POU5F1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.RBPJ.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.SIX1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.SIX2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.TEAD2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ZNF24.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ZNF384.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.ZNF282.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ZSCAN4.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.RORB.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.RORC.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.TCF7L1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HINFP1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ZNF238.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ZNF306.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ZNF524.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group55.done'
    shell:
        'touch {output}'
rule rawTF_group56:
    input:
        '{path}footprints/operations/{mergedsample}.ZNF75A.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ZNF784.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.HSFY2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.NFATC1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.POU2F3.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.POU5F1P1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.BHLHB2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.BHLHB3.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.CART1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.HOXA1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.HOXB5.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.HOXD8.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.IRX5.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.PHOX2B.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.RAXL1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.ESRRG.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.THRB.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.Trp53.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.Trp73.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB49.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group56.done'
    shell:
        'touch {output}'
rule rawTF_group57:
    input:
        '{path}footprints/operations/{mergedsample}.ZNF232.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ZNF435.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ZNF713.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ARID5A.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.BARHL1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.BBX.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.BCL3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.CHD1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.CHD2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.CREB3L2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.DBX2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.DMC1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.EBF3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.EP300.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.EZH2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.FOXJ1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.FOXN1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.GMEB1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.GTF2F1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.GTF2I.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group57.done'
    shell:
        'touch {output}'
rule rawTF_group58:
    input:
        '{path}footprints/operations/{mergedsample}.GZF1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.HCFC1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.HDX.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.HIVEP1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.HLX.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.HOXA11.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.HOXA3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.HOXA4.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.HOXA6.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.HOXA7.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.HOXA9.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.HOXB1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.HOXB4.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.HOXB6.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.HOXB7.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.HOXB8.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.HOXC4.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.HOXC5.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.HOXC6.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.HOXC8.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group58.done'
    shell:
        'touch {output}'
rule rawTF_group59:
    input:
        '{path}footprints/operations/{mergedsample}.HOXC9.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.HOXD10.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.HOXD1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.HOXD4.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.HOXD9.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.IKZF2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.IRX4.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.IRX6.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.KLF7.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.LHX1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.LHX5.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.MECOM.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.MTA3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.OSR1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.OSR2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.OTP.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.PATZ1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.PGR.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.PML.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.PRDM14.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group59.done'
    shell:
        'touch {output}'
rule rawTF_group60:
    input:
        '{path}footprints/operations/{mergedsample}.RAD21.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.RCOR1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.RFX7.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.RHOXF2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.SIN3A.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.SIX3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.SIX4.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.SIX5.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.SIX6.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.SMARCC1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.SMARCC2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.SMC3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.SOX12.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.SOX30.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.SOX6.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.SP100.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.STAT5A.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.STAT5B.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.TAF1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.TBL1XR1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group60.done'
    shell:
        'touch {output}'
rule rawTF_group61:
    input:
        '{path}footprints/operations/{mergedsample}.TCF21.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.TFAP2E.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.TFCP2L1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.TLX2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.UBP1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/{mergedsample}.WRNIP1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/{mergedsample}.YBX1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB14.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB16.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/{mergedsample}.ZBTB3.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/{mergedsample}.ZKSCAN1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/{mergedsample}.ZKSCAN3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/{mergedsample}.ZNF148.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/{mergedsample}.ZNF219.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/{mergedsample}.ZNF274.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/{mergedsample}.ZNF281.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/{mergedsample}.ZNF333.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/{mergedsample}.ZNF350.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/{mergedsample}.ZNF35.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/{mergedsample}.ZNF423.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group61.done'
    shell:
        'touch {output}'
rule rawTF_group62:
    input:
        '{path}footprints/operations/{mergedsample}.ZNF652.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/{mergedsample}.ZNF691.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/{mergedsample}.ZNF711.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/{mergedsample}.ZNF8.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/{mergedsample}.Sox4.rawFPanalysis.bamcopy5.done'
    output:
        '{path}footprints/operations/{mergedsample}.rawTF.group62.done'
    shell:
        'touch {output}'



########################################################################################################################################
#### CREATE LOCAL PWM SCAN DATABASE ####################################################################################################
########################################################################################################################################

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

rule PWMscan_group_aggregator:
    input:
        'sites/operations/groups/PWMscan.group1.done',
        'sites/operations/groups/PWMscan.group2.done',
        'sites/operations/groups/PWMscan.group3.done',
        'sites/operations/groups/PWMscan.group4.done',
        'sites/operations/groups/PWMscan.group5.done',
        'sites/operations/groups/PWMscan.group6.done',
        'sites/operations/groups/PWMscan.group7.done',
        'sites/operations/groups/PWMscan.group8.done',
        'sites/operations/groups/PWMscan.group9.done',
        'sites/operations/groups/PWMscan.group10.done',
        'sites/operations/groups/PWMscan.group11.done',
        'sites/operations/groups/PWMscan.group12.done',
        'sites/operations/groups/PWMscan.group13.done',
        'sites/operations/groups/PWMscan.group14.done',
        'sites/operations/groups/PWMscan.group15.done',
        'sites/operations/groups/PWMscan.group16.done',
        'sites/operations/groups/PWMscan.group17.done',
        'sites/operations/groups/PWMscan.group18.done',
        'sites/operations/groups/PWMscan.group19.done',
        'sites/operations/groups/PWMscan.group20.done',
        'sites/operations/groups/PWMscan.group21.done',
        'sites/operations/groups/PWMscan.group22.done',
        'sites/operations/groups/PWMscan.group23.done',
        'sites/operations/groups/PWMscan.group24.done',
        'sites/operations/groups/PWMscan.group25.done',
        'sites/operations/groups/PWMscan.group26.done',
        'sites/operations/groups/PWMscan.group27.done',
        'sites/operations/groups/PWMscan.group28.done',
        'sites/operations/groups/PWMscan.group29.done',
        'sites/operations/groups/PWMscan.group30.done',
        'sites/operations/groups/PWMscan.group31.done',
        'sites/operations/groups/PWMscan.group32.done',
        'sites/operations/groups/PWMscan.group33.done',
        'sites/operations/groups/PWMscan.group34.done',
        'sites/operations/groups/PWMscan.group35.done',
        'sites/operations/groups/PWMscan.group36.done',
        'sites/operations/groups/PWMscan.group37.done',
        'sites/operations/groups/PWMscan.group38.done',
        'sites/operations/groups/PWMscan.group39.done',
        'sites/operations/groups/PWMscan.group40.done',
        'sites/operations/groups/PWMscan.group41.done',
        'sites/operations/groups/PWMscan.group42.done',
        'sites/operations/groups/PWMscan.group43.done',
        'sites/operations/groups/PWMscan.group44.done',
        'sites/operations/groups/PWMscan.group45.done',
        'sites/operations/groups/PWMscan.group46.done',
        'sites/operations/groups/PWMscan.group47.done',
        'sites/operations/groups/PWMscan.group48.done',
        'sites/operations/groups/PWMscan.group49.done',
        'sites/operations/groups/PWMscan.group50.done',
        'sites/operations/groups/PWMscan.group51.done',
        'sites/operations/groups/PWMscan.group52.done',
        'sites/operations/groups/PWMscan.group53.done',
        'sites/operations/groups/PWMscan.group54.done',
        'sites/operations/groups/PWMscan.group55.done',
        'sites/operations/groups/PWMscan.group56.done',
        'sites/operations/groups/PWMscan.group57.done',
        'sites/operations/groups/PWMscan.group58.done',
        'sites/operations/groups/PWMscan.group59.done',
        'sites/operations/groups/PWMscan.group60.done',
        'sites/operations/groups/PWMscan.group61.done',
        'sites/operations/groups/PWMscan.group62.done'
    output:
        "sites/operations/groups/PWMscan.allgroups.done"
    shell:
        "touch {output}"
rule PWMscan_group1:
    input:
        'sites/operations/scans/TFAP2A.PWMscan.done', 
        'sites/operations/scans/NFIL3.PWMscan.done', 
        'sites/operations/scans/HLF.PWMscan.done', 
        'sites/operations/scans/NHLH1.PWMscan.done', 
        'sites/operations/scans/MAX.PWMscan.done', 
        'sites/operations/scans/USF1.PWMscan.done', 
        'sites/operations/scans/CEBPA.PWMscan.done', 
        'sites/operations/scans/EBF1.PWMscan.done', 
        'sites/operations/scans/CEBPB.PWMscan.done', 
        'sites/operations/scans/FOS.PWMscan.done', 
        'sites/operations/scans/FOSL1.PWMscan.done', 
        'sites/operations/scans/FOSL2.PWMscan.done', 
        'sites/operations/scans/JUN.PWMscan.done', 
        'sites/operations/scans/JUNB.PWMscan.done', 
        'sites/operations/scans/JUND.PWMscan.done', 
        'sites/operations/scans/MAFF.PWMscan.done', 
        'sites/operations/scans/MAFK.PWMscan.done', 
        'sites/operations/scans/TFAP2C.PWMscan.done', 
        'sites/operations/scans/USF2.PWMscan.done', 
        'sites/operations/scans/SREBF1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group1.done'
    shell:
        'touch {output}'
rule PWMscan_group2:
    input:
        'sites/operations/scans/SREBF2.PWMscan.done', 
        'sites/operations/scans/AHR.PWMscan.done', 
        'sites/operations/scans/TFAP4.PWMscan.done', 
        'sites/operations/scans/ARNT.PWMscan.done', 
        'sites/operations/scans/ATF6.PWMscan.done', 
        'sites/operations/scans/BACH1.PWMscan.done', 
        'sites/operations/scans/BACH2.PWMscan.done', 
        'sites/operations/scans/CREB1.PWMscan.done', 
        'sites/operations/scans/ATF2.PWMscan.done', 
        'sites/operations/scans/TCF3.PWMscan.done', 
        'sites/operations/scans/XBP1.PWMscan.done', 
        'sites/operations/scans/ARID5B.PWMscan.done', 
        'sites/operations/scans/MYOD1.PWMscan.done', 
        'sites/operations/scans/NFE2.PWMscan.done', 
        'sites/operations/scans/MYCN.PWMscan.done', 
        'sites/operations/scans/NFE2L1.PWMscan.done', 
        'sites/operations/scans/TEF.PWMscan.done', 
        'sites/operations/scans/ATF3.PWMscan.done', 
        'sites/operations/scans/BATF.PWMscan.done', 
        'sites/operations/scans/TCF12.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group2.done'
    shell:
        'touch {output}'
rule PWMscan_group3:
    input:
        'sites/operations/scans/MYC.PWMscan.done', 
        'sites/operations/scans/MXI1.PWMscan.done', 
        'sites/operations/scans/BHLHE40.PWMscan.done', 
        'sites/operations/scans/ARNTL.PWMscan.done', 
        'sites/operations/scans/ATF4.PWMscan.done', 
        'sites/operations/scans/ATF7.PWMscan.done', 
        'sites/operations/scans/BATF3.PWMscan.done', 
        'sites/operations/scans/BHLHA15.PWMscan.done', 
        'sites/operations/scans/BHLHE41.PWMscan.done', 
        'sites/operations/scans/BHLHE22.PWMscan.done', 
        'sites/operations/scans/BHLHE23.PWMscan.done', 
        'sites/operations/scans/CEBPD.PWMscan.done', 
        'sites/operations/scans/CEBPE.PWMscan.done', 
        'sites/operations/scans/CEBPG.PWMscan.done', 
        'sites/operations/scans/CLOCK.PWMscan.done', 
        'sites/operations/scans/CREB3.PWMscan.done', 
        'sites/operations/scans/CREB3L1.PWMscan.done', 
        'sites/operations/scans/DBP.PWMscan.done', 
        'sites/operations/scans/FIGLA.PWMscan.done', 
        'sites/operations/scans/HES5.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group3.done'
    shell:
        'touch {output}'
rule PWMscan_group4:
    input:
        'sites/operations/scans/HES7.PWMscan.done', 
        'sites/operations/scans/HEY1.PWMscan.done', 
        'sites/operations/scans/HEY2.PWMscan.done', 
        'sites/operations/scans/ID4.PWMscan.done', 
        'sites/operations/scans/JDP2.PWMscan.done', 
        'sites/operations/scans/MAFG.PWMscan.done', 
        'sites/operations/scans/MESP1.PWMscan.done', 
        'sites/operations/scans/MGA.PWMscan.done', 
        'sites/operations/scans/MLX.PWMscan.done', 
        'sites/operations/scans/MLXIPL.PWMscan.done', 
        'sites/operations/scans/MNT.PWMscan.done', 
        'sites/operations/scans/MSC.PWMscan.done', 
        'sites/operations/scans/MYF6.PWMscan.done', 
        'sites/operations/scans/NEUROD2.PWMscan.done', 
        'sites/operations/scans/NEUROG2.PWMscan.done', 
        'sites/operations/scans/NRL.PWMscan.done', 
        'sites/operations/scans/OLIG1.PWMscan.done', 
        'sites/operations/scans/OLIG2.PWMscan.done', 
        'sites/operations/scans/OLIG3.PWMscan.done', 
        'sites/operations/scans/TCF4.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group4.done'
    shell:
        'touch {output}'
rule PWMscan_group5:
    input:
        'sites/operations/scans/TFAP2B.PWMscan.done', 
        'sites/operations/scans/TFE3.PWMscan.done', 
        'sites/operations/scans/TFEB.PWMscan.done', 
        'sites/operations/scans/TFEC.PWMscan.done', 
        'sites/operations/scans/TFAP2D.PWMscan.done', 
        'sites/operations/scans/ARID3A.PWMscan.done', 
        'sites/operations/scans/ARNT2.PWMscan.done', 
        'sites/operations/scans/ATF1.PWMscan.done', 
        'sites/operations/scans/ATF5.PWMscan.done', 
        'sites/operations/scans/CREM.PWMscan.done', 
        'sites/operations/scans/DDIT3.PWMscan.done', 
        'sites/operations/scans/EPAS1.PWMscan.done', 
        'sites/operations/scans/FOSB.PWMscan.done', 
        'sites/operations/scans/HAND1.PWMscan.done', 
        'sites/operations/scans/HES1.PWMscan.done', 
        'sites/operations/scans/HIF1A.PWMscan.done', 
        'sites/operations/scans/HMGA1.PWMscan.done', 
        'sites/operations/scans/HMGA2.PWMscan.done', 
        'sites/operations/scans/MAFA.PWMscan.done', 
        'sites/operations/scans/MAFB.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group5.done'
    shell:
        'touch {output}'
rule PWMscan_group6:
    input:
        'sites/operations/scans/MAF.PWMscan.done', 
        'sites/operations/scans/MITF.PWMscan.done', 
        'sites/operations/scans/MYOG.PWMscan.done', 
        'sites/operations/scans/NEUROD1.PWMscan.done', 
        'sites/operations/scans/NFE2L2.PWMscan.done', 
        'sites/operations/scans/PTF1A.PWMscan.done', 
        'sites/operations/scans/TAL1.PWMscan.done', 
        'sites/operations/scans/TWIST1.PWMscan.done', 
        'sites/operations/scans/AIRE.PWMscan.done', 
        'sites/operations/scans/ALX1.PWMscan.done', 
        'sites/operations/scans/ALX3.PWMscan.done', 
        'sites/operations/scans/ALX4.PWMscan.done', 
        'sites/operations/scans/ANDR.PWMscan.done', 
        'sites/operations/scans/AP2A.PWMscan.done', 
        'sites/operations/scans/AP2B.PWMscan.done', 
        'sites/operations/scans/AP2C.PWMscan.done', 
        'sites/operations/scans/AP2D.PWMscan.done', 
        'sites/operations/scans/ARI3A.PWMscan.done', 
        'sites/operations/scans/ARI5B.PWMscan.done', 
        'sites/operations/scans/ARX.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group6.done'
    shell:
        'touch {output}'
rule PWMscan_group7:
    input:
        'sites/operations/scans/ASCL2.PWMscan.done', 
        'sites/operations/scans/ATF6A.PWMscan.done', 
        'sites/operations/scans/ATOH1.PWMscan.done', 
        'sites/operations/scans/BARH1.PWMscan.done', 
        'sites/operations/scans/BARH2.PWMscan.done', 
        'sites/operations/scans/BARX1.PWMscan.done', 
        'sites/operations/scans/BARX2.PWMscan.done', 
        'sites/operations/scans/BC11A.PWMscan.done', 
        'sites/operations/scans/BCL6B.PWMscan.done', 
        'sites/operations/scans/BCL6.PWMscan.done', 
        'sites/operations/scans/BHA15.PWMscan.done', 
        'sites/operations/scans/BHE22.PWMscan.done', 
        'sites/operations/scans/BHE23.PWMscan.done', 
        'sites/operations/scans/BHE40.PWMscan.done', 
        'sites/operations/scans/BHE41.PWMscan.done', 
        'sites/operations/scans/BMAL1.PWMscan.done', 
        'sites/operations/scans/BPTF.PWMscan.done', 
        'sites/operations/scans/BRAC.PWMscan.done', 
        'sites/operations/scans/BRCA1.PWMscan.done', 
        'sites/operations/scans/BSH.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group7.done'
    shell:
        'touch {output}'
rule PWMscan_group8:
    input:
        'sites/operations/scans/CDC5L.PWMscan.done', 
        'sites/operations/scans/CDX1.PWMscan.done', 
        'sites/operations/scans/CDX2.PWMscan.done', 
        'sites/operations/scans/CEBPZ.PWMscan.done', 
        'sites/operations/scans/CENPB.PWMscan.done', 
        'sites/operations/scans/COE1.PWMscan.done', 
        'sites/operations/scans/COT1.PWMscan.done', 
        'sites/operations/scans/COT2.PWMscan.done', 
        'sites/operations/scans/CPEB1.PWMscan.done', 
        'sites/operations/scans/CR3L1.PWMscan.done', 
        'sites/operations/scans/CR3L2.PWMscan.done', 
        'sites/operations/scans/CREB5.PWMscan.done', 
        'sites/operations/scans/CRX.PWMscan.done', 
        'sites/operations/scans/CTCFL.PWMscan.done', 
        'sites/operations/scans/CTCF.PWMscan.done', 
        'sites/operations/scans/CUX1.PWMscan.done', 
        'sites/operations/scans/CUX2.PWMscan.done', 
        'sites/operations/scans/CXXC1.PWMscan.done', 
        'sites/operations/scans/DLX1.PWMscan.done', 
        'sites/operations/scans/DLX2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group8.done'
    shell:
        'touch {output}'
rule PWMscan_group9:
    input:
        'sites/operations/scans/DLX3.PWMscan.done', 
        'sites/operations/scans/DLX4.PWMscan.done', 
        'sites/operations/scans/DLX5.PWMscan.done', 
        'sites/operations/scans/DLX6.PWMscan.done', 
        'sites/operations/scans/DMBX1.PWMscan.done', 
        'sites/operations/scans/DPRX.PWMscan.done', 
        'sites/operations/scans/DRGX.PWMscan.done', 
        'sites/operations/scans/DUXA.PWMscan.done', 
        'sites/operations/scans/E2F1.PWMscan.done', 
        'sites/operations/scans/E2F2.PWMscan.done', 
        'sites/operations/scans/E2F3.PWMscan.done', 
        'sites/operations/scans/E2F4.PWMscan.done', 
        'sites/operations/scans/E2F5.PWMscan.done', 
        'sites/operations/scans/E2F6.PWMscan.done', 
        'sites/operations/scans/E2F7.PWMscan.done', 
        'sites/operations/scans/E2F8.PWMscan.done', 
        'sites/operations/scans/E4F1.PWMscan.done', 
        'sites/operations/scans/EGR1.PWMscan.done', 
        'sites/operations/scans/EGR2.PWMscan.done', 
        'sites/operations/scans/EGR3.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group9.done'
    shell:
        'touch {output}'
rule PWMscan_group10:
    input:
        'sites/operations/scans/EGR4.PWMscan.done', 
        'sites/operations/scans/EHF.PWMscan.done', 
        'sites/operations/scans/ELF1.PWMscan.done', 
        'sites/operations/scans/ELF2.PWMscan.done', 
        'sites/operations/scans/ELF3.PWMscan.done', 
        'sites/operations/scans/ELF5.PWMscan.done', 
        'sites/operations/scans/ELK1.PWMscan.done', 
        'sites/operations/scans/ELK3.PWMscan.done', 
        'sites/operations/scans/ELK4.PWMscan.done', 
        'sites/operations/scans/EMX1.PWMscan.done', 
        'sites/operations/scans/EMX2.PWMscan.done', 
        'sites/operations/scans/EOMES.PWMscan.done', 
        'sites/operations/scans/ERF.PWMscan.done', 
        'sites/operations/scans/ERG.PWMscan.done', 
        'sites/operations/scans/ERR1.PWMscan.done', 
        'sites/operations/scans/ERR2.PWMscan.done', 
        'sites/operations/scans/ERR3.PWMscan.done', 
        'sites/operations/scans/ESR1.PWMscan.done', 
        'sites/operations/scans/ESR2.PWMscan.done', 
        'sites/operations/scans/ESX1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group10.done'
    shell:
        'touch {output}'
rule PWMscan_group11:
    input:
        'sites/operations/scans/ETS1.PWMscan.done', 
        'sites/operations/scans/ETS2.PWMscan.done', 
        'sites/operations/scans/ETV1.PWMscan.done', 
        'sites/operations/scans/ETV2.PWMscan.done', 
        'sites/operations/scans/ETV3.PWMscan.done', 
        'sites/operations/scans/ETV4.PWMscan.done', 
        'sites/operations/scans/ETV5.PWMscan.done', 
        'sites/operations/scans/ETV6.PWMscan.done', 
        'sites/operations/scans/ETV7.PWMscan.done', 
        'sites/operations/scans/EVI1.PWMscan.done', 
        'sites/operations/scans/EVX1.PWMscan.done', 
        'sites/operations/scans/EVX2.PWMscan.done', 
        'sites/operations/scans/FEV.PWMscan.done', 
        'sites/operations/scans/FLI1.PWMscan.done', 
        'sites/operations/scans/FOXA1.PWMscan.done', 
        'sites/operations/scans/FOXA2.PWMscan.done', 
        'sites/operations/scans/FOXA3.PWMscan.done', 
        'sites/operations/scans/FOXB1.PWMscan.done', 
        'sites/operations/scans/FOXC1.PWMscan.done', 
        'sites/operations/scans/FOXC2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group11.done'
    shell:
        'touch {output}'
rule PWMscan_group12:
    input:
        'sites/operations/scans/FOXD1.PWMscan.done', 
        'sites/operations/scans/FOXD2.PWMscan.done', 
        'sites/operations/scans/FOXD3.PWMscan.done', 
        'sites/operations/scans/FOXF1.PWMscan.done', 
        'sites/operations/scans/FOXF2.PWMscan.done', 
        'sites/operations/scans/FOXG1.PWMscan.done', 
        'sites/operations/scans/FOXH1.PWMscan.done', 
        'sites/operations/scans/FOXI1.PWMscan.done', 
        'sites/operations/scans/FOXJ2.PWMscan.done', 
        'sites/operations/scans/FOXJ3.PWMscan.done', 
        'sites/operations/scans/FOXK1.PWMscan.done', 
        'sites/operations/scans/FOXL1.PWMscan.done', 
        'sites/operations/scans/FOXM1.PWMscan.done', 
        'sites/operations/scans/FOXO1.PWMscan.done', 
        'sites/operations/scans/FOXO3.PWMscan.done', 
        'sites/operations/scans/FOXO4.PWMscan.done', 
        'sites/operations/scans/FOXO6.PWMscan.done', 
        'sites/operations/scans/FOXP2.PWMscan.done', 
        'sites/operations/scans/FOXP3.PWMscan.done', 
        'sites/operations/scans/FOXQ1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group12.done'
    shell:
        'touch {output}'
rule PWMscan_group13:
    input:
        'sites/operations/scans/FUBP1.PWMscan.done', 
        'sites/operations/scans/GABP1.PWMscan.done', 
        'sites/operations/scans/GABPA.PWMscan.done', 
        'sites/operations/scans/GATA1.PWMscan.done', 
        'sites/operations/scans/GATA2.PWMscan.done', 
        'sites/operations/scans/GATA3.PWMscan.done', 
        'sites/operations/scans/GATA4.PWMscan.done', 
        'sites/operations/scans/GATA5.PWMscan.done', 
        'sites/operations/scans/GATA6.PWMscan.done', 
        'sites/operations/scans/GBX1.PWMscan.done', 
        'sites/operations/scans/GBX2.PWMscan.done', 
        'sites/operations/scans/GCM1.PWMscan.done', 
        'sites/operations/scans/GCM2.PWMscan.done', 
        'sites/operations/scans/GCR.PWMscan.done', 
        'sites/operations/scans/GFI1B.PWMscan.done', 
        'sites/operations/scans/GFI1.PWMscan.done', 
        'sites/operations/scans/GLI1.PWMscan.done', 
        'sites/operations/scans/GLI2.PWMscan.done', 
        'sites/operations/scans/GLI3.PWMscan.done', 
        'sites/operations/scans/GLIS1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group13.done'
    shell:
        'touch {output}'
rule PWMscan_group14:
    input:
        'sites/operations/scans/GLIS2.PWMscan.done', 
        'sites/operations/scans/GLIS3.PWMscan.done', 
        'sites/operations/scans/GMEB2.PWMscan.done', 
        'sites/operations/scans/GRHL1.PWMscan.done', 
        'sites/operations/scans/GSC2.PWMscan.done', 
        'sites/operations/scans/GSC.PWMscan.done', 
        'sites/operations/scans/GSX1.PWMscan.done', 
        'sites/operations/scans/GSX2.PWMscan.done', 
        'sites/operations/scans/HBP1.PWMscan.done', 
        'sites/operations/scans/HEN1.PWMscan.done', 
        'sites/operations/scans/HESX1.PWMscan.done', 
        'sites/operations/scans/HIC1.PWMscan.done', 
        'sites/operations/scans/HIC2.PWMscan.done', 
        'sites/operations/scans/HINFP.PWMscan.done', 
        'sites/operations/scans/HLTF.PWMscan.done', 
        'sites/operations/scans/HMBX1.PWMscan.done', 
        'sites/operations/scans/HME1.PWMscan.done', 
        'sites/operations/scans/HME2.PWMscan.done', 
        'sites/operations/scans/HMX1.PWMscan.done', 
        'sites/operations/scans/HMX2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group14.done'
    shell:
        'touch {output}'
rule PWMscan_group15:
    input:
        'sites/operations/scans/HMX3.PWMscan.done', 
        'sites/operations/scans/HNF1A.PWMscan.done', 
        'sites/operations/scans/HNF1B.PWMscan.done', 
        'sites/operations/scans/HNF4A.PWMscan.done', 
        'sites/operations/scans/HNF4G.PWMscan.done', 
        'sites/operations/scans/HNF6.PWMscan.done', 
        'sites/operations/scans/HOMEZ.PWMscan.done', 
        'sites/operations/scans/HSF1.PWMscan.done', 
        'sites/operations/scans/HSF2.PWMscan.done', 
        'sites/operations/scans/HSF4.PWMscan.done', 
        'sites/operations/scans/HSFY1.PWMscan.done', 
        'sites/operations/scans/HTF4.PWMscan.done', 
        'sites/operations/scans/HXA10.PWMscan.done', 
        'sites/operations/scans/HXA11.PWMscan.done', 
        'sites/operations/scans/HXA13.PWMscan.done', 
        'sites/operations/scans/HXA1.PWMscan.done', 
        'sites/operations/scans/HXA2.PWMscan.done', 
        'sites/operations/scans/HXA5.PWMscan.done', 
        'sites/operations/scans/HXA7.PWMscan.done', 
        'sites/operations/scans/HXA9.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group15.done'
    shell:
        'touch {output}'
rule PWMscan_group16:
    input:
        'sites/operations/scans/HXB13.PWMscan.done', 
        'sites/operations/scans/HXB1.PWMscan.done', 
        'sites/operations/scans/HXB2.PWMscan.done', 
        'sites/operations/scans/HXB3.PWMscan.done', 
        'sites/operations/scans/HXB6.PWMscan.done', 
        'sites/operations/scans/HXB7.PWMscan.done', 
        'sites/operations/scans/HXB8.PWMscan.done', 
        'sites/operations/scans/HXC10.PWMscan.done', 
        'sites/operations/scans/HXC11.PWMscan.done', 
        'sites/operations/scans/HXC12.PWMscan.done', 
        'sites/operations/scans/HXC13.PWMscan.done', 
        'sites/operations/scans/HXC6.PWMscan.done', 
        'sites/operations/scans/HXC8.PWMscan.done', 
        'sites/operations/scans/HXD10.PWMscan.done', 
        'sites/operations/scans/HXD11.PWMscan.done', 
        'sites/operations/scans/HXD12.PWMscan.done', 
        'sites/operations/scans/HXD13.PWMscan.done', 
        'sites/operations/scans/HXD3.PWMscan.done', 
        'sites/operations/scans/HXD4.PWMscan.done', 
        'sites/operations/scans/HXD8.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group16.done'
    shell:
        'touch {output}'
rule PWMscan_group17:
    input:
        'sites/operations/scans/HXD9.PWMscan.done', 
        'sites/operations/scans/IKZF1.PWMscan.done', 
        'sites/operations/scans/INSM1.PWMscan.done', 
        'sites/operations/scans/IRF1.PWMscan.done', 
        'sites/operations/scans/IRF2.PWMscan.done', 
        'sites/operations/scans/IRF3.PWMscan.done', 
        'sites/operations/scans/IRF4.PWMscan.done', 
        'sites/operations/scans/IRF5.PWMscan.done', 
        'sites/operations/scans/IRF7.PWMscan.done', 
        'sites/operations/scans/IRF8.PWMscan.done', 
        'sites/operations/scans/IRF9.PWMscan.done', 
        'sites/operations/scans/IRX2.PWMscan.done', 
        'sites/operations/scans/IRX3.PWMscan.done', 
        'sites/operations/scans/ISL1.PWMscan.done', 
        'sites/operations/scans/ISL2.PWMscan.done', 
        'sites/operations/scans/ISX.PWMscan.done', 
        'sites/operations/scans/ITF2.PWMscan.done', 
        'sites/operations/scans/KAISO.PWMscan.done', 
        'sites/operations/scans/KLF12.PWMscan.done', 
        'sites/operations/scans/KLF13.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group17.done'
    shell:
        'touch {output}'
rule PWMscan_group18:
    input:
        'sites/operations/scans/KLF14.PWMscan.done', 
        'sites/operations/scans/KLF15.PWMscan.done', 
        'sites/operations/scans/KLF16.PWMscan.done', 
        'sites/operations/scans/KLF1.PWMscan.done', 
        'sites/operations/scans/KLF3.PWMscan.done', 
        'sites/operations/scans/KLF4.PWMscan.done', 
        'sites/operations/scans/KLF6.PWMscan.done', 
        'sites/operations/scans/KLF8.PWMscan.done', 
        'sites/operations/scans/LBX2.PWMscan.done', 
        'sites/operations/scans/LEF1.PWMscan.done', 
        'sites/operations/scans/LHX2.PWMscan.done', 
        'sites/operations/scans/LHX3.PWMscan.done', 
        'sites/operations/scans/LHX4.PWMscan.done', 
        'sites/operations/scans/LHX6.PWMscan.done', 
        'sites/operations/scans/LHX8.PWMscan.done', 
        'sites/operations/scans/LHX9.PWMscan.done', 
        'sites/operations/scans/LMX1A.PWMscan.done', 
        'sites/operations/scans/LMX1B.PWMscan.done', 
        'sites/operations/scans/MAZ.PWMscan.done', 
        'sites/operations/scans/MBD2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group18.done'
    shell:
        'touch {output}'
rule PWMscan_group19:
    input:
        'sites/operations/scans/MCR.PWMscan.done', 
        'sites/operations/scans/MECP2.PWMscan.done', 
        'sites/operations/scans/MEF2A.PWMscan.done', 
        'sites/operations/scans/MEF2B.PWMscan.done', 
        'sites/operations/scans/MEF2C.PWMscan.done', 
        'sites/operations/scans/MEF2D.PWMscan.done', 
        'sites/operations/scans/MEIS1.PWMscan.done', 
        'sites/operations/scans/MEIS2.PWMscan.done', 
        'sites/operations/scans/MEIS3.PWMscan.done', 
        'sites/operations/scans/MEOX1.PWMscan.done', 
        'sites/operations/scans/MEOX2.PWMscan.done', 
        'sites/operations/scans/MGAP.PWMscan.done', 
        'sites/operations/scans/MIXL1.PWMscan.done', 
        'sites/operations/scans/MLXPL.PWMscan.done', 
        'sites/operations/scans/MNX1.PWMscan.done', 
        'sites/operations/scans/MSX1.PWMscan.done', 
        'sites/operations/scans/MSX2.PWMscan.done', 
        'sites/operations/scans/MTF1.PWMscan.done', 
        'sites/operations/scans/MUSC.PWMscan.done', 
        'sites/operations/scans/MYBA.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group19.done'
    shell:
        'touch {output}'
rule PWMscan_group20:
    input:
        'sites/operations/scans/MYBB.PWMscan.done', 
        'sites/operations/scans/MYB.PWMscan.done', 
        'sites/operations/scans/MZF1.PWMscan.done', 
        'sites/operations/scans/NANOG.PWMscan.done', 
        'sites/operations/scans/NDF1.PWMscan.done', 
        'sites/operations/scans/NDF2.PWMscan.done', 
        'sites/operations/scans/NF2L1.PWMscan.done', 
        'sites/operations/scans/NF2L2.PWMscan.done', 
        'sites/operations/scans/NFAC1.PWMscan.done', 
        'sites/operations/scans/NFAC2.PWMscan.done', 
        'sites/operations/scans/NFAC3.PWMscan.done', 
        'sites/operations/scans/NFAC4.PWMscan.done', 
        'sites/operations/scans/NFAT5.PWMscan.done', 
        'sites/operations/scans/NFIA.PWMscan.done', 
        'sites/operations/scans/NFIC.PWMscan.done', 
        'sites/operations/scans/NFKB1.PWMscan.done', 
        'sites/operations/scans/NFKB2.PWMscan.done', 
        'sites/operations/scans/NFYA.PWMscan.done', 
        'sites/operations/scans/NFYB.PWMscan.done', 
        'sites/operations/scans/NFYC.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group20.done'
    shell:
        'touch {output}'
rule PWMscan_group21:
    input:
        'sites/operations/scans/NGN2.PWMscan.done', 
        'sites/operations/scans/NKX21.PWMscan.done', 
        'sites/operations/scans/NKX22.PWMscan.done', 
        'sites/operations/scans/NKX23.PWMscan.done', 
        'sites/operations/scans/NKX25.PWMscan.done', 
        'sites/operations/scans/NKX28.PWMscan.done', 
        'sites/operations/scans/NKX31.PWMscan.done', 
        'sites/operations/scans/NKX32.PWMscan.done', 
        'sites/operations/scans/NKX61.PWMscan.done', 
        'sites/operations/scans/NKX62.PWMscan.done', 
        'sites/operations/scans/NOBOX.PWMscan.done', 
        'sites/operations/scans/NOTO.PWMscan.done', 
        'sites/operations/scans/NR0B1.PWMscan.done', 
        'sites/operations/scans/NR1D1.PWMscan.done', 
        'sites/operations/scans/NR1H2.PWMscan.done', 
        'sites/operations/scans/NR1H4.PWMscan.done', 
        'sites/operations/scans/NR1I2.PWMscan.done', 
        'sites/operations/scans/NR1I3.PWMscan.done', 
        'sites/operations/scans/NR2C1.PWMscan.done', 
        'sites/operations/scans/NR2C2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group21.done'
    shell:
        'touch {output}'
rule PWMscan_group22:
    input:
        'sites/operations/scans/NR2E1.PWMscan.done', 
        'sites/operations/scans/NR2E3.PWMscan.done', 
        'sites/operations/scans/NR2F6.PWMscan.done', 
        'sites/operations/scans/NR4A1.PWMscan.done', 
        'sites/operations/scans/NR4A2.PWMscan.done', 
        'sites/operations/scans/NR4A3.PWMscan.done', 
        'sites/operations/scans/NR5A2.PWMscan.done', 
        'sites/operations/scans/NR6A1.PWMscan.done', 
        'sites/operations/scans/NRF1.PWMscan.done', 
        'sites/operations/scans/ONEC2.PWMscan.done', 
        'sites/operations/scans/ONEC3.PWMscan.done', 
        'sites/operations/scans/OTX1.PWMscan.done', 
        'sites/operations/scans/OTX2.PWMscan.done', 
        'sites/operations/scans/OVOL1.PWMscan.done', 
        'sites/operations/scans/P53.PWMscan.done', 
        'sites/operations/scans/P5F1B.PWMscan.done', 
        'sites/operations/scans/P63.PWMscan.done', 
        'sites/operations/scans/P73.PWMscan.done', 
        'sites/operations/scans/PAX1.PWMscan.done', 
        'sites/operations/scans/PAX2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group22.done'
    shell:
        'touch {output}'
rule PWMscan_group23:
    input:
        'sites/operations/scans/PAX3.PWMscan.done', 
        'sites/operations/scans/PAX4.PWMscan.done', 
        'sites/operations/scans/PAX5.PWMscan.done', 
        'sites/operations/scans/PAX6.PWMscan.done', 
        'sites/operations/scans/PAX7.PWMscan.done', 
        'sites/operations/scans/PAX8.PWMscan.done', 
        'sites/operations/scans/PBX1.PWMscan.done', 
        'sites/operations/scans/PBX2.PWMscan.done', 
        'sites/operations/scans/PBX3.PWMscan.done', 
        'sites/operations/scans/PDX1.PWMscan.done', 
        'sites/operations/scans/PEBB.PWMscan.done', 
        'sites/operations/scans/PHX2A.PWMscan.done', 
        'sites/operations/scans/PHX2B.PWMscan.done', 
        'sites/operations/scans/PIT1.PWMscan.done', 
        'sites/operations/scans/PITX1.PWMscan.done', 
        'sites/operations/scans/PITX2.PWMscan.done', 
        'sites/operations/scans/PITX3.PWMscan.done', 
        'sites/operations/scans/PKNX1.PWMscan.done', 
        'sites/operations/scans/PKNX2.PWMscan.done', 
        'sites/operations/scans/PLAG1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group23.done'
    shell:
        'touch {output}'
rule PWMscan_group24:
    input:
        'sites/operations/scans/PLAL1.PWMscan.done', 
        'sites/operations/scans/PO2F1.PWMscan.done', 
        'sites/operations/scans/PO2F2.PWMscan.done', 
        'sites/operations/scans/PO2F3.PWMscan.done', 
        'sites/operations/scans/PO3F1.PWMscan.done', 
        'sites/operations/scans/PO3F2.PWMscan.done', 
        'sites/operations/scans/PO3F3.PWMscan.done', 
        'sites/operations/scans/PO3F4.PWMscan.done', 
        'sites/operations/scans/PO4F1.PWMscan.done', 
        'sites/operations/scans/PO4F2.PWMscan.done', 
        'sites/operations/scans/PO4F3.PWMscan.done', 
        'sites/operations/scans/PO5F1.PWMscan.done', 
        'sites/operations/scans/PO6F1.PWMscan.done', 
        'sites/operations/scans/PO6F2.PWMscan.done', 
        'sites/operations/scans/PPARA.PWMscan.done', 
        'sites/operations/scans/PPARD.PWMscan.done', 
        'sites/operations/scans/PPARG.PWMscan.done', 
        'sites/operations/scans/PRD14.PWMscan.done', 
        'sites/operations/scans/PRDM1.PWMscan.done', 
        'sites/operations/scans/PRDM4.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group24.done'
    shell:
        'touch {output}'
rule PWMscan_group25:
    input:
        'sites/operations/scans/PRGR.PWMscan.done', 
        'sites/operations/scans/PROP1.PWMscan.done', 
        'sites/operations/scans/PROX1.PWMscan.done', 
        'sites/operations/scans/PRRX1.PWMscan.done', 
        'sites/operations/scans/PRRX2.PWMscan.done', 
        'sites/operations/scans/PURA.PWMscan.done', 
        'sites/operations/scans/RARA.PWMscan.done', 
        'sites/operations/scans/RARB.PWMscan.done', 
        'sites/operations/scans/RARG.PWMscan.done', 
        'sites/operations/scans/RAX2.PWMscan.done', 
        'sites/operations/scans/RELB.PWMscan.done', 
        'sites/operations/scans/REL.PWMscan.done', 
        'sites/operations/scans/REST.PWMscan.done', 
        'sites/operations/scans/RFX1.PWMscan.done', 
        'sites/operations/scans/RFX2.PWMscan.done', 
        'sites/operations/scans/RFX3.PWMscan.done', 
        'sites/operations/scans/RFX4.PWMscan.done', 
        'sites/operations/scans/RFX5.PWMscan.done', 
        'sites/operations/scans/RHXF1.PWMscan.done', 
        'sites/operations/scans/RORA.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group25.done'
    shell:
        'touch {output}'
rule PWMscan_group26:
    input:
        'sites/operations/scans/RORG.PWMscan.done', 
        'sites/operations/scans/RREB1.PWMscan.done', 
        'sites/operations/scans/RUNX1.PWMscan.done', 
        'sites/operations/scans/RUNX2.PWMscan.done', 
        'sites/operations/scans/RUNX3.PWMscan.done', 
        'sites/operations/scans/RXRA.PWMscan.done', 
        'sites/operations/scans/RXRB.PWMscan.done', 
        'sites/operations/scans/RXRG.PWMscan.done', 
        'sites/operations/scans/RX.PWMscan.done', 
        'sites/operations/scans/SCRT1.PWMscan.done', 
        'sites/operations/scans/SCRT2.PWMscan.done', 
        'sites/operations/scans/SHOX2.PWMscan.done', 
        'sites/operations/scans/SHOX.PWMscan.done', 
        'sites/operations/scans/SMAD1.PWMscan.done', 
        'sites/operations/scans/SMAD2.PWMscan.done', 
        'sites/operations/scans/SMAD3.PWMscan.done', 
        'sites/operations/scans/SMAD4.PWMscan.done', 
        'sites/operations/scans/SMRC1.PWMscan.done', 
        'sites/operations/scans/SNAI1.PWMscan.done', 
        'sites/operations/scans/SNAI2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group26.done'
    shell:
        'touch {output}'
rule PWMscan_group27:
    input:
        'sites/operations/scans/SOX10.PWMscan.done', 
        'sites/operations/scans/SOX11.PWMscan.done', 
        'sites/operations/scans/SOX13.PWMscan.done', 
        'sites/operations/scans/SOX15.PWMscan.done', 
        'sites/operations/scans/SOX17.PWMscan.done', 
        'sites/operations/scans/SOX18.PWMscan.done', 
        'sites/operations/scans/SOX1.PWMscan.done', 
        'sites/operations/scans/SOX21.PWMscan.done', 
        'sites/operations/scans/SOX2.PWMscan.done', 
        'sites/operations/scans/SOX3.PWMscan.done', 
        'sites/operations/scans/SOX4.PWMscan.done', 
        'sites/operations/scans/SOX5.PWMscan.done', 
        'sites/operations/scans/SOX7.PWMscan.done', 
        'sites/operations/scans/SOX8.PWMscan.done', 
        'sites/operations/scans/SOX9.PWMscan.done', 
        'sites/operations/scans/SP1.PWMscan.done', 
        'sites/operations/scans/SP2.PWMscan.done', 
        'sites/operations/scans/SP3.PWMscan.done', 
        'sites/operations/scans/SP4.PWMscan.done', 
        'sites/operations/scans/SPDEF.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group27.done'
    shell:
        'touch {output}'
rule PWMscan_group28:
    input:
        'sites/operations/scans/SPI1.PWMscan.done', 
        'sites/operations/scans/SPIB.PWMscan.done', 
        'sites/operations/scans/SPIC.PWMscan.done', 
        'sites/operations/scans/SPZ1.PWMscan.done', 
        'sites/operations/scans/SRBP1.PWMscan.done', 
        'sites/operations/scans/SRBP2.PWMscan.done', 
        'sites/operations/scans/SRF.PWMscan.done', 
        'sites/operations/scans/SRY.PWMscan.done', 
        'sites/operations/scans/STA5A.PWMscan.done', 
        'sites/operations/scans/STA5B.PWMscan.done', 
        'sites/operations/scans/STAT1.PWMscan.done', 
        'sites/operations/scans/STAT2.PWMscan.done', 
        'sites/operations/scans/STAT3.PWMscan.done', 
        'sites/operations/scans/STAT4.PWMscan.done', 
        'sites/operations/scans/STAT6.PWMscan.done', 
        'sites/operations/scans/STF1.PWMscan.done', 
        'sites/operations/scans/SUH.PWMscan.done', 
        'sites/operations/scans/TBP.PWMscan.done', 
        'sites/operations/scans/TBR1.PWMscan.done', 
        'sites/operations/scans/TBX15.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group28.done'
    shell:
        'touch {output}'
rule PWMscan_group29:
    input:
        'sites/operations/scans/TBX19.PWMscan.done', 
        'sites/operations/scans/TBX1.PWMscan.done', 
        'sites/operations/scans/TBX20.PWMscan.done', 
        'sites/operations/scans/TBX21.PWMscan.done', 
        'sites/operations/scans/TBX2.PWMscan.done', 
        'sites/operations/scans/TBX3.PWMscan.done', 
        'sites/operations/scans/TBX4.PWMscan.done', 
        'sites/operations/scans/TBX5.PWMscan.done', 
        'sites/operations/scans/TCF7.PWMscan.done', 
        'sites/operations/scans/TEAD1.PWMscan.done', 
        'sites/operations/scans/TEAD3.PWMscan.done', 
        'sites/operations/scans/TEAD4.PWMscan.done', 
        'sites/operations/scans/TF2LX.PWMscan.done', 
        'sites/operations/scans/TF65.PWMscan.done', 
        'sites/operations/scans/TF7L1.PWMscan.done', 
        'sites/operations/scans/TF7L2.PWMscan.done', 
        'sites/operations/scans/TFCP2.PWMscan.done', 
        'sites/operations/scans/TFDP1.PWMscan.done', 
        'sites/operations/scans/TFE2.PWMscan.done', 
        'sites/operations/scans/TGIF1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group29.done'
    shell:
        'touch {output}'
rule PWMscan_group30:
    input:
        'sites/operations/scans/TGIF2.PWMscan.done', 
        'sites/operations/scans/THAP1.PWMscan.done', 
        'sites/operations/scans/THA.PWMscan.done', 
        'sites/operations/scans/THB.PWMscan.done', 
        'sites/operations/scans/TLX1.PWMscan.done', 
        'sites/operations/scans/TWST1.PWMscan.done', 
        'sites/operations/scans/TYY1.PWMscan.done', 
        'sites/operations/scans/TYY2.PWMscan.done', 
        'sites/operations/scans/UBIP1.PWMscan.done', 
        'sites/operations/scans/UNC4.PWMscan.done', 
        'sites/operations/scans/VAX1.PWMscan.done', 
        'sites/operations/scans/VAX2.PWMscan.done', 
        'sites/operations/scans/VDR.PWMscan.done', 
        'sites/operations/scans/VENTX.PWMscan.done', 
        'sites/operations/scans/VSX1.PWMscan.done', 
        'sites/operations/scans/VSX2.PWMscan.done', 
        'sites/operations/scans/WT1.PWMscan.done', 
        'sites/operations/scans/YBOX1.PWMscan.done', 
        'sites/operations/scans/ZBED1.PWMscan.done', 
        'sites/operations/scans/ZBT18.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group30.done'
    shell:
        'touch {output}'
rule PWMscan_group31:
    input:
        'sites/operations/scans/ZBT49.PWMscan.done', 
        'sites/operations/scans/ZBT7A.PWMscan.done', 
        'sites/operations/scans/ZBT7B.PWMscan.done', 
        'sites/operations/scans/ZBTB4.PWMscan.done', 
        'sites/operations/scans/ZBTB6.PWMscan.done', 
        'sites/operations/scans/ZEB1.PWMscan.done', 
        'sites/operations/scans/ZEP1.PWMscan.done', 
        'sites/operations/scans/ZEP2.PWMscan.done', 
        'sites/operations/scans/ZFHX3.PWMscan.done', 
        'sites/operations/scans/ZFX.PWMscan.done', 
        'sites/operations/scans/ZIC1.PWMscan.done', 
        'sites/operations/scans/ZIC2.PWMscan.done', 
        'sites/operations/scans/ZIC3.PWMscan.done', 
        'sites/operations/scans/ZIC4.PWMscan.done', 
        'sites/operations/scans/ZKSC1.PWMscan.done', 
        'sites/operations/scans/ZKSC3.PWMscan.done', 
        'sites/operations/scans/ZN143.PWMscan.done', 
        'sites/operations/scans/ZN148.PWMscan.done', 
        'sites/operations/scans/ZN219.PWMscan.done', 
        'sites/operations/scans/ZN232.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group31.done'
    shell:
        'touch {output}'
rule PWMscan_group32:
    input:
        'sites/operations/scans/ZN282.PWMscan.done', 
        'sites/operations/scans/ZN333.PWMscan.done', 
        'sites/operations/scans/ZN350.PWMscan.done', 
        'sites/operations/scans/ZN384.PWMscan.done', 
        'sites/operations/scans/ZN410.PWMscan.done', 
        'sites/operations/scans/ZN423.PWMscan.done', 
        'sites/operations/scans/ZN524.PWMscan.done', 
        'sites/operations/scans/ZN589.PWMscan.done', 
        'sites/operations/scans/ZN639.PWMscan.done', 
        'sites/operations/scans/ZN652.PWMscan.done', 
        'sites/operations/scans/ZN713.PWMscan.done', 
        'sites/operations/scans/ZN740.PWMscan.done', 
        'sites/operations/scans/ZN784.PWMscan.done', 
        'sites/operations/scans/ZSC16.PWMscan.done', 
        'sites/operations/scans/ZSCA4.PWMscan.done', 
        'sites/operations/scans/ABCF2.PWMscan.done', 
        'sites/operations/scans/A1CF.PWMscan.done', 
        'sites/operations/scans/ACO1.PWMscan.done', 
        'sites/operations/scans/ADARB1.PWMscan.done', 
        'sites/operations/scans/AFF4.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group32.done'
    shell:
        'touch {output}'
rule PWMscan_group33:
    input:
        'sites/operations/scans/AGGF1.PWMscan.done', 
        'sites/operations/scans/AKR1A1.PWMscan.done', 
        'sites/operations/scans/ANXA1.PWMscan.done', 
        'sites/operations/scans/ANXA11.PWMscan.done', 
        'sites/operations/scans/APEX2.PWMscan.done', 
        'sites/operations/scans/ARFGAP1.PWMscan.done', 
        'sites/operations/scans/ASCC1.PWMscan.done', 
        'sites/operations/scans/ASPSCR1.PWMscan.done', 
        'sites/operations/scans/AVEN.PWMscan.done', 
        'sites/operations/scans/BAD.PWMscan.done', 
        'sites/operations/scans/GPANK1.PWMscan.done', 
        'sites/operations/scans/BAX.PWMscan.done', 
        'sites/operations/scans/BCL11A.PWMscan.done', 
        'sites/operations/scans/BOLL.PWMscan.done', 
        'sites/operations/scans/CELF4.PWMscan.done', 
        'sites/operations/scans/CELF5.PWMscan.done', 
        'sites/operations/scans/CELF6.PWMscan.done', 
        'sites/operations/scans/C19orf25.PWMscan.done', 
        'sites/operations/scans/C19orf40.PWMscan.done', 
        'sites/operations/scans/EXO5.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group33.done'
    shell:
        'touch {output}'
rule PWMscan_group34:
    input:
        'sites/operations/scans/LINC00471.PWMscan.done', 
        'sites/operations/scans/C9orf156.PWMscan.done', 
        'sites/operations/scans/CANX.PWMscan.done', 
        'sites/operations/scans/CAT.PWMscan.done', 
        'sites/operations/scans/CBFA2T2.PWMscan.done', 
        'sites/operations/scans/CBFB.PWMscan.done', 
        'sites/operations/scans/CBX7.PWMscan.done', 
        'sites/operations/scans/ZNF830.PWMscan.done', 
        'sites/operations/scans/CCDC25.PWMscan.done', 
        'sites/operations/scans/CD59.PWMscan.done', 
        'sites/operations/scans/CDK2AP1.PWMscan.done', 
        'sites/operations/scans/AGAP2.PWMscan.done', 
        'sites/operations/scans/CFL2.PWMscan.done', 
        'sites/operations/scans/FOXN3.PWMscan.done', 
        'sites/operations/scans/CKMT1B.PWMscan.done', 
        'sites/operations/scans/CLK1.PWMscan.done', 
        'sites/operations/scans/CNOT6.PWMscan.done', 
        'sites/operations/scans/NELFB.PWMscan.done', 
        'sites/operations/scans/CPSF4.PWMscan.done', 
        'sites/operations/scans/CSNK2B.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group34.done'
    shell:
        'touch {output}'
rule PWMscan_group35:
    input:
        'sites/operations/scans/CSTF2.PWMscan.done', 
        'sites/operations/scans/CYB5R1.PWMscan.done', 
        'sites/operations/scans/CYCS.PWMscan.done', 
        'sites/operations/scans/DAB2.PWMscan.done', 
        'sites/operations/scans/DAZAP1.PWMscan.done', 
        'sites/operations/scans/ASAP3.PWMscan.done', 
        'sites/operations/scans/DDX20.PWMscan.done', 
        'sites/operations/scans/DDX4.PWMscan.done', 
        'sites/operations/scans/DDX43.PWMscan.done', 
        'sites/operations/scans/DDX53.PWMscan.done', 
        'sites/operations/scans/DGCR8.PWMscan.done', 
        'sites/operations/scans/DHX36.PWMscan.done', 
        'sites/operations/scans/DIABLO.PWMscan.done', 
        'sites/operations/scans/DIS3.PWMscan.done', 
        'sites/operations/scans/DNMT3A.PWMscan.done', 
        'sites/operations/scans/DTL.PWMscan.done', 
        'sites/operations/scans/DUS3L.PWMscan.done', 
        'sites/operations/scans/DUSP22.PWMscan.done', 
        'sites/operations/scans/DUSP26.PWMscan.done', 
        'sites/operations/scans/ECSIT.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group35.done'
    shell:
        'touch {output}'
rule PWMscan_group36:
    input:
        'sites/operations/scans/EDN1.PWMscan.done', 
        'sites/operations/scans/EEF1D.PWMscan.done', 
        'sites/operations/scans/EIF5A2.PWMscan.done', 
        'sites/operations/scans/ENO1.PWMscan.done', 
        'sites/operations/scans/ESRRA.PWMscan.done', 
        'sites/operations/scans/ETFB.PWMscan.done', 
        'sites/operations/scans/EWSR1.PWMscan.done', 
        'sites/operations/scans/EXOSC3.PWMscan.done', 
        'sites/operations/scans/METTL21B.PWMscan.done', 
        'sites/operations/scans/FAM127B.PWMscan.done', 
        'sites/operations/scans/FEZ1.PWMscan.done', 
        'sites/operations/scans/FEZF2.PWMscan.done', 
        'sites/operations/scans/FGF19.PWMscan.done', 
        'sites/operations/scans/FHL2.PWMscan.done', 
        'sites/operations/scans/FIP1L1.PWMscan.done', 
        'sites/operations/scans/SRRM3.PWMscan.done', 
        'sites/operations/scans/FOXP4.PWMscan.done', 
        'sites/operations/scans/GADD45A.PWMscan.done', 
        'sites/operations/scans/GIT2.PWMscan.done', 
        'sites/operations/scans/GLYCTK.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group36.done'
    shell:
        'touch {output}'
rule PWMscan_group37:
    input:
        'sites/operations/scans/GOT1.PWMscan.done', 
        'sites/operations/scans/GPAM.PWMscan.done', 
        'sites/operations/scans/GPD1.PWMscan.done', 
        'sites/operations/scans/GRHPR.PWMscan.done', 
        'sites/operations/scans/GTF2B.PWMscan.done', 
        'sites/operations/scans/GTF2H3.PWMscan.done', 
        'sites/operations/scans/GTF3C2.PWMscan.done', 
        'sites/operations/scans/GTF3C5.PWMscan.done', 
        'sites/operations/scans/GTPBP1.PWMscan.done', 
        'sites/operations/scans/GTPBP6.PWMscan.done', 
        'sites/operations/scans/H1FX.PWMscan.done', 
        'sites/operations/scans/H2AFY.PWMscan.done', 
        'sites/operations/scans/H2AFZ.PWMscan.done', 
        'sites/operations/scans/HCFC2.PWMscan.done', 
        'sites/operations/scans/HCLS1.PWMscan.done', 
        'sites/operations/scans/HDAC8.PWMscan.done', 
        'sites/operations/scans/HHAT.PWMscan.done', 
        'sites/operations/scans/HHEX.PWMscan.done', 
        'sites/operations/scans/UBE2K.PWMscan.done', 
        'sites/operations/scans/HIRIP3.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group37.done'
    shell:
        'touch {output}'
rule PWMscan_group38:
    input:
        'sites/operations/scans/HIST1H2BN.PWMscan.done', 
        'sites/operations/scans/HIST2H2AB.PWMscan.done', 
        'sites/operations/scans/HIST2H2BE.PWMscan.done', 
        'sites/operations/scans/HLCS.PWMscan.done', 
        'sites/operations/scans/HMG20A.PWMscan.done', 
        'sites/operations/scans/HNRNPA0.PWMscan.done', 
        'sites/operations/scans/HNRNPA1.PWMscan.done', 
        'sites/operations/scans/HNRNPC.PWMscan.done', 
        'sites/operations/scans/HNRNPH3.PWMscan.done', 
        'sites/operations/scans/HNRNPLL.PWMscan.done', 
        'sites/operations/scans/HOXB13.PWMscan.done', 
        'sites/operations/scans/HOXB9.PWMscan.done', 
        'sites/operations/scans/HOXD3.PWMscan.done', 
        'sites/operations/scans/HP1BP3.PWMscan.done', 
        'sites/operations/scans/HSPA1L.PWMscan.done', 
        'sites/operations/scans/HSPA5.PWMscan.done', 
        'sites/operations/scans/HTATIP2.PWMscan.done', 
        'sites/operations/scans/ID2.PWMscan.done', 
        'sites/operations/scans/IL24.PWMscan.done', 
        'sites/operations/scans/ING3.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group38.done'
    shell:
        'touch {output}'
rule PWMscan_group39:
    input:
        'sites/operations/scans/IRF6.PWMscan.done', 
        'sites/operations/scans/IVD.PWMscan.done', 
        'sites/operations/scans/KDM5A.PWMscan.done', 
        'sites/operations/scans/KDM5D.PWMscan.done', 
        'sites/operations/scans/KCNIP1.PWMscan.done', 
        'sites/operations/scans/KIAA0907.PWMscan.done', 
        'sites/operations/scans/KIF22.PWMscan.done', 
        'sites/operations/scans/LARP1.PWMscan.done', 
        'sites/operations/scans/LARP4.PWMscan.done', 
        'sites/operations/scans/LAS1L.PWMscan.done', 
        'sites/operations/scans/CERS4.PWMscan.done', 
        'sites/operations/scans/UBXN1.PWMscan.done', 
        'sites/operations/scans/CBX3.PWMscan.done', 
        'sites/operations/scans/LRRFIP1.PWMscan.done', 
        'sites/operations/scans/LSM6.PWMscan.done', 
        'sites/operations/scans/LUZP1.PWMscan.done', 
        'sites/operations/scans/LUZP2.PWMscan.done', 
        'sites/operations/scans/MAGEA8.PWMscan.done', 
        'sites/operations/scans/MAGED4B.PWMscan.done', 
        'sites/operations/scans/MAGEF1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group39.done'
    shell:
        'touch {output}'
rule PWMscan_group40:
    input:
        'sites/operations/scans/MAGOH.PWMscan.done', 
        'sites/operations/scans/MAP4K2.PWMscan.done', 
        'sites/operations/scans/MAPK1.PWMscan.done', 
        'sites/operations/scans/MBTPS2.PWMscan.done', 
        'sites/operations/scans/MCTP2.PWMscan.done', 
        'sites/operations/scans/MDM2.PWMscan.done', 
        'sites/operations/scans/GLTPD1.PWMscan.done', 
        'sites/operations/scans/RBM42.PWMscan.done', 
        'sites/operations/scans/WDR83.PWMscan.done', 
        'sites/operations/scans/MORN1.PWMscan.done', 
        'sites/operations/scans/MRPL1.PWMscan.done', 
        'sites/operations/scans/MRPL2.PWMscan.done', 
        'sites/operations/scans/MRPS25.PWMscan.done', 
        'sites/operations/scans/MSI1.PWMscan.done', 
        'sites/operations/scans/MSI2.PWMscan.done', 
        'sites/operations/scans/MSRA.PWMscan.done', 
        'sites/operations/scans/MSRB3.PWMscan.done', 
        'sites/operations/scans/MTHFD1.PWMscan.done', 
        'sites/operations/scans/MXD4.PWMscan.done', 
        'sites/operations/scans/MYEF2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group40.done'
    shell:
        'touch {output}'
rule PWMscan_group41:
    input:
        'sites/operations/scans/MYLK.PWMscan.done', 
        'sites/operations/scans/NANOS1.PWMscan.done', 
        'sites/operations/scans/NAP1L1.PWMscan.done', 
        'sites/operations/scans/NCALD.PWMscan.done', 
        'sites/operations/scans/NCBP2.PWMscan.done', 
        'sites/operations/scans/NFATC3.PWMscan.done', 
        'sites/operations/scans/NFATC4.PWMscan.done', 
        'sites/operations/scans/NFIB.PWMscan.done', 
        'sites/operations/scans/NFIX.PWMscan.done', 
        'sites/operations/scans/NME1.PWMscan.done', 
        'sites/operations/scans/NMI.PWMscan.done', 
        'sites/operations/scans/NMRAL1.PWMscan.done', 
        'sites/operations/scans/NNT.PWMscan.done', 
        'sites/operations/scans/NOC2L.PWMscan.done', 
        'sites/operations/scans/GAR1.PWMscan.done', 
        'sites/operations/scans/NONO.PWMscan.done', 
        'sites/operations/scans/NR2F1.PWMscan.done', 
        'sites/operations/scans/NUCB1.PWMscan.done', 
        'sites/operations/scans/NUP107.PWMscan.done', 
        'sites/operations/scans/NUP133.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group41.done'
    shell:
        'touch {output}'
rule PWMscan_group42:
    input:
        'sites/operations/scans/NXPH3.PWMscan.done', 
        'sites/operations/scans/ODC1.PWMscan.done', 
        'sites/operations/scans/OTUD4.PWMscan.done', 
        'sites/operations/scans/P4HB.PWMscan.done', 
        'sites/operations/scans/PAXIP1.PWMscan.done', 
        'sites/operations/scans/PCK2.PWMscan.done', 
        'sites/operations/scans/PDCD11.PWMscan.done', 
        'sites/operations/scans/PDE6H.PWMscan.done', 
        'sites/operations/scans/PDLIM5.PWMscan.done', 
        'sites/operations/scans/PGAM2.PWMscan.done', 
        'sites/operations/scans/PHLDA2.PWMscan.done', 
        'sites/operations/scans/PHOX2A.PWMscan.done', 
        'sites/operations/scans/PHTF1.PWMscan.done', 
        'sites/operations/scans/PICK1.PWMscan.done', 
        'sites/operations/scans/PIK3C3.PWMscan.done', 
        'sites/operations/scans/PIR.PWMscan.done', 
        'sites/operations/scans/PKM.PWMscan.done', 
        'sites/operations/scans/PKNOX2.PWMscan.done', 
        'sites/operations/scans/PLAGL1.PWMscan.done', 
        'sites/operations/scans/PLG.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group42.done'
    shell:
        'touch {output}'
rule PWMscan_group43:
    input:
        'sites/operations/scans/POLE3.PWMscan.done', 
        'sites/operations/scans/POLI.PWMscan.done', 
        'sites/operations/scans/POU3F2.PWMscan.done', 
        'sites/operations/scans/POU4F3.PWMscan.done', 
        'sites/operations/scans/PPP1R10.PWMscan.done', 
        'sites/operations/scans/PPP2R3B.PWMscan.done', 
        'sites/operations/scans/PPP5C.PWMscan.done', 
        'sites/operations/scans/PQBP1.PWMscan.done', 
        'sites/operations/scans/PRDX5.PWMscan.done', 
        'sites/operations/scans/PRKRIR.PWMscan.done', 
        'sites/operations/scans/PRNP.PWMscan.done', 
        'sites/operations/scans/PSMA6.PWMscan.done', 
        'sites/operations/scans/PSMC2.PWMscan.done', 
        'sites/operations/scans/PTCD1.PWMscan.done', 
        'sites/operations/scans/PTPMT1.PWMscan.done', 
        'sites/operations/scans/PURG.PWMscan.done', 
        'sites/operations/scans/R3HDM2.PWMscan.done', 
        'sites/operations/scans/RAB14.PWMscan.done', 
        'sites/operations/scans/RAB18.PWMscan.done', 
        'sites/operations/scans/RAB2A.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group43.done'
    shell:
        'touch {output}'
rule PWMscan_group44:
    input:
        'sites/operations/scans/RAB7A.PWMscan.done', 
        'sites/operations/scans/RAN.PWMscan.done', 
        'sites/operations/scans/RAX.PWMscan.done', 
        'sites/operations/scans/RBBP5.PWMscan.done', 
        'sites/operations/scans/RBBP9.PWMscan.done', 
        'sites/operations/scans/RBM17.PWMscan.done', 
        'sites/operations/scans/RBM22.PWMscan.done', 
        'sites/operations/scans/RBM3.PWMscan.done', 
        'sites/operations/scans/ESRP1.PWMscan.done', 
        'sites/operations/scans/ESRP2.PWMscan.done', 
        'sites/operations/scans/RBM7.PWMscan.done', 
        'sites/operations/scans/RBM8A.PWMscan.done', 
        'sites/operations/scans/RBFOX2.PWMscan.done', 
        'sites/operations/scans/RBMS1.PWMscan.done', 
        'sites/operations/scans/RFC2.PWMscan.done', 
        'sites/operations/scans/RFC3.PWMscan.done', 
        'sites/operations/scans/RFXANK.PWMscan.done', 
        'sites/operations/scans/RIOK2.PWMscan.done', 
        'sites/operations/scans/MEX3C.PWMscan.done', 
        'sites/operations/scans/RNASEH2C.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group44.done'
    shell:
        'touch {output}'
rule PWMscan_group45:
    input:
        'sites/operations/scans/RNF138.PWMscan.done', 
        'sites/operations/scans/RPL35.PWMscan.done', 
        'sites/operations/scans/RPL6.PWMscan.done', 
        'sites/operations/scans/RPP25.PWMscan.done', 
        'sites/operations/scans/RPS10.PWMscan.done', 
        'sites/operations/scans/RPS4X.PWMscan.done', 
        'sites/operations/scans/RPS6KA5.PWMscan.done', 
        'sites/operations/scans/RUFY3.PWMscan.done', 
        'sites/operations/scans/RUVBL1.PWMscan.done', 
        'sites/operations/scans/SCAND2P.PWMscan.done', 
        'sites/operations/scans/PDS5A.PWMscan.done', 
        'sites/operations/scans/SCMH1.PWMscan.done', 
        'sites/operations/scans/SEMA4A.PWMscan.done', 
        'sites/operations/scans/SF1.PWMscan.done', 
        'sites/operations/scans/SF3B1.PWMscan.done', 
        'sites/operations/scans/SFT2D1.PWMscan.done', 
        'sites/operations/scans/SLC18A1.PWMscan.done', 
        'sites/operations/scans/SMAP2.PWMscan.done', 
        'sites/operations/scans/SMCR7L.PWMscan.done', 
        'sites/operations/scans/SMPX.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group45.done'
    shell:
        'touch {output}'
rule PWMscan_group46:
    input:
        'sites/operations/scans/SMUG1.PWMscan.done', 
        'sites/operations/scans/SNAPC4.PWMscan.done', 
        'sites/operations/scans/SNAPC5.PWMscan.done', 
        'sites/operations/scans/SND1.PWMscan.done', 
        'sites/operations/scans/SNRNP70.PWMscan.done', 
        'sites/operations/scans/SNRPB2.PWMscan.done', 
        'sites/operations/scans/SOCS4.PWMscan.done', 
        'sites/operations/scans/SOD1.PWMscan.done', 
        'sites/operations/scans/SOX14.PWMscan.done', 
        'sites/operations/scans/SPAG7.PWMscan.done', 
        'sites/operations/scans/SPATS2.PWMscan.done', 
        'sites/operations/scans/SPR.PWMscan.done', 
        'sites/operations/scans/SRBD1.PWMscan.done', 
        'sites/operations/scans/SRP9.PWMscan.done', 
        'sites/operations/scans/SSBP3.PWMscan.done', 
        'sites/operations/scans/SSX2.PWMscan.done', 
        'sites/operations/scans/SSX3.PWMscan.done', 
        'sites/operations/scans/STAU2.PWMscan.done', 
        'sites/operations/scans/STUB1.PWMscan.done', 
        'sites/operations/scans/SUCLG1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group46.done'
    shell:
        'touch {output}'
rule PWMscan_group47:
    input:
        'sites/operations/scans/TAF1A.PWMscan.done', 
        'sites/operations/scans/TAF9.PWMscan.done', 
        'sites/operations/scans/TAGLN2.PWMscan.done', 
        'sites/operations/scans/TBPL1.PWMscan.done', 
        'sites/operations/scans/TCEAL2.PWMscan.done', 
        'sites/operations/scans/TCEAL6.PWMscan.done', 
        'sites/operations/scans/TFAM.PWMscan.done', 
        'sites/operations/scans/TGIF2LX.PWMscan.done', 
        'sites/operations/scans/THAP5.PWMscan.done', 
        'sites/operations/scans/THRA.PWMscan.done', 
        'sites/operations/scans/MED30.PWMscan.done', 
        'sites/operations/scans/TIA1.PWMscan.done', 
        'sites/operations/scans/TIMELESS.PWMscan.done', 
        'sites/operations/scans/TIMM44.PWMscan.done', 
        'sites/operations/scans/TIMM8A.PWMscan.done', 
        'sites/operations/scans/TMSB4XP8.PWMscan.done', 
        'sites/operations/scans/TOB2.PWMscan.done', 
        'sites/operations/scans/TP73.PWMscan.done', 
        'sites/operations/scans/TPI1.PWMscan.done', 
        'sites/operations/scans/TPPP.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group47.done'
    shell:
        'touch {output}'
rule PWMscan_group48:
    input:
        'sites/operations/scans/TRIM21.PWMscan.done', 
        'sites/operations/scans/TRIM69.PWMscan.done', 
        'sites/operations/scans/TRIP10.PWMscan.done', 
        'sites/operations/scans/TRMT1.PWMscan.done', 
        'sites/operations/scans/TROVE2.PWMscan.done', 
        'sites/operations/scans/TSC22D4.PWMscan.done', 
        'sites/operations/scans/TSN.PWMscan.done', 
        'sites/operations/scans/TSNAX.PWMscan.done', 
        'sites/operations/scans/TULP1.PWMscan.done', 
        'sites/operations/scans/U2AF1.PWMscan.done', 
        'sites/operations/scans/UBB.PWMscan.done', 
        'sites/operations/scans/UBE2V1.PWMscan.done', 
        'sites/operations/scans/UGP2.PWMscan.done', 
        'sites/operations/scans/UQCRB.PWMscan.done', 
        'sites/operations/scans/USP39.PWMscan.done', 
        'sites/operations/scans/UTP18.PWMscan.done', 
        'sites/operations/scans/VAMP3.PWMscan.done', 
        'sites/operations/scans/EZR.PWMscan.done', 
        'sites/operations/scans/VPS4B.PWMscan.done', 
        'sites/operations/scans/NELFA.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group48.done'
    shell:
        'touch {output}'
rule PWMscan_group49:
    input:
        'sites/operations/scans/WISP2.PWMscan.done', 
        'sites/operations/scans/XG.PWMscan.done', 
        'sites/operations/scans/XRCC1.PWMscan.done', 
        'sites/operations/scans/YEATS4.PWMscan.done', 
        'sites/operations/scans/YWHAE.PWMscan.done', 
        'sites/operations/scans/YWHAZ.PWMscan.done', 
        'sites/operations/scans/ZBTB12.PWMscan.done', 
        'sites/operations/scans/ZBTB25.PWMscan.done', 
        'sites/operations/scans/ZBTB43.PWMscan.done', 
        'sites/operations/scans/ZBTB46.PWMscan.done', 
        'sites/operations/scans/ZC3H7A.PWMscan.done', 
        'sites/operations/scans/ZCCHC14.PWMscan.done', 
        'sites/operations/scans/ZCCHC17.PWMscan.done', 
        'sites/operations/scans/ZDHHC15.PWMscan.done', 
        'sites/operations/scans/ZDHHC5.PWMscan.done', 
        'sites/operations/scans/ZFP3.PWMscan.done', 
        'sites/operations/scans/ZHX3.PWMscan.done', 
        'sites/operations/scans/ZMAT2.PWMscan.done', 
        'sites/operations/scans/ZMAT4.PWMscan.done', 
        'sites/operations/scans/ZNF124.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group49.done'
    shell:
        'touch {output}'
rule PWMscan_group50:
    input:
        'sites/operations/scans/ZNF131.PWMscan.done', 
        'sites/operations/scans/ZNF160.PWMscan.done', 
        'sites/operations/scans/ZKSCAN8.PWMscan.done', 
        'sites/operations/scans/ZSCAN9.PWMscan.done', 
        'sites/operations/scans/ZNF205.PWMscan.done', 
        'sites/operations/scans/ZNF207.PWMscan.done', 
        'sites/operations/scans/ZBTB18.PWMscan.done', 
        'sites/operations/scans/ZNF250.PWMscan.done', 
        'sites/operations/scans/ZNF26.PWMscan.done', 
        'sites/operations/scans/ZNF3.PWMscan.done', 
        'sites/operations/scans/ZNF304.PWMscan.done', 
        'sites/operations/scans/RNF114.PWMscan.done', 
        'sites/operations/scans/ZSCAN31.PWMscan.done', 
        'sites/operations/scans/ZNF326.PWMscan.done', 
        'sites/operations/scans/ZNF385A.PWMscan.done', 
        'sites/operations/scans/ZNF503.PWMscan.done', 
        'sites/operations/scans/ZNF510.PWMscan.done', 
        'sites/operations/scans/ZNF655.PWMscan.done', 
        'sites/operations/scans/ZNF671.PWMscan.done', 
        'sites/operations/scans/ZNF695.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group50.done'
    shell:
        'touch {output}'
rule PWMscan_group51:
    input:
        'sites/operations/scans/ZNF706.PWMscan.done', 
        'sites/operations/scans/ZNF71.PWMscan.done', 
        'sites/operations/scans/ZNF720.PWMscan.done', 
        'sites/operations/scans/ZNF76.PWMscan.done', 
        'sites/operations/scans/ZNF766.PWMscan.done', 
        'sites/operations/scans/ZRSR2.PWMscan.done', 
        'sites/operations/scans/ZSWIM1.PWMscan.done', 
        'sites/operations/scans/Myf.PWMscan.done', 
        'sites/operations/scans/Pax6.PWMscan.done', 
        'sites/operations/scans/RORA_1.PWMscan.done', 
        'sites/operations/scans/RORA_2.PWMscan.done', 
        'sites/operations/scans/YY1.PWMscan.done', 
        'sites/operations/scans/TP53.PWMscan.done', 
        'sites/operations/scans/RELA.PWMscan.done', 
        'sites/operations/scans/ZNF354C.PWMscan.done', 
        'sites/operations/scans/MIZF.PWMscan.done', 
        'sites/operations/scans/AP1.PWMscan.done', 
        'sites/operations/scans/DUX4.PWMscan.done', 
        'sites/operations/scans/FOXP1.PWMscan.done', 
        'sites/operations/scans/POU2F2.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group51.done'
    shell:
        'touch {output}'
rule PWMscan_group52:
    input:
        'sites/operations/scans/TCF7L2.PWMscan.done', 
        'sites/operations/scans/TP63.PWMscan.done', 
        'sites/operations/scans/ZBTB33.PWMscan.done', 
        'sites/operations/scans/ZNF263.PWMscan.done', 
        'sites/operations/scans/AR.PWMscan.done', 
        'sites/operations/scans/KLF5.PWMscan.done', 
        'sites/operations/scans/T.PWMscan.done', 
        'sites/operations/scans/EN1.PWMscan.done', 
        'sites/operations/scans/ZNF143.PWMscan.done', 
        'sites/operations/scans/NR3C1.PWMscan.done', 
        'sites/operations/scans/ESRRB.PWMscan.done', 
        'sites/operations/scans/HOXA5.PWMscan.done', 
        'sites/operations/scans/DMRT3.PWMscan.done', 
        'sites/operations/scans/LBX1.PWMscan.done', 
        'sites/operations/scans/POU6F1.PWMscan.done', 
        'sites/operations/scans/BARHL2.PWMscan.done', 
        'sites/operations/scans/ELF4.PWMscan.done', 
        'sites/operations/scans/EN2.PWMscan.done', 
        'sites/operations/scans/HOXA13.PWMscan.done', 
        'sites/operations/scans/HOXC11.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group52.done'
    shell:
        'touch {output}'
rule PWMscan_group53:
    input:
        'sites/operations/scans/ONECUT1.PWMscan.done', 
        'sites/operations/scans/POU4F2.PWMscan.done', 
        'sites/operations/scans/ZBTB7B.PWMscan.done', 
        'sites/operations/scans/ZBTB7C.PWMscan.done', 
        'sites/operations/scans/RHOXF1.PWMscan.done', 
        'sites/operations/scans/UNCX.PWMscan.done', 
        'sites/operations/scans/NR3C2.PWMscan.done', 
        'sites/operations/scans/SP8.PWMscan.done', 
        'sites/operations/scans/YY2.PWMscan.done', 
        'sites/operations/scans/ZBTB7A.PWMscan.done', 
        'sites/operations/scans/ZNF410.PWMscan.done', 
        'sites/operations/scans/ZNF740.PWMscan.done', 
        'sites/operations/scans/ONECUT2.PWMscan.done', 
        'sites/operations/scans/ONECUT3.PWMscan.done', 
        'sites/operations/scans/MYBL1.PWMscan.done', 
        'sites/operations/scans/MYBL2.PWMscan.done', 
        'sites/operations/scans/PAX9.PWMscan.done', 
        'sites/operations/scans/PKNOX1.PWMscan.done', 
        'sites/operations/scans/POU1F1.PWMscan.done', 
        'sites/operations/scans/POU2F1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group53.done'
    shell:
        'touch {output}'
rule PWMscan_group54:
    input:
        'sites/operations/scans/POU3F1.PWMscan.done', 
        'sites/operations/scans/POU3F3.PWMscan.done', 
        'sites/operations/scans/POU3F4.PWMscan.done', 
        'sites/operations/scans/POU4F1.PWMscan.done', 
        'sites/operations/scans/POU5F1B.PWMscan.done', 
        'sites/operations/scans/POU6F2.PWMscan.done', 
        'sites/operations/scans/HOXD12.PWMscan.done', 
        'sites/operations/scans/BSX.PWMscan.done', 
        'sites/operations/scans/HMBOX1.PWMscan.done', 
        'sites/operations/scans/HOXA10.PWMscan.done', 
        'sites/operations/scans/HOXA2.PWMscan.done', 
        'sites/operations/scans/HOXB2.PWMscan.done', 
        'sites/operations/scans/HOXB3.PWMscan.done', 
        'sites/operations/scans/HOXC10.PWMscan.done', 
        'sites/operations/scans/HOXC12.PWMscan.done', 
        'sites/operations/scans/HOXC13.PWMscan.done', 
        'sites/operations/scans/HOXD11.PWMscan.done', 
        'sites/operations/scans/HOXD13.PWMscan.done', 
        'sites/operations/scans/NFATC2.PWMscan.done', 
        'sites/operations/scans/ASCL1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group54.done'
    shell:
        'touch {output}'
rule PWMscan_group55:
    input:
        'sites/operations/scans/FOXK2.PWMscan.done', 
        'sites/operations/scans/GRHL2.PWMscan.done', 
        'sites/operations/scans/KLF9.PWMscan.done', 
        'sites/operations/scans/NR2F2.PWMscan.done', 
        'sites/operations/scans/POU5F1.PWMscan.done', 
        'sites/operations/scans/RBPJ.PWMscan.done', 
        'sites/operations/scans/SIX1.PWMscan.done', 
        'sites/operations/scans/SIX2.PWMscan.done', 
        'sites/operations/scans/TEAD2.PWMscan.done', 
        'sites/operations/scans/ZNF24.PWMscan.done', 
        'sites/operations/scans/ZNF384.PWMscan.done', 
        'sites/operations/scans/ZNF282.PWMscan.done', 
        'sites/operations/scans/ZSCAN4.PWMscan.done', 
        'sites/operations/scans/RORB.PWMscan.done', 
        'sites/operations/scans/RORC.PWMscan.done', 
        'sites/operations/scans/TCF7L1.PWMscan.done', 
        'sites/operations/scans/HINFP1.PWMscan.done', 
        'sites/operations/scans/ZNF238.PWMscan.done', 
        'sites/operations/scans/ZNF306.PWMscan.done', 
        'sites/operations/scans/ZNF524.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group55.done'
    shell:
        'touch {output}'
rule PWMscan_group56:
    input:
        'sites/operations/scans/ZNF75A.PWMscan.done', 
        'sites/operations/scans/ZNF784.PWMscan.done', 
        'sites/operations/scans/HSFY2.PWMscan.done', 
        'sites/operations/scans/NFATC1.PWMscan.done', 
        'sites/operations/scans/POU2F3.PWMscan.done', 
        'sites/operations/scans/POU5F1P1.PWMscan.done', 
        'sites/operations/scans/BHLHB2.PWMscan.done', 
        'sites/operations/scans/BHLHB3.PWMscan.done', 
        'sites/operations/scans/CART1.PWMscan.done', 
        'sites/operations/scans/HOXA1.PWMscan.done', 
        'sites/operations/scans/HOXB5.PWMscan.done', 
        'sites/operations/scans/HOXD8.PWMscan.done', 
        'sites/operations/scans/IRX5.PWMscan.done', 
        'sites/operations/scans/PHOX2B.PWMscan.done', 
        'sites/operations/scans/RAXL1.PWMscan.done', 
        'sites/operations/scans/ESRRG.PWMscan.done', 
        'sites/operations/scans/THRB.PWMscan.done', 
        'sites/operations/scans/Trp53.PWMscan.done', 
        'sites/operations/scans/Trp73.PWMscan.done', 
        'sites/operations/scans/ZBTB49.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group56.done'
    shell:
        'touch {output}'
rule PWMscan_group57:
    input:
        'sites/operations/scans/ZNF232.PWMscan.done', 
        'sites/operations/scans/ZNF435.PWMscan.done', 
        'sites/operations/scans/ZNF713.PWMscan.done', 
        'sites/operations/scans/ARID5A.PWMscan.done', 
        'sites/operations/scans/BARHL1.PWMscan.done', 
        'sites/operations/scans/BBX.PWMscan.done', 
        'sites/operations/scans/BCL3.PWMscan.done', 
        'sites/operations/scans/CHD1.PWMscan.done', 
        'sites/operations/scans/CHD2.PWMscan.done', 
        'sites/operations/scans/CREB3L2.PWMscan.done', 
        'sites/operations/scans/DBX2.PWMscan.done', 
        'sites/operations/scans/DMC1.PWMscan.done', 
        'sites/operations/scans/EBF3.PWMscan.done', 
        'sites/operations/scans/EP300.PWMscan.done', 
        'sites/operations/scans/EZH2.PWMscan.done', 
        'sites/operations/scans/FOXJ1.PWMscan.done', 
        'sites/operations/scans/FOXN1.PWMscan.done', 
        'sites/operations/scans/GMEB1.PWMscan.done', 
        'sites/operations/scans/GTF2F1.PWMscan.done', 
        'sites/operations/scans/GTF2I.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group57.done'
    shell:
        'touch {output}'
rule PWMscan_group58:
    input:
        'sites/operations/scans/GZF1.PWMscan.done', 
        'sites/operations/scans/HCFC1.PWMscan.done', 
        'sites/operations/scans/HDX.PWMscan.done', 
        'sites/operations/scans/HIVEP1.PWMscan.done', 
        'sites/operations/scans/HLX.PWMscan.done', 
        'sites/operations/scans/HOXA11.PWMscan.done', 
        'sites/operations/scans/HOXA3.PWMscan.done', 
        'sites/operations/scans/HOXA4.PWMscan.done', 
        'sites/operations/scans/HOXA6.PWMscan.done', 
        'sites/operations/scans/HOXA7.PWMscan.done', 
        'sites/operations/scans/HOXA9.PWMscan.done', 
        'sites/operations/scans/HOXB1.PWMscan.done', 
        'sites/operations/scans/HOXB4.PWMscan.done', 
        'sites/operations/scans/HOXB6.PWMscan.done', 
        'sites/operations/scans/HOXB7.PWMscan.done', 
        'sites/operations/scans/HOXB8.PWMscan.done', 
        'sites/operations/scans/HOXC4.PWMscan.done', 
        'sites/operations/scans/HOXC5.PWMscan.done', 
        'sites/operations/scans/HOXC6.PWMscan.done', 
        'sites/operations/scans/HOXC8.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group58.done'
    shell:
        'touch {output}'
rule PWMscan_group59:
    input:
        'sites/operations/scans/HOXC9.PWMscan.done', 
        'sites/operations/scans/HOXD10.PWMscan.done', 
        'sites/operations/scans/HOXD1.PWMscan.done', 
        'sites/operations/scans/HOXD4.PWMscan.done', 
        'sites/operations/scans/HOXD9.PWMscan.done', 
        'sites/operations/scans/IKZF2.PWMscan.done', 
        'sites/operations/scans/IRX4.PWMscan.done', 
        'sites/operations/scans/IRX6.PWMscan.done', 
        'sites/operations/scans/KLF7.PWMscan.done', 
        'sites/operations/scans/LHX1.PWMscan.done', 
        'sites/operations/scans/LHX5.PWMscan.done', 
        'sites/operations/scans/MECOM.PWMscan.done', 
        'sites/operations/scans/MTA3.PWMscan.done', 
        'sites/operations/scans/OSR1.PWMscan.done', 
        'sites/operations/scans/OSR2.PWMscan.done', 
        'sites/operations/scans/OTP.PWMscan.done', 
        'sites/operations/scans/PATZ1.PWMscan.done', 
        'sites/operations/scans/PGR.PWMscan.done', 
        'sites/operations/scans/PML.PWMscan.done', 
        'sites/operations/scans/PRDM14.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group59.done'
    shell:
        'touch {output}'
rule PWMscan_group60:
    input:
        'sites/operations/scans/RAD21.PWMscan.done', 
        'sites/operations/scans/RCOR1.PWMscan.done', 
        'sites/operations/scans/RFX7.PWMscan.done', 
        'sites/operations/scans/RHOXF2.PWMscan.done', 
        'sites/operations/scans/SIN3A.PWMscan.done', 
        'sites/operations/scans/SIX3.PWMscan.done', 
        'sites/operations/scans/SIX4.PWMscan.done', 
        'sites/operations/scans/SIX5.PWMscan.done', 
        'sites/operations/scans/SIX6.PWMscan.done', 
        'sites/operations/scans/SMARCC1.PWMscan.done', 
        'sites/operations/scans/SMARCC2.PWMscan.done', 
        'sites/operations/scans/SMC3.PWMscan.done', 
        'sites/operations/scans/SOX12.PWMscan.done', 
        'sites/operations/scans/SOX30.PWMscan.done', 
        'sites/operations/scans/SOX6.PWMscan.done', 
        'sites/operations/scans/SP100.PWMscan.done', 
        'sites/operations/scans/STAT5A.PWMscan.done', 
        'sites/operations/scans/STAT5B.PWMscan.done', 
        'sites/operations/scans/TAF1.PWMscan.done', 
        'sites/operations/scans/TBL1XR1.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group60.done'
    shell:
        'touch {output}'
rule PWMscan_group61:
    input:
        'sites/operations/scans/TCF21.PWMscan.done', 
        'sites/operations/scans/TFAP2E.PWMscan.done', 
        'sites/operations/scans/TFCP2L1.PWMscan.done', 
        'sites/operations/scans/TLX2.PWMscan.done', 
        'sites/operations/scans/UBP1.PWMscan.done', 
        'sites/operations/scans/WRNIP1.PWMscan.done', 
        'sites/operations/scans/YBX1.PWMscan.done', 
        'sites/operations/scans/ZBTB14.PWMscan.done', 
        'sites/operations/scans/ZBTB16.PWMscan.done', 
        'sites/operations/scans/ZBTB3.PWMscan.done', 
        'sites/operations/scans/ZKSCAN1.PWMscan.done', 
        'sites/operations/scans/ZKSCAN3.PWMscan.done', 
        'sites/operations/scans/ZNF148.PWMscan.done', 
        'sites/operations/scans/ZNF219.PWMscan.done', 
        'sites/operations/scans/ZNF274.PWMscan.done', 
        'sites/operations/scans/ZNF281.PWMscan.done', 
        'sites/operations/scans/ZNF333.PWMscan.done', 
        'sites/operations/scans/ZNF350.PWMscan.done', 
        'sites/operations/scans/ZNF35.PWMscan.done', 
        'sites/operations/scans/ZNF423.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group61.done'
    shell:
        'touch {output}'
rule PWMscan_group62:
    input:
        'sites/operations/scans/ZNF652.PWMscan.done', 
        'sites/operations/scans/ZNF691.PWMscan.done', 
        'sites/operations/scans/ZNF711.PWMscan.done', 
        'sites/operations/scans/ZNF8.PWMscan.done', 
        'sites/operations/scans/Sox4.PWMscan.done'
    output:
        'sites/operations/groups/PWMscan.group62.done'
    shell:
        'touch {output}'
