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
#### SPOOL INDIVIDUAL OPERATIONS #######################################################################################################
########################################################################################################################################

rule snu61wt01_pantf_analysis:
    input:
        "snu61/wt01/pantf/operations/SNU61-WT-01.alltf.analysis.done.txt"

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
#### FOOTPRINT ANALYSIS RULES ##########################################################################################################
########################################################################################################################################

rule STEP35_make_footprint_signal_by_chr:
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        "sites/{gene}.sites.Rdata",
        "{path}operations/{mergedsample}-downsample.final.txt",
        "{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt"
    output:
        "{path}footprints/operations/{mergedsample}.{gene}.{chr}.done.bychr.txt"
    resources:
        fp_by_chr=1
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
    resources:
        raw_fp_graph=1
    script:
        "scripts/snakeGenerateMergedFPGraph.R"

rule STEP38_parse_footprint_signals_and_generate_graphs:
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        "sites/{gene}.sites.Rdata",
        "{path}footprints/operations/{mergedsample}.{gene}.merged.done.txt",
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak",
        "{path}footprints/operations/{mergedsample}.{gene}.graphs.done.txt"
    output:
        "{path}footprints/operations/{mergedsample}.{gene}.parsed.done.txt"
    resources:
        parse_fp=1
    script:
        "scripts/snakeParseFP.R"

# rule STEP39_merge_footprint_signal_motifs:
#     input:
#         "{path}footprints/operations/{mergedsample}.{gene}.parsed.done.txt"
#     output:
#         "{path}footprints/data/motifmerge/{mergedsample}.{gene}.{nummotif}.mergedmotif.Rdata"
#     script:
#         "scripts/snakeMergeMotifs.R"

rule AGGREGATOR_COADMR_footprinting:
    input:
        "{path}footprints/operations/{mergedsample}.CDX2.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.TCF7.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.HOXA3.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.MNX1.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.POU5F1B.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.OVOL1.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.ESRRA.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.ASCL2.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.HNF4A.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.GMEB2.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.ZSWIM1.parsed.done.txt",
        "{path}footprints/operations/{mergedsample}.CBFA2T2.parsed.done.txt"
    output:
        "{path}footprints/operations/{mergedsample}.footprints.coadmr.done.txt"
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
#### PAN-TF ANALYSIS RULES #############################################################################################################
########################################################################################################################################

rule PANTF_make_footprint_signal_by_chr:
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        "sites/{gene}.sites.Rdata",
        "{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt"
    output:
        "{path}pantf/operations/{mergedsample}.{gene}.{chr}.pantf.done.bychr.txt"
    resources:
        fp_by_chr=1
    script:
        "scripts/panTF/snakeMakeFPbyChrPanTF.R"

rule PANTF_merge_footprint_signal_by_chr:
    input:
        "sites/{gene}.sites.Rdata",
        "{path}pantf/operations/{mergedsample}.{gene}.chr1.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr2.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr3.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr4.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr5.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr6.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr7.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr8.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr9.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr10.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr11.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr12.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr13.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr14.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr15.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr16.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr17.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr18.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr19.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr20.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr21.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chr22.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chrX.pantf.done.bychr.txt",
        "{path}pantf/operations/{mergedsample}.{gene}.chrY.pantf.done.bychr.txt"
    output:
        "{path}pantf/operations/{mergedsample}.{gene}.merged.pantf.done.txt"
    script:
        "scripts/panTF/snakeMergeFPbyChrPanTF.R"

rule PANTF_parse_footprint_signals:
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        "sites/{gene}.sites.Rdata",
        "{path}pantf/operations/{mergedsample}.{gene}.merged.pantf.done.txt",
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
    output:
        "{path}pantf/operations/{mergedsample}.{gene}.parsed.pantf.done.txt"
    resources:
        parse_fp=1
    script:
        "scripts/panTF/snakeParseFPPanTF.R"

rule PANTF_TFgroup_aggregator:
    input:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
        # '{path}pantf/operations/{mergedsample}.pantf.group2.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group3.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group4.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group5.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group6.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group7.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group8.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group9.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group10.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group11.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group12.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group13.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group14.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group15.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group16.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group17.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group18.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group19.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group20.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group21.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group22.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group23.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group24.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group25.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group26.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group27.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group28.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group29.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group30.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group31.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group32.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group33.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group34.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group35.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group36.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group37.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group38.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group39.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group40.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group41.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group42.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group43.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group44.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group45.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group46.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group47.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group48.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group49.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group50.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group51.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group52.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group53.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group54.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group55.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group56.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group57.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group58.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group59.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group60.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group61.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group62.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group63.done.txt',
        # '{path}pantf/operations/{mergedsample}.pantf.group64.done.txt'
    output:
        "{path}pantf/operations/{mergedsample}.alltf.analysis.done.txt"
    shell:
        "touch {output}"

########################################################################################################################################
#### PAN-TF GROUPS #####################################################################################################################
########################################################################################################################################

rule PANTF_group1:
    input:
        '{path}pantf/operations/{mergedsample}.TFAP2A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFIL3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HLF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NHLH1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.USF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CEBPA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EBF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CEBPB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOS.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUNB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUND.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAFF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAFK.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP2C.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.USF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SREBF1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group2:
    input:
        '{path}pantf/operations/{mergedsample}.SREBF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AHR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARNT.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATF6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BACH1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BACH2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CREB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TCF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.XBP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARID5B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYOD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFE2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYCN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFE2L1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TEF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BATF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TCF12.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group3:
    input:
        '{path}pantf/operations/{mergedsample}.MYC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MXI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHLHE40.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARNTL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATF4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATF7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BATF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHLHA15.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHLHE41.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHLHE22.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHLHE23.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CEBPD.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CEBPE.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CEBPG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CLOCK.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CREB3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CREB3L1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DBP.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FIGLA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HES5.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group4:
    input:
        '{path}pantf/operations/{mergedsample}.HES7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HEY1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HEY2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ID4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JDP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAFG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MESP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MGA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MLX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MLXIPL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MNT.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MSC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYF6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NEUROD2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NEUROG2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NRL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OLIG1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OLIG2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OLIG3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TCF4.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group5:
    input:
        '{path}pantf/operations/{mergedsample}.TFAP2B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFE3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFEB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFEC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP2D.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARID3A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARNT2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATF5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CREM.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DDIT3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EPAS1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HAND1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HES1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HIF1A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HMGA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HMGA2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAFA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAFB.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group6:
    input:
        '{path}pantf/operations/{mergedsample}.MAF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MITF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYOG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NEUROD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFE2L2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PTF1A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TAL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TWIST1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AIRE.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ALX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ALX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ALX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ANDR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AP2A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AP2B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AP2C.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AP2D.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARI3A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARI5B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARX.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group7:
    input:
        '{path}pantf/operations/{mergedsample}.ASCL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATF6A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ATOH1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BARH1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BARH2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BARX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BARX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BC11A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BCL6B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BCL6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHA15.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHE22.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHE23.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHE40.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHE41.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BMAL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BPTF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BRAC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BRCA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BSH.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group8:
    input:
        '{path}pantf/operations/{mergedsample}.CDC5L.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CDX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CDX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CEBPZ.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CENPB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.COE1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.COT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.COT2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CPEB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CR3L1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CR3L2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CREB5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CRX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CTCFL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CTCF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CUX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CUX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CXXC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DLX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DLX2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group9:
    input:
        '{path}pantf/operations/{mergedsample}.DLX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DLX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DLX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DLX6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DMBX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DPRX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DRGX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DUXA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E2F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E2F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E2F3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E2F4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E2F5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E2F6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E2F7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E2F8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.E4F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EGR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EGR2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EGR3.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group10:
    input:
        '{path}pantf/operations/{mergedsample}.EGR4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EHF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ELF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ELF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ELF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ELF5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ELK1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ELK3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ELK4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EMX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EMX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EOMES.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ERF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ERG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ERR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ERR2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ERR3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ESR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ESR2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ESX1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group11:
    input:
        '{path}pantf/operations/{mergedsample}.ETS1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETS2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETV1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETV2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETV3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETV4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETV5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETV6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETV7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EVI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EVX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EVX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FEV.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FLI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXA2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXA3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXC2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group12:
    input:
        '{path}pantf/operations/{mergedsample}.FOXD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXD2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXD3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXG1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXH1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXJ2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXJ3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXK1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXM1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXO1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXO3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXO4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXO6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXP3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXQ1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group13:
    input:
        '{path}pantf/operations/{mergedsample}.FUBP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GABP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GABPA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GATA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GATA2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GATA3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GATA4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GATA5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GATA6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GBX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GBX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GCM1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GCM2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GCR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GFI1B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GFI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GLI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GLI2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GLI3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GLIS1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group14:
    input:
        '{path}pantf/operations/{mergedsample}.GLIS2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GLIS3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GMEB2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GRHL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GSC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GSC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GSX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GSX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HBP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HEN1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HESX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HIC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HIC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HINFP.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HLTF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HMBX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HME1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HME2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HMX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HMX2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group15:
    input:
        '{path}pantf/operations/{mergedsample}.HMX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNF1A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNF1B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNF4A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNF4G.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNF6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOMEZ.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HSF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HSF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HSF4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HSFY1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HTF4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXA10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXA11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXA13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXA2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXA5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXA7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXA9.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group16:
    input:
        '{path}pantf/operations/{mergedsample}.HXB13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXB2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXB3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXB6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXB7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXB8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXC10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXC11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXC12.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXC13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXC6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXC8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXD10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXD11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXD12.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXD13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXD3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXD4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HXD8.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group17:
    input:
        '{path}pantf/operations/{mergedsample}.HXD9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IKZF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.INSM1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ISL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ISL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ISX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ITF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KAISO.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF12.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF13.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group18:
    input:
        '{path}pantf/operations/{mergedsample}.KLF14.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF15.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF16.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LBX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LEF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LHX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LHX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LHX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LHX6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LHX8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LHX9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LMX1A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LMX1B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAZ.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MBD2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group19:
    input:
        '{path}pantf/operations/{mergedsample}.MCR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MECP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEF2A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEF2B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEF2C.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEF2D.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEIS1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEIS2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEIS3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEOX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEOX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MGAP.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MIXL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MLXPL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MNX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MSX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MSX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MTF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MUSC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYBA.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group20:
    input:
        '{path}pantf/operations/{mergedsample}.MYBB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MZF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NANOG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NDF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NDF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NF2L1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NF2L2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFAC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFAC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFAC3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFAC4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFAT5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFIA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFIC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFKB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFKB2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFYA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFYB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFYC.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group21:
    input:
        '{path}pantf/operations/{mergedsample}.NGN2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX21.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX22.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX23.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX25.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX28.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX31.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX32.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX61.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX62.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NOBOX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NOTO.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR0B1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR1D1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR1H2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR1H4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR1I2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR1I3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR2C1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR2C2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group22:
    input:
        '{path}pantf/operations/{mergedsample}.NR2E1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR2E3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR2F6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR4A1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR4A2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR4A3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR5A2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR6A1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NRF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ONEC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ONEC3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OTX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OTX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OVOL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.P53.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.P5F1B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.P63.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.P73.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAX2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group23:
    input:
        '{path}pantf/operations/{mergedsample}.PAX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAX6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAX7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAX8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PBX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PBX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PBX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PDX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PEBB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PHX2A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PHX2B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PIT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PITX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PITX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PITX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PKNX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PKNX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PLAG1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group24:
    input:
        '{path}pantf/operations/{mergedsample}.PLAL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO2F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO2F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO2F3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO3F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO3F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO3F3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO3F4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO4F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO4F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO4F3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO5F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO6F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PO6F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PPARA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PPARD.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PPARG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRD14.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRDM1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRDM4.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group25:
    input:
        '{path}pantf/operations/{mergedsample}.PRGR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PROP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PROX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRRX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRRX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PURA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RARA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RARB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RARG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RAX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RELB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.REL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.REST.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RHXF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RORA.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group26:
    input:
        '{path}pantf/operations/{mergedsample}.RORG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RREB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RUNX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RUNX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RUNX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RXRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RXRB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RXRG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SCRT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SCRT2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SHOX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SHOX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMAD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMAD2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMAD3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMAD4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMRC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SNAI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SNAI2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group27:
    input:
        '{path}pantf/operations/{mergedsample}.SOX10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX15.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX17.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX18.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX21.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SP3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SP4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SPDEF.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group28:
    input:
        '{path}pantf/operations/{mergedsample}.SPI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SPIB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SPIC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SPZ1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SRBP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SRBP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SRF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SRY.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STA5A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STA5B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SUH.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBP.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBX15.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group29:
    input:
        '{path}pantf/operations/{mergedsample}.TBX19.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBX20.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBX21.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TCF7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TEAD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TEAD3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TEAD4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TF2LX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TF65.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TF7L1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TF7L2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFCP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFDP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFE2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TGIF1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group30:
    input:
        '{path}pantf/operations/{mergedsample}.TGIF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.THAP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.THA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.THB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TLX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TWST1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TYY1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TYY2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UBIP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UNC4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.VAX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.VAX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.VDR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.VENTX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.VSX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.VSX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.WT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.YBOX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBED1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBT18.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group31:
    input:
        '{path}pantf/operations/{mergedsample}.ZBT49.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBT7A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBT7B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZEB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZEP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZEP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZFHX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZFX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZIC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZIC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZIC3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZIC4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZKSC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZKSC3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN143.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN148.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN219.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN232.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group32:
    input:
        '{path}pantf/operations/{mergedsample}.ZN282.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN333.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN350.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN384.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN410.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN423.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN524.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN589.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN639.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN652.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN713.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN740.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZN784.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZSC16.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZSCA4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ABCF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.A1CF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ACO1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ADARB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AFF4.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group33:
    input:
        '{path}pantf/operations/{mergedsample}.AGGF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AKR1A1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ANXA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ANXA11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.APEX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARFGAP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ASCC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ASPSCR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AVEN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BAD.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GPANK1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BAX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BCL11A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BOLL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CELF4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CELF5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CELF6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.C19orf25.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.C19orf40.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EXO5.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group34:
    input:
        '{path}pantf/operations/{mergedsample}.LINC00471.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.C9orf156.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CANX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CAT.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CBFA2T2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CBFB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CBX7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF830.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CCDC25.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CD59.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CDK2AP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AGAP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CFL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXN3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CKMT1B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CLK1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CNOT6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NELFB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CPSF4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CSNK2B.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group35:
    input:
        '{path}pantf/operations/{mergedsample}.CSTF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CYB5R1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CYCS.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DAB2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DAZAP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ASAP3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DDX20.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DDX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DDX43.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DDX53.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DGCR8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DHX36.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DIABLO.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DIS3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DNMT3A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DTL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DUS3L.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DUSP22.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DUSP26.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group36:
    input:
        '{path}pantf/operations/{mergedsample}.ECSIT.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EDN1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EEF1D.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EIF5A2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ENO1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ESRRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ETFB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EWSR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EXOSC3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.METTL21B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FAM127B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FEZ1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FEZF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FGF19.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FHL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FIP1L1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SRRM3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXP4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GADD45A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GIT2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group37:
    input:
        '{path}pantf/operations/{mergedsample}.GLYCTK.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GOT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GPAM.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GPD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GRHPR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GTF2B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GTF2H3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GTF3C2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GTF3C5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GTPBP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GTPBP6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.H1FX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.H2AFY.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.H2AFZ.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HCFC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HCLS1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HDAC8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HHAT.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HHEX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UBE2K.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group38:
    input:
        '{path}pantf/operations/{mergedsample}.HIRIP3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HIST1H2BN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HIST2H2AB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HIST2H2BE.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HLCS.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HMG20A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNRNPA0.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNRNPA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNRNPC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNRNPH3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HNRNPLL.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HP1BP3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HSPA1L.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HSPA5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HTATIP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ID2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IL24.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group39:
    input:
        '{path}pantf/operations/{mergedsample}.ING3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRF6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IVD.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KDM5A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KDM5D.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KCNIP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KIAA0907.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KIF22.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LARP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LARP4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LAS1L.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CERS4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UBXN1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CBX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LRRFIP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LSM6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LUZP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LUZP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAGEA8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAGED4B.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group40:
    input:
        '{path}pantf/operations/{mergedsample}.MAGEF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAGOH.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAP4K2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAPK1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MBTPS2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MCTP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MDM2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEF2BNB-MEF2B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GLTPD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBM42.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.WDR83.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MORN1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MRPL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MRPL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MRPS25.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MSI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MSI2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MSRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MSRB3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MTHFD1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group41:
    input:
        '{path}pantf/operations/{mergedsample}.MXD4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYEF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYLK.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NANOS1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NAP1L1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NCALD.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NCBP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFATC3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFATC4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFIB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFIX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX2-3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NME1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NMI.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NMRAL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NNT.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NOC2L.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GAR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NONO.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR2F1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group42:
    input:
        '{path}pantf/operations/{mergedsample}.NUCB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NUP107.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NUP133.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NXPH3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ODC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OTUD4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.P4HB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAXIP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PCK2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PDCD11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PDE6H.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PDLIM5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PGAM2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PHLDA2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PHOX2A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PHTF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PICK1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PIK3C3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PIR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PKM.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group43:
    input:
        '{path}pantf/operations/{mergedsample}.PKNOX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PLAGL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PLG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POLE3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POLI.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU3F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU4F3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PPP1R10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PPP2R3B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PPP5C.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PQBP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRDX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRKRIR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRNP.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PSMA6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PSMC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PTCD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PTPMT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PURG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.R3HDM2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group44:
    input:
        '{path}pantf/operations/{mergedsample}.RAB14.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RAB18.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RAB2A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RAB7A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RAN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RAX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBBP5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBBP9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBM17.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBM22.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBM3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ESRP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ESRP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBM7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBM8A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBFOX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBMS1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFC3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFXANK.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group45:
    input:
        '{path}pantf/operations/{mergedsample}.RIOK2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MEX3C.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RNASEH2C.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RNF138.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RPL35.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RPL6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RPP25.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RPS10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RPS4X.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RPS6KA5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RUFY3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RUVBL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SCAND2P.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PDS5A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SCMH1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SEMA4A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SF3B1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SFT2D1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SLC18A1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group46:
    input:
        '{path}pantf/operations/{mergedsample}.SMAP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMCR7L.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMPX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMUG1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SNAPC4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SNAPC5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SND1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SNRNP70.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SNRPB2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOCS4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX14.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SPAG7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SPATS2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SPR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SRBD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SRP9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SSBP3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SSX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SSX3.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group47:
    input:
        '{path}pantf/operations/{mergedsample}.STAU2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STUB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SUCLG1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TAF1A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TAF9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TAGLN2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBPL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TCEAL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TCEAL6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAM.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TGIF2LX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.THAP5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.THRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MED30.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TIA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TIMELESS.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TIMM44.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TIMM8A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TMSB4XP8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TOB2.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group48:
    input:
        '{path}pantf/operations/{mergedsample}.TP73.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TPI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TPPP.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TRIM21.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TRIM69.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TRIP10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TRMT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TROVE2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TSC22D4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TSN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TSNAX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TULP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.U2AF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UBB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UBE2V1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UGP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UQCRB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.USP39.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UTP18.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.VAMP3.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group49:
    input:
        '{path}pantf/operations/{mergedsample}.EZR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.VPS4B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NELFA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.WISP2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.XG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.XRCC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.YEATS4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.YWHAE.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.YWHAZ.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB12.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB25.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB43.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB46.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZC3H7A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZCCHC14.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZCCHC17.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZDHHC15.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZDHHC5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZFP3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZHX3.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group50:
    input:
        '{path}pantf/operations/{mergedsample}.ZMAT2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZMAT4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF124.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF131.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF160.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZKSCAN8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZSCAN9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF205.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF207.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB18.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF250.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF26.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF304.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RNF114.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZSCAN31.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF326.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF385A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF503.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF510.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group51:
    input:
        '{path}pantf/operations/{mergedsample}.ZNF655.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF671.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF695.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF706.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF71.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF720.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF76.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF766.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZRSR2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZSWIM1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.Myf.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MZF1_1-4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MZF1_5-13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYC::MAX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.Pax6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RORA_1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RORA_2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RXRA::VDR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TAL1::TCF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.YY1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group52:
    input:
        '{path}pantf/operations/{mergedsample}.TP53.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RELA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR1H2::RXRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TLX1::NFIC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX3-1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF354C.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MIZF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EWSR1-FLI1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RXR::RAR_DR5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HIF1A::ARNT.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BATF::JUN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DUX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUN (var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUND (var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFE2::MAF.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU2F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMAD2::SMAD3::SMAD4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT2::STAT1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group53:
    input:
        '{path}pantf/operations/{mergedsample}.TCF7L2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TP63.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB33.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF263.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.AR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TAL1::GATA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.T.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EN1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MZF1(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAX::MYC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PPARG::RXRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RORA(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF143.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR3C1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFIC::TLX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX3-2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GATA1::TAL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ESRRB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA5.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group54:
    input:
        '{path}pantf/operations/{mergedsample}.RARA::RXRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUN(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUND(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MAF::NFE2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT1::STAT2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DMRT3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LBX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU6F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BARHL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ELF4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EN2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXC11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JDP2(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX2-8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX6-1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX6-2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ONECUT1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU4F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB7B.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group55:
    input:
        '{path}pantf/operations/{mergedsample}.ZBTB7C.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RHOXF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.UNCX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR3C2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RARA(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SP8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.YY2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB7A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF410.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF740.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ONECUT2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ONECUT3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYBL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MYBL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PAX9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PKNOX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU1F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU2F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU3F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU3F3.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group56:
    input:
        '{path}pantf/operations/{mergedsample}.POU3F4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU4F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU5F1B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU6F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP2A(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP2B(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP2B(var.3).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP2C(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP2C(var.3).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SREBF2(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TFAP2A(var.3).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD12.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BSX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HMBOX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXC10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXC12.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group57:
    input:
        '{path}pantf/operations/{mergedsample}.HOXC13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD13.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFATC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOS::JUN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARNT::HIF1A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ASCL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXK2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GRHL2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.KLF9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR2F2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU5F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RBPJ.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SIX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SIX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TEAD2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF24.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF384.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF282.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZSCAN4.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group58:
    input:
        '{path}pantf/operations/{mergedsample}.FOS::JUN(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSB::JUN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL1::JUN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL1::JUN(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL2::JUN.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL2::JUN(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUN::JUNB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUN::JUNB(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOS::JUNB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSB::JUNB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSB::JUNB(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL1::JUNB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL2::JUNB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL2::JUNB(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.JUNB(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOS::JUND.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL1::JUND.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL1::JUND(var.2).parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL2::JUND.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOSL2::JUND(var.2).parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group59:
    input:
        '{path}pantf/operations/{mergedsample}.NR1A4::RXRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NR4A2::RXRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PPARA::RXRA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RARA::RXRG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RORB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RORC.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TCF7L1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HINFP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF238.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF306.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF524.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF75A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF784.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HSFY2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NFATC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU2F3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.POU5F1P1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHLHB2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BHLHB3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CART1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group60:
    input:
        '{path}pantf/operations/{mergedsample}.HOXA1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PHOX2B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RAXL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ESRRG.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.THRB.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.Trp53.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.Trp73.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZBTB49.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF232.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF435.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ZNF713.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NA.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.ARID5A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BARHL1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BBX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.BCL3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CHD1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group61:
    input:
        '{path}pantf/operations/{mergedsample}.CHD2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.CREB3L2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DBX2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.DMC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EBF3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EP300.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.EZH2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXJ1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.FOXN1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GMEB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GTF2F1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GTF2I.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.GZF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HCFC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HDX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HIVEP1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HLX.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA11.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA4.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group62:
    input:
        '{path}pantf/operations/{mergedsample}.HOXA6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXA9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXB8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXC4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXC5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXC6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXC8.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXC9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD10.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.HOXD9.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IKZF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.IRX6.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group63:
    input:
        '{path}pantf/operations/{mergedsample}.KLF7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LHX1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.LHX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MECOM.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.MTA3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX1-1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX1-2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX2-1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX2-2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX2-4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX2-5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX2-6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.NKX6-3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OSR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OSR2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.OTP.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PATZ1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PGR.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PML.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.PRDM14.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'

rule PANTF_group64:
    input:
        '{path}pantf/operations/{mergedsample}.RAD21.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RCOR1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RFX7.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.RHOXF2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SIN3A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SIX3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SIX4.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SIX5.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SIX6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMARCC1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMARCC2.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SMC3.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX12.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX30.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SOX6.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.SP100.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT5A.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.STAT5B.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TAF1.parsed.pantf.done.txt', 
        '{path}pantf/operations/{mergedsample}.TBL1XR1.parsed.pantf.done.txt', 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
    shell:
        'touch {output}'