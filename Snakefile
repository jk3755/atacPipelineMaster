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
        "sites/operations/PWMscan.allgroups.done"

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
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt',
        '{path}pantf/operations/{mergedsample}.pantf.group2.done.txt'
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
        '{path}pantf/operations/{mergedsample}.FOS.parsed.pantf.done.txt' 
    output:
        '{path}pantf/operations/{mergedsample}.pantf.group1.done.txt'
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
        'sites/operations/{gene}.PWMscan.done'
    resources:
        scanPWM=1
    script:
        'scripts/scanPWM/snakeScanPWM.R'

rule PWMscan_group_aggregator:
    input:
        'sites/operations/PWMscan.group1.done',
        'sites/operations/PWMscan.group2.done',
        'sites/operations/PWMscan.group3.done',
        'sites/operations/PWMscan.group4.done',
        'sites/operations/PWMscan.group5.done',
        'sites/operations/PWMscan.group6.done',
        'sites/operations/PWMscan.group7.done',
        'sites/operations/PWMscan.group8.done',
        'sites/operations/PWMscan.group9.done',
        'sites/operations/PWMscan.group10.done',
        'sites/operations/PWMscan.group11.done',
        'sites/operations/PWMscan.group12.done',
        'sites/operations/PWMscan.group13.done',
        'sites/operations/PWMscan.group14.done',
        'sites/operations/PWMscan.group15.done',
        'sites/operations/PWMscan.group16.done',
        'sites/operations/PWMscan.group17.done',
        'sites/operations/PWMscan.group18.done',
        'sites/operations/PWMscan.group19.done',
        'sites/operations/PWMscan.group20.done',
        'sites/operations/PWMscan.group21.done',
        'sites/operations/PWMscan.group22.done',
        'sites/operations/PWMscan.group23.done',
        'sites/operations/PWMscan.group24.done',
        'sites/operations/PWMscan.group25.done',
        'sites/operations/PWMscan.group26.done',
        'sites/operations/PWMscan.group27.done',
        'sites/operations/PWMscan.group28.done',
        'sites/operations/PWMscan.group29.done',
        'sites/operations/PWMscan.group30.done',
        'sites/operations/PWMscan.group31.done',
        'sites/operations/PWMscan.group32.done',
        'sites/operations/PWMscan.group33.done',
        'sites/operations/PWMscan.group34.done',
        'sites/operations/PWMscan.group35.done',
        'sites/operations/PWMscan.group36.done',
        'sites/operations/PWMscan.group37.done',
        'sites/operations/PWMscan.group38.done',
        'sites/operations/PWMscan.group39.done',
        'sites/operations/PWMscan.group40.done',
        'sites/operations/PWMscan.group41.done',
        'sites/operations/PWMscan.group42.done',
        'sites/operations/PWMscan.group43.done',
        'sites/operations/PWMscan.group44.done',
        'sites/operations/PWMscan.group45.done',
        'sites/operations/PWMscan.group46.done',
        'sites/operations/PWMscan.group47.done',
        'sites/operations/PWMscan.group48.done',
        'sites/operations/PWMscan.group49.done',
        'sites/operations/PWMscan.group50.done',
        'sites/operations/PWMscan.group51.done',
        'sites/operations/PWMscan.group52.done',
        'sites/operations/PWMscan.group53.done',
        'sites/operations/PWMscan.group54.done',
        'sites/operations/PWMscan.group55.done',
        'sites/operations/PWMscan.group56.done',
        'sites/operations/PWMscan.group57.done',
        'sites/operations/PWMscan.group58.done',
        'sites/operations/PWMscan.group59.done',
        'sites/operations/PWMscan.group60.done',
        'sites/operations/PWMscan.group61.done',
        'sites/operations/PWMscan.group62.done'
    output:
        "sites/operations/PWMscan.allgroups.done"
    shell:
        "touch {output}"

rule PWMscan_group1:
    input:
        'sites/operations/TFAP2A.PWMscan.done', 
        'sites/operations/NFIL3.PWMscan.done', 
        'sites/operations/HLF.PWMscan.done', 
        'sites/operations/NHLH1.PWMscan.done', 
        'sites/operations/MAX.PWMscan.done', 
        'sites/operations/USF1.PWMscan.done', 
        'sites/operations/CEBPA.PWMscan.done', 
        'sites/operations/EBF1.PWMscan.done', 
        'sites/operations/CEBPB.PWMscan.done', 
        'sites/operations/FOS.PWMscan.done', 
        'sites/operations/FOSL1.PWMscan.done', 
        'sites/operations/FOSL2.PWMscan.done', 
        'sites/operations/JUN.PWMscan.done', 
        'sites/operations/JUNB.PWMscan.done', 
        'sites/operations/JUND.PWMscan.done', 
        'sites/operations/MAFF.PWMscan.done', 
        'sites/operations/MAFK.PWMscan.done', 
        'sites/operations/TFAP2C.PWMscan.done', 
        'sites/operations/USF2.PWMscan.done', 
        'sites/operations/SREBF1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group1.done'
    shell:
        'touch {output}'

rule PWMscan_group2:
    input:
        'sites/operations/SREBF2.PWMscan.done', 
        'sites/operations/AHR.PWMscan.done', 
        'sites/operations/TFAP4.PWMscan.done', 
        'sites/operations/ARNT.PWMscan.done', 
        'sites/operations/ATF6.PWMscan.done', 
        'sites/operations/BACH1.PWMscan.done', 
        'sites/operations/BACH2.PWMscan.done', 
        'sites/operations/CREB1.PWMscan.done', 
        'sites/operations/ATF2.PWMscan.done', 
        'sites/operations/TCF3.PWMscan.done', 
        'sites/operations/XBP1.PWMscan.done', 
        'sites/operations/ARID5B.PWMscan.done', 
        'sites/operations/MYOD1.PWMscan.done', 
        'sites/operations/NFE2.PWMscan.done', 
        'sites/operations/MYCN.PWMscan.done', 
        'sites/operations/NFE2L1.PWMscan.done', 
        'sites/operations/TEF.PWMscan.done', 
        'sites/operations/ATF3.PWMscan.done', 
        'sites/operations/BATF.PWMscan.done', 
        'sites/operations/TCF12.PWMscan.done'
    output:
        'sites/operations/PWMscan.group2.done'
    shell:
        'touch {output}'

rule PWMscan_group3:
    input:
        'sites/operations/MYC.PWMscan.done', 
        'sites/operations/MXI1.PWMscan.done', 
        'sites/operations/BHLHE40.PWMscan.done', 
        'sites/operations/ARNTL.PWMscan.done', 
        'sites/operations/ATF4.PWMscan.done', 
        'sites/operations/ATF7.PWMscan.done', 
        'sites/operations/BATF3.PWMscan.done', 
        'sites/operations/BHLHA15.PWMscan.done', 
        'sites/operations/BHLHE41.PWMscan.done', 
        'sites/operations/BHLHE22.PWMscan.done', 
        'sites/operations/BHLHE23.PWMscan.done', 
        'sites/operations/CEBPD.PWMscan.done', 
        'sites/operations/CEBPE.PWMscan.done', 
        'sites/operations/CEBPG.PWMscan.done', 
        'sites/operations/CLOCK.PWMscan.done', 
        'sites/operations/CREB3.PWMscan.done', 
        'sites/operations/CREB3L1.PWMscan.done', 
        'sites/operations/DBP.PWMscan.done', 
        'sites/operations/FIGLA.PWMscan.done', 
        'sites/operations/HES5.PWMscan.done'
    output:
        'sites/operations/PWMscan.group3.done'
    shell:
        'touch {output}'

rule PWMscan_group4:
    input:
        'sites/operations/HES7.PWMscan.done', 
        'sites/operations/HEY1.PWMscan.done', 
        'sites/operations/HEY2.PWMscan.done', 
        'sites/operations/ID4.PWMscan.done', 
        'sites/operations/JDP2.PWMscan.done', 
        'sites/operations/MAFG.PWMscan.done', 
        'sites/operations/MESP1.PWMscan.done', 
        'sites/operations/MGA.PWMscan.done', 
        'sites/operations/MLX.PWMscan.done', 
        'sites/operations/MLXIPL.PWMscan.done', 
        'sites/operations/MNT.PWMscan.done', 
        'sites/operations/MSC.PWMscan.done', 
        'sites/operations/MYF6.PWMscan.done', 
        'sites/operations/NEUROD2.PWMscan.done', 
        'sites/operations/NEUROG2.PWMscan.done', 
        'sites/operations/NRL.PWMscan.done', 
        'sites/operations/OLIG1.PWMscan.done', 
        'sites/operations/OLIG2.PWMscan.done', 
        'sites/operations/OLIG3.PWMscan.done', 
        'sites/operations/TCF4.PWMscan.done'
    output:
        'sites/operations/PWMscan.group4.done'
    shell:
        'touch {output}'

rule PWMscan_group5:
    input:
        'sites/operations/TFAP2B.PWMscan.done', 
        'sites/operations/TFE3.PWMscan.done', 
        'sites/operations/TFEB.PWMscan.done', 
        'sites/operations/TFEC.PWMscan.done', 
        'sites/operations/TFAP2D.PWMscan.done', 
        'sites/operations/ARID3A.PWMscan.done', 
        'sites/operations/ARNT2.PWMscan.done', 
        'sites/operations/ATF1.PWMscan.done', 
        'sites/operations/ATF5.PWMscan.done', 
        'sites/operations/CREM.PWMscan.done', 
        'sites/operations/DDIT3.PWMscan.done', 
        'sites/operations/EPAS1.PWMscan.done', 
        'sites/operations/FOSB.PWMscan.done', 
        'sites/operations/HAND1.PWMscan.done', 
        'sites/operations/HES1.PWMscan.done', 
        'sites/operations/HIF1A.PWMscan.done', 
        'sites/operations/HMGA1.PWMscan.done', 
        'sites/operations/HMGA2.PWMscan.done', 
        'sites/operations/MAFA.PWMscan.done', 
        'sites/operations/MAFB.PWMscan.done'
    output:
        'sites/operations/PWMscan.group5.done'
    shell:
        'touch {output}'

rule PWMscan_group6:
    input:
        'sites/operations/MAF.PWMscan.done', 
        'sites/operations/MITF.PWMscan.done', 
        'sites/operations/MYOG.PWMscan.done', 
        'sites/operations/NEUROD1.PWMscan.done', 
        'sites/operations/NFE2L2.PWMscan.done', 
        'sites/operations/PTF1A.PWMscan.done', 
        'sites/operations/TAL1.PWMscan.done', 
        'sites/operations/TWIST1.PWMscan.done', 
        'sites/operations/AIRE.PWMscan.done', 
        'sites/operations/ALX1.PWMscan.done', 
        'sites/operations/ALX3.PWMscan.done', 
        'sites/operations/ALX4.PWMscan.done', 
        'sites/operations/ANDR.PWMscan.done', 
        'sites/operations/AP2A.PWMscan.done', 
        'sites/operations/AP2B.PWMscan.done', 
        'sites/operations/AP2C.PWMscan.done', 
        'sites/operations/AP2D.PWMscan.done', 
        'sites/operations/ARI3A.PWMscan.done', 
        'sites/operations/ARI5B.PWMscan.done', 
        'sites/operations/ARX.PWMscan.done'
    output:
        'sites/operations/PWMscan.group6.done'
    shell:
        'touch {output}'

rule PWMscan_group7:
    input:
        'sites/operations/ASCL2.PWMscan.done', 
        'sites/operations/ATF6A.PWMscan.done', 
        'sites/operations/ATOH1.PWMscan.done', 
        'sites/operations/BARH1.PWMscan.done', 
        'sites/operations/BARH2.PWMscan.done', 
        'sites/operations/BARX1.PWMscan.done', 
        'sites/operations/BARX2.PWMscan.done', 
        'sites/operations/BC11A.PWMscan.done', 
        'sites/operations/BCL6B.PWMscan.done', 
        'sites/operations/BCL6.PWMscan.done', 
        'sites/operations/BHA15.PWMscan.done', 
        'sites/operations/BHE22.PWMscan.done', 
        'sites/operations/BHE23.PWMscan.done', 
        'sites/operations/BHE40.PWMscan.done', 
        'sites/operations/BHE41.PWMscan.done', 
        'sites/operations/BMAL1.PWMscan.done', 
        'sites/operations/BPTF.PWMscan.done', 
        'sites/operations/BRAC.PWMscan.done', 
        'sites/operations/BRCA1.PWMscan.done', 
        'sites/operations/BSH.PWMscan.done'
    output:
        'sites/operations/PWMscan.group7.done'
    shell:
        'touch {output}'

rule PWMscan_group8:
    input:
        'sites/operations/CDC5L.PWMscan.done', 
        'sites/operations/CDX1.PWMscan.done', 
        'sites/operations/CDX2.PWMscan.done', 
        'sites/operations/CEBPZ.PWMscan.done', 
        'sites/operations/CENPB.PWMscan.done', 
        'sites/operations/COE1.PWMscan.done', 
        'sites/operations/COT1.PWMscan.done', 
        'sites/operations/COT2.PWMscan.done', 
        'sites/operations/CPEB1.PWMscan.done', 
        'sites/operations/CR3L1.PWMscan.done', 
        'sites/operations/CR3L2.PWMscan.done', 
        'sites/operations/CREB5.PWMscan.done', 
        'sites/operations/CRX.PWMscan.done', 
        'sites/operations/CTCFL.PWMscan.done', 
        'sites/operations/CTCF.PWMscan.done', 
        'sites/operations/CUX1.PWMscan.done', 
        'sites/operations/CUX2.PWMscan.done', 
        'sites/operations/CXXC1.PWMscan.done', 
        'sites/operations/DLX1.PWMscan.done', 
        'sites/operations/DLX2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group8.done'
    shell:
        'touch {output}'

rule PWMscan_group9:
    input:
        'sites/operations/DLX3.PWMscan.done', 
        'sites/operations/DLX4.PWMscan.done', 
        'sites/operations/DLX5.PWMscan.done', 
        'sites/operations/DLX6.PWMscan.done', 
        'sites/operations/DMBX1.PWMscan.done', 
        'sites/operations/DPRX.PWMscan.done', 
        'sites/operations/DRGX.PWMscan.done', 
        'sites/operations/DUXA.PWMscan.done', 
        'sites/operations/E2F1.PWMscan.done', 
        'sites/operations/E2F2.PWMscan.done', 
        'sites/operations/E2F3.PWMscan.done', 
        'sites/operations/E2F4.PWMscan.done', 
        'sites/operations/E2F5.PWMscan.done', 
        'sites/operations/E2F6.PWMscan.done', 
        'sites/operations/E2F7.PWMscan.done', 
        'sites/operations/E2F8.PWMscan.done', 
        'sites/operations/E4F1.PWMscan.done', 
        'sites/operations/EGR1.PWMscan.done', 
        'sites/operations/EGR2.PWMscan.done', 
        'sites/operations/EGR3.PWMscan.done'
    output:
        'sites/operations/PWMscan.group9.done'
    shell:
        'touch {output}'

rule PWMscan_group10:
    input:
        'sites/operations/EGR4.PWMscan.done', 
        'sites/operations/EHF.PWMscan.done', 
        'sites/operations/ELF1.PWMscan.done', 
        'sites/operations/ELF2.PWMscan.done', 
        'sites/operations/ELF3.PWMscan.done', 
        'sites/operations/ELF5.PWMscan.done', 
        'sites/operations/ELK1.PWMscan.done', 
        'sites/operations/ELK3.PWMscan.done', 
        'sites/operations/ELK4.PWMscan.done', 
        'sites/operations/EMX1.PWMscan.done', 
        'sites/operations/EMX2.PWMscan.done', 
        'sites/operations/EOMES.PWMscan.done', 
        'sites/operations/ERF.PWMscan.done', 
        'sites/operations/ERG.PWMscan.done', 
        'sites/operations/ERR1.PWMscan.done', 
        'sites/operations/ERR2.PWMscan.done', 
        'sites/operations/ERR3.PWMscan.done', 
        'sites/operations/ESR1.PWMscan.done', 
        'sites/operations/ESR2.PWMscan.done', 
        'sites/operations/ESX1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group10.done'
    shell:
        'touch {output}'

rule PWMscan_group11:
    input:
        'sites/operations/ETS1.PWMscan.done', 
        'sites/operations/ETS2.PWMscan.done', 
        'sites/operations/ETV1.PWMscan.done', 
        'sites/operations/ETV2.PWMscan.done', 
        'sites/operations/ETV3.PWMscan.done', 
        'sites/operations/ETV4.PWMscan.done', 
        'sites/operations/ETV5.PWMscan.done', 
        'sites/operations/ETV6.PWMscan.done', 
        'sites/operations/ETV7.PWMscan.done', 
        'sites/operations/EVI1.PWMscan.done', 
        'sites/operations/EVX1.PWMscan.done', 
        'sites/operations/EVX2.PWMscan.done', 
        'sites/operations/FEV.PWMscan.done', 
        'sites/operations/FLI1.PWMscan.done', 
        'sites/operations/FOXA1.PWMscan.done', 
        'sites/operations/FOXA2.PWMscan.done', 
        'sites/operations/FOXA3.PWMscan.done', 
        'sites/operations/FOXB1.PWMscan.done', 
        'sites/operations/FOXC1.PWMscan.done', 
        'sites/operations/FOXC2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group11.done'
    shell:
        'touch {output}'

rule PWMscan_group12:
    input:
        'sites/operations/FOXD1.PWMscan.done', 
        'sites/operations/FOXD2.PWMscan.done', 
        'sites/operations/FOXD3.PWMscan.done', 
        'sites/operations/FOXF1.PWMscan.done', 
        'sites/operations/FOXF2.PWMscan.done', 
        'sites/operations/FOXG1.PWMscan.done', 
        'sites/operations/FOXH1.PWMscan.done', 
        'sites/operations/FOXI1.PWMscan.done', 
        'sites/operations/FOXJ2.PWMscan.done', 
        'sites/operations/FOXJ3.PWMscan.done', 
        'sites/operations/FOXK1.PWMscan.done', 
        'sites/operations/FOXL1.PWMscan.done', 
        'sites/operations/FOXM1.PWMscan.done', 
        'sites/operations/FOXO1.PWMscan.done', 
        'sites/operations/FOXO3.PWMscan.done', 
        'sites/operations/FOXO4.PWMscan.done', 
        'sites/operations/FOXO6.PWMscan.done', 
        'sites/operations/FOXP2.PWMscan.done', 
        'sites/operations/FOXP3.PWMscan.done', 
        'sites/operations/FOXQ1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group12.done'
    shell:
        'touch {output}'

rule PWMscan_group13:
    input:
        'sites/operations/FUBP1.PWMscan.done', 
        'sites/operations/GABP1.PWMscan.done', 
        'sites/operations/GABPA.PWMscan.done', 
        'sites/operations/GATA1.PWMscan.done', 
        'sites/operations/GATA2.PWMscan.done', 
        'sites/operations/GATA3.PWMscan.done', 
        'sites/operations/GATA4.PWMscan.done', 
        'sites/operations/GATA5.PWMscan.done', 
        'sites/operations/GATA6.PWMscan.done', 
        'sites/operations/GBX1.PWMscan.done', 
        'sites/operations/GBX2.PWMscan.done', 
        'sites/operations/GCM1.PWMscan.done', 
        'sites/operations/GCM2.PWMscan.done', 
        'sites/operations/GCR.PWMscan.done', 
        'sites/operations/GFI1B.PWMscan.done', 
        'sites/operations/GFI1.PWMscan.done', 
        'sites/operations/GLI1.PWMscan.done', 
        'sites/operations/GLI2.PWMscan.done', 
        'sites/operations/GLI3.PWMscan.done', 
        'sites/operations/GLIS1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group13.done'
    shell:
        'touch {output}'

rule PWMscan_group14:
    input:
        'sites/operations/GLIS2.PWMscan.done', 
        'sites/operations/GLIS3.PWMscan.done', 
        'sites/operations/GMEB2.PWMscan.done', 
        'sites/operations/GRHL1.PWMscan.done', 
        'sites/operations/GSC2.PWMscan.done', 
        'sites/operations/GSC.PWMscan.done', 
        'sites/operations/GSX1.PWMscan.done', 
        'sites/operations/GSX2.PWMscan.done', 
        'sites/operations/HBP1.PWMscan.done', 
        'sites/operations/HEN1.PWMscan.done', 
        'sites/operations/HESX1.PWMscan.done', 
        'sites/operations/HIC1.PWMscan.done', 
        'sites/operations/HIC2.PWMscan.done', 
        'sites/operations/HINFP.PWMscan.done', 
        'sites/operations/HLTF.PWMscan.done', 
        'sites/operations/HMBX1.PWMscan.done', 
        'sites/operations/HME1.PWMscan.done', 
        'sites/operations/HME2.PWMscan.done', 
        'sites/operations/HMX1.PWMscan.done', 
        'sites/operations/HMX2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group14.done'
    shell:
        'touch {output}'

rule PWMscan_group15:
    input:
        'sites/operations/HMX3.PWMscan.done', 
        'sites/operations/HNF1A.PWMscan.done', 
        'sites/operations/HNF1B.PWMscan.done', 
        'sites/operations/HNF4A.PWMscan.done', 
        'sites/operations/HNF4G.PWMscan.done', 
        'sites/operations/HNF6.PWMscan.done', 
        'sites/operations/HOMEZ.PWMscan.done', 
        'sites/operations/HSF1.PWMscan.done', 
        'sites/operations/HSF2.PWMscan.done', 
        'sites/operations/HSF4.PWMscan.done', 
        'sites/operations/HSFY1.PWMscan.done', 
        'sites/operations/HTF4.PWMscan.done', 
        'sites/operations/HXA10.PWMscan.done', 
        'sites/operations/HXA11.PWMscan.done', 
        'sites/operations/HXA13.PWMscan.done', 
        'sites/operations/HXA1.PWMscan.done', 
        'sites/operations/HXA2.PWMscan.done', 
        'sites/operations/HXA5.PWMscan.done', 
        'sites/operations/HXA7.PWMscan.done', 
        'sites/operations/HXA9.PWMscan.done'
    output:
        'sites/operations/PWMscan.group15.done'
    shell:
        'touch {output}'

rule PWMscan_group16:
    input:
        'sites/operations/HXB13.PWMscan.done', 
        'sites/operations/HXB1.PWMscan.done', 
        'sites/operations/HXB2.PWMscan.done', 
        'sites/operations/HXB3.PWMscan.done', 
        'sites/operations/HXB6.PWMscan.done', 
        'sites/operations/HXB7.PWMscan.done', 
        'sites/operations/HXB8.PWMscan.done', 
        'sites/operations/HXC10.PWMscan.done', 
        'sites/operations/HXC11.PWMscan.done', 
        'sites/operations/HXC12.PWMscan.done', 
        'sites/operations/HXC13.PWMscan.done', 
        'sites/operations/HXC6.PWMscan.done', 
        'sites/operations/HXC8.PWMscan.done', 
        'sites/operations/HXD10.PWMscan.done', 
        'sites/operations/HXD11.PWMscan.done', 
        'sites/operations/HXD12.PWMscan.done', 
        'sites/operations/HXD13.PWMscan.done', 
        'sites/operations/HXD3.PWMscan.done', 
        'sites/operations/HXD4.PWMscan.done', 
        'sites/operations/HXD8.PWMscan.done'
    output:
        'sites/operations/PWMscan.group16.done'
    shell:
        'touch {output}'

rule PWMscan_group17:
    input:
        'sites/operations/HXD9.PWMscan.done', 
        'sites/operations/IKZF1.PWMscan.done', 
        'sites/operations/INSM1.PWMscan.done', 
        'sites/operations/IRF1.PWMscan.done', 
        'sites/operations/IRF2.PWMscan.done', 
        'sites/operations/IRF3.PWMscan.done', 
        'sites/operations/IRF4.PWMscan.done', 
        'sites/operations/IRF5.PWMscan.done', 
        'sites/operations/IRF7.PWMscan.done', 
        'sites/operations/IRF8.PWMscan.done', 
        'sites/operations/IRF9.PWMscan.done', 
        'sites/operations/IRX2.PWMscan.done', 
        'sites/operations/IRX3.PWMscan.done', 
        'sites/operations/ISL1.PWMscan.done', 
        'sites/operations/ISL2.PWMscan.done', 
        'sites/operations/ISX.PWMscan.done', 
        'sites/operations/ITF2.PWMscan.done', 
        'sites/operations/KAISO.PWMscan.done', 
        'sites/operations/KLF12.PWMscan.done', 
        'sites/operations/KLF13.PWMscan.done'
    output:
        'sites/operations/PWMscan.group17.done'
    shell:
        'touch {output}'

rule PWMscan_group18:
    input:
        'sites/operations/KLF14.PWMscan.done', 
        'sites/operations/KLF15.PWMscan.done', 
        'sites/operations/KLF16.PWMscan.done', 
        'sites/operations/KLF1.PWMscan.done', 
        'sites/operations/KLF3.PWMscan.done', 
        'sites/operations/KLF4.PWMscan.done', 
        'sites/operations/KLF6.PWMscan.done', 
        'sites/operations/KLF8.PWMscan.done', 
        'sites/operations/LBX2.PWMscan.done', 
        'sites/operations/LEF1.PWMscan.done', 
        'sites/operations/LHX2.PWMscan.done', 
        'sites/operations/LHX3.PWMscan.done', 
        'sites/operations/LHX4.PWMscan.done', 
        'sites/operations/LHX6.PWMscan.done', 
        'sites/operations/LHX8.PWMscan.done', 
        'sites/operations/LHX9.PWMscan.done', 
        'sites/operations/LMX1A.PWMscan.done', 
        'sites/operations/LMX1B.PWMscan.done', 
        'sites/operations/MAZ.PWMscan.done', 
        'sites/operations/MBD2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group18.done'
    shell:
        'touch {output}'

rule PWMscan_group19:
    input:
        'sites/operations/MCR.PWMscan.done', 
        'sites/operations/MECP2.PWMscan.done', 
        'sites/operations/MEF2A.PWMscan.done', 
        'sites/operations/MEF2B.PWMscan.done', 
        'sites/operations/MEF2C.PWMscan.done', 
        'sites/operations/MEF2D.PWMscan.done', 
        'sites/operations/MEIS1.PWMscan.done', 
        'sites/operations/MEIS2.PWMscan.done', 
        'sites/operations/MEIS3.PWMscan.done', 
        'sites/operations/MEOX1.PWMscan.done', 
        'sites/operations/MEOX2.PWMscan.done', 
        'sites/operations/MGAP.PWMscan.done', 
        'sites/operations/MIXL1.PWMscan.done', 
        'sites/operations/MLXPL.PWMscan.done', 
        'sites/operations/MNX1.PWMscan.done', 
        'sites/operations/MSX1.PWMscan.done', 
        'sites/operations/MSX2.PWMscan.done', 
        'sites/operations/MTF1.PWMscan.done', 
        'sites/operations/MUSC.PWMscan.done', 
        'sites/operations/MYBA.PWMscan.done'
    output:
        'sites/operations/PWMscan.group19.done'
    shell:
        'touch {output}'

rule PWMscan_group20:
    input:
        'sites/operations/MYBB.PWMscan.done', 
        'sites/operations/MYB.PWMscan.done', 
        'sites/operations/MZF1.PWMscan.done', 
        'sites/operations/NANOG.PWMscan.done', 
        'sites/operations/NDF1.PWMscan.done', 
        'sites/operations/NDF2.PWMscan.done', 
        'sites/operations/NF2L1.PWMscan.done', 
        'sites/operations/NF2L2.PWMscan.done', 
        'sites/operations/NFAC1.PWMscan.done', 
        'sites/operations/NFAC2.PWMscan.done', 
        'sites/operations/NFAC3.PWMscan.done', 
        'sites/operations/NFAC4.PWMscan.done', 
        'sites/operations/NFAT5.PWMscan.done', 
        'sites/operations/NFIA.PWMscan.done', 
        'sites/operations/NFIC.PWMscan.done', 
        'sites/operations/NFKB1.PWMscan.done', 
        'sites/operations/NFKB2.PWMscan.done', 
        'sites/operations/NFYA.PWMscan.done', 
        'sites/operations/NFYB.PWMscan.done', 
        'sites/operations/NFYC.PWMscan.done'
    output:
        'sites/operations/PWMscan.group20.done'
    shell:
        'touch {output}'

rule PWMscan_group21:
    input:
        'sites/operations/NGN2.PWMscan.done', 
        'sites/operations/NKX21.PWMscan.done', 
        'sites/operations/NKX22.PWMscan.done', 
        'sites/operations/NKX23.PWMscan.done', 
        'sites/operations/NKX25.PWMscan.done', 
        'sites/operations/NKX28.PWMscan.done', 
        'sites/operations/NKX31.PWMscan.done', 
        'sites/operations/NKX32.PWMscan.done', 
        'sites/operations/NKX61.PWMscan.done', 
        'sites/operations/NKX62.PWMscan.done', 
        'sites/operations/NOBOX.PWMscan.done', 
        'sites/operations/NOTO.PWMscan.done', 
        'sites/operations/NR0B1.PWMscan.done', 
        'sites/operations/NR1D1.PWMscan.done', 
        'sites/operations/NR1H2.PWMscan.done', 
        'sites/operations/NR1H4.PWMscan.done', 
        'sites/operations/NR1I2.PWMscan.done', 
        'sites/operations/NR1I3.PWMscan.done', 
        'sites/operations/NR2C1.PWMscan.done', 
        'sites/operations/NR2C2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group21.done'
    shell:
        'touch {output}'

rule PWMscan_group22:
    input:
        'sites/operations/NR2E1.PWMscan.done', 
        'sites/operations/NR2E3.PWMscan.done', 
        'sites/operations/NR2F6.PWMscan.done', 
        'sites/operations/NR4A1.PWMscan.done', 
        'sites/operations/NR4A2.PWMscan.done', 
        'sites/operations/NR4A3.PWMscan.done', 
        'sites/operations/NR5A2.PWMscan.done', 
        'sites/operations/NR6A1.PWMscan.done', 
        'sites/operations/NRF1.PWMscan.done', 
        'sites/operations/ONEC2.PWMscan.done', 
        'sites/operations/ONEC3.PWMscan.done', 
        'sites/operations/OTX1.PWMscan.done', 
        'sites/operations/OTX2.PWMscan.done', 
        'sites/operations/OVOL1.PWMscan.done', 
        'sites/operations/P53.PWMscan.done', 
        'sites/operations/P5F1B.PWMscan.done', 
        'sites/operations/P63.PWMscan.done', 
        'sites/operations/P73.PWMscan.done', 
        'sites/operations/PAX1.PWMscan.done', 
        'sites/operations/PAX2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group22.done'
    shell:
        'touch {output}'

rule PWMscan_group23:
    input:
        'sites/operations/PAX3.PWMscan.done', 
        'sites/operations/PAX4.PWMscan.done', 
        'sites/operations/PAX5.PWMscan.done', 
        'sites/operations/PAX6.PWMscan.done', 
        'sites/operations/PAX7.PWMscan.done', 
        'sites/operations/PAX8.PWMscan.done', 
        'sites/operations/PBX1.PWMscan.done', 
        'sites/operations/PBX2.PWMscan.done', 
        'sites/operations/PBX3.PWMscan.done', 
        'sites/operations/PDX1.PWMscan.done', 
        'sites/operations/PEBB.PWMscan.done', 
        'sites/operations/PHX2A.PWMscan.done', 
        'sites/operations/PHX2B.PWMscan.done', 
        'sites/operations/PIT1.PWMscan.done', 
        'sites/operations/PITX1.PWMscan.done', 
        'sites/operations/PITX2.PWMscan.done', 
        'sites/operations/PITX3.PWMscan.done', 
        'sites/operations/PKNX1.PWMscan.done', 
        'sites/operations/PKNX2.PWMscan.done', 
        'sites/operations/PLAG1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group23.done'
    shell:
        'touch {output}'

rule PWMscan_group24:
    input:
        'sites/operations/PLAL1.PWMscan.done', 
        'sites/operations/PO2F1.PWMscan.done', 
        'sites/operations/PO2F2.PWMscan.done', 
        'sites/operations/PO2F3.PWMscan.done', 
        'sites/operations/PO3F1.PWMscan.done', 
        'sites/operations/PO3F2.PWMscan.done', 
        'sites/operations/PO3F3.PWMscan.done', 
        'sites/operations/PO3F4.PWMscan.done', 
        'sites/operations/PO4F1.PWMscan.done', 
        'sites/operations/PO4F2.PWMscan.done', 
        'sites/operations/PO4F3.PWMscan.done', 
        'sites/operations/PO5F1.PWMscan.done', 
        'sites/operations/PO6F1.PWMscan.done', 
        'sites/operations/PO6F2.PWMscan.done', 
        'sites/operations/PPARA.PWMscan.done', 
        'sites/operations/PPARD.PWMscan.done', 
        'sites/operations/PPARG.PWMscan.done', 
        'sites/operations/PRD14.PWMscan.done', 
        'sites/operations/PRDM1.PWMscan.done', 
        'sites/operations/PRDM4.PWMscan.done'
    output:
        'sites/operations/PWMscan.group24.done'
    shell:
        'touch {output}'

rule PWMscan_group25:
    input:
        'sites/operations/PRGR.PWMscan.done', 
        'sites/operations/PROP1.PWMscan.done', 
        'sites/operations/PROX1.PWMscan.done', 
        'sites/operations/PRRX1.PWMscan.done', 
        'sites/operations/PRRX2.PWMscan.done', 
        'sites/operations/PURA.PWMscan.done', 
        'sites/operations/RARA.PWMscan.done', 
        'sites/operations/RARB.PWMscan.done', 
        'sites/operations/RARG.PWMscan.done', 
        'sites/operations/RAX2.PWMscan.done', 
        'sites/operations/RELB.PWMscan.done', 
        'sites/operations/REL.PWMscan.done', 
        'sites/operations/REST.PWMscan.done', 
        'sites/operations/RFX1.PWMscan.done', 
        'sites/operations/RFX2.PWMscan.done', 
        'sites/operations/RFX3.PWMscan.done', 
        'sites/operations/RFX4.PWMscan.done', 
        'sites/operations/RFX5.PWMscan.done', 
        'sites/operations/RHXF1.PWMscan.done', 
        'sites/operations/RORA.PWMscan.done'
    output:
        'sites/operations/PWMscan.group25.done'
    shell:
        'touch {output}'

rule PWMscan_group26:
    input:
        'sites/operations/RORG.PWMscan.done', 
        'sites/operations/RREB1.PWMscan.done', 
        'sites/operations/RUNX1.PWMscan.done', 
        'sites/operations/RUNX2.PWMscan.done', 
        'sites/operations/RUNX3.PWMscan.done', 
        'sites/operations/RXRA.PWMscan.done', 
        'sites/operations/RXRB.PWMscan.done', 
        'sites/operations/RXRG.PWMscan.done', 
        'sites/operations/RX.PWMscan.done', 
        'sites/operations/SCRT1.PWMscan.done', 
        'sites/operations/SCRT2.PWMscan.done', 
        'sites/operations/SHOX2.PWMscan.done', 
        'sites/operations/SHOX.PWMscan.done', 
        'sites/operations/SMAD1.PWMscan.done', 
        'sites/operations/SMAD2.PWMscan.done', 
        'sites/operations/SMAD3.PWMscan.done', 
        'sites/operations/SMAD4.PWMscan.done', 
        'sites/operations/SMRC1.PWMscan.done', 
        'sites/operations/SNAI1.PWMscan.done', 
        'sites/operations/SNAI2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group26.done'
    shell:
        'touch {output}'

rule PWMscan_group27:
    input:
        'sites/operations/SOX10.PWMscan.done', 
        'sites/operations/SOX11.PWMscan.done', 
        'sites/operations/SOX13.PWMscan.done', 
        'sites/operations/SOX15.PWMscan.done', 
        'sites/operations/SOX17.PWMscan.done', 
        'sites/operations/SOX18.PWMscan.done', 
        'sites/operations/SOX1.PWMscan.done', 
        'sites/operations/SOX21.PWMscan.done', 
        'sites/operations/SOX2.PWMscan.done', 
        'sites/operations/SOX3.PWMscan.done', 
        'sites/operations/SOX4.PWMscan.done', 
        'sites/operations/SOX5.PWMscan.done', 
        'sites/operations/SOX7.PWMscan.done', 
        'sites/operations/SOX8.PWMscan.done', 
        'sites/operations/SOX9.PWMscan.done', 
        'sites/operations/SP1.PWMscan.done', 
        'sites/operations/SP2.PWMscan.done', 
        'sites/operations/SP3.PWMscan.done', 
        'sites/operations/SP4.PWMscan.done', 
        'sites/operations/SPDEF.PWMscan.done'
    output:
        'sites/operations/PWMscan.group27.done'
    shell:
        'touch {output}'

rule PWMscan_group28:
    input:
        'sites/operations/SPI1.PWMscan.done', 
        'sites/operations/SPIB.PWMscan.done', 
        'sites/operations/SPIC.PWMscan.done', 
        'sites/operations/SPZ1.PWMscan.done', 
        'sites/operations/SRBP1.PWMscan.done', 
        'sites/operations/SRBP2.PWMscan.done', 
        'sites/operations/SRF.PWMscan.done', 
        'sites/operations/SRY.PWMscan.done', 
        'sites/operations/STA5A.PWMscan.done', 
        'sites/operations/STA5B.PWMscan.done', 
        'sites/operations/STAT1.PWMscan.done', 
        'sites/operations/STAT2.PWMscan.done', 
        'sites/operations/STAT3.PWMscan.done', 
        'sites/operations/STAT4.PWMscan.done', 
        'sites/operations/STAT6.PWMscan.done', 
        'sites/operations/STF1.PWMscan.done', 
        'sites/operations/SUH.PWMscan.done', 
        'sites/operations/TBP.PWMscan.done', 
        'sites/operations/TBR1.PWMscan.done', 
        'sites/operations/TBX15.PWMscan.done'
    output:
        'sites/operations/PWMscan.group28.done'
    shell:
        'touch {output}'

rule PWMscan_group29:
    input:
        'sites/operations/TBX19.PWMscan.done', 
        'sites/operations/TBX1.PWMscan.done', 
        'sites/operations/TBX20.PWMscan.done', 
        'sites/operations/TBX21.PWMscan.done', 
        'sites/operations/TBX2.PWMscan.done', 
        'sites/operations/TBX3.PWMscan.done', 
        'sites/operations/TBX4.PWMscan.done', 
        'sites/operations/TBX5.PWMscan.done', 
        'sites/operations/TCF7.PWMscan.done', 
        'sites/operations/TEAD1.PWMscan.done', 
        'sites/operations/TEAD3.PWMscan.done', 
        'sites/operations/TEAD4.PWMscan.done', 
        'sites/operations/TF2LX.PWMscan.done', 
        'sites/operations/TF65.PWMscan.done', 
        'sites/operations/TF7L1.PWMscan.done', 
        'sites/operations/TF7L2.PWMscan.done', 
        'sites/operations/TFCP2.PWMscan.done', 
        'sites/operations/TFDP1.PWMscan.done', 
        'sites/operations/TFE2.PWMscan.done', 
        'sites/operations/TGIF1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group29.done'
    shell:
        'touch {output}'

rule PWMscan_group30:
    input:
        'sites/operations/TGIF2.PWMscan.done', 
        'sites/operations/THAP1.PWMscan.done', 
        'sites/operations/THA.PWMscan.done', 
        'sites/operations/THB.PWMscan.done', 
        'sites/operations/TLX1.PWMscan.done', 
        'sites/operations/TWST1.PWMscan.done', 
        'sites/operations/TYY1.PWMscan.done', 
        'sites/operations/TYY2.PWMscan.done', 
        'sites/operations/UBIP1.PWMscan.done', 
        'sites/operations/UNC4.PWMscan.done', 
        'sites/operations/VAX1.PWMscan.done', 
        'sites/operations/VAX2.PWMscan.done', 
        'sites/operations/VDR.PWMscan.done', 
        'sites/operations/VENTX.PWMscan.done', 
        'sites/operations/VSX1.PWMscan.done', 
        'sites/operations/VSX2.PWMscan.done', 
        'sites/operations/WT1.PWMscan.done', 
        'sites/operations/YBOX1.PWMscan.done', 
        'sites/operations/ZBED1.PWMscan.done', 
        'sites/operations/ZBT18.PWMscan.done'
    output:
        'sites/operations/PWMscan.group30.done'
    shell:
        'touch {output}'

rule PWMscan_group31:
    input:
        'sites/operations/ZBT49.PWMscan.done', 
        'sites/operations/ZBT7A.PWMscan.done', 
        'sites/operations/ZBT7B.PWMscan.done', 
        'sites/operations/ZBTB4.PWMscan.done', 
        'sites/operations/ZBTB6.PWMscan.done', 
        'sites/operations/ZEB1.PWMscan.done', 
        'sites/operations/ZEP1.PWMscan.done', 
        'sites/operations/ZEP2.PWMscan.done', 
        'sites/operations/ZFHX3.PWMscan.done', 
        'sites/operations/ZFX.PWMscan.done', 
        'sites/operations/ZIC1.PWMscan.done', 
        'sites/operations/ZIC2.PWMscan.done', 
        'sites/operations/ZIC3.PWMscan.done', 
        'sites/operations/ZIC4.PWMscan.done', 
        'sites/operations/ZKSC1.PWMscan.done', 
        'sites/operations/ZKSC3.PWMscan.done', 
        'sites/operations/ZN143.PWMscan.done', 
        'sites/operations/ZN148.PWMscan.done', 
        'sites/operations/ZN219.PWMscan.done', 
        'sites/operations/ZN232.PWMscan.done'
    output:
        'sites/operations/PWMscan.group31.done'
    shell:
        'touch {output}'

rule PWMscan_group32:
    input:
        'sites/operations/ZN282.PWMscan.done', 
        'sites/operations/ZN333.PWMscan.done', 
        'sites/operations/ZN350.PWMscan.done', 
        'sites/operations/ZN384.PWMscan.done', 
        'sites/operations/ZN410.PWMscan.done', 
        'sites/operations/ZN423.PWMscan.done', 
        'sites/operations/ZN524.PWMscan.done', 
        'sites/operations/ZN589.PWMscan.done', 
        'sites/operations/ZN639.PWMscan.done', 
        'sites/operations/ZN652.PWMscan.done', 
        'sites/operations/ZN713.PWMscan.done', 
        'sites/operations/ZN740.PWMscan.done', 
        'sites/operations/ZN784.PWMscan.done', 
        'sites/operations/ZSC16.PWMscan.done', 
        'sites/operations/ZSCA4.PWMscan.done', 
        'sites/operations/ABCF2.PWMscan.done', 
        'sites/operations/A1CF.PWMscan.done', 
        'sites/operations/ACO1.PWMscan.done', 
        'sites/operations/ADARB1.PWMscan.done', 
        'sites/operations/AFF4.PWMscan.done'
    output:
        'sites/operations/PWMscan.group32.done'
    shell:
        'touch {output}'

rule PWMscan_group33:
    input:
        'sites/operations/AGGF1.PWMscan.done', 
        'sites/operations/AKR1A1.PWMscan.done', 
        'sites/operations/ANXA1.PWMscan.done', 
        'sites/operations/ANXA11.PWMscan.done', 
        'sites/operations/APEX2.PWMscan.done', 
        'sites/operations/ARFGAP1.PWMscan.done', 
        'sites/operations/ASCC1.PWMscan.done', 
        'sites/operations/ASPSCR1.PWMscan.done', 
        'sites/operations/AVEN.PWMscan.done', 
        'sites/operations/BAD.PWMscan.done', 
        'sites/operations/GPANK1.PWMscan.done', 
        'sites/operations/BAX.PWMscan.done', 
        'sites/operations/BCL11A.PWMscan.done', 
        'sites/operations/BOLL.PWMscan.done', 
        'sites/operations/CELF4.PWMscan.done', 
        'sites/operations/CELF5.PWMscan.done', 
        'sites/operations/CELF6.PWMscan.done', 
        'sites/operations/C19orf25.PWMscan.done', 
        'sites/operations/C19orf40.PWMscan.done', 
        'sites/operations/EXO5.PWMscan.done'
    output:
        'sites/operations/PWMscan.group33.done'
    shell:
        'touch {output}'

rule PWMscan_group34:
    input:
        'sites/operations/LINC00471.PWMscan.done', 
        'sites/operations/C9orf156.PWMscan.done', 
        'sites/operations/CANX.PWMscan.done', 
        'sites/operations/CAT.PWMscan.done', 
        'sites/operations/CBFA2T2.PWMscan.done', 
        'sites/operations/CBFB.PWMscan.done', 
        'sites/operations/CBX7.PWMscan.done', 
        'sites/operations/ZNF830.PWMscan.done', 
        'sites/operations/CCDC25.PWMscan.done', 
        'sites/operations/CD59.PWMscan.done', 
        'sites/operations/CDK2AP1.PWMscan.done', 
        'sites/operations/AGAP2.PWMscan.done', 
        'sites/operations/CFL2.PWMscan.done', 
        'sites/operations/FOXN3.PWMscan.done', 
        'sites/operations/CKMT1B.PWMscan.done', 
        'sites/operations/CLK1.PWMscan.done', 
        'sites/operations/CNOT6.PWMscan.done', 
        'sites/operations/NELFB.PWMscan.done', 
        'sites/operations/CPSF4.PWMscan.done', 
        'sites/operations/CSNK2B.PWMscan.done'
    output:
        'sites/operations/PWMscan.group34.done'
    shell:
        'touch {output}'

rule PWMscan_group35:
    input:
        'sites/operations/CSTF2.PWMscan.done', 
        'sites/operations/CYB5R1.PWMscan.done', 
        'sites/operations/CYCS.PWMscan.done', 
        'sites/operations/DAB2.PWMscan.done', 
        'sites/operations/DAZAP1.PWMscan.done', 
        'sites/operations/ASAP3.PWMscan.done', 
        'sites/operations/DDX20.PWMscan.done', 
        'sites/operations/DDX4.PWMscan.done', 
        'sites/operations/DDX43.PWMscan.done', 
        'sites/operations/DDX53.PWMscan.done', 
        'sites/operations/DGCR8.PWMscan.done', 
        'sites/operations/DHX36.PWMscan.done', 
        'sites/operations/DIABLO.PWMscan.done', 
        'sites/operations/DIS3.PWMscan.done', 
        'sites/operations/DNMT3A.PWMscan.done', 
        'sites/operations/DTL.PWMscan.done', 
        'sites/operations/DUS3L.PWMscan.done', 
        'sites/operations/DUSP22.PWMscan.done', 
        'sites/operations/DUSP26.PWMscan.done', 
        'sites/operations/ECSIT.PWMscan.done'
    output:
        'sites/operations/PWMscan.group35.done'
    shell:
        'touch {output}'

rule PWMscan_group36:
    input:
        'sites/operations/EDN1.PWMscan.done', 
        'sites/operations/EEF1D.PWMscan.done', 
        'sites/operations/EIF5A2.PWMscan.done', 
        'sites/operations/ENO1.PWMscan.done', 
        'sites/operations/ESRRA.PWMscan.done', 
        'sites/operations/ETFB.PWMscan.done', 
        'sites/operations/EWSR1.PWMscan.done', 
        'sites/operations/EXOSC3.PWMscan.done', 
        'sites/operations/METTL21B.PWMscan.done', 
        'sites/operations/FAM127B.PWMscan.done', 
        'sites/operations/FEZ1.PWMscan.done', 
        'sites/operations/FEZF2.PWMscan.done', 
        'sites/operations/FGF19.PWMscan.done', 
        'sites/operations/FHL2.PWMscan.done', 
        'sites/operations/FIP1L1.PWMscan.done', 
        'sites/operations/SRRM3.PWMscan.done', 
        'sites/operations/FOXP4.PWMscan.done', 
        'sites/operations/GADD45A.PWMscan.done', 
        'sites/operations/GIT2.PWMscan.done', 
        'sites/operations/GLYCTK.PWMscan.done'
    output:
        'sites/operations/PWMscan.group36.done'
    shell:
        'touch {output}'

rule PWMscan_group37:
    input:
        'sites/operations/GOT1.PWMscan.done', 
        'sites/operations/GPAM.PWMscan.done', 
        'sites/operations/GPD1.PWMscan.done', 
        'sites/operations/GRHPR.PWMscan.done', 
        'sites/operations/GTF2B.PWMscan.done', 
        'sites/operations/GTF2H3.PWMscan.done', 
        'sites/operations/GTF3C2.PWMscan.done', 
        'sites/operations/GTF3C5.PWMscan.done', 
        'sites/operations/GTPBP1.PWMscan.done', 
        'sites/operations/GTPBP6.PWMscan.done', 
        'sites/operations/H1FX.PWMscan.done', 
        'sites/operations/H2AFY.PWMscan.done', 
        'sites/operations/H2AFZ.PWMscan.done', 
        'sites/operations/HCFC2.PWMscan.done', 
        'sites/operations/HCLS1.PWMscan.done', 
        'sites/operations/HDAC8.PWMscan.done', 
        'sites/operations/HHAT.PWMscan.done', 
        'sites/operations/HHEX.PWMscan.done', 
        'sites/operations/UBE2K.PWMscan.done', 
        'sites/operations/HIRIP3.PWMscan.done'
    output:
        'sites/operations/PWMscan.group37.done'
    shell:
        'touch {output}'

rule PWMscan_group38:
    input:
        'sites/operations/HIST1H2BN.PWMscan.done', 
        'sites/operations/HIST2H2AB.PWMscan.done', 
        'sites/operations/HIST2H2BE.PWMscan.done', 
        'sites/operations/HLCS.PWMscan.done', 
        'sites/operations/HMG20A.PWMscan.done', 
        'sites/operations/HNRNPA0.PWMscan.done', 
        'sites/operations/HNRNPA1.PWMscan.done', 
        'sites/operations/HNRNPC.PWMscan.done', 
        'sites/operations/HNRNPH3.PWMscan.done', 
        'sites/operations/HNRNPLL.PWMscan.done', 
        'sites/operations/HOXB13.PWMscan.done', 
        'sites/operations/HOXB9.PWMscan.done', 
        'sites/operations/HOXD3.PWMscan.done', 
        'sites/operations/HP1BP3.PWMscan.done', 
        'sites/operations/HSPA1L.PWMscan.done', 
        'sites/operations/HSPA5.PWMscan.done', 
        'sites/operations/HTATIP2.PWMscan.done', 
        'sites/operations/ID2.PWMscan.done', 
        'sites/operations/IL24.PWMscan.done', 
        'sites/operations/ING3.PWMscan.done'
    output:
        'sites/operations/PWMscan.group38.done'
    shell:
        'touch {output}'

rule PWMscan_group39:
    input:
        'sites/operations/IRF6.PWMscan.done', 
        'sites/operations/IVD.PWMscan.done', 
        'sites/operations/KDM5A.PWMscan.done', 
        'sites/operations/KDM5D.PWMscan.done', 
        'sites/operations/KCNIP1.PWMscan.done', 
        'sites/operations/KIAA0907.PWMscan.done', 
        'sites/operations/KIF22.PWMscan.done', 
        'sites/operations/LARP1.PWMscan.done', 
        'sites/operations/LARP4.PWMscan.done', 
        'sites/operations/LAS1L.PWMscan.done', 
        'sites/operations/CERS4.PWMscan.done', 
        'sites/operations/UBXN1.PWMscan.done', 
        'sites/operations/CBX3.PWMscan.done', 
        'sites/operations/LRRFIP1.PWMscan.done', 
        'sites/operations/LSM6.PWMscan.done', 
        'sites/operations/LUZP1.PWMscan.done', 
        'sites/operations/LUZP2.PWMscan.done', 
        'sites/operations/MAGEA8.PWMscan.done', 
        'sites/operations/MAGED4B.PWMscan.done', 
        'sites/operations/MAGEF1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group39.done'
    shell:
        'touch {output}'

rule PWMscan_group40:
    input:
        'sites/operations/MAGOH.PWMscan.done', 
        'sites/operations/MAP4K2.PWMscan.done', 
        'sites/operations/MAPK1.PWMscan.done', 
        'sites/operations/MBTPS2.PWMscan.done', 
        'sites/operations/MCTP2.PWMscan.done', 
        'sites/operations/MDM2.PWMscan.done', 
        'sites/operations/GLTPD1.PWMscan.done', 
        'sites/operations/RBM42.PWMscan.done', 
        'sites/operations/WDR83.PWMscan.done', 
        'sites/operations/MORN1.PWMscan.done', 
        'sites/operations/MRPL1.PWMscan.done', 
        'sites/operations/MRPL2.PWMscan.done', 
        'sites/operations/MRPS25.PWMscan.done', 
        'sites/operations/MSI1.PWMscan.done', 
        'sites/operations/MSI2.PWMscan.done', 
        'sites/operations/MSRA.PWMscan.done', 
        'sites/operations/MSRB3.PWMscan.done', 
        'sites/operations/MTHFD1.PWMscan.done', 
        'sites/operations/MXD4.PWMscan.done', 
        'sites/operations/MYEF2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group40.done'
    shell:
        'touch {output}'

rule PWMscan_group41:
    input:
        'sites/operations/MYLK.PWMscan.done', 
        'sites/operations/NANOS1.PWMscan.done', 
        'sites/operations/NAP1L1.PWMscan.done', 
        'sites/operations/NCALD.PWMscan.done', 
        'sites/operations/NCBP2.PWMscan.done', 
        'sites/operations/NFATC3.PWMscan.done', 
        'sites/operations/NFATC4.PWMscan.done', 
        'sites/operations/NFIB.PWMscan.done', 
        'sites/operations/NFIX.PWMscan.done', 
        'sites/operations/NME1.PWMscan.done', 
        'sites/operations/NMI.PWMscan.done', 
        'sites/operations/NMRAL1.PWMscan.done', 
        'sites/operations/NNT.PWMscan.done', 
        'sites/operations/NOC2L.PWMscan.done', 
        'sites/operations/GAR1.PWMscan.done', 
        'sites/operations/NONO.PWMscan.done', 
        'sites/operations/NR2F1.PWMscan.done', 
        'sites/operations/NUCB1.PWMscan.done', 
        'sites/operations/NUP107.PWMscan.done', 
        'sites/operations/NUP133.PWMscan.done'
    output:
        'sites/operations/PWMscan.group41.done'
    shell:
        'touch {output}'

rule PWMscan_group42:
    input:
        'sites/operations/NXPH3.PWMscan.done', 
        'sites/operations/ODC1.PWMscan.done', 
        'sites/operations/OTUD4.PWMscan.done', 
        'sites/operations/P4HB.PWMscan.done', 
        'sites/operations/PAXIP1.PWMscan.done', 
        'sites/operations/PCK2.PWMscan.done', 
        'sites/operations/PDCD11.PWMscan.done', 
        'sites/operations/PDE6H.PWMscan.done', 
        'sites/operations/PDLIM5.PWMscan.done', 
        'sites/operations/PGAM2.PWMscan.done', 
        'sites/operations/PHLDA2.PWMscan.done', 
        'sites/operations/PHOX2A.PWMscan.done', 
        'sites/operations/PHTF1.PWMscan.done', 
        'sites/operations/PICK1.PWMscan.done', 
        'sites/operations/PIK3C3.PWMscan.done', 
        'sites/operations/PIR.PWMscan.done', 
        'sites/operations/PKM.PWMscan.done', 
        'sites/operations/PKNOX2.PWMscan.done', 
        'sites/operations/PLAGL1.PWMscan.done', 
        'sites/operations/PLG.PWMscan.done'
    output:
        'sites/operations/PWMscan.group42.done'
    shell:
        'touch {output}'

rule PWMscan_group43:
    input:
        'sites/operations/POLE3.PWMscan.done', 
        'sites/operations/POLI.PWMscan.done', 
        'sites/operations/POU3F2.PWMscan.done', 
        'sites/operations/POU4F3.PWMscan.done', 
        'sites/operations/PPP1R10.PWMscan.done', 
        'sites/operations/PPP2R3B.PWMscan.done', 
        'sites/operations/PPP5C.PWMscan.done', 
        'sites/operations/PQBP1.PWMscan.done', 
        'sites/operations/PRDX5.PWMscan.done', 
        'sites/operations/PRKRIR.PWMscan.done', 
        'sites/operations/PRNP.PWMscan.done', 
        'sites/operations/PSMA6.PWMscan.done', 
        'sites/operations/PSMC2.PWMscan.done', 
        'sites/operations/PTCD1.PWMscan.done', 
        'sites/operations/PTPMT1.PWMscan.done', 
        'sites/operations/PURG.PWMscan.done', 
        'sites/operations/R3HDM2.PWMscan.done', 
        'sites/operations/RAB14.PWMscan.done', 
        'sites/operations/RAB18.PWMscan.done', 
        'sites/operations/RAB2A.PWMscan.done'
    output:
        'sites/operations/PWMscan.group43.done'
    shell:
        'touch {output}'

rule PWMscan_group44:
    input:
        'sites/operations/RAB7A.PWMscan.done', 
        'sites/operations/RAN.PWMscan.done', 
        'sites/operations/RAX.PWMscan.done', 
        'sites/operations/RBBP5.PWMscan.done', 
        'sites/operations/RBBP9.PWMscan.done', 
        'sites/operations/RBM17.PWMscan.done', 
        'sites/operations/RBM22.PWMscan.done', 
        'sites/operations/RBM3.PWMscan.done', 
        'sites/operations/ESRP1.PWMscan.done', 
        'sites/operations/ESRP2.PWMscan.done', 
        'sites/operations/RBM7.PWMscan.done', 
        'sites/operations/RBM8A.PWMscan.done', 
        'sites/operations/RBFOX2.PWMscan.done', 
        'sites/operations/RBMS1.PWMscan.done', 
        'sites/operations/RFC2.PWMscan.done', 
        'sites/operations/RFC3.PWMscan.done', 
        'sites/operations/RFXANK.PWMscan.done', 
        'sites/operations/RIOK2.PWMscan.done', 
        'sites/operations/MEX3C.PWMscan.done', 
        'sites/operations/RNASEH2C.PWMscan.done'
    output:
        'sites/operations/PWMscan.group44.done'
    shell:
        'touch {output}'

rule PWMscan_group45:
    input:
        'sites/operations/RNF138.PWMscan.done', 
        'sites/operations/RPL35.PWMscan.done', 
        'sites/operations/RPL6.PWMscan.done', 
        'sites/operations/RPP25.PWMscan.done', 
        'sites/operations/RPS10.PWMscan.done', 
        'sites/operations/RPS4X.PWMscan.done', 
        'sites/operations/RPS6KA5.PWMscan.done', 
        'sites/operations/RUFY3.PWMscan.done', 
        'sites/operations/RUVBL1.PWMscan.done', 
        'sites/operations/SCAND2P.PWMscan.done', 
        'sites/operations/PDS5A.PWMscan.done', 
        'sites/operations/SCMH1.PWMscan.done', 
        'sites/operations/SEMA4A.PWMscan.done', 
        'sites/operations/SF1.PWMscan.done', 
        'sites/operations/SF3B1.PWMscan.done', 
        'sites/operations/SFT2D1.PWMscan.done', 
        'sites/operations/SLC18A1.PWMscan.done', 
        'sites/operations/SMAP2.PWMscan.done', 
        'sites/operations/SMCR7L.PWMscan.done', 
        'sites/operations/SMPX.PWMscan.done'
    output:
        'sites/operations/PWMscan.group45.done'
    shell:
        'touch {output}'

rule PWMscan_group46:
    input:
        'sites/operations/SMUG1.PWMscan.done', 
        'sites/operations/SNAPC4.PWMscan.done', 
        'sites/operations/SNAPC5.PWMscan.done', 
        'sites/operations/SND1.PWMscan.done', 
        'sites/operations/SNRNP70.PWMscan.done', 
        'sites/operations/SNRPB2.PWMscan.done', 
        'sites/operations/SOCS4.PWMscan.done', 
        'sites/operations/SOD1.PWMscan.done', 
        'sites/operations/SOX14.PWMscan.done', 
        'sites/operations/SPAG7.PWMscan.done', 
        'sites/operations/SPATS2.PWMscan.done', 
        'sites/operations/SPR.PWMscan.done', 
        'sites/operations/SRBD1.PWMscan.done', 
        'sites/operations/SRP9.PWMscan.done', 
        'sites/operations/SSBP3.PWMscan.done', 
        'sites/operations/SSX2.PWMscan.done', 
        'sites/operations/SSX3.PWMscan.done', 
        'sites/operations/STAU2.PWMscan.done', 
        'sites/operations/STUB1.PWMscan.done', 
        'sites/operations/SUCLG1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group46.done'
    shell:
        'touch {output}'

rule PWMscan_group47:
    input:
        'sites/operations/TAF1A.PWMscan.done', 
        'sites/operations/TAF9.PWMscan.done', 
        'sites/operations/TAGLN2.PWMscan.done', 
        'sites/operations/TBPL1.PWMscan.done', 
        'sites/operations/TCEAL2.PWMscan.done', 
        'sites/operations/TCEAL6.PWMscan.done', 
        'sites/operations/TFAM.PWMscan.done', 
        'sites/operations/TGIF2LX.PWMscan.done', 
        'sites/operations/THAP5.PWMscan.done', 
        'sites/operations/THRA.PWMscan.done', 
        'sites/operations/MED30.PWMscan.done', 
        'sites/operations/TIA1.PWMscan.done', 
        'sites/operations/TIMELESS.PWMscan.done', 
        'sites/operations/TIMM44.PWMscan.done', 
        'sites/operations/TIMM8A.PWMscan.done', 
        'sites/operations/TMSB4XP8.PWMscan.done', 
        'sites/operations/TOB2.PWMscan.done', 
        'sites/operations/TP73.PWMscan.done', 
        'sites/operations/TPI1.PWMscan.done', 
        'sites/operations/TPPP.PWMscan.done'
    output:
        'sites/operations/PWMscan.group47.done'
    shell:
        'touch {output}'

rule PWMscan_group48:
    input:
        'sites/operations/TRIM21.PWMscan.done', 
        'sites/operations/TRIM69.PWMscan.done', 
        'sites/operations/TRIP10.PWMscan.done', 
        'sites/operations/TRMT1.PWMscan.done', 
        'sites/operations/TROVE2.PWMscan.done', 
        'sites/operations/TSC22D4.PWMscan.done', 
        'sites/operations/TSN.PWMscan.done', 
        'sites/operations/TSNAX.PWMscan.done', 
        'sites/operations/TULP1.PWMscan.done', 
        'sites/operations/U2AF1.PWMscan.done', 
        'sites/operations/UBB.PWMscan.done', 
        'sites/operations/UBE2V1.PWMscan.done', 
        'sites/operations/UGP2.PWMscan.done', 
        'sites/operations/UQCRB.PWMscan.done', 
        'sites/operations/USP39.PWMscan.done', 
        'sites/operations/UTP18.PWMscan.done', 
        'sites/operations/VAMP3.PWMscan.done', 
        'sites/operations/EZR.PWMscan.done', 
        'sites/operations/VPS4B.PWMscan.done', 
        'sites/operations/NELFA.PWMscan.done'
    output:
        'sites/operations/PWMscan.group48.done'
    shell:
        'touch {output}'

rule PWMscan_group49:
    input:
        'sites/operations/WISP2.PWMscan.done', 
        'sites/operations/XG.PWMscan.done', 
        'sites/operations/XRCC1.PWMscan.done', 
        'sites/operations/YEATS4.PWMscan.done', 
        'sites/operations/YWHAE.PWMscan.done', 
        'sites/operations/YWHAZ.PWMscan.done', 
        'sites/operations/ZBTB12.PWMscan.done', 
        'sites/operations/ZBTB25.PWMscan.done', 
        'sites/operations/ZBTB43.PWMscan.done', 
        'sites/operations/ZBTB46.PWMscan.done', 
        'sites/operations/ZC3H7A.PWMscan.done', 
        'sites/operations/ZCCHC14.PWMscan.done', 
        'sites/operations/ZCCHC17.PWMscan.done', 
        'sites/operations/ZDHHC15.PWMscan.done', 
        'sites/operations/ZDHHC5.PWMscan.done', 
        'sites/operations/ZFP3.PWMscan.done', 
        'sites/operations/ZHX3.PWMscan.done', 
        'sites/operations/ZMAT2.PWMscan.done', 
        'sites/operations/ZMAT4.PWMscan.done', 
        'sites/operations/ZNF124.PWMscan.done'
    output:
        'sites/operations/PWMscan.group49.done'
    shell:
        'touch {output}'

rule PWMscan_group50:
    input:
        'sites/operations/ZNF131.PWMscan.done', 
        'sites/operations/ZNF160.PWMscan.done', 
        'sites/operations/ZKSCAN8.PWMscan.done', 
        'sites/operations/ZSCAN9.PWMscan.done', 
        'sites/operations/ZNF205.PWMscan.done', 
        'sites/operations/ZNF207.PWMscan.done', 
        'sites/operations/ZBTB18.PWMscan.done', 
        'sites/operations/ZNF250.PWMscan.done', 
        'sites/operations/ZNF26.PWMscan.done', 
        'sites/operations/ZNF3.PWMscan.done', 
        'sites/operations/ZNF304.PWMscan.done', 
        'sites/operations/RNF114.PWMscan.done', 
        'sites/operations/ZSCAN31.PWMscan.done', 
        'sites/operations/ZNF326.PWMscan.done', 
        'sites/operations/ZNF385A.PWMscan.done', 
        'sites/operations/ZNF503.PWMscan.done', 
        'sites/operations/ZNF510.PWMscan.done', 
        'sites/operations/ZNF655.PWMscan.done', 
        'sites/operations/ZNF671.PWMscan.done', 
        'sites/operations/ZNF695.PWMscan.done'
    output:
        'sites/operations/PWMscan.group50.done'
    shell:
        'touch {output}'

rule PWMscan_group51:
    input:
        'sites/operations/ZNF706.PWMscan.done', 
        'sites/operations/ZNF71.PWMscan.done', 
        'sites/operations/ZNF720.PWMscan.done', 
        'sites/operations/ZNF76.PWMscan.done', 
        'sites/operations/ZNF766.PWMscan.done', 
        'sites/operations/ZRSR2.PWMscan.done', 
        'sites/operations/ZSWIM1.PWMscan.done', 
        'sites/operations/Myf.PWMscan.done', 
        'sites/operations/Pax6.PWMscan.done', 
        'sites/operations/RORA_1.PWMscan.done', 
        'sites/operations/RORA_2.PWMscan.done', 
        'sites/operations/YY1.PWMscan.done', 
        'sites/operations/TP53.PWMscan.done', 
        'sites/operations/RELA.PWMscan.done', 
        'sites/operations/ZNF354C.PWMscan.done', 
        'sites/operations/MIZF.PWMscan.done', 
        'sites/operations/AP1.PWMscan.done', 
        'sites/operations/DUX4.PWMscan.done', 
        'sites/operations/FOXP1.PWMscan.done', 
        'sites/operations/POU2F2.PWMscan.done'
    output:
        'sites/operations/PWMscan.group51.done'
    shell:
        'touch {output}'

rule PWMscan_group52:
    input:
        'sites/operations/TCF7L2.PWMscan.done', 
        'sites/operations/TP63.PWMscan.done', 
        'sites/operations/ZBTB33.PWMscan.done', 
        'sites/operations/ZNF263.PWMscan.done', 
        'sites/operations/AR.PWMscan.done', 
        'sites/operations/KLF5.PWMscan.done', 
        'sites/operations/T.PWMscan.done', 
        'sites/operations/EN1.PWMscan.done', 
        'sites/operations/ZNF143.PWMscan.done', 
        'sites/operations/NR3C1.PWMscan.done', 
        'sites/operations/ESRRB.PWMscan.done', 
        'sites/operations/HOXA5.PWMscan.done', 
        'sites/operations/DMRT3.PWMscan.done', 
        'sites/operations/LBX1.PWMscan.done', 
        'sites/operations/POU6F1.PWMscan.done', 
        'sites/operations/BARHL2.PWMscan.done', 
        'sites/operations/ELF4.PWMscan.done', 
        'sites/operations/EN2.PWMscan.done', 
        'sites/operations/HOXA13.PWMscan.done', 
        'sites/operations/HOXC11.PWMscan.done'
    output:
        'sites/operations/PWMscan.group52.done'
    shell:
        'touch {output}'

rule PWMscan_group53:
    input:
        'sites/operations/ONECUT1.PWMscan.done', 
        'sites/operations/POU4F2.PWMscan.done', 
        'sites/operations/ZBTB7B.PWMscan.done', 
        'sites/operations/ZBTB7C.PWMscan.done', 
        'sites/operations/RHOXF1.PWMscan.done', 
        'sites/operations/UNCX.PWMscan.done', 
        'sites/operations/NR3C2.PWMscan.done', 
        'sites/operations/SP8.PWMscan.done', 
        'sites/operations/YY2.PWMscan.done', 
        'sites/operations/ZBTB7A.PWMscan.done', 
        'sites/operations/ZNF410.PWMscan.done', 
        'sites/operations/ZNF740.PWMscan.done', 
        'sites/operations/ONECUT2.PWMscan.done', 
        'sites/operations/ONECUT3.PWMscan.done', 
        'sites/operations/MYBL1.PWMscan.done', 
        'sites/operations/MYBL2.PWMscan.done', 
        'sites/operations/PAX9.PWMscan.done', 
        'sites/operations/PKNOX1.PWMscan.done', 
        'sites/operations/POU1F1.PWMscan.done', 
        'sites/operations/POU2F1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group53.done'
    shell:
        'touch {output}'

rule PWMscan_group54:
    input:
        'sites/operations/POU3F1.PWMscan.done', 
        'sites/operations/POU3F3.PWMscan.done', 
        'sites/operations/POU3F4.PWMscan.done', 
        'sites/operations/POU4F1.PWMscan.done', 
        'sites/operations/POU5F1B.PWMscan.done', 
        'sites/operations/POU6F2.PWMscan.done', 
        'sites/operations/HOXD12.PWMscan.done', 
        'sites/operations/BSX.PWMscan.done', 
        'sites/operations/HMBOX1.PWMscan.done', 
        'sites/operations/HOXA10.PWMscan.done', 
        'sites/operations/HOXA2.PWMscan.done', 
        'sites/operations/HOXB2.PWMscan.done', 
        'sites/operations/HOXB3.PWMscan.done', 
        'sites/operations/HOXC10.PWMscan.done', 
        'sites/operations/HOXC12.PWMscan.done', 
        'sites/operations/HOXC13.PWMscan.done', 
        'sites/operations/HOXD11.PWMscan.done', 
        'sites/operations/HOXD13.PWMscan.done', 
        'sites/operations/NFATC2.PWMscan.done', 
        'sites/operations/ASCL1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group54.done'
    shell:
        'touch {output}'

rule PWMscan_group55:
    input:
        'sites/operations/FOXK2.PWMscan.done', 
        'sites/operations/GRHL2.PWMscan.done', 
        'sites/operations/KLF9.PWMscan.done', 
        'sites/operations/NR2F2.PWMscan.done', 
        'sites/operations/POU5F1.PWMscan.done', 
        'sites/operations/RBPJ.PWMscan.done', 
        'sites/operations/SIX1.PWMscan.done', 
        'sites/operations/SIX2.PWMscan.done', 
        'sites/operations/TEAD2.PWMscan.done', 
        'sites/operations/ZNF24.PWMscan.done', 
        'sites/operations/ZNF384.PWMscan.done', 
        'sites/operations/ZNF282.PWMscan.done', 
        'sites/operations/ZSCAN4.PWMscan.done', 
        'sites/operations/RORB.PWMscan.done', 
        'sites/operations/RORC.PWMscan.done', 
        'sites/operations/TCF7L1.PWMscan.done', 
        'sites/operations/HINFP1.PWMscan.done', 
        'sites/operations/ZNF238.PWMscan.done', 
        'sites/operations/ZNF306.PWMscan.done', 
        'sites/operations/ZNF524.PWMscan.done'
    output:
        'sites/operations/PWMscan.group55.done'
    shell:
        'touch {output}'

rule PWMscan_group56:
    input:
        'sites/operations/ZNF75A.PWMscan.done', 
        'sites/operations/ZNF784.PWMscan.done', 
        'sites/operations/HSFY2.PWMscan.done', 
        'sites/operations/NFATC1.PWMscan.done', 
        'sites/operations/POU2F3.PWMscan.done', 
        'sites/operations/POU5F1P1.PWMscan.done', 
        'sites/operations/BHLHB2.PWMscan.done', 
        'sites/operations/BHLHB3.PWMscan.done', 
        'sites/operations/CART1.PWMscan.done', 
        'sites/operations/HOXA1.PWMscan.done', 
        'sites/operations/HOXB5.PWMscan.done', 
        'sites/operations/HOXD8.PWMscan.done', 
        'sites/operations/IRX5.PWMscan.done', 
        'sites/operations/PHOX2B.PWMscan.done', 
        'sites/operations/RAXL1.PWMscan.done', 
        'sites/operations/ESRRG.PWMscan.done', 
        'sites/operations/THRB.PWMscan.done', 
        'sites/operations/Trp53.PWMscan.done', 
        'sites/operations/Trp73.PWMscan.done', 
        'sites/operations/ZBTB49.PWMscan.done'
    output:
        'sites/operations/PWMscan.group56.done'
    shell:
        'touch {output}'

rule PWMscan_group57:
    input:
        'sites/operations/ZNF232.PWMscan.done', 
        'sites/operations/ZNF435.PWMscan.done', 
        'sites/operations/ZNF713.PWMscan.done', 
        'sites/operations/ARID5A.PWMscan.done', 
        'sites/operations/BARHL1.PWMscan.done', 
        'sites/operations/BBX.PWMscan.done', 
        'sites/operations/BCL3.PWMscan.done', 
        'sites/operations/CHD1.PWMscan.done', 
        'sites/operations/CHD2.PWMscan.done', 
        'sites/operations/CREB3L2.PWMscan.done', 
        'sites/operations/DBX2.PWMscan.done', 
        'sites/operations/DMC1.PWMscan.done', 
        'sites/operations/EBF3.PWMscan.done', 
        'sites/operations/EP300.PWMscan.done', 
        'sites/operations/EZH2.PWMscan.done', 
        'sites/operations/FOXJ1.PWMscan.done', 
        'sites/operations/FOXN1.PWMscan.done', 
        'sites/operations/GMEB1.PWMscan.done', 
        'sites/operations/GTF2F1.PWMscan.done', 
        'sites/operations/GTF2I.PWMscan.done'
    output:
        'sites/operations/PWMscan.group57.done'
    shell:
        'touch {output}'

rule PWMscan_group58:
    input:
        'sites/operations/GZF1.PWMscan.done', 
        'sites/operations/HCFC1.PWMscan.done', 
        'sites/operations/HDX.PWMscan.done', 
        'sites/operations/HIVEP1.PWMscan.done', 
        'sites/operations/HLX.PWMscan.done', 
        'sites/operations/HOXA11.PWMscan.done', 
        'sites/operations/HOXA3.PWMscan.done', 
        'sites/operations/HOXA4.PWMscan.done', 
        'sites/operations/HOXA6.PWMscan.done', 
        'sites/operations/HOXA7.PWMscan.done', 
        'sites/operations/HOXA9.PWMscan.done', 
        'sites/operations/HOXB1.PWMscan.done', 
        'sites/operations/HOXB4.PWMscan.done', 
        'sites/operations/HOXB6.PWMscan.done', 
        'sites/operations/HOXB7.PWMscan.done', 
        'sites/operations/HOXB8.PWMscan.done', 
        'sites/operations/HOXC4.PWMscan.done', 
        'sites/operations/HOXC5.PWMscan.done', 
        'sites/operations/HOXC6.PWMscan.done', 
        'sites/operations/HOXC8.PWMscan.done'
    output:
        'sites/operations/PWMscan.group58.done'
    shell:
        'touch {output}'

rule PWMscan_group59:
    input:
        'sites/operations/HOXC9.PWMscan.done', 
        'sites/operations/HOXD10.PWMscan.done', 
        'sites/operations/HOXD1.PWMscan.done', 
        'sites/operations/HOXD4.PWMscan.done', 
        'sites/operations/HOXD9.PWMscan.done', 
        'sites/operations/IKZF2.PWMscan.done', 
        'sites/operations/IRX4.PWMscan.done', 
        'sites/operations/IRX6.PWMscan.done', 
        'sites/operations/KLF7.PWMscan.done', 
        'sites/operations/LHX1.PWMscan.done', 
        'sites/operations/LHX5.PWMscan.done', 
        'sites/operations/MECOM.PWMscan.done', 
        'sites/operations/MTA3.PWMscan.done', 
        'sites/operations/OSR1.PWMscan.done', 
        'sites/operations/OSR2.PWMscan.done', 
        'sites/operations/OTP.PWMscan.done', 
        'sites/operations/PATZ1.PWMscan.done', 
        'sites/operations/PGR.PWMscan.done', 
        'sites/operations/PML.PWMscan.done', 
        'sites/operations/PRDM14.PWMscan.done'
    output:
        'sites/operations/PWMscan.group59.done'
    shell:
        'touch {output}'

rule PWMscan_group60:
    input:
        'sites/operations/RAD21.PWMscan.done', 
        'sites/operations/RCOR1.PWMscan.done', 
        'sites/operations/RFX7.PWMscan.done', 
        'sites/operations/RHOXF2.PWMscan.done', 
        'sites/operations/SIN3A.PWMscan.done', 
        'sites/operations/SIX3.PWMscan.done', 
        'sites/operations/SIX4.PWMscan.done', 
        'sites/operations/SIX5.PWMscan.done', 
        'sites/operations/SIX6.PWMscan.done', 
        'sites/operations/SMARCC1.PWMscan.done', 
        'sites/operations/SMARCC2.PWMscan.done', 
        'sites/operations/SMC3.PWMscan.done', 
        'sites/operations/SOX12.PWMscan.done', 
        'sites/operations/SOX30.PWMscan.done', 
        'sites/operations/SOX6.PWMscan.done', 
        'sites/operations/SP100.PWMscan.done', 
        'sites/operations/STAT5A.PWMscan.done', 
        'sites/operations/STAT5B.PWMscan.done', 
        'sites/operations/TAF1.PWMscan.done', 
        'sites/operations/TBL1XR1.PWMscan.done'
    output:
        'sites/operations/PWMscan.group60.done'
    shell:
        'touch {output}'

rule PWMscan_group61:
    input:
        'sites/operations/TCF21.PWMscan.done', 
        'sites/operations/TFAP2E.PWMscan.done', 
        'sites/operations/TFCP2L1.PWMscan.done', 
        'sites/operations/TLX2.PWMscan.done', 
        'sites/operations/UBP1.PWMscan.done', 
        'sites/operations/WRNIP1.PWMscan.done', 
        'sites/operations/YBX1.PWMscan.done', 
        'sites/operations/ZBTB14.PWMscan.done', 
        'sites/operations/ZBTB16.PWMscan.done', 
        'sites/operations/ZBTB3.PWMscan.done', 
        'sites/operations/ZKSCAN1.PWMscan.done', 
        'sites/operations/ZKSCAN3.PWMscan.done', 
        'sites/operations/ZNF148.PWMscan.done', 
        'sites/operations/ZNF219.PWMscan.done', 
        'sites/operations/ZNF274.PWMscan.done', 
        'sites/operations/ZNF281.PWMscan.done', 
        'sites/operations/ZNF333.PWMscan.done', 
        'sites/operations/ZNF350.PWMscan.done', 
        'sites/operations/ZNF35.PWMscan.done', 
        'sites/operations/ZNF423.PWMscan.done'
    output:
        'sites/operations/PWMscan.group61.done'
    shell:
        'touch {output}'

rule PWMscan_group62:
    input:
        'sites/operations/ZNF652.PWMscan.done', 
        'sites/operations/ZNF691.PWMscan.done', 
        'sites/operations/ZNF711.PWMscan.done', 
        'sites/operations/ZNF8.PWMscan.done', 
        'sites/operations/Sox4.PWMscan.done'
    output:
        'sites/operations/PWMscan.group62.done'
    shell:
        'touch {output}'