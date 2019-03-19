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

# rule run_h508wt01:
#     input:
#         "h508/wt01/preprocessing/logs/H508-WT-01.preprocessing.cleaning.done.txt"
# rule run_ls1034wt01:
#     input:
#         "ls1034/wt01/preprocessing/logs/LS1034-WT-01.preprocessing.cleaning.done.txt"
# rule run_snu61wt01:
#     input:
#         "snu61/wt01/preprocessing/logs/SNU61-WT-01.preprocessing.cleaning.done.txt"

# rule run_mdst8wt01:
#   input:
#       "mdst8/wt01/preprocessing/logs/MDST8-WT-01.preprocessing.cleaning.done.txt"

########################################################################################################################################
#### SPOOL FOOTPRINTING ################################################################################################################
########################################################################################################################################

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
#### RESET PIPELINE ####################################################################################################################
########################################################################################################################################
# Note - use these commands if you want to clear all the files associated with a specific processing run, if you plan to rerun it

########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################

rule test:
    input:
        "test01/preprocessing/12bigwig/test-repmerged.bw",
        "test01/preprocessing/12bigwig/test-REP1of2.bw",
        "test01/preprocessing/12bigwig/test-REP2of2.bw"

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
        mkdir -p -v {wildcards.path}saturation/complexity {wildcards.path}saturation/footprints {wildcards.path}saturation/peaks
        mkdir -p -v {wildcards.path}saturation/footprints/data
        mkdir -p -v {wildcards.path}saturation/footprints/data/merged {wildcards.path}saturation/footprints/data/motifmerge {wildcards.path}saturation/footprints/data/parsed {wildcards.path}saturation/footprints/data/bychr
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}footprints/graphs {wildcards.path}footprints/heatmaps {wildcards.path}footprints/data
        mkdir -p -v {wildcards.path}footprints/data/merged {wildcards.path}footprints/data/motifmerge {wildcards.path}footprints/parsed {wildcards.path}footprints/bychr
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

rule STEP17_makebigwig_bamcov_merged:
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
        b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai"
    output:
        "{path}preprocessing/12bigwig/{mergedsample}-repmerged.bw"
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"


# rule STEP16_MACS2_peaks_individual_globalnormilization:
#     # notes:
#     # because we are going to use the TCGA data downstream likely as a reference point,
#     # we will need to call the peaks in the exact same way as they did in this paper:
#     # http://science.sciencemag.org/content/sci/suppl/2018/10/24/362.6413.eaav1898.DC1/aav1898_Corces_SM.pdf
#     # which is "macs2 callpeak --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
#     ## params:
#     # -t input bam file (treatment)
#     # -n base name for output files
#     # --outdir output directory
#     # --shift find all tags in the bam, and shift them by 75 bp
#     # --extsize extend all shifted tags by 150 bp (should be roughly equal to avg frag size in lib)
#     # --nomodel do not use the macs2 function to determine shifting model
#     # --call-summits call the peak summits, detect subpeaks within a peaks
#     # --nolambda do not use local bias correction, use background nolambda
#     # --keep-dup all keep all duplicate reads (bam should be purged of PCR duplicates at this point)
#     # -p set the p-value cutoff for peak calling
#     input:
#         a="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam"
#         b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai"
#     output:
#         "{path}preprocessing/11peaks/{sample}-{REP}_peaks.xls"
#     shell:
#         "macs2 callpeak -t {input.a} -n {wildcards.sample}-{wildcards.REP} --outdir {wildcards.path}preprocessing/11peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# rule STEP17_callpeaksmacs2merged:
#     # notes:
#     # because we are going to use the TCGA data downstream likely as a reference point,
#     # we will need to call the peaks in the exact same way as they did in this paper:
#     # http://science.sciencemag.org/content/sci/suppl/2018/10/24/362.6413.eaav1898.DC1/aav1898_Corces_SM.pdf
#     # which is "macs2 callpeak --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
#     # params:
#     # -t input bam file (treatment)
#     # -n base name for output files
#     # --outdir output directory
#     # --shift find all tags in the bam, and shift them by 75 bp
#     # --extsize extend all shifted tags by 150 bp (should be roughly equal to avg frag size in lib)
#     # --nomodel do not use the macs2 function to determine shifting model
#     # --call-summits call the peak summits, detect subpeaks within a peaks
#     # --nolambda do not use local bias correction, use background nolambda
#     # --keep-dup all keep all duplicate reads (bam should be purged of PCR duplicates at this point)
#     # -p set the p-value cutoff for peak calling
#     input:
#         a="{path}preprocessing/12all/{sample}.all.bam",
#         b="{path}preprocessing/12all/{sample}.all.bai"
#     output:
#         "{path}preprocessing/13allpeaks/{sample}.all_peaks.xls"
#     shell:
#         "macs2 callpeak -t {input.a} -n {wildcards.sample}.all --outdir {wildcards.path}preprocessing/13allpeaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# rule STEP18_callpeaksmacs2merged_localnorm:
#     # notes:
#     # same as above, but with local background correction enabled
#     input:
#         a="{path}preprocessing/12all/{sample}.all.bam",
#         b="{path}preprocessing/12all/{sample}.all.bai"
#     output:
#         "{path}preprocessing/13allpeaks/{sample}.localnorm.all_peaks.xls"
#     shell:
#         "macs2 callpeak -t {input.a} -n {wildcards.sample}.localnorm.all --outdir {wildcards.path}preprocessing/13allpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

# rule STEP19_plotcorrspearman:
#     # parameters:
#     # -b input bam files
#     # -o output file name
#     # -bs set the bin size used for comparison, default is 10000 bp
#     # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
#     # -p set the number of computing processors to use
#     # -v verbose mode
#     input:
#         a="{path}preprocessing/10unique/{sample}-REP1.u.bam",
#         b="{path}preprocessing/10unique/{sample}-REP2.u.bam",
#         c="{path}preprocessing/10unique/{sample}-REP3.u.bam",
#         d="{path}preprocessing/10unique/{sample}-REP1.u.bai",
#         e="{path}preprocessing/10unique/{sample}-REP2.u.bai",
#         f="{path}preprocessing/10unique/{sample}-REP3.u.bai"
#     output:
#         "{path}preprocessing/14qcplots/{sample}.spearman.corrTest"
#     shell:
#         "multiBamSummary bins -b {input.a} {input.b} {input.c} -o {output} -bs 10000 -p 20 -v"

# rule STEP20_makecorrheatmap:
#     input:
#         "{path}preprocessing/14qcplots/{sample}.spearman.corrTest"
#     output:
#         "{path}preprocessing/14qcplots/{sample}.spearman.heatmap.svg"
#     shell:
#         "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"



# rule STEP23_downsamplebam:
#     # params:
#     # -Xmx50g set java mem limit to X gb
#     input:
#         "{path}preprocessing/8merged/{sample}.m.bam"
#     output:
#         "{path}preprocessing/15downsample/raw/{sample}.{prob}.bam"
#     shell:
#         "java -jar programs/picard/picard.jar DownsampleSam \
#         I={input} \
#         O={output} \
#         PROBABILITY=0.{wildcards.prob}"

# rule STEP24_sortdownsampled:
#     # params:
#     # -Xmx50g set java mem limit to X gb
#     input:
#         "{path}preprocessing/15downsample/raw/{sample}.{prob}.bam"
#     output:
#         "{path}preprocessing/15downsample/raw/{sample}.{prob}.cs.bam"
#     shell:
#         "java -jar programs/picard/picard.jar SortSam \
#         I={input} \
#         O={output} \
#         SORT_ORDER=coordinate"

# rule STEP25_markdupdownsampled:
#     # params:
#     # -Xmx50g set java mem limit to X gb
#     input:
#         "{path}preprocessing/15downsample/raw/{sample}.{prob}.cs.bam"
#     output:
#         a="{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bam",
#         b="{path}preprocessing/15downsample/complexity/{sample}.{prob}.dupmetrics.txt",
#     shell:
#         "java -Xmx5g -jar programs/picard/picard.jar MarkDuplicates \
#         I={input} \
#         O={output.a} \
#         M={output.b} \
#         REMOVE_DUPLICATES=true \
#         ASSUME_SORTED=true"

# rule STEP26_indexdownsampled:
#     input:
#         "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bam"
#     output:
#         "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bai"
#     shell:
#         "java -jar programs/picard/picard.jar BuildBamIndex \
#         I={input} \
#         O={output}"

# rule STEP27_callpeaksmacs2downsampled:
#     # notes:
#     # because we are going to use the TCGA data downstream likely as a reference point,
#     # we will need to call the peaks in the exact same way as they did in this paper:
#     # http://science.sciencemag.org/content/sci/suppl/2018/10/24/362.6413.eaav1898.DC1/aav1898_Corces_SM.pdf
#     # which is "macs2 callpeak --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
#     # params:
#     # -t input bam file (treatment)
#     # -n base name for output files
#     # --outdir output directory
#     # --shift find all tags in the bam, and shift them by 75 bp
#     # --extsize extend all shifted tags by 150 bp (should be roughly equal to avg frag size in lib)
#     # --nomodel do not use the macs2 function to determine shifting model
#     # --call-summits call the peak summits, detect subpeaks within a peaks
#     # --nolambda do not use local bias correction, use background nolambda
#     # --keep-dup all keep all duplicate reads (bam should be purged of PCR duplicates at this point)
#     # -p set the p-value cutoff for peak calling
#     input:
#         "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bam"
#     output:
#         "{path}preprocessing/15downsample/peaks/{sample}.{prob}_peaks.xls"
#     shell:
#         "macs2 callpeak -t {input} -n {wildcards.sample}.{wildcards.prob} --outdir preprocessing/15downsample/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# rule AGGREGATE_preprocessing:
#     input:
#         "{path}preprocessing/10unique/{sample}-REP1.u.bai",
#         "{path}preprocessing/10unique/{sample}-REP2.u.bai",
#         "{path}preprocessing/10unique/{sample}-REP3.u.bai",
#         "{path}preprocessing/12all/{sample}.all.bai",
#         "{path}preprocessing/11peaks/{sample}-REP1_peaks.xls",
#         "{path}preprocessing/11peaks/{sample}-REP2_peaks.xls",
#         "{path}preprocessing/11peaks/{sample}-REP3_peaks.xls",
#         "{path}preprocessing/13allpeaks/{sample}.all_peaks.xls",
#         "{path}preprocessing/13allpeaks/{sample}.localnorm.all_peaks.xls",
#         "{path}preprocessing/14qcplots/{sample}.spearman.heatmap.svg",
#         "{path}preprocessing/16bigwig/{sample}-REP1.u.bw",
#         "{path}preprocessing/16bigwig/{sample}-REP2.u.bw",
#         "{path}preprocessing/16bigwig/{sample}-REP3.u.bw",
#         "{path}preprocessing/16bigwig/{sample}.all.bw",
#         "{path}preprocessing/15downsample/complexity/{sample}.9.md.bai",
#         "{path}preprocessing/15downsample/complexity/{sample}.8.md.bai",
#         "{path}preprocessing/15downsample/complexity/{sample}.7.md.bai",
#         "{path}preprocessing/15downsample/complexity/{sample}.6.md.bai",
#         "{path}preprocessing/15downsample/complexity/{sample}.5.md.bai",
#         "{path}preprocessing/15downsample/complexity/{sample}.4.md.bai",
#         "{path}preprocessing/15downsample/complexity/{sample}.3.md.bai",
#         "{path}preprocessing/15downsample/complexity/{sample}.2.md.bai",
#         "{path}preprocessing/15downsample/complexity/{sample}.1.md.bai",
#         "{path}preprocessing/15downsample/peaks/{sample}.9_peaks.xls",
#         "{path}preprocessing/15downsample/peaks/{sample}.8_peaks.xls",
#         "{path}preprocessing/15downsample/peaks/{sample}.7_peaks.xls",
#         "{path}preprocessing/15downsample/peaks/{sample}.6_peaks.xls",
#         "{path}preprocessing/15downsample/peaks/{sample}.5_peaks.xls",
#         "{path}preprocessing/15downsample/peaks/{sample}.4_peaks.xls",
#         "{path}preprocessing/15downsample/peaks/{sample}.3_peaks.xls",
#         "{path}preprocessing/15downsample/peaks/{sample}.2_peaks.xls",
#         "{path}preprocessing/15downsample/peaks/{sample}.1_peaks.xls"
#     output:
#         "{path}preprocessing/logs/{sample}.preprocessing.done.txt"
#     shell:
#         "touch {wildcards.path}preprocessing/logs/{wildcards.sample}.preprocessing.done.txt"

# rule gather metrics (write me):
#               c="{path}preprocessing/6rawbam/mitochondrial_reads.txt"
      #  samtools view -c {input} chrM >> {output.c}
      #  samtools view -c {input} >> {output.c}

# rule CLEAN_preprocessing:
#     input:
#         "{path}preprocessing/logs/{sample}.preprocessing.done.txt"
#     output:
#         "{path}preprocessing/logs/{sample}.preprocessing.cleaning.done.txt"
#     shell:
#         """
#         rm {wildcards.path}preprocessing/2fastq/*.fastq
#         rm {wildcards.path}preprocessing/3goodfastq/*.fq
#         rm {wildcards.path}preprocessing/4mycoalign/*.sam
#         rm {wildcards.path}preprocessing/5hg38align/*.sam
#         rm {wildcards.path}preprocessing/6rawbam/*.bam
#         rm {wildcards.path}preprocessing/7rgsort/*.bam
#         rm {wildcards.path}preprocessing/8merged/*.bam
#         rm {wildcards.path}preprocessing/9dedup/*.bam
#         touch {output}
#         """

# ########################################################################################################################################
# #### Library Complexity Saturation Analysis Rules ######################################################################################
# ########################################################################################################################################
# rule analyzecomplexitysaturation:
#     input:
#         a="{path}preprocessing/15downsample/complexity/{sample}.9.dupmetrics.txt",
#         b="{path}preprocessing/15downsample/complexity/{sample}.8.dupmetrics.txt",
#         c="{path}preprocessing/15downsample/complexity/{sample}.7.dupmetrics.txt",
#         d="{path}preprocessing/15downsample/complexity/{sample}.6.dupmetrics.txt",
#         e="{path}preprocessing/15downsample/complexity/{sample}.5.dupmetrics.txt",
#         f="{path}preprocessing/15downsample/complexity/{sample}.4.dupmetrics.txt",
#         g="{path}preprocessing/15downsample/complexity/{sample}.3.dupmetrics.txt",
#         h="{path}preprocessing/15downsample/complexity/{sample}.2.dupmetrics.txt",
#         i="{path}preprocessing/15downsample/complexity/{sample}.1.dupmetrics.txt"
#     output:
#         "{path}saturation/{sample}.downsampled_lib_sizes.txt"
#     shell:
#         """
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.a} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.b} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.c} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.d} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.e} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.f} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.g} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.h} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.i} >> {output}
#         """
# ########################################################################################################################################
# #### Peaks Saturation Analysis Rules ###################################################################################################
# ########################################################################################################################################
# rule analyzepeaksaturation:
#     input:
#         "{path}15downsample/peaks/{mergedsample}.9_peaks.xls",
#         "{path}15downsample/peaks/{mergedsample}.8_peaks.xls",
#         "{path}15downsample/peaks/{mergedsample}.7_peaks.xls",
#         "{path}15downsample/peaks/{mergedsample}.6_peaks.xls",
#         "{path}15downsample/peaks/{mergedsample}.5_peaks.xls",
#         "{path}15downsample/peaks/{mergedsample}.4_peaks.xls",
#         "{path}15downsample/peaks/{mergedsample}.3_peaks.xls",
#         "{path}15downsample/peaks/{mergedsample}.2_peaks.xls",
#         "{path}15downsample/peaks/{mergedsample}.1_peaks.xls"
#     output:
#         "{path}15downsample/peaks/{mergedsample}.downsampled_numpeaks.txt"
#     shell:
#         "wl -l < {input} >> {output}"

# rule AGGREGATE_saturationanalysis:
#     input:
#         "{path}logs/{mergedsample}.preprocessing.done.txt",
#         "{path}15downsample/complexity/{mergedsample}.downsampled_lib_sizes.txt",
#         "{path}15downsample/peaks/{mergedsample}.downsampled_numpeaks.txt"
#     output:
#         "{path}logs/{mergedsample}.saturation_analysis.done.txt"
#     shell:
#         "touch {output}"

# ########################################################################################################################################
# #### Footprints Saturation Analysis Rules ##############################################################################################
# ########################################################################################################################################
# rule makefpbychr_downsampled:
#     input:
#         "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bam",
#         "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bai",
#         "sites/{gene}.sites.Rdata"
#     output:
#         "{path}preprocessing/15downsample/footprints/chr/{mergedsample}.{prob}.{gene}.{chr}.done.txt"
#     script:
#         "scripts/snakeMakeFPbyChrDownsampled.R"

# rule mergefpchr_downsampled:
#     input:
#         "sites/{gene}.sites.Rdata",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr1.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr2.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr3.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr4.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr5.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr6.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr7.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr8.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr9.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr10.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr11.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr12.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr13.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr14.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr15.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr16.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr17.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr18.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr19.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr20.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr21.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr22.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chrY.done.txt",
#         "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chrX.done.txt"
#     output:
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.{prob}.{gene}.done.merged.txt"
#     script:
#         "scripts/snakeMergeFPbyChrDownsampled.R"

# rule allprob_aggregator:
#     input:
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.9.{gene}.done.merged.txt",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.8.{gene}.done.merged.txt",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.7.{gene}.done.merged.txt",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.6.{gene}.done.merged.txt",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.5.{gene}.done.merged.txt",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.4.{gene}.done.merged.txt",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.3.{gene}.done.merged.txt",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.2.{gene}.done.merged.txt",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.1.{gene}.done.merged.txt"
#     output:
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.allprob.{gene}.done.txt"
#     shell:
#         "touch {output}"

# rule makefpgraph_downsampled:
#     input:
#         "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bam",
#         "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bai",
#         "sites/{gene}.sites.Rdata",
#         "{path}preprocessing/15downsample/footprints/merged/{sample}.{prob}.{gene}.done.merged.txt"
#     output:
#         "{path}saturation/footprints/{sample}.{prob}.{gene}.graphs.done.txt"
#     script:
#         "scripts/snakeGenerateMergedFPGraphDownsampled.R"

# rule allgraph_aggregator:
#     input:
#         "{path}saturation/footprints/{sample}.9.{gene}.graphs.done.txt",
#         "{path}saturation/footprints/{sample}.8.{gene}.graphs.done.txt",
#         "{path}saturation/footprints/{sample}.7.{gene}.graphs.done.txt",
#         "{path}saturation/footprints/{sample}.6.{gene}.graphs.done.txt",
#         "{path}saturation/footprints/{sample}.5.{gene}.graphs.done.txt",
#         "{path}saturation/footprints/{sample}.4.{gene}.graphs.done.txt",
#         "{path}saturation/footprints/{sample}.3.{gene}.graphs.done.txt",
#         "{path}saturation/footprints/{sample}.2.{gene}.graphs.done.txt",
#         "{path}saturation/footprints/{sample}.1.{gene}.graphs.done.txt"
#     output:
#         "{path}preprocessing/15downsample/footprints/graphs/{sample}.allprob.{gene}.done.graphs.txt"
#     shell:
#         "touch {output}"

# ########################################################################################################################################
# #### Footprint Analysis Rules ##########################################################################################################
# ########################################################################################################################################
# #rule run_fp_coadmr_mdst8_wt01:
# #    input:
# #        "mdst8/wt01/footprints/parsed/MDST8-WT-01.coadmr.parsed.done.txt"

# rule aggregate_coadmr_fp:
#   input:
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.CDX2.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.TCF7.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.HOXA3.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.MNX1.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.POU5F1B.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.OVOL1.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.ESRRA.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.ASCL2.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.HNF4A.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.GMEB2.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.ZSWIM1.parsed.done.txt",
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.CBFA2T2.parsed.done.txt"
#   output:
#       "mdst8/wt01/footprints/parsed/MDST8-WT-01.coadmr.parsed.done.txt"
#   shell:
#       "touch {output}"

# rule make_heatmaps:
#   input:
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif3.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.CBFA2T2.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ZSWIM1.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.MNX1.motif2.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.POU5F1B.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ASCL2.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.OVOL1.motif2.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.GMEB2.motif3.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.TCF7.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif9.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.OVOL1.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif6.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.TCF7.motif5.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.TCF7.motif6.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.GMEB2.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.TCF7.motif4.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif2.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.TCF7.motif2.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.HOXA3.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.GMEB2.motif2.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ASCL2.motif2.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.HNF4A.motif2.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.MNX1.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.TCF7.motif3.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.HNF4A.motif1.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif7.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif8.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif5.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.ESRRA.motif4.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.GMEB2.motif4.heatmap.svg",
#       "mdst8/wt01/footprints/heatmaps/MDST8-WT-01.CDX2.motif1.heatmap.svg"

# rule makefp_by_chr:
#     input:
#         "{path}preprocessing/12all/{mergedsample}.all.bam",
#         "{path}preprocessing/12all/{mergedsample}.all.bai",
#         "sites/{gene}.sites.Rdata"
#     output:
#         "{path}footprints/temp/{mergedsample}.{gene}.{chr}.done.txt"
#     script:
#         "scripts/snakeMakeFPbyChr.R"

# rule merge_chr:
#     input:
#         "sites/{gene}.sites.Rdata",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr1.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr2.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr3.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr4.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr5.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr6.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr7.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr8.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr9.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr10.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr11.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr12.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr13.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr14.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr15.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr16.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr17.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr18.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr19.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr20.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr21.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chr22.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chrX.done.txt",
#         "{path}footprints/temp/{mergedsample}.{gene}.chrY.done.txt"
#     output:
#         "{path}footprints/merged/{mergedsample}.{gene}.merged.done.txt"
#     script:
#         "scripts/snakeMergeFPbyChr.R"

# rule make_graphs:
#     input:
#         "{path}preprocessing/12all/{mergedsample}.all.bam",
#         "{path}preprocessing/12all/{mergedsample}.all.bai",
#         "sites/{gene}.sites.Rdata",
#         "{path}footprints/merged/{mergedsample}.{gene}.merged.done.txt"
#     output:
#         "{path}footprints/graphs/{mergedsample}.{gene}.graphs.done.txt"
#     script:
#         "scripts/snakeGenerateMergedFPGraph.R"

# rule parse_footprints:
#     input:
#         "{path}preprocessing/12all/{mergedsample}.all.bam",
#         "{path}preprocessing/12all/{mergedsample}.all.bai",
#         "sites/{gene}.sites.Rdata",
#         "{path}footprints/merged/{mergedsample}.{gene}.merged.done.txt",
#         "{path}preprocessing/13allpeaks/{mergedsample}.all_peaks.narrowPeak"
#     output:
#         "{path}footprints/parsed/{mergedsample}.{gene}.parsed.done.txt"
#     script:
#         "scripts/snakeParseFP.R"

# rule make_parsed_heatmaps:
#     input:
#         "{path}footprints/parsed/{mergedsample}.{gene}.motif{motif}.info.Rdata",
#     output:
#         "{path}footprints/heatmaps/{mergedsample}.{gene}.motif{motif}.heatmap.svg"
#     script:
#         "scripts/snakeFootprintHeatmaps.R"

# rule make_merged_motifs:
#     input:
#         "{path}parsed/{mergedsample}.{gene}.parsed.done.txt"
#     output:
#         "{path}merged_motifs/{mergedsample}.{gene}.{nummotif}.mergedmotif.Rdata"
#     script:
#         "scripts/snakeMergeMotifs.R"

# rule make_aracne_overlap:
#     input:
#         "{path}merged_motifs/{mergedsample}.{gene}.{nummotif}.mergedmotif.Rdata"
#     output:
#         "{path}aracne/{mergedsample}.{gene}.{nummotif}.{entrez}.aracne.Rdata"
#     script:
#         "scripts/snakeFindARACNeFootprintOverlap.R"

# ########################################################################################################################################
# #### Analysis Rules ####################################################################################################################
# ########################################################################################################################################
# rule threesample_plotcorrspearman:
#     # parameters:
#     # -b input bam files
#     # -o output file name
#     # -bs set the bin size used for comparison, default is 10000 bp
#     # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
#     # -p set the number of computing processors to use
#     # -v verbose mode
#     input:
#         a="{sample1}/{wt1}{num1}/preprocessing/12all/{s1}.bam",
#         b="{sample1}/{wt1}{num1}/preprocessing/12all/{s2}.bam",
#         c="{sample1}/{wt1}{num1}/preprocessing/12all/{s3}.bam"
#     output:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.corrTest"
#     shell:
#         "multiBamSummary bins -b {input.a} {input.b} {input.c} -o {output} -bs 10000 -p 20 -v"

# rule threesample_makecorrheatmap:
#     input:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.corrTest"
#     output:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.heatmap.svg"
#     shell:
#         "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"

# rule ninesample_plotcorrspearman:
#     # parameters:
#     # -b input bam files
#     # -o output file name
#     # -bs set the bin size used for comparison, default is 10000 bp
#     # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
#     # -p set the number of computing processors to use
#     # -v verbose mode
#     # For processing nine sample, restrict analysis to chr1, or computation will take forever
#     input:
#         a="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
#         b="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
#         c="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam",
#         d="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
#         e="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
#         f="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam",
#         g="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
#         h="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
#         i="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam"
#     output:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.REPS.spearman.corrTest"
#     shell:
#         "multiBamSummary bins -b {input.a} {input.b} {input.c} {input.d} {input.e} {input.f} {input.g} {input.h} {input.i} -o {output} -bs 10000 -p 20 -v -r chr1"

# rule ninesample_makecorrheatmap:
#     input:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.REPS.spearman.corrTest"
#     output:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.heatmap.svg"
#     shell:
#         "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"