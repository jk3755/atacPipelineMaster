###########################################################################################################################################
################################ GENERAL INFO #############################################################################################
###########################################################################################################################################
## Snakemake execution guide
# A dry run of the pipeline can be run with:
# snakemake -np h508go
# On one of the virtualization servers, run the pipeline with the following to allocate 20 threads and 90 gb max memory (to avoid crashing the process)
# snakemake -j 20 h508go --resources mem_gb=90
#
## Raw file info
# H508-1_S3_L001_R1_001.fastq.gz - Sample 1
# H508-2_S2_L001_R1_001.fastq.gz - Sample 2
# H508-3_S1_L001_R1_001.fastq.gz - Sample 3
## Before running pipeline, rename these to:
# H508-WT-01_REP1_L1_R1.fastq.gz
# H508-WT-01_REP2_L1_R1.fastq.gz
# H508-WT-01_REP3_L1_R1.fastq.gz
#
#######################################################################################################################
#### PIPELINE SPOOLING COMMANDS #######################################################################################
#######################################################################################################################
rule run_h508wt01:
    input:
        "h508/wt01/preprocessing/logs/H508-WT-01.preprocessing.cleaning.done.txt"
rule run_ls1034wt01:
    input:
        "ls1034/wt01/preprocessing/logs/LS1034-WT-01.preprocessing.cleaning.done.txt"
rule run_snu61wt01:
    input:
        "snu61/wt01/preprocessing/logs/SNU61-WT-01.preprocessing.cleaning.done.txt"
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
        "{path}preprocessing/logs/dirtree.built.done"
    shell:
        """
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}preprocessing/2fastq {wildcards.path}preprocessing/3goodfastq {wildcards.path}preprocessing/4mycoalign {wildcards.path}preprocessing/5hg38align
        mkdir -p -v {wildcards.path}preprocessing/6rawbam {wildcards.path}preprocessing/7rgsort {wildcards.path}preprocessing/8merged {wildcards.path}preprocessing/9dedup
        mkdir -p -v {wildcards.path}preprocessing/10unique {wildcards.path}preprocessing/11peaks {wildcards.path}preprocessing/12all {wildcards.path}preprocessing/13allpeaks
        mkdir -p -v {wildcards.path}preprocessing/14qcplots {wildcards.path}preprocessing/15downsample {wildcards.path}preprocessing/16bigwig {wildcards.path}preprocessing/logs
        mkdir -p -v {wildcards.path}preprocessing/15downsample/complexity {wildcards.path}preprocessing/15downsample/footprints {wildcards.path}preprocessing/15downsample/peaks {wildcards.path}preprocessing/15downsample/raw
        mkdir -p -v {wildcards.path}preprocessing/15downsample/footprints/chr {wildcards.path}preprocessing/15downsample/footprints/merged
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}footprints/graphs {wildcards.path}footprints/merged {wildcards.path}footprints/merged_motifs {wildcards.path}footprints/parsed {wildcards.path}footprints/temp
        mkdir -p -v {wildcards.path}saturation
        mkdir -p -v {wildcards.path}saturation/complexity {wildcards.path}saturation/footprints {wildcards.path}saturation/peaks
        mkdir -p -v {wildcards.path}correlation
        touch {output}
        """
rule STEP1_simplifynames_gunzip:
    # params: -k keep original files, -c write to standard output
    input: 
        a="{path}data/{sample}_L00{lane}_R{read}_001.fastq.gz",
        b="{path}preprocessing/logs/dirtree.built.done"
    output: 
        "{path}preprocessing/2fastq/{sample}_L{lane}_R{read}.fastq"
    log: 
        "{path}preprocessing/logs/{sample}.L{lane}.R{read}.simplifynames_gunzip.txt"
    shell: 
        "gunzip -k -c {input.a} > {output}"
rule STEP2_afterqc_fastqfiltering:
    # params: -s is the shortest trimmed read length allowed past QC filter
    input:
        a="{path}preprocessing/2fastq/{sample}_R1.fastq",
        b="{path}preprocessing/2fastq/{sample}_R2.fastq"
    output:
        c="{path}preprocessing/3goodfastq/{sample}_R1.good.fq",
        d="{path}preprocessing/3goodfastq/{sample}_R2.good.fq"
    log:
        "{path}preprocessing/logs/{sample}.afterqc_fastqfiltering.txt"
    shell:
        "after.py -1 {input.a} -2 {input.b} -g {wildcards.path}preprocessing/3goodfastq -b {wildcards.path}preprocessing/3goodfastq -s 15"
rule STEP3_mycoalign:
    # params:
    # -q fastq input
    # -p num threads
    # -X1000 align to a maximum of 1000 bp frag length
    # -1/2 inputs
    # -S output
    input:
        a="{path}preprocessing/3goodfastq/{sample}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}_R2.good.fq"
    output:
        "{path}preprocessing/4mycoalign/{sample}.myco.sam"
    log:
        "{path}preprocessing/logs/{sample}.mycoalign.txt"
    shell:
        "bowtie2 -q -p 20 -X1000 -x /home/ubuntu2/genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}preprocessing/4mycoalign/{wildcards.sample}alignment_metrics.txt"
rule STEP4_hg38align:
    # params:
    # -q fastq input
    # -p num threads
    # -X1000 align to a maximum of 1000 bp frag length
    # -1/2 inputs
    # -S output
    input:
        a="{path}preprocessing/3goodfastq/{sample}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}_R2.good.fq",
        c="{path}preprocessing/4mycoalign/{sample}.myco.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}.hg38.sam"
    log:
        "{path}preprocessing/logs/{sample}.hg38align.txt"
    shell:
        "bowtie2 -q -p 20 -X1000 -x /home/ubuntu2/genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}preprocessing/5hg38align/{wildcards.sample}alignment_metrics.txt"
rule STEP5_coordsort_sam:
    # coordinate sort the sam files to prepare for blacklist filtering
    input:
        "{path}preprocessing/5hg38align/{sample}.hg38.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}.hg38.cs.sam"
    shell:
        "samtools sort {input} -o {output} -O sam"
rule STEP6_blacklistfilter_bamconversion:
    # remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
    input:
        "{path}preprocessing/5hg38align/{sample}.hg38.cs.sam"
    output:
        a="{path}preprocessing/6rawbam/{sample}.blacklist.bam",
        b="{path}preprocessing/6rawbam/{sample}.bam"
    shell:
        "samtools view -b -h -o {output.a} -L /home/ubuntu2/genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 10 {input}"
rule STEP7_addrgandcsbam:
    # note - proper specification of RG tags is critical
    # see: https://software.broadinstitute.org/gatk/documentation/article.php?id=6472
    # Required @RG parameter specifications:
    # RGID (read group ID) - this must be a globally unique string. for illumina data, use flowcell + lane
    # RGLB (read group library) - This is used by MarkDuplicates to collect reads from the same library on different lanes, so it must be common to all files from the same library
    # RGPL (read group platform) - ILLUMINA
    # RGPU (read group platform unit) - The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
    # RGSM (read group sample name) - the name of the sample sequenced in this file. should be consistent across different files from different lanes
    input:
        "{path}preprocessing/6rawbam/{sample}_L{lane}.bam"
    output:
        "{path}preprocessing/7rgsort/{sample}_L{lane}.rg.cs.bam"
    log:
        "{path}preprocessing/logs/{sample}.L{lane}.addrgandcsbam.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID=H5YHHBGX3.{wildcards.lane} \
        RGLB={wildcards.sample} \
        RGPL=ILLUMINA \
        RGPU=H5YHHBGX3.{wildcards.lane}.{wildcards.sample} \
        RGSM={wildcards.sample}"
rule STEP8_cleanbam:
    # params:
    # -Xmx50g set java mem limit to X gb
    input:
        "{path}preprocessing/7rgsort/{sample}_L{lane}.rg.cs.bam"
    output:
        "{path}preprocessing/7rgsort/{sample}_L{lane}.clean.bam"
    log:
        "{path}preprocessing/logs/{sample}.L{lane}.cleanbam.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar CleanSam \
        I={input} \
        O={output}"
rule STEP9_mergelanes:
    input:
        a="{path}preprocessing/7rgsort/{sample}_L1.clean.bam",
        b="{path}preprocessing/7rgsort/{sample}_L2.clean.bam",
        c="{path}preprocessing/7rgsort/{sample}_L3.clean.bam",
        d="{path}preprocessing/7rgsort/{sample}_L4.clean.bam"
    output:
        "{path}preprocessing/8merged/{sample}.m.bam"
    log:
        "{path}preprocessing/logs/{sample}.mergelanes.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar MergeSamFiles \
        I={input.a} \
        I={input.b} \
        I={input.c} \
        I={input.d} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"
rule STEP10_purgeduplicates:
    # params:
    # -Xmx50g set java mem limit to X gb
    input:
        "{path}preprocessing/8merged/{sample}.m.bam"
    output:
        a="{path}preprocessing/9dedup/{sample}.dp.bam",
        b="{path}preprocessing/9dedup/{sample}.metrics.txt"
    log:
        "{path}preprocessing/logs/{sample}.purgeduplicates.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar MarkDuplicates \
        I={input} \
        O={output.a} \
        M={output.b} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"
rule STEP11_mapqfilter:
    # STEP 10 - REMOVE MULTI MAPPING READS WITH SAMTOOLS
    # Notes:
    # for an explanation of how bowtie2 calculates mapq scores:
    # http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
    # for bowtie2, mapq higher than 2 is a uniquely mapped read
    # params:
    # -h include the header in the output
    # -q only include reads with mapping quality X or higher
    # -b output as a bam file
    input:
        "{path}preprocessing/9dedup/{sample}.dp.bam"
    output:
        "{path}preprocessing/10unique/{sample}.u.bam"
    log:
        "{path}preprocessing/logs/{sample}.mapqfilter.txt"
    shell:
        "samtools view -h -q 2 -b {input} > {output}"
rule STEP12_buildindex:
    # params:
    # -XmxXg set java mem limit to X gb
    input:
        "{path}preprocessing/10unique/{sample}.u.bam"
    output:
        "{path}preprocessing/10unique/{sample}.u.bai"
    log:
        "{path}preprocessing/logs/{sample}.buildindex.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"
rule STEP13_mergereplicates:
    # params:
    # -Xmx50g set java mem limit to X gb
    input:
        a="{path}preprocessing/10unique/{sample}-REP1.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP2.u.bam",
        c="{path}preprocessing/10unique/{sample}-REP3.u.bam"
    output:
        "{path}preprocessing/12all/{sample}.all.bam"
    log:
        "{path}preprocessing/logs/{sample}.mergereplicates.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar MergeSamFiles \
        I={input.a} \
        I={input.b} \
        I={input.c} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"
rule STEP14_indexmerged:
    # params:
    # -Xmx50g set java mem limit to X gb
    input:
        "{path}preprocessing/12all/{mergedsample}.all.bam"
    output:
        "{path}preprocessing/12all/{mergedsample}.all.bai"
    log:
        "{path}preprocessing/logs/{mergedsample}.indexmerged.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"
rule STEP15_callpeaksmac2replicates:
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
        a="{path}preprocessing/10unique/{sample}-{REP}.u.bam",
        b="{path}preprocessing/10unique/{sample}-{REP}.u.bai"
    output:
        "{path}preprocessing/11peaks/{sample}-{REP}_peaks.xls"
    log:
        "{path}preprocessing/logs/{sample}-{REP}.callpeaksmac2replicates.txt"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-{wildcards.REP} --outdir {wildcards.path}preprocessing/11peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
rule STEP16_callpeaksmacs2merged:
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
        a="{path}preprocessing/12all/{sample}.all.bam",
        b="{path}preprocessing/12all/{sample}.all.bai"
    output:
        "{path}preprocessing/13allpeaks/{sample}.all_peaks.xls"
    log:
        "{path}preprocessing/logs/{sample}.callpeaksmac2merged.txt"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.all --outdir {wildcards.path}preprocessing/13allpeaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
rule STEP17_callpeaksmacs2merged_localnorm:
    # notes:
    # same as above, but with local background correction enabled
    input:
        a="{path}preprocessing/12all/{sample}.all.bam",
        b="{path}preprocessing/12all/{sample}.all.bai"
    output:
        "{path}preprocessing/13allpeaks/{sample}.localnorm.all_peaks.xls"
    log:
        "{path}preprocessing/logs/{sample}.callpeaksmac2merged.localnorm.txt"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.localnorm.all --outdir {wildcards.path}preprocessing/13allpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"
rule STEP18_plotcorrspearman:
    # parameters:
    # -b input bam files
    # -o output file name
    # -bs set the bin size used for comparison, default is 10000 bp
    # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
    # -p set the number of computing processors to use
    # -v verbose mode
    input:
        a="{path}preprocessing/10unique/{sample}-REP1.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP2.u.bam",
        c="{path}preprocessing/10unique/{sample}-REP3.u.bam",
        d="{path}preprocessing/10unique/{sample}-REP1.u.bai",
        e="{path}preprocessing/10unique/{sample}-REP2.u.bai",
        f="{path}preprocessing/10unique/{sample}-REP3.u.bai"
    output:
        "{path}preprocessing/14qcplots/{sample}.spearman.corrTest"
    log:
        "{path}preprocessing/logs/{sample}.plotcorrspearman.txt"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} {input.c} -o {output} -bs 10000 -p 20 -v"
rule STEP19_makecorrheatmap:
    input:
        "{path}preprocessing/14qcplots/{sample}.spearman.corrTest"
    output:
        "{path}preprocessing/14qcplots/{sample}.spearman.heatmap.svg"
    log:
        "{path}preprocessing/logs/{sample}.makecorrheatmap.txt"
    shell:
        "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"
rule STEP20_makebigwig_bamcov_individual:
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
    input:
        a="{path}preprocessing/10unique/{sample}.u.bam",
        b="{path}preprocessing/10unique/{sample}.u.bai"
    output:
        "{path}preprocessing/16bigwig/{sample}.u.bw"
    log:
        "{path}preprocessing/logs/{sample}.makebigwig_bamcov_individual.txt"
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"
rule STEP21_makebigwig_bamcov_merged:
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
    input:
        a="{path}preprocessing/12all/{mergedsample}.all.bam",
        b="{path}preprocessing/12all/{mergedsample}.all.bai"
    output:
        "{path}preprocessing/16bigwig/{mergedsample}.all.bw"
    log:
        "{path}preprocessing/logs/{mergedsample}.makebigwig_bamcov_merged.txt"
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"
rule STEP22_downsamplebam:
    # params:
    # -Xmx50g set java mem limit to X gb
    input:
        "{path}preprocessing/8merged/{sample}.m.bam"
    output:
        "{path}preprocessing/15downsample/raw/{sample}.{prob}.bam"
    log:
        "{path}preprocessing/logs/{sample}.{prob}.downsamplebam.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar DownsampleSam \
        I={input} \
        O={output} \
        PROBABILITY=0.{wildcards.prob}"
rule STEP23_sortdownsampled:
    # params:
    # -Xmx50g set java mem limit to X gb
    input:
        "{path}preprocessing/15downsample/raw/{sample}.{prob}.bam"
    output:
        "{path}preprocessing/15downsample/raw/{sample}.{prob}.cs.bam"
    log:
        "{path}preprocessing/logs/{sample}.{prob}.sortdownampled.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar SortSam \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate"
rule STEP24_markdupdownsampled:
    # params:
    # -Xmx50g set java mem limit to X gb
    input:
        "{path}preprocessing/15downsample/raw/{sample}.{prob}.cs.bam"
    output:
        a="{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bam",
        b="{path}preprocessing/15downsample/complexity/{sample}.{prob}.dupmetrics.txt",
    log:
        "{path}preprocessing/logs/{sample}.{prob}.markdupdownsampled.txt"
    shell:
        "java -Xmx5g -jar /home/ubuntu2/programs/picard/picard.jar MarkDuplicates \
        I={input} \
        O={output.a} \
        M={output.b} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"
rule STEP25_indexdownsampled:
    input:
        "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bam"
    output:
        "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bai"
    log:
        "{path}preprocessing/logs/{sample}.{prob}.indexdownsampled.txt"
    shell:
        "java -jar /home/ubuntu2/programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"
rule STEP26_callpeaksmacs2downsampled:
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
        "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bam"
    output:
        "{path}preprocessing/15downsample/peaks/{sample}.{prob}_peaks.xls"
    shell:
        "macs2 callpeak -t {input} -n {wildcards.sample}.{wildcards.prob} --outdir preprocessing/15downsample/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
rule AGGREGATE_preprocessing:
    input:
        "{path}preprocessing/10unique/{sample}-REP1.u.bai",
        "{path}preprocessing/10unique/{sample}-REP2.u.bai",
        "{path}preprocessing/10unique/{sample}-REP3.u.bai",
        "{path}preprocessing/12all/{sample}.all.bai",
        "{path}preprocessing/11peaks/{sample}-REP1_peaks.xls",
        "{path}preprocessing/11peaks/{sample}-REP2_peaks.xls",
        "{path}preprocessing/11peaks/{sample}-REP3_peaks.xls",
        "{path}preprocessing/13allpeaks/{sample}.all_peaks.xls",
        "{path}preprocessing/13allpeaks/{sample}.localnorm.all_peaks.xls",
        "{path}preprocessing/14qcplots/{sample}.spearman.heatmap.svg",
        "{path}preprocessing/16bigwig/{sample}-REP1.u.bw",
        "{path}preprocessing/16bigwig/{sample}-REP2.u.bw",
        "{path}preprocessing/16bigwig/{sample}-REP3.u.bw",
        "{path}preprocessing/16bigwig/{sample}.all.bw",
        "{path}preprocessing/15downsample/complexity/{sample}.9.md.bai",
        "{path}preprocessing/15downsample/complexity/{sample}.8.md.bai",
        "{path}preprocessing/15downsample/complexity/{sample}.7.md.bai",
        "{path}preprocessing/15downsample/complexity/{sample}.6.md.bai",
        "{path}preprocessing/15downsample/complexity/{sample}.5.md.bai",
        "{path}preprocessing/15downsample/complexity/{sample}.4.md.bai",
        "{path}preprocessing/15downsample/complexity/{sample}.3.md.bai",
        "{path}preprocessing/15downsample/complexity/{sample}.2.md.bai",
        "{path}preprocessing/15downsample/complexity/{sample}.1.md.bai",
        "{path}preprocessing/15downsample/peaks/{sample}.9_peaks.xls",
        "{path}preprocessing/15downsample/peaks/{sample}.8_peaks.xls",
        "{path}preprocessing/15downsample/peaks/{sample}.7_peaks.xls",
        "{path}preprocessing/15downsample/peaks/{sample}.6_peaks.xls",
        "{path}preprocessing/15downsample/peaks/{sample}.5_peaks.xls",
        "{path}preprocessing/15downsample/peaks/{sample}.4_peaks.xls",
        "{path}preprocessing/15downsample/peaks/{sample}.3_peaks.xls",
        "{path}preprocessing/15downsample/peaks/{sample}.2_peaks.xls",
        "{path}preprocessing/15downsample/peaks/{sample}.1_peaks.xls"
    output:
        "{path}preprocessing/logs/{sample}.preprocessing.done.txt"
    shell:
        "touch {wildcards.path}preprocessing/logs/{wildcards.sample}.preprocessing.done.txt"
rule CLEAN_preprocessing:
    input:
        "{path}preprocessing/logs/{sample}.preprocessing.done.txt"
    output:
        "{path}preprocessing/logs/{sample}.preprocessing.cleaning.done.txt"
    shell:
        """
        rm {wildcards.path}preprocessing/2fastq/*.fastq
        rm {wildcards.path}preprocessing/3goodfastq/*.fq
        rm {wildcards.path}preprocessing/4mycoalign/*.sam
        rm {wildcards.path}preprocessing/5hg38align/*.sam
        rm {wildcards.path}preprocessing/6rawbam/*.bam
        rm {wildcards.path}preprocessing/7rgsort/*.bam
        rm {wildcards.path}preprocessing/8merged/*.bam
        rm {wildcards.path}preprocessing/9dedup/*.bam
        touch {output}
        """

rule test:
    input:
        "ls1034/wt01/saturation/LS1034-WT-01.downsampled_lib_sizes.txt"

########################################################################################################################################
#### Library Complexity Saturation Analysis Rules ######################################################################################
########################################################################################################################################
rule analyzecomplexitysaturation:
    input:
        a="{path}preprocessing/15downsample/complexity/{sample}.9.dupmetrics.txt",
        b="{path}preprocessing/15downsample/complexity/{sample}.8.dupmetrics.txt",
        c="{path}preprocessing/15downsample/complexity/{sample}.7.dupmetrics.txt",
        d="{path}preprocessing/15downsample/complexity/{sample}.6.dupmetrics.txt",
        e="{path}preprocessing/15downsample/complexity/{sample}.5.dupmetrics.txt",
        f="{path}preprocessing/15downsample/complexity/{sample}.4.dupmetrics.txt",
        g="{path}preprocessing/15downsample/complexity/{sample}.3.dupmetrics.txt",
        h="{path}preprocessing/15downsample/complexity/{sample}.2.dupmetrics.txt",
        i="{path}preprocessing/15downsample/complexity/{sample}.1.dupmetrics.txt"
    output:
        "{path}saturation/{sample}.downsampled_lib_sizes.txt"
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
########################################################################################################################################
#### Peaks Saturation Analysis Rules ###################################################################################################
########################################################################################################################################
rule analyzepeaksaturation:
    input:
        "{path}15downsample/peaks/{mergedsample}.9_peaks.xls",
        "{path}15downsample/peaks/{mergedsample}.8_peaks.xls",
        "{path}15downsample/peaks/{mergedsample}.7_peaks.xls",
        "{path}15downsample/peaks/{mergedsample}.6_peaks.xls",
        "{path}15downsample/peaks/{mergedsample}.5_peaks.xls",
        "{path}15downsample/peaks/{mergedsample}.4_peaks.xls",
        "{path}15downsample/peaks/{mergedsample}.3_peaks.xls",
        "{path}15downsample/peaks/{mergedsample}.2_peaks.xls",
        "{path}15downsample/peaks/{mergedsample}.1_peaks.xls"
    output:
        "{path}15downsample/peaks/{mergedsample}.downsampled_numpeaks.txt"
    shell:
        "wl -l < {input} >> {output}"
rule AGGREGATE_saturationanalysis:
    input:
        "{path}logs/{mergedsample}.preprocessing.done.txt",
        "{path}15downsample/complexity/{mergedsample}.downsampled_lib_sizes.txt",
        "{path}15downsample/peaks/{mergedsample}.downsampled_numpeaks.txt"
    output:
        "{path}logs/{mergedsample}.saturation_analysis.done.txt"
    shell:
        "touch {output}"
########################################################################################################################################
#### Footprints Saturation Analysis Rules ##############################################################################################
########################################################################################################################################
rule makefpbychr_downsampled:
    input:
        "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bam",
        "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bai",
        "sites/{gene}.sites.Rdata"
    output:
        "{path}preprocessing/15downsample/footprints/chr/{mergedsample}.{prob}.{gene}.{chr}.done.txt"
    script:
        "scripts/snakeMakeFPbyChrDownsampled.R"
rule mergefpchr_downsampled:
    input:
        "sites/{gene}.sites.Rdata",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr1.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr2.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr3.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr4.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr5.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr6.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr7.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr8.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr9.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr10.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr11.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr12.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr13.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr14.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr15.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr16.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr17.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr18.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr19.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr20.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr21.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chr22.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chrY.done.txt",
        "{path}preprocessing/15downsample/footprints/chr/{sample}.{prob}.{gene}.chrX.done.txt"
    output:
        "{path}preprocessing/15downsample/footprints/merged/{sample}.{prob}.{gene}.done.merged.txt"
    script:
        "scripts/snakeMergeFPbyChrDownsampled.R"
rule allprob_aggregator:
    input:
        "{path}preprocessing/15downsample/footprints/merged/{sample}.9.{gene}.done.merged.txt",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.8.{gene}.done.merged.txt",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.7.{gene}.done.merged.txt",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.6.{gene}.done.merged.txt",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.5.{gene}.done.merged.txt",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.4.{gene}.done.merged.txt",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.3.{gene}.done.merged.txt",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.2.{gene}.done.merged.txt",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.1.{gene}.done.merged.txt"
    output:
        "{path}preprocessing/15downsample/footprints/merged/{sample}.allprob.{gene}.done.txt"
    shell:
        "touch {output}"
rule makefpgraph_downsampled:
    input:
        "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bam",
        "{path}preprocessing/15downsample/complexity/{sample}.{prob}.md.bai",
        "sites/{gene}.sites.Rdata",
        "{path}preprocessing/15downsample/footprints/merged/{sample}.{prob}.{gene}.done.merged.txt"
    output:
        "{path}saturation/footprints/{sample}.{prob}.{gene}.graphs.done.txt"
    script:
        "scripts/snakeGenerateMergedFPGraphDownsampled.R"
rule allgraph_aggregator:
    input:
        "{path}saturation/footprints/{sample}.9.{gene}.graphs.done.txt",
        "{path}saturation/footprints/{sample}.8.{gene}.graphs.done.txt",
        "{path}saturation/footprints/{sample}.7.{gene}.graphs.done.txt",
        "{path}saturation/footprints/{sample}.6.{gene}.graphs.done.txt",
        "{path}saturation/footprints/{sample}.5.{gene}.graphs.done.txt",
        "{path}saturation/footprints/{sample}.4.{gene}.graphs.done.txt",
        "{path}saturation/footprints/{sample}.3.{gene}.graphs.done.txt",
        "{path}saturation/footprints/{sample}.2.{gene}.graphs.done.txt",
        "{path}saturation/footprints/{sample}.1.{gene}.graphs.done.txt"
    output:
        "{path}preprocessing/15downsample/footprints/graphs/{sample}.allprob.{gene}.done.graphs.txt"
    shell:
        "touch {output}"
########################################################################################################################################
#### Footprint Analysis Rules ##########################################################################################################
########################################################################################################################################
rule makefp_by_chr:
    input:
        "{path}preprocessing/12all/{mergedsample}.all.bam",
        "{path}preprocessing/12all/{mergedsample}.all.bai",
        "sites/{gene}.sites.Rdata"
    output:
        "{path}footprints/temp/{mergedsample}.{gene}.{chr}.done.txt"
    script:
        "scripts/snakeMakeFPbyChr.R"
rule merge_chr:
    input:
        "sites/{gene}.sites.Rdata",
        "{path}footprints/temp/{mergedsample}.{gene}.chr1.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr2.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr3.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr4.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr5.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr6.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr7.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr8.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr9.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr10.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr11.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr12.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr13.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr14.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr15.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr16.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr17.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr18.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr19.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr20.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr21.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chr22.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chrX.done.txt",
        "{path}footprints/temp/{mergedsample}.{gene}.chrY.done.txt"
    output:
        "{path}footprints/merged/{mergedsample}.{gene}.merged.done.txt"
    script:
        "scripts/snakeMergeFPbyChr.R"
rule make_graphs:
    input:
        "{path}preprocessing/12all/{mergedsample}.all.bam",
        "{path}preprocessing/12all/{mergedsample}.all.bai",
        "sites/{gene}.sites.Rdata",
        "{path}footprints/merged/{mergedsample}.{gene}.merged.done.txt"
    output:
        "{path}footprints/graphs/{mergedsample}.{gene}.graphs.done.txt"
    script:
        "scripts/snakeGenerateMergedFPGraph.R"
rule parse_footprints:
    input:
        "{path}preprocessing/12all/{mergedsample}.all.bam",
        "{path}preprocessing/12all/{mergedsample}.all.bai",
        "sites/{gene}.sites.Rdata",
        "{path}footprints/merged/{mergedsample}.{gene}.merged.done.txt",
        "{path}preprocessing/13allpeaks/{mergedsample}.all_peaks.narrowPeak"
    output:
        "{path}footprints/parsed/{mergedsample}.{gene}.parsed.done.txt"
    script:
        "scripts/snakeParseFP.R"
rule make_parsed_heatmaps:
    input:
        "{path}footprints/parsed/{mergedsample}.{gene}.motif{motif}.info.Rdata",
    output:
        "{path}footprints/heatmaps/{mergedsample}.{gene}.motif{motif}.heatmap.svg"
    script:
        "scripts/snakeFootprintHeatmaps.R"
rule make_merged_motifs:
    input:
        "{path}parsed/{mergedsample}.{gene}.parsed.done.txt"
    output:
        "{path}merged_motifs/{mergedsample}.{gene}.{nummotif}.mergedmotif.Rdata"
    script:
        "scripts/snakeMergeMotifs.R"
rule make_aracne_overlap:
    input:
        "{path}merged_motifs/{mergedsample}.{gene}.{nummotif}.mergedmotif.Rdata"
    output:
        "{path}aracne/{mergedsample}.{gene}.{nummotif}.{entrez}.aracne.Rdata"
    script:
        "scripts/snakeFindARACNeFootprintOverlap.R"

########################################################################################################################################
#### Analysis Rules ####################################################################################################################
########################################################################################################################################
rule threesample_plotcorrspearman:
    # parameters:
    # -b input bam files
    # -o output file name
    # -bs set the bin size used for comparison, default is 10000 bp
    # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
    # -p set the number of computing processors to use
    # -v verbose mode
    input:
        a="{sample1}/{wt1}{num1}/preprocessing/12all/{s1}.bam",
        b="{sample1}/{wt1}{num1}/preprocessing/12all/{s2}.bam",
        c="{sample1}/{wt1}{num1}/preprocessing/12all/{s3}.bam"
    output:
        "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.corrTest"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} {input.c} -o {output} -bs 10000 -p 20 -v"
rule threesample_makecorrheatmap:
    input:
        "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.corrTest"
    output:
        "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.heatmap.svg"
    shell:
        "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"
rule ninesample_plotcorrspearman:
    # parameters:
    # -b input bam files
    # -o output file name
    # -bs set the bin size used for comparison, default is 10000 bp
    # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
    # -p set the number of computing processors to use
    # -v verbose mode
    # For processing nine sample, restrict analysis to chr1, or computation will take forever
    input:
        a="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
        b="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
        c="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam",
        d="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
        e="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
        f="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam",
        g="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
        h="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
        i="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam"
    output:
        "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.REPS.spearman.corrTest"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} {input.c} {input.d} {input.e} {input.f} {input.g} {input.h} {input.i} -o {output} -bs 10000 -p 20 -v -r chr1"
rule ninesample_makecorrheatmap:
    input:
        "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.REPS.spearman.corrTest"
    output:
        "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.heatmap.svg"
    shell:
        "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"