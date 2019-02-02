####################################################################################################################################################################
################################ File Targets ######################################################################################################################
####################################################################################################################################################################
rule snu61_init:
        input:
            "snu61/wt01/preprocessing/2fastq/SNU61-WT-01-S3_S1_L1_R1.fastq",
            "snu61/wt01/preprocessing/2fastq/SNU61-WT-01-S3_S1_L1_R2.fastq",
            "snu61/wt01/preprocessing/2fastq/SNU61-WT-01-S3_S1_L2_R1.fastq",
            "snu61/wt01/preprocessing/2fastq/SNU61-WT-01-S3_S1_L2_R2.fastq",
            "snu61/wt01/preprocessing/2fastq/SNU61-WT-01-S3_S1_L3_R1.fastq",
            "snu61/wt01/preprocessing/2fastq/SNU61-WT-01-S3_S1_L3_R2.fastq",
            "snu61/wt01/preprocessing/2fastq/SNU61-WT-01-S3_S1_L4_R1.fastq",
            "snu61/wt01/preprocessing/2fastq/SNU61-WT-01-S3_S1_L4_R2.fastq"

## Run the pipeline all the way up to indexed, deduplicated, uniquely mapped bam files
rule snu61_bai:
        input:
            "snu61/wt01/preprocessing/10unique/SNU61-WT-01-S3_S1.u.bai"

##
rule snu61_downsample:
        input:
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.09.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.08.bam"



## Generate peak calling FILES
rule peaks:
        input:
            "11peaks/LS1034-WT-01-S1_S5_peaks.xls",
            "11peaks/LS1034-WT-01-S2_S3_peaks.xls",
            "11peaks/LS1034-WT-01-S3_S2_peaks.xls"
## Run the mycoplasma alignment test
rule myco:
        input:
            "4myco/LS1034-WT-01-S1_S5_L1.myco.sam",
            "4myco/LS1034-WT-01-S1_S5_L2.myco.sam",
            "4myco/LS1034-WT-01-S1_S5_L3.myco.sam",
            "4myco/LS1034-WT-01-S1_S5_L4.myco.sam",
            #
            "4myco/LS1034-WT-01-S2_S3_L1.myco.sam",
            "4myco/LS1034-WT-01-S2_S3_L2.myco.sam",
            "4myco/LS1034-WT-01-S2_S3_L3.myco.sam",
            "4myco/LS1034-WT-01-S2_S3_L4.myco.sam",
            #
            "4myco/LS1034-WT-01-S3_S2_L1.myco.sam",
            "4myco/LS1034-WT-01-S3_S2_L2.myco.sam",
            "4myco/LS1034-WT-01-S3_S2_L3.myco.sam",
            "4myco/LS1034-WT-01-S3_S2_L4.myco.sam"
## Run pairwise correlation test spearman model
rule corr_spearman:
        input:
            "QC/LS1034-WT-01_spearman_corrTest"
## Make heatmaps
rule corr_heatmap:
        input:
            "QC/LS1034-WT-01_spearman_heatmap.svg"
## Perform sequencing saturation analysis
rule sat_frags:
        input:
            "sat/LS1034-WT-01-S1_S5.9.bam",
            "sat/LS1034-WT-01-S2_S3.9.bam",
            "sat/LS1034-WT-01-S3_S2.9.bam",
            "sat/LS1034-WT-01-S1_S5.8.bam",
            "sat/LS1034-WT-01-S2_S3.8.bam",
            "sat/LS1034-WT-01-S3_S2.8.bam",
            "sat/LS1034-WT-01-S1_S5.7.bam",
            "sat/LS1034-WT-01-S2_S3.7.bam",
            "sat/LS1034-WT-01-S3_S2.7.bam",
            "sat/LS1034-WT-01-S1_S5.6.bam",
            "sat/LS1034-WT-01-S2_S3.6.bam",
            "sat/LS1034-WT-01-S3_S2.6.bam",
            "sat/LS1034-WT-01-S1_S5.5.bam",
            "sat/LS1034-WT-01-S2_S3.5.bam",
            "sat/LS1034-WT-01-S3_S2.5.bam",
            "sat/LS1034-WT-01-S1_S5.4.bam",
            "sat/LS1034-WT-01-S2_S3.4.bam",
            "sat/LS1034-WT-01-S3_S2.4.bam",
            "sat/LS1034-WT-01-S1_S5.3.bam",
            "sat/LS1034-WT-01-S2_S3.3.bam",
            "sat/LS1034-WT-01-S3_S2.3.bam",
            "sat/LS1034-WT-01-S1_S5.2.bam",
            "sat/LS1034-WT-01-S2_S3.2.bam",
            "sat/LS1034-WT-01-S3_S2.2.bam",
            "sat/LS1034-WT-01-S1_S5.1.bam",
            "sat/LS1034-WT-01-S2_S3.1.bam",
            "sat/LS1034-WT-01-S3_S2.1.bam"
##
rule sat_frags_sort:
        input:
            "sat/LS1034-WT-01-S1_S5.9.cs.bam",
            "sat/LS1034-WT-01-S2_S3.9.cs.bam",
            "sat/LS1034-WT-01-S3_S2.9.cs.bam",
            "sat/LS1034-WT-01-S1_S5.8.cs.bam",
            "sat/LS1034-WT-01-S2_S3.8.cs.bam",
            "sat/LS1034-WT-01-S3_S2.8.cs.bam",
            "sat/LS1034-WT-01-S1_S5.7.cs.bam",
            "sat/LS1034-WT-01-S2_S3.7.cs.bam",
            "sat/LS1034-WT-01-S3_S2.7.cs.bam",
            "sat/LS1034-WT-01-S1_S5.6.cs.bam",
            "sat/LS1034-WT-01-S2_S3.6.cs.bam",
            "sat/LS1034-WT-01-S3_S2.6.cs.bam",
            "sat/LS1034-WT-01-S1_S5.5.cs.bam",
            "sat/LS1034-WT-01-S2_S3.5.cs.bam",
            "sat/LS1034-WT-01-S3_S2.5.cs.bam",
            "sat/LS1034-WT-01-S1_S5.4.cs.bam",
            "sat/LS1034-WT-01-S2_S3.4.cs.bam",
            "sat/LS1034-WT-01-S3_S2.4.cs.bam",
            "sat/LS1034-WT-01-S1_S5.3.cs.bam",
            "sat/LS1034-WT-01-S2_S3.3.cs.bam",
            "sat/LS1034-WT-01-S3_S2.3.cs.bam",
            "sat/LS1034-WT-01-S1_S5.2.cs.bam",
            "sat/LS1034-WT-01-S2_S3.2.cs.bam",
            "sat/LS1034-WT-01-S3_S2.2.cs.bam",
            "sat/LS1034-WT-01-S1_S5.1.cs.bam",
            "sat/LS1034-WT-01-S2_S3.1.cs.bam",
            "sat/LS1034-WT-01-S3_S2.1.cs.bam"
##
rule sat_frags_duplicate:
        input:
            "sat/LS1034-WT-01-S1_S5.9.md.bam",
            "sat/LS1034-WT-01-S2_S3.9.md.bam",
            "sat/LS1034-WT-01-S3_S2.9.md.bam",
            "sat/LS1034-WT-01-S1_S5.8.md.bam",
            "sat/LS1034-WT-01-S2_S3.8.md.bam",
            "sat/LS1034-WT-01-S3_S2.8.md.bam",
            "sat/LS1034-WT-01-S1_S5.7.md.bam",
            "sat/LS1034-WT-01-S2_S3.7.md.bam",
            "sat/LS1034-WT-01-S3_S2.7.md.bam",
            "sat/LS1034-WT-01-S1_S5.6.md.bam",
            "sat/LS1034-WT-01-S2_S3.6.md.bam",
            "sat/LS1034-WT-01-S3_S2.6.md.bam",
            "sat/LS1034-WT-01-S1_S5.5.md.bam",
            "sat/LS1034-WT-01-S2_S3.5.md.bam",
            "sat/LS1034-WT-01-S3_S2.5.md.bam",
            "sat/LS1034-WT-01-S1_S5.4.md.bam",
            "sat/LS1034-WT-01-S2_S3.4.md.bam",
            "sat/LS1034-WT-01-S3_S2.4.md.bam",
            "sat/LS1034-WT-01-S1_S5.3.md.bam",
            "sat/LS1034-WT-01-S2_S3.3.md.bam",
            "sat/LS1034-WT-01-S3_S2.3.md.bam",
            "sat/LS1034-WT-01-S1_S5.2.md.bam",
            "sat/LS1034-WT-01-S2_S3.2.md.bam",
            "sat/LS1034-WT-01-S3_S2.2.md.bam",
            "sat/LS1034-WT-01-S1_S5.1.md.bam",
            "sat/LS1034-WT-01-S2_S3.1.md.bam",
            "sat/LS1034-WT-01-S3_S2.1.md.bam"
##
rule sat_peaks:
    input:
        "peaksat/LS1034-WT-01-S1_S5.9_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.9_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.9_peaks.xls",
        "peaksat/LS1034-WT-01-S1_S5.8_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.8_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.8_peaks.xls",
        "peaksat/LS1034-WT-01-S1_S5.7_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.7_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.7_peaks.xls",
        "peaksat/LS1034-WT-01-S1_S5.6_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.6_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.6_peaks.xls",
        "peaksat/LS1034-WT-01-S1_S5.5_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.5_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.5_peaks.xls",
        "peaksat/LS1034-WT-01-S1_S5.4_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.4_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.4_peaks.xls",
        "peaksat/LS1034-WT-01-S1_S5.3_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.3_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.3_peaks.xls",
        "peaksat/LS1034-WT-01-S1_S5.2_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.2_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.2_peaks.xls",
        "peaksat/LS1034-WT-01-S1_S5.1_peaks.xls",
        "peaksat/LS1034-WT-01-S2_S3.1_peaks.xls",
        "peaksat/LS1034-WT-01-S3_S2.1_peaks.xls"

####################################################################################################################################################################
################################ Processing Rules ##################################################################################################################
####################################################################################################################################################################

# STEP 1 - GUNZIP FASTQ FILES
# params: -k keep original files, -c write to standard output
rule gunzip_namechange:
        input:
            "{path}1gz/{sample}_L00{lane}_R{read}_001.fastq.gz"
        output:
            "{path}2fastq/{sample}_L{lane}_R{read}.fastq"
        shell:
            "gunzip -k -c {input} > {output}"

# STEP 2 - FASTQ QUALITY FILTERING WITH AFTERQC
# param -s is the shortest trimmed read length allowed past QC filter
rule afterqc_qc:
        input:
            a="{path}2fastq/{sample}_R1.fastq",
            b="{path}2fastq/{sample}_R2.fastq"
        output:
            c="{path}3goodfastq/{sample}_R1.good.fq",
            d="{path}3goodfastq/{sample}_R2.good.fq"
        shell:
            "after.py -1 {input.a} -2 {input.b} -g 3goodfastq -b 3goodfastq -s 15"

# STEP 3 - ALIGN TO MYCO WITH BOWTIE2
rule myco_align:
        input:
            a="{path}3goodfastq/{sample}_R1.good.fq",
            b="{path}3goodfastq/{sample}_R2.good.fq"
        output:
            "{path}4mycoalign/{sample}.myco.sam"
        shell:
            "bowtie2 -q -p 20 -X1000 -x /home/ubuntu2/genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} > {path}4mycoalign/{sample}alignment_metrics.txt"

# STEP 4 - ALIGN TO HG38 WITH BOWTIE2
## params:
# -q fastq input
# -p num threads
# -X1000 align to a maximum of 1000 bp frag length
# -1/2 inputs
# -S output
rule hg38_align:
        input:
            a="{path}3goodfastq/{sample}_R1.good.fq",
            b="{path}3goodfastq/{sample}_R2.good.fq"
        output:
            "{path}5hg38align/{sample}.hg38.sam"
        shell:
            "bowtie2 -q -p 20 -X1000 -x /home/ubuntu2/genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} > {path}5hg38align/{sample}alignment_metrics.txt"

# STEP 5 - CONVERT SAM TO BAM
## params:
# -Xmx50g set java mem limit to 50 gb
rule sam_to_bam:
        input:
            "{path}5hg38align/{sample}.hg38.sam"
        output:
            "{path}6rawbam/{sample}.bam"
        shell:
            "java -Xmx50g -jar /home/ubuntu2/programs/picard/picard.jar SamFormatConverter \
            I={input} \
            O={output}"

# STEP 6
# note - proper specification of RG tags is critical
# see: https://software.broadinstitute.org/gatk/documentation/article.php?id=6472
## Required @RG parameter specifications:
# RGID (read group ID) - this must be a globally unique string. for illumina data, use flowcell + lane
# RGLB (read group library) - This is used by MarkDuplicates to collect reads from the same library on different lanes, so it must be common to all files from the same library
# RGPL (read group platform) - ILLUMINA
# RGPU (read group platform unit) - The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
# RGSM (read group sample name) - the name of the sample sequenced in this file. should be consistent across different files from different lanes
rule add_rg_and_cs_bam:
        input:
            "{path}6rawbam/{sample}_L{lane}.bam"
        output:
            "{path}7rgsort/{sample}_L{lane}.rg.cs.bam"
        shell:
            "java -Xmx50g -jar /home/ubuntu2/programs/picard/picard.jar AddOrReplaceReadGroups \
            I={input} \
            O={output} \
            SORT_ORDER=coordinate \
            RGID=H5YHHBGX3.{wildcards.lane} \
            RGLB={wildcards.sample} \
            RGPL=ILLUMINA \
            RGPU=H5YHHBGX3.{wildcards.lane}.{wildcards.sample} \
            RGSM={wildcards.sample}"

# STEP 7 - CLEAN BAM FILES
rule clean_bam:
        input:
            "{path}7rgsort/{sample}_L{lane}.rg.cs.bam"
        output:
            "{path}7rgsort/{sample}_L{lane}.clean.bam"
        shell:
            "java -Xmx50g -jar /home/ubuntu2/programs/picard/picard.jar CleanSam \
            I={input} \
            O={output}"

# STEP 8 - MERGE LANES
rule merge_lanes:
        input:
            a="{path}7rgsort/{sample}_L1.clean.bam",
            b="{path}7rgsort/{sample}_L2.clean.bam",
            c="{path}7rgsort/{sample}_L3.clean.bam",
            d="{path}7rgsort/{sample}_L4.clean.bam"
        output:
            "{path}8merged/{sample}.m.bam"
        shell:
            "java -Xmx80g -jar /home/ubuntu2/programs/picard/picard.jar MergeSamFiles \
            I={input.a} \
            I={input.b} \
            I={input.c} \
            I={input.d} \
            O={output} \
            SORT_ORDER=coordinate \
            ASSUME_SORTED=true \
            MERGE_SEQUENCE_DICTIONARIES=true \
            USE_THREADING=true"

# STEP 9 - PURGE PCR DUPLICATES
rule purge_duplicates:
        input:
            "{path}8merged/{sample}.m.bam"
        output:
            a="{path}9dedup/{sample}.dp.bam",
            b="{path}9dedup/{sample}.metrics.txt"
        shell:
            "java -Xmx50g -jar /home/ubuntu2/programs/picard/picard.jar MarkDuplicates \
            I={input} \
            O={output.a} \
            M={output.b} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true"

# STEP 10 - REMOVE MULTI MAPPING READS WITH SAMTOOLS
## Notes:
# for an explanation of how bowtie2 calculates mapq scores:
# http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
# for bowtie2, mapq higher than 2 is a uniquely mapped read
## params:
# -h include the header in the output
# -q only include reads with mapping quality X or higher
# -b output as a bam file
rule mapq_filter:
        input:
            "{path}9dedup/{sample}.dp.bam"
        output:
            a="{path}10unique/{sample}.u.bam",
        shell:
            "samtools view -h -q 2 -b {input} > {output}"

# STEP 11 - BUILD AN INDEX OF THE FINAL BAM FILE
rule build_index:
        input:
            "{path}10unique/{sample}.u.bam"
        output:
            "{path}10unique/{sample}.u.bai"
        shell:
            "java -Xmx50g -jar /home/ubuntu2/programs/picard/picard.jar BuildBamIndex \
            I={input} \
            O={output}"

# STEP 12 - CALL PEAKS WITH MACS2
## notes:
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
rule peaks_macs2:
        input:
            "{path}10unique/{sample}.u.bam"
        output:
            "{path}11peaks/{sample}.peaks.xls"
        shell:
            "macs2 callpeak -t {input} -n {wildcards.sample} --outdir 11peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# STEP 13 - PLOT REPLICATE CORRELATION
rule plot_corr_spearman:
        input:
            a="{path}10unique/{sample}-{rep1}.u.bam",
            b="{path}10unique/{sample}-{rep2}.u.bam",
            c="{path}10unique/{sample}-{rep3}.u.bam",
        output:
            "{path}12qcplots/{sample}.spearman.corrTest"
        shell:
            "multiBamSummary bins --bamfiles {input.a} {input.b} {input.c} --outFileName {output}"

# STEP 14 - MAKE CORRELATION HEATMAP
rule make_corr_heatmap:
        input:
            "{path}12qcplots/{sample}.spearman.corrTest"
        output:
            "{path}12qcplots/{sample}.spearman.heatmap.svg"
        shell:
            "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"

# STEP 15 - DOWNSAMPLE FOR SATURATION ANALYSIS
rule downsample_bam:
        input:
            "{path}8merged/{sample}.m.bam"
        output:
            "{path}13downsample/{sample}.{prob}.bam"
        shell:
            "java -Xmx5g -jar /home/ubuntu2/programs/picard/picard.jar DownsampleSam \
             I={input} \
             O={output.a} \
             PROBABILITY={wildcards.prob}"
##
rule downsample_frags8:
        input:
            a="sat/{sample}.9.bam",
            b="8merged/{sample}.m.bam"
        output:
            "sat/{sample}.8.bam"
        shell:
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input.b} \
             O={output} \
             PROBABILITY=0.8"
##
rule downsample_frags7:
        input:
            a="sat/{sample}.8.bam",
            b="8merged/{sample}.m.bam"
        output:
            "sat/{sample}.7.bam"
        shell:
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input.b} \
             O={output} \
             PROBABILITY=0.7"
##
rule downsample_frags6:
        input:
            a="sat/{sample}.7.bam",
            b="8merged/{sample}.m.bam"
        output:
            "sat/{sample}.6.bam"
        shell:
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input.b} \
             O={output} \
             PROBABILITY=0.6"
##
rule downsample_frags5:
        input:
            a="sat/{sample}.6.bam",
            b="8merged/{sample}.m.bam"
        output:
            "sat/{sample}.5.bam"
        shell:
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input.b} \
             O={output} \
             PROBABILITY=0.5"
##
rule downsample_frags4:
        input:
            a="sat/{sample}.5.bam",
            b="8merged/{sample}.m.bam"
        output:
            "sat/{sample}.4.bam"
        shell:
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input.b} \
             O={output} \
             PROBABILITY=0.4"
##
rule downsample_frags3:
        input:
            a="sat/{sample}.4.bam",
            b="8merged/{sample}.m.bam"
        output:
            "sat/{sample}.3.bam"
        shell:
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input.b} \
             O={output} \
             PROBABILITY=0.3"
##
rule downsample_frags2:
        input:
            a="sat/{sample}.3.bam",
            b="8merged/{sample}.m.bam"
        output:
            "sat/{sample}.2.bam"
        shell:
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input.b} \
             O={output} \
             PROBABILITY=0.2"
##
rule downsample_frags1:
        input:
            a="sat/{sample}.2.bam",
            b="8merged/{sample}.m.bam"
        output:
            "sat/{sample}.1.bam"
        shell:
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input.b} \
             O={output} \
             PROBABILITY=0.1"
#######################################################################################################

#######################################################################################################
#### STEP 16 - CORD SORT DOWNSAMPLED ####
#########################################
## notes:
rule sort_downsampled:
        input:
            "sat/{sample}.{num}.bam"
        output:
            "sat/{sample}.{num}.cs.bam"
        shell:
            "java -Xmx5g -jar /home/ubuntu1/programs/picard/picard.jar SortSam \
             I={input} \
             O={output} \
             SORT_ORDER=coordinate"
##
#######################################################################################################
#### STEP 17 - MARK DUPLICATES DOWNSAMPLED ####
###############################################
## notes:
rule markdup_downsampled:
        input:
            "sat/{sample}.{num}.cs.bam"
        output:
            "sat/{sample}.{num}.md.bam"
        shell:
            "java -Xmx5g -jar /home/ubuntu1/programs/picard/picard.jar MarkDuplicates \
             I={input} \
             O={output} \
             M={wildcards.sample}.{wildcards.num}.dupmetrics.txt"
##
#######################################################################################################
#### STEP 17 - PEAK SATURATION ANALYSIS ####
############################################
rule peak_sat_macs2:
        input:
            "sat/{sample}.{num}.cs.bam"
        output:
            "peaksat/{sample}.{num}_peaks.xls"
        shell:
            "macs2 callpeak -t {input} -n {wildcards.sample}.{wildcards.num} --outdir peaksat --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

#######################################################################################################
#### STEP 17 - MERGE REPLICATES ####
####################################
rule merge_replicates:
        input:
            a="10unique/{sample}-S1_S5.u.bam",
            b="10unique/{sample}-S2_S3.u.bam",
            c="10unique/{sample}-S3_S2.u.bam"
        output:
            "12all/{sample}.all.bam"
        shell:
            "java -Xmx80g -jar /home/ubuntu1/programs/picard/picard.jar MergeSamFiles \
            I={input.a} \
            I={input.b} \
            I={input.c} \
            O={output} \
            SORT_ORDER=coordinate \
            ASSUME_SORTED=true \
            MERGE_SEQUENCE_DICTIONARIES=true \
            USE_THREADING=true"

#######################################################################################################
#### STEP 18 - INDEX MERGED ####
################################
rule index_merged:
        input:
            "12all/{sample}.all.bam"
        output:
            "12all/{sample}.all.bai"
        shell:
            "java -Xmx50g -jar /home/ubuntu1/programs/picard/picard.jar BuildBamIndex \
            I={input} \
            O={output}"

#######################################################################################################
#### STEP 17 - FOOTPRINTS ####
##############################

rule generate_coad_mr_footprints:
    input:
        "12all/{sample}.all.bam",
        "12all/{sample}.all.bai",
        "/home/ubuntu1/atac/coad_sites/{gene}.sites.Rdata"
    output:
        "coad_footprints/{sample}.{gene}.done.txt"
    script:
        "scripts/snake_FP_COADMR_WG.R"
#######################################################################################################

## Generate footprints for different TFs
rule group1:
        input:
            "fp/SNU61-WT-01-S3_S1.TFAP2A.Rdata",
            "fp/SNU61-WT-01-S3_S1.NFIL3.Rdata",
            "fp/SNU61-WT-01-S3_S1.HLF.Rdata",
            "fp/SNU61-WT-01-S3_S1.NHLH1.Rdata",
            "fp/SNU61-WT-01-S3_S1.MAX.Rdata",
            "fp/SNU61-WT-01-S3_S1.USF1.Rdata",
            "fp/SNU61-WT-01-S3_S1.CEBPA.Rdata",
            "fp/SNU61-WT-01-S3_S1.EBF1.Rdata",
            "fp/SNU61-WT-01-S3_S1.CEBPB.Rdata",
            "fp/SNU61-WT-01-S3_S1.FOS.Rdata",
            "fp/SNU61-WT-01-S3_S1.FOSL1.Rdata",
            "fp/SNU61-WT-01-S3_S1.FOSL2.Rdata",
            "fp/SNU61-WT-01-S3_S1.JUN.Rdata",
            "fp/SNU61-WT-01-S3_S1.JUNB.Rdata",
            "fp/SNU61-WT-01-S3_S1.JUND.Rdata",
            "fp/SNU61-WT-01-S3_S1.MAFF.Rdata",
            "fp/SNU61-WT-01-S3_S1.MAFK.Rdata",
            "fp/SNU61-WT-01-S3_S1.TFAP2C.Rdata",
            "fp/SNU61-WT-01-S3_S1.USF2.Rdata",
            "fp/SNU61-WT-01-S3_S1.SREBF1.Rdata"

####################################################################################

###################################################################################################################
#### STEP 1 - Get FOOTPRINT with R ####
#####################################
rule make_fp:
        input:
            "bam/{sample}.u.bam",
            "bam/{sample}.u.bai"
        output:
            "fp/{sample}.{gene}.Rdata"
        script:
            "scripts/snake_FP.R"
##################################################################################################################
