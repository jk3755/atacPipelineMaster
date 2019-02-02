####################################################################################################################################################################
################################ File Targets ######################################################################################################################
####################################################################################################################################################################

## Run the pipeline all the way up to indexed, deduplicated, uniquely mapped bam files
rule snu61_bai:
        input:
            "snu61/wt01/preprocessing/10unique/SNU61-WT-01-REP1.u.bai",
            "snu61/wt01/preprocessing/10unique/SNU61-WT-01-REP2.u.bai",
            "snu61/wt01/preprocessing/10unique/SNU61-WT-01-REP3.u.bai"

rule snu61_downsample:
        input:
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.09.md.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.08.md.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.07.md.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.06.md.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.05.md.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.04.md.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.03.md.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.02.md.bam",
            "snu61/wt01/preprocessing/13downsample/SNU61-WT-01.01.md.bam"


## Used to generate ATAC-seq footprints with ATACseqQC package coad_sites
## This method will generate the footprint signals by chromosome and then merge the results
## Can be safely run with 10 provided cores on server without exceeding 100 gb mem
rule snu61_coad_footprints:
        input:
            "snu61/wt01/footprints/graphs/SNU61-WT-01.ASCL2.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.TCF7.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.POU5F1B.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.HNF4A.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.OVOL1.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.CBFA2T2.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.HOXA3.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.MNX1.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.ZSWIM1.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.CDX2.graphs.done.txt"

rule snu61_parsed_coad_footprints:
    input:
        "snu61/wt01/footprints/parsed/SNU61-WT-01.ASCL2.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.ESRRA.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.TCF7.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.POU5F1B.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.HNF4A.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.OVOL1.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.GMEB2.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.CBFA2T2.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.HOXA3.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.MNX1.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.ZSWIM1.parsed.done.txt",
        "snu61/wt01/footprints/parsed/SNU61-WT-01.CDX2.parsed.done.txt"


## Generate peak calling FILES
rule peaks:
        input:
            "11peaks/LS1034-WT-01-S1_S5_peaks.xls",
            "11peaks/LS1034-WT-01-S2_S3_peaks.xls",
            "11peaks/LS1034-WT-01-S3_S2_peaks.xls"

## Run pairwise correlation test spearman model
rule corr_spearman:
        input:
            "QC/LS1034-WT-01_spearman_corrTest"

## Make heatmaps
rule corr_heatmap:
        input:
            "QC/LS1034-WT-01_spearman_heatmap.svg"

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


####################################################################################################################################################################
################################ Preprocessing Rules ###############################################################################################################
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
            "after.py -1 {input.a} -2 {input.b} -g {wildcards.path}3goodfastq -b {wildcards.path}3goodfastq -s 15"

# STEP 3 - ALIGN TO MYCO WITH BOWTIE2
rule myco_align:
        input:
            a="{path}3goodfastq/{sample}_R1.good.fq",
            b="{path}3goodfastq/{sample}_R2.good.fq"
        output:
            "{path}4mycoalign/{sample}.myco.sam"
        shell:
            "bowtie2 -q -p 20 -X1000 -x /home/ubuntu1/genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}4mycoalign/{wildcards.sample}alignment_metrics.txt"

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
            b="{path}3goodfastq/{sample}_R2.good.fq",
            c="{path}4mycoalign/{sample}.myco.sam"
        output:
            "{path}5hg38align/{sample}.hg38.sam"
        shell:
            "bowtie2 -q -p 20 -X1000 -x /home/ubuntu1/genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}5hg38align/{wildcards.sample}alignment_metrics.txt"

# STEP 5 - CONVERT SAM TO BAM
## params:
# -Xmx50g set java mem limit to 50 gb
rule sam_to_bam:
        input:
            "{path}5hg38align/{sample}.hg38.sam"
        output:
            "{path}6rawbam/{sample}.bam"
        shell:
            "java -Xmx50g -jar /home/ubuntu1/programs/picard/picard.jar SamFormatConverter \
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
            "java -Xmx50g -jar /home/ubuntu1/programs/picard/picard.jar AddOrReplaceReadGroups \
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
            "java -Xmx50g -jar /home/ubuntu1/programs/picard/picard.jar CleanSam \
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
            "java -Xmx80g -jar /home/ubuntu1/programs/picard/picard.jar MergeSamFiles \
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
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar MarkDuplicates \
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

# STEP 11 - BUILD AN INDEX OF THE FINAL BAM FILES
rule build_index:
        input:
            "{path}10unique/{sample}.u.bam"
        output:
            "{path}10unique/{sample}.u.bai"
        shell:
            "java -Xmx50g -jar /home/ubuntu1/programs/picard/picard.jar BuildBamIndex \
            I={input} \
            O={output}"

# STEP 12 - MERGE REPLICATES
rule merge_replicates:
        input:
            a="{path}10unique/{sample}-REP1.u.bam",
            b="{path}10unique/{sample}-REP2.u.bam",
            c="{path}10unique/{sample}-REP3.u.bam",
            d="{path}10unique/{sample}-REP1.u.bai",
            e="{path}10unique/{sample}-REP2.u.bai",
            f="{path}10unique/{sample}-REP3.u.bai"
        output:
            "{path}12all/{sample}.all.bam"
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

# STEP 13 - INDEX BAM MERGED REPLICATES
rule index_merged:
        input:
            "{path}12all/{sample}.all.bam"
        output:
            "{path}12all/{sample}.all.bai"
        shell:
            "java -Xmx50g -jar /home/ubuntu1/programs/picard/picard.jar BuildBamIndex \
            I={input} \
            O={output}"

# STEP 14 - IND CALL PEAKS WITH MACS2
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
rule peaks_macs2_ind:
        input:
            "{path}10unique/{sample}.u.bam"
        output:
            "{path}11peaks/{sample}.peaks.xls"
        shell:
            "macs2 callpeak -t {input} -n {wildcards.sample} --outdir 11peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# STEP 14 - MERGED CALL PEAKS WITH MACS2
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
rule peaks_macs2_merged:
        input:
            a="{path}12all/{sample}.all.bam",
            b="{path}11peaks/{sample}-REP1.peaks.xls",
            c="{path}11peaks/{sample}-REP2.peaks.xls",
            d="{path}11peaks/{sample}-REP3.peaks.xls"
        output:
            "{path}11peaks/{sample}.all.peaks.xls"
        shell:
            "macs2 callpeak -t {input.a} -n {wildcards.sample} --outdir 11peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# STEP 15 - PLOT REPLICATE CORRELATION
rule plot_corr_spearman:
        input:
            a="{path}10unique/{sample}-REP1.u.bam",
            b="{path}10unique/{sample}-REP2.u.bam",
            c="{path}10unique/{sample}-REP3.u.bam",
            d="{path}11peaks/{sample}.all.peaks.xls"
        output:
            "{path}13qcplots/{sample}.spearman.corrTest"
        shell:
            "multiBamSummary bins --bamfiles {input.a} {input.b} {input.c} --outFileName {output}"

# STEP 16 - MAKE CORRELATION HEATMAP
rule make_corr_heatmap:
        input:
            "{path}13qcplots/{sample}.spearman.corrTest"
        output:
            "{path}13qcplots/{sample}.spearman.heatmap.svg"
        shell:
            "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"

# STEP 17 - DOWNSAMPLE FOR SATURATION ANALYSIS
rule downsample_bam:
        input:
            "{path}8merged/{sample}.m.bam"
        output:
            "{path}14downsample/{sample}.{prob}.bam"
        shell:
            "java -Xmx5g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input} \
             O={output.a} \
             PROBABILITY={wildcards.prob}"

# STEP 18 - COORDINATE SORT DOWNSAMPLED
rule sort_downsampled:
        input:
            "{path}14downsample/{sample}.{prob}.bam"
        output:
            "{path}14downsample/{sample}.{prob}.cs.bam"
        shell:
            "java -Xmx5g -jar /home/ubuntu1/programs/picard/picard.jar SortSam \
             I={input} \
             O={output} \
             SORT_ORDER=coordinate"

# STEP 19 - MARK DUPLICATES DOWNSAMPLED AND LIBRARY COMPLEXITY SATURATION ANALYSIS
rule markdup_downsampled:
        input:
            "{path}14downsample/{sample}.{prob}.cs.bam"
        output:
            "{path}14downsample/complexity/{sample}.{prob}.md.bam"
        shell:
            "java -Xmx5g -jar /home/ubuntu1/programs/picard/picard.jar MarkDuplicates \
             I={input} \
             O={output} \
             M={wildcards.sample}.{wildcards.prob}.dupmetrics.txt"

# STEP 20 - PEAK SATURATION ANALYSIS
rule peak_saturation_macs2:
        input:
            "{path}14downsample/complexity/{sample}.{prob}.cs.bam"
        output:
            "{path}14downsample/peaks/{sample}.{prob}.peaks.xls"
        shell:
            "macs2 callpeak -t {input} -n {wildcards.sample}.{wildcards.num} --outdir {path}14downsample/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# STEP 21 - FOOTPRINT SATURATION analysis
# must use a few, include CTCF, CDX2, etc
rule footprint_saturation:
        input:
            "{path}14downsample/complexity/{sample}.{prob}.cs.bam"
        output:
            "{path}14downsample/footprints/{sample}.{prob}.peaks.xls"
        shell:
            "scripts/snakeFootprintSaturation.R"

####################################################################################################################################################################
################################ Footprint Analysis Rules ##########################################################################################################
####################################################################################################################################################################
rule makefp_by_chr:
    input:
        "{path}bam/{sample}.all.bam",
        "{path}bam/{sample}.all.bai",
        "sites/{gene}.sites.Rdata"
    output:
        "{path}footprints/{sample}.{gene}.{chr}.done.txt"
    script:
        "scripts/snakeMakeFPbyChr.R"

rule merge_chr:
    input:
        "sites/{gene}.sites.Rdata",
        "{path}footprints/temp/{sample}.{gene}.chr1.done.txt",
        "{path}footprints/temp/{sample}.{gene}.chr2.done.txt",
        "{path}footprints/temp/{sample}.{gene}.chr3.done.txt",
        "{path}footprints/temp/{sample}.{gene}.chr4.done.txt",
        "{path}footprints/temp/{sample}.{gene}.chr5.done.txt",
        "{path}footprints/temp/{sample}.{gene}.chr6.done.txt",
        "{path}footprints/temp/{sample}.{gene}.chr7.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr8.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr9.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr10.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr11.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr12.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr13.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr14.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr15.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr16.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr17.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr18.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr19.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr20.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr21.done.txt",
        "{path}footprints/temp{sample}.{gene}.chr22.done.txt",
        "{path}footprints/temp{sample}.{gene}.chrX.done.txt",
        "{path}footprints/temp{sample}.{gene}.chrY.done.txt"
    output:
        "{path}footprints/merged/{sample}.{gene}.merged.done.txt"
    script:
        "scripts/snakeMergeFPbyChr.R"

rule make_graphs:
    input:
        "{path}preprocessing/12all/{sample}.all.bam",
        "{path}preprocessing/12all/{sample}.all.bai",
        "sites/{gene}.sites.Rdata",
        "{path}footprints/merged/{sample}.{gene}.merged.done.txt"
    output:
        "{path}footprints/graphs/{sample}.{gene}.graphs.done.txt"
    script:
        "scripts/snakeGenerateMergedFPGraph.R"

rule parse_footprints:
    input:
        "ls1034/wt01/bam/{sample}.all.bam",
        "ls1034/wt01/bam/{sample}.all.bai",
        "sites/{gene}.sites.Rdata",
        "ls1034/wt01/merged/{sample}.{gene}.merged.done.txt",
        "ls1034/wt01/peaks/{sample}_peaks.narrowPeak"
    output:
        "ls1034/wt01/parsed/{sample}.{gene}.parsed.done.txt"
    script:
        "ls1034/wt01/scripts/snakeParseFP.R"
