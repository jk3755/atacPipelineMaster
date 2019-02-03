####################################################################################################################################################################
################################ File Targets ######################################################################################################################
####################################################################################################################################################################

########################
##### SNU61 WT 01 ######
########################
## Run the pipeline all the way up to indexed, deduplicated, uniquely mapped bam files
rule snu61_bai:
        input:
            "snu61/wt01/preprocessing/10unique/SNU61-WT-01-REP1.u.bai",
            "snu61/wt01/preprocessing/10unique/SNU61-WT-01-REP2.u.bai",
            "snu61/wt01/preprocessing/10unique/SNU61-WT-01-REP3.u.bai"
rule snu61_peaks_ind:
        input:
            "snu61/wt01/preprocessing/11peaks/SNU61-WT-01-REP1_peaks.xls",
            "snu61/wt01/preprocessing/11peaks/SNU61-WT-01-REP2_peaks.xls",
            "snu61/wt01/preprocessing/11peaks/SNU61-WT-01-REP3_peaks.xls"
rule snu61_bai_all:
        input:
            "snu61/wt01/preprocessing/12all/SNU61-WT-01.all.bai"
rule snu61_peaks_all:
        input:
            "snu61/wt01/preprocessing/13allpeaks/SNU61-WT-01.all_peaks.xls"
rule snu61_corr_heatmap:
        input:
            "snu61/wt01/preprocessing/14qcplots/SNU61-WT-01.spearman.heatmap.svg"
rule snu61_index_downsampled:
        input:
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.9.md.bai",
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.8.md.bai",
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.7.md.bai",
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.6.md.bai",
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.5.md.bai",
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.4.md.bai",
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.3.md.bai",
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.2.md.bai",
            "snu61/wt01/preprocessing/15downsample/complexity/SNU61-WT-01.1.md.bai"
## The below rule bw_cov can run the whole pipeline from start to finish
rule snu61_bwcov:
        input:
            "snu61/wt01/preprocessing/16bigwig/SNU61-WT-01.all.bw"

### Use these for footprinting saturation analysis
rule snu61_footprint_saturation:
        input:
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.9.CTCF.done.txt",
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.8.CTCF.done.txt",
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.7.CTCF.done.txt",
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.6.CTCF.done.txt",
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.5.CTCF.done.txt",
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.4.CTCF.done.txt",
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.3.CTCF.done.txt",
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.2.CTCF.done.txt",
            "snu61/wt01/preprocessing/15downsample/footprints/parsed/SNU61-WT-01.1.CTCF.done.txt"


## Used to generate ATAC-seq footprints with ATACseqQC package coad_sites
## This method will generate the footprint signals by chromosome and then merge the results
## NOTE it makes more sense to run the first operations, calc by chr and merge separately, as these can handle 20 threads
## To parse, use only 4 threads or it may crash from memory overload
rule snu61_make_fp_bychr_10threads:
        input:
            "snu61/wt01/footprints/graphs/SNU61-WT-01.ASCL2.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.ESRRA.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.TCF7.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.POU5F1B.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.HNF4A.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.OVOL1.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.GMEB2.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.CBFA2T2.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.HOXA3.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.MNX1.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.ZSWIM1.graphs.done.txt",
            "snu61/wt01/footprints/graphs/SNU61-WT-01.CDX2.graphs.done.txt"

## Run me with 4 threads maximum
rule snu61_coad_footprints_parsed_4threads:
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

########################
##### LS1034 WT 01 #####
########################
rule ls1034_bwcov:
        input:
            "ls1034/wt01/preprocessing/16bigwig/LS1034-WT-01.all.bw"

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
            "java -Xmx30g -jar /home/ubuntu1/programs/picard/picard.jar BuildBamIndex \
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
            "{path}12all/{mergedsample}.all.bam"
        output:
            "{path}12all/{mergedsample}.all.bai"
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
            "{path}10unique/{sample}-{REP}.u.bam"
        output:
            "{path}11peaks/{sample}-{REP}_peaks.xls"
        shell:
            "macs2 callpeak -t {input} -n {wildcards.sample}-{wildcards.REP} --outdir {wildcards.path}11peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
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
            b="{path}12all/{sample}.all.bai",
            c="{path}11peaks/{sample}-REP1_peaks.xls",
            d="{path}11peaks/{sample}-REP2_peaks.xls",
            e="{path}11peaks/{sample}-REP3_peaks.xls"
        output:
            "{path}13allpeaks/{sample}.all_peaks.xls"
        shell:
            "macs2 callpeak -t {input.a} -n {wildcards.sample}.all --outdir {wildcards.path}13allpeaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
# STEP 15 - PLOT REPLICATE CORRELATION
rule plot_corr_spearman:
        input:
            a="{path}10unique/{sample}-REP1.u.bam",
            b="{path}10unique/{sample}-REP2.u.bam",
            c="{path}10unique/{sample}-REP3.u.bam",
            d="{path}13allpeaks/{sample}.all_peaks.xls"
        output:
            "{path}14qcplots/{sample}.spearman.corrTest"
        shell:
            "multiBamSummary bins --bamfiles {input.a} {input.b} {input.c} --outFileName {output}"
# STEP 16 - MAKE CORRELATION HEATMAP
rule make_corr_heatmap:
        input:
            "{path}14qcplots/{sample}.spearman.corrTest"
        output:
            "{path}14qcplots/{sample}.spearman.heatmap.svg"
        shell:
            "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"
# STEP 17 - DOWNSAMPLE FOR SATURATION ANALYSIS
rule downsample_bam:
        input:
            "{path}12all/{mergedsample}.all.bam"
        output:
            "{path}15downsample/{mergedsample}.{prob}.bam"
        shell:
            "java -Xmx9g -jar /home/ubuntu1/programs/picard/picard.jar DownsampleSam \
             I={input} \
             O={output} \
             PROBABILITY=0.{wildcards.prob}"
# STEP 18 - COORDINATE SORT DOWNSAMPLED
rule sort_downsampled:
        input:
            "{path}15downsample/{mergedsample}.{prob}.bam"
        output:
            "{path}15downsample/{mergedsample}.{prob}.cs.bam"
        shell:
            "java -Xmx9g -jar /home/ubuntu1/programs/picard/picard.jar SortSam \
             I={input} \
             O={output} \
             SORT_ORDER=coordinate"
# STEP 19 - MARK DUPLICATES DOWNSAMPLED AND LIBRARY COMPLEXITY SATURATION ANALYSIS
rule markdup_downsampled:
        input:
            "{path}15downsample/{mergedsample}.{prob}.cs.bam"
        output:
            "{path}15downsample/complexity/{mergedsample}.{prob}.md.bam"
        shell:
            "java -Xmx5g -jar /home/ubuntu1/programs/picard/picard.jar MarkDuplicates \
             I={input} \
             O={output} \
             M={wildcards.path}15downsample/complexity/{wildcards.mergedsample}.{wildcards.prob}.dupmetrics.txt \
             REMOVE_DUPLICATES=true \
             ASSUME_SORTED=true"
# STEP 20 - INDEX DUPLICATE PURGED DOWNSAMPLES BAM
rule index_downsampled:
        input:
            "{path}15downsample/complexity/{mergedsample}.{prob}.md.bam"
        output:
            "{path}15downsample/complexity/{mergedsample}.{prob}.md.bai"
        shell:
            "java -Xmx5g -jar /home/ubuntu1/programs/picard/picard.jar BuildBamIndex \
            I={input} \
            O={output}"
# STEP 21 - Generate bigwig files with deeptools/bamCoverage/bamCompate for genome browser viewing
# parameters: -b bam input, -o output file, -of output format, -bs bin size in bp, -p number of processors to use, -v verbose mode, --normalizeUsing RPKM (reads per kilobase per million mapped)
rule make_bigwig_bamcov:
        input:
            a="{path}12all/{mergedsample}.all.bam",
            b="{path}12all/{mergedsample}.all.bai",
        output:
            "{path}16bigwig/{mergedsample}.all.bw"
        shell:
            "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v --normalizeUsing RPKM"
# STEP 21 - PEAK SATURATION ANALYSIS
rule peak_saturation_macs2:
        input:
            "{path}14downsample/complexity/{sample}.{prob}.cs.bam"
        output:
            "{path}14downsample/peaks/{sample}.{prob}.peaks.xls"
        shell:
            "macs2 callpeak -t {input} -n {wildcards.sample}.{wildcards.num} --outdir {path}14downsample/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
# STEP 22 - FOOTPRINT SATURATION Analysis
rule saturation_makefp_by_chr:
    input:
        "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bam",
        "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bai",
        "sites/{gene}.sites.Rdata"
    output:
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.{chr}.done.txt",
    script:
        "scripts/snakeMakeFPbyChrDownsampled.R"
#
rule saturation_merge_chr:
    input:
        "sites/{gene}.sites.Rdata",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr1.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr2.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr3.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr4.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr5.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr6.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr7.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr8.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr9.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr10.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr11.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr12.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr13.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr14.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr15.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr16.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr17.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr18.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr19.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr20.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr21.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chr22.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chrY.done.txt",
        "{path}preprocessing/15downsample/footprints/temp/{mergedsample}.{prob}.{gene}.chrX.done.txt"
    output:
        "{path}preprocessing/15downsample/footprints/merged/{mergedsample}.{prob}.{gene}.merged.done.txt"
    script:
        "scripts/snakeMergeFPbyChrDownsampled.R"
#
rule saturation_make_graphs:
    input:
        "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bam",
        "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bai",
        "sites/{gene}.sites.Rdata",
        "{path}preprocessing/15downsample/footprints/merged/{mergedsample}.{prob}.{gene}.merged.done.txt"
    output:
        "{path}preprocessing/15downsample/footprints/graphs/{mergedsample}.{prob}.{gene}.graphs.done.txt"
    script:
        "scripts/snakeGenerateMergedFPGraphDownsample.R"
#
rule saturation_parse_footprints:
    input:
        "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bam",
        "{path}preprocessing/15downsample/complexity/{mergedsample}.{prob}.md.bai",
        "sites/{gene}.sites.Rdata",
        "{path}preprocessing/15downsample/footprints/merged/{mergedsample}.{prob}.{gene}.merged.done.txt",
        "{path}preprocessing/13allpeaks/{mergedsample}.all_peaks.narrowPeak"
    output:
        "{path}preprocessing/15downsample/footprints/parsed/{mergedsample}.{prob}.{gene}.parsed.done.txt"
    script:
        "scripts/snakeParseFPDownsample.R"
#
rule makefp_by_chr_downsampled:
    input:
        "{path}complexity/{sample}.{prob}.cs.bam",
        "{path}complexity/{sample}.{prob}.cs.bai",
        "sites/{gene}.sites.Rdata"
    output:
        "{path}footprints/temp/{sample}.{gene}.{prob}.{chr}.done.txt"
    script:
        "scripts/snakeMakeFPbyChrDownsampled.R"
#
rule merge_chr_downsampled:
    input:
        "sites/{gene}.sites.Rdata",
        "{path}temp/{sample}.{gene}.{prob}.chr1.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr2.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr3.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr4.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr5.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr6.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr7.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr8.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr9.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr10.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr11.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr12.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr13.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr14.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr15.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr16.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr17.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr18.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr19.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr20.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr21.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chr22.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chrX.done.txt",
        "{path}temp/{sample}.{gene}.{prob}.chrY.done.txt"
    output:
        "{path}merged/{sample}.{gene}.{prob}.merged.done.txt"
    script:
        "scripts/snakeMergeFPbyChrDownsampled.R"
#
# STEP 21 - FOOTPRINT SATURATION analysis
# must use a few, include CTCF, CDX2, etc
rule footprint_saturation:
        input:
            "{path}14downsample/complexity/{sample}.09.cs.bam",
            "{path}14downsample/complexity/{sample}.09.cs.bai",
            "{path}14downsample/complexity/{sample}.08.cs.bam",
            "{path}14downsample/complexity/{sample}.08.cs.bai",
            "{path}14downsample/complexity/{sample}.07.cs.bam",
            "{path}14downsample/complexity/{sample}.07.cs.bai",
            "{path}14downsample/complexity/{sample}.06.cs.bam",
            "{path}14downsample/complexity/{sample}.06.cs.bai",
            "{path}14downsample/complexity/{sample}.05.cs.bam",
            "{path}14downsample/complexity/{sample}.05.cs.bai",
            "{path}14downsample/complexity/{sample}.04.cs.bam",
            "{path}14downsample/complexity/{sample}.04.cs.bai",
            "{path}14downsample/complexity/{sample}.03.cs.bam",
            "{path}14downsample/complexity/{sample}.03.cs.bai",
            "{path}14downsample/complexity/{sample}.02.cs.bam",
            "{path}14downsample/complexity/{sample}.02.cs.bai",
            "{path}14downsample/complexity/{sample}.01.cs.bam",
            "{path}14downsample/complexity/{sample}.01.cs.bai"
        output:
            "{path}14downsample/footprints/{sample}.{gene}.done.txt"
        shell:
            "scripts/snakeFootprintSaturation.R"
####################################################################################################################################################################
################################ Footprint Analysis Rules ##########################################################################################################
####################################################################################################################################################################
rule makefp_by_chr:
    input:
        "{path}preprocessing/12all/{mergedsample}.all.bam",
        "{path}preprocessing/12all/{mergedsample}.all.bai",
        "sites/{gene}.sites.Rdata"
    output:
        "{path}footprints/temp/{mergedsample}.{gene}.{chr}.done.txt"
    script:
        "scripts/snakeMakeFPbyChr.R"
#
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
#
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
#
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
#
rule make_parsed_heatmaps:
    input:
        "{path}footprints/parsed/{mergedsample}.{gene}.motif{motif}.info.Rdata",
    output:
        "{path}footprints/heatmaps/{mergedsample}.{gene}.motif{motif}.heatmap.svg"
    script:
        "scripts/snakeFootprintHeatmaps.R"
########################################################################
