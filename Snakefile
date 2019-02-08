####################################################################################################################################################################
################################ H508 WT 01 File Targets ###########################################################################################################
####################################################################################################################################################################
## Raw file info
# H508-1_S3_L001_R1_001.fastq.gz - Sample 1
# H508-2_S2_L001_R1_001.fastq.gz - Sample 2
# H508-3_S1_L001_R1_001.fastq.gz - Sample 3
#
## Before running pipeline, rename these to:
# H508-WT-01_REP1_L1_R1.fastq.gz
# H508-WT-01_REP2_L1_R1.fastq.gz
# H508-WT-01_REP3_L1_R1.fastq.gz

rule h508go:
    input:
        "h508/wt01/preprocessing/logs/H508-WT-01.preprocessing.done.txt"


########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################
rule STEP1_simplifynames_gunzip:
        # params: -k keep original files, -c write to standard output
        input: "{path}1gz/{sample}_L00{lane}_R{read}_001.fastq.gz"
        output: "{path}2fastq/{sample}_L{lane}_R{read}.fastq"
        log: "{path}logs/{sample}.L{lane}.R{read}.simplifynames_gunzip.txt"
        shell: "gunzip -k -c {input} > {output}"
rule STEP2_afterqc_fastqfiltering:
        # params: -s is the shortest trimmed read length allowed past QC filter
        input:
            a="{path}2fastq/{sample}_R1.fastq",
            b="{path}2fastq/{sample}_R2.fastq"
        output:
            c="{path}3goodfastq/{sample}_R1.good.fq",
            d="{path}3goodfastq/{sample}_R2.good.fq"
        log:
            "{path}logs/{sample}.afterqc_fastqfiltering.txt"
        shell:
            "after.py -1 {input.a} -2 {input.b} -g {wildcards.path}3goodfastq -b {wildcards.path}3goodfastq -s 15"
rule STEP3_mycoalign:
        input:
            a="{path}3goodfastq/{sample}_R1.good.fq",
            b="{path}3goodfastq/{sample}_R2.good.fq"
        output:
            "{path}4mycoalign/{sample}.myco.sam"
        log:
            "{path}logs/{sample}.mycoalign.txt"
        shell:
            "bowtie2 -q -p 20 -X1000 -x /home/ubuntu2/genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}4mycoalign/{wildcards.sample}alignment_metrics.txt"
rule STEP4_hg38align:
        # params:
        # -q fastq input
        # -p num threads
        # -X1000 align to a maximum of 1000 bp frag length
        # -1/2 inputs
        # -S output
        input:
            a="{path}3goodfastq/{sample}_R1.good.fq",
            b="{path}3goodfastq/{sample}_R2.good.fq",
            c="{path}4mycoalign/{sample}.myco.sam"
        output:
            "{path}5hg38align/{sample}.hg38.sam"
        log:
            "{path}logs/{sample}.hg38align.txt"
        shell:
            "bowtie2 -q -p 20 -X1000 -x /home/ubuntu2/genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}5hg38align/{wildcards.sample}alignment_metrics.txt"
rule STEP5_samtobam:
        # params:
        # -Xmx50g set java mem limit to X gb
        input:
            "{path}5hg38align/{sample}.hg38.sam"
        output:
            "{path}6rawbam/{sample}.bam"
        log:
            "{path}logs/{sample}.samtobam.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar SamFormatConverter \
            I={input} \
            O={output}"
rule STEP6_addrgandcsbam:
        # note - proper specification of RG tags is critical
        # see: https://software.broadinstitute.org/gatk/documentation/article.php?id=6472
        # Required @RG parameter specifications:
        # RGID (read group ID) - this must be a globally unique string. for illumina data, use flowcell + lane
        # RGLB (read group library) - This is used by MarkDuplicates to collect reads from the same library on different lanes, so it must be common to all files from the same library
        # RGPL (read group platform) - ILLUMINA
        # RGPU (read group platform unit) - The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
        # RGSM (read group sample name) - the name of the sample sequenced in this file. should be consistent across different files from different lanes
        input:
            "{path}6rawbam/{sample}_L{lane}.bam"
        output:
            "{path}7rgsort/{sample}_L{lane}.rg.cs.bam"
        log:
            "{path}logs/{sample}.L{lane}.addrgandcsbam.txt"
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
rule STEP7_cleanbam:
        # params:
        # -Xmx50g set java mem limit to X gb
        input:
            "{path}7rgsort/{sample}_L{lane}.rg.cs.bam"
        output:
            "{path}7rgsort/{sample}_L{lane}.clean.bam"
        log:
            "{path}logs/{sample}.L{lane}.cleanbam.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar CleanSam \
            I={input} \
            O={output}"
rule STEP8_mergelanes:
        input:
            a="{path}7rgsort/{sample}_L1.clean.bam",
            b="{path}7rgsort/{sample}_L2.clean.bam",
            c="{path}7rgsort/{sample}_L3.clean.bam",
            d="{path}7rgsort/{sample}_L4.clean.bam"
        output:
            "{path}8merged/{sample}.m.bam"
        log:
            "{path}logs/{sample}.mergelanes.txt"
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
rule STEP9_purgeduplicates:
        # params:
        # -Xmx50g set java mem limit to X gb
        input:
            "{path}8merged/{sample}.m.bam"
        output:
            a="{path}9dedup/{sample}.dp.bam",
            b="{path}9dedup/{sample}.metrics.txt"
        log:
            "{path}logs/{sample}.purgeduplicates.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar MarkDuplicates \
            I={input} \
            O={output.a} \
            M={output.b} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true"
rule STEP10_mapqfilter:
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
            "{path}9dedup/{sample}.dp.bam"
        output:
            "{path}10unique/{sample}.u.bam"
        log:
            "{path}logs/{sample}.mapqfilter.txt"
        shell:
            "samtools view -h -q 2 -b {input} > {output}"
rule STEP11_buildindex:
        # params:
        # -XmxXg set java mem limit to X gb
        input:
            "{path}10unique/{sample}.u.bam"
        output:
            "{path}10unique/{sample}.u.bai"
        log:
            "{path}logs/{sample}.buildindex.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar BuildBamIndex \
            I={input} \
            O={output}"
rule STEP12_mergereplicates:
        # params:
        # -Xmx50g set java mem limit to X gb
        input:
            a="{path}10unique/{sample}-REP1.u.bam",
            b="{path}10unique/{sample}-REP2.u.bam",
            c="{path}10unique/{sample}-REP3.u.bam"
        output:
            "{path}12all/{sample}.all.bam"
        log:
            "{path}logs/{sample}.mergereplicates.txt"
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
rule STEP13_indexmerged:
        # params:
        # -Xmx50g set java mem limit to X gb
        input:
            "{path}12all/{mergedsample}.all.bam"
        output:
            "{path}12all/{mergedsample}.all.bai"
        log:
            "{path}logs/{mergedsample}.indexmerged.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar BuildBamIndex \
            I={input} \
            O={output}"
rule STEP14_callpeaksmac2replicates:
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
            "{path}10unique/{sample}-{REP}.u.bam"
        output:
            "{path}11peaks/{sample}-{REP}_peaks.xls"
        log:
            "{path}logs/{sample}-{REP}.callpeaksmac2replicates.txt"
        shell:
            "macs2 callpeak -t {input} -n {wildcards.sample}-{wildcards.REP} --outdir {wildcards.path}11peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
rule STEP15_callpeaksmacs2merged:
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
            "{path}12all/{sample}.all.bam"
        output:
            "{path}13allpeaks/{sample}.all_peaks.xls"
        log:
            "{path}logs/{sample}.callpeaksmac2merged.txt"
        shell:
            "macs2 callpeak -t {input.a} -n {wildcards.sample}.all --outdir {wildcards.path}13allpeaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
rule STEP16_plotcorrspearman:
        input:
            a="{path}10unique/{sample}-REP1.u.bam",
            b="{path}10unique/{sample}-REP2.u.bam",
            c="{path}10unique/{sample}-REP3.u.bam"
        output:
            "{path}14qcplots/{sample}.spearman.corrTest"
        log:
            "{path}logs/{sample}.plotcorrspearman.txt"
        shell:
            "multiBamSummary bins --bamfiles {input.a} {input.b} {input.c} --outFileName {output}"
rule STEP17_makecorrheatmap:
        input:
            "{path}14qcplots/{sample}.spearman.corrTest"
        output:
            "{path}14qcplots/{sample}.spearman.heatmap.svg"
        log:
            "{path}logs/{sample}.makecorrheatmap.txt"
        shell:
            "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"

rule STEP18_makebigwig_bamcov_individual:
        # params:
        # -b bam input
        # -o output file
        # -of output format
        # -bs binsize in bp
        # -p number of processors to use
        # -v verbose mode
        # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
        input:
            a="{path}10unique/{sample}.u.bam"
        output:
            "{path}16bigwig/{sample}.u.bw"
        log:
            "{path}logs/{sample}.makebigwig_bamcov_individual.txt"
        shell:
            "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

rule STEP19_makebigwig_bamcov_merged:
        # params:
        # -b bam input
        # -o output file
        # -of output format
        # -bs binsize in bp
        # -p number of processors to use
        # -v verbose mode
        # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
        input:
            "{path}12all/{mergedsample}.all.bam"
        output:
            "{path}16bigwig/{mergedsample}.all.bw"
        log:
            "{path}logs/{mergedsample}.makebigwig_bamcov_merged.txt"
        shell:
            "bamCoverage -b {input} -o {output} -of bigwig -bs 1 -p 20 -v"

rule STEP20_downsamplebam:
        # params:
        # -Xmx50g set java mem limit to X gb
        input:
            "{path}12all/{mergedsample}.all.bam"
        output:
            "{path}15downsample/{mergedsample}.{prob}.bam"
        log:
            "{path}logs/{mergedsample}.{prob}.downsamplebam.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar DownsampleSam \
             I={input} \
             O={output} \
             PROBABILITY=0.{wildcards.prob}"

rule STEP21_sortdownsampled:
        # params:
        # -Xmx50g set java mem limit to X gb
        input:
            "{path}15downsample/{mergedsample}.{prob}.bam"
        output:
            "{path}15downsample/{mergedsample}.{prob}.cs.bam"
        log:
            "{path}logs/{mergedsample}.{prob}.sortdownampled.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar SortSam \
             I={input} \
             O={output} \
             SORT_ORDER=coordinate"

rule STEP22_markdupdownsampled:
        # params:
        # -Xmx50g set java mem limit to X gb
        input:
            "{path}15downsample/{mergedsample}.{prob}.cs.bam"
        output:
            "{path}15downsample/complexity/{mergedsample}.{prob}.md.bam"
        log:
            "{path}logs/{mergedsample}.{prob}.markdupdownsampled.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar MarkDuplicates \
             I={input} \
             O={output} \
             M={wildcards.path}15downsample/complexity/{wildcards.mergedsample}.{wildcards.prob}.dupmetrics.txt \
             REMOVE_DUPLICATES=true \
             ASSUME_SORTED=true"

rule AGGREGATE_preprocessing:
 		input:
 			"{path}10unique/{sample}-REP1.u.bai",
 			"{path}10unique/{sample}-REP2.u.bai",
 			"{path}10unique/{sample}-REP3.u.bai",
 			"{path}12all/{sample}.all.bai",
 			"{path}11peaks/{sample}-REP1_peaks.xls",
 			"{path}11peaks/{sample}-REP2_peaks.xls",
 			"{path}11peaks/{sample}-REP3_peaks.xls",
 			"{path}13allpeaks/{sample}.all_peaks.xls",
 			"{path}14qcplots/{sample}.spearman.heatmap.svg",
 			"{path}16bigwig/{sample}-REP1.u.bw",
 			"{path}16bigwig/{sample}-REP2.u.bw",
 			"{path}16bigwig/{sample}-REP3.u.bw",
 			"{path}16bigwig/{sample}.all.bw"
 		output:
 			"{path}logs/{sample}.preprocessing.done.txt"
 		shell:
 			"touch {path}logs/{sample}.preprocessing.done.txt"



rule STEP23_indexdownsampled:
        input:
            "{path}15downsample/complexity/{mergedsample}.{prob}.md.bam"
        output:
            "{path}15downsample/complexity/{mergedsample}.{prob}.md.bai"
        log:
            "{path}logs/{mergedsample}.{prob}.indexdownsampled.txt"
        shell:
            "java -jar /home/ubuntu2/programs/picard/picard.jar BuildBamIndex \
            I={input} \
            O={output}"

########################################################################################################################################
#### Saturation Analysis Rules #########################################################################################################
########################################################################################################################################

#rule STEP23_analyze_saturation_complexity
#rule STEP24_saturation_libcomplexity

#rule STEP25_saturation_peaks
#input:
#    "{path}14downsample/complexity/{sample}.{prob}.cs.bam"
#output:
#    "{path}14downsample/peaks/{sample}.{prob}.peaks.xls"
#log:
#    "{path}logs/{sample}.peaksaturation.txt"
#shell:
#    "macs2 callpeak -t {input} -n {wildcards.sample}.{wildcards.num} --outdir {path}14downsample/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"


#rule STEP26_saturation_footprints

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
#
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