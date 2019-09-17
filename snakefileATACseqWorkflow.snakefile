########################################################################################################################################
#### IMPORT MODULES AND CONFIG #########################################################################################################
########################################################################################################################################
configfile: "snakeResources/config/config.yaml"
include: "snakeResources/modules/spoolPreprocessing.snakefile"
include: "snakeResources/modules/spoolFootprinting.snakefile"
include: "snakeResources/modules/spoolFullAnalysis.snakefile"
include: "snakeResources/modules/spoolSampleCorrelation.snakefile"

########################################################################################################################################
#### AGGREGATOR RULES ##################################################################################################################
########################################################################################################################################

## This rule determines what is run in the full analysis spooling option
rule AGGREGATOR_full_analysis:
    input:
        "{path}operations/modules/{sample}.rep{repnum}.ref{refgenome}.preprocessing.complete",
        "{path}operations/modules/{sample}.rep{repnum}.ref{refgenome}.generate_figures.complete",
        "{path}operations/modules/{sample}.rep{repnum}.ref{refgenome}.footprinting_analysis.complete"
    output:
        "{path}operations/modules/{sample}.rep{repnum}.ref{refgenome}.full_analysis.finished"
    shell:
        "touch {output}"

########################################################################################################################################
#### PREPROCESSING AGGREGATORS #########################################################################################################
########################################################################################################################################

## Clean up all the intermediate files just before touching the final flag file
rule AGGREGATOR_preprocessing_clean_intermediate_data:
    input:
        "{path}operations/preprocessing/{sample}.rep{repnum}.ref{refgenome}.preprocessing_preclean.complete"  
    output:
        "{path}operations/modules/{sample}.rep{repnum}.ref{refgenome}.preprocessing.complete"
    shell:
        """
        rm -rf {wildcards.path}preprocessing/2fastq
        rm -rf {wildcards.path}preprocessing/3goodfastq
        rm -rf {wildcards.path}preprocessing/4mycoalign
        rm -rf {wildcards.path}preprocessing/5align
        rm -rf {wildcards.path}preprocessing/6raw
        rm -rf {wildcards.path}preprocessing/7rgsort
        rm -rf {wildcards.path}preprocessing/8merged
        rm -rf {wildcards.path}preprocessing/9dedup
        rm -rf {wildcards.path}preprocessing/10unique
        rm -rf {wildcards.path}preprocessing/11bigwig
        rm -rf {wildcards.path}preprocessing/saturation
        touch {output}
        """

## This rule determines what is run in the preprocessing spooling option
rule AGGREGATOR_preprocessing:
    input:
        "{path}operations/dir/all.built",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}bigwig/{sample}.rep{repnum}.ref{refgenome}.bw",
        "{path}operations/preprocessing/{sample}.rep{repnum}.ref{refgenome}.peak_calling.complete",
        "{path}operations/preprocessing/{sample}.rep{repnum}.ref{refgenome}.preprocessing_metrics.complete",
        "{path}peaks/samplemerged/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "{path}operations/saturation/{sample}.rep{repnum}.ref{refgenome}.saturation_analysis.complete"
    output:
        "{path}operations/preprocessing/{sample}.rep{repnum}.ref{refgenome}.preprocessing_preclean.complete"
    priority:
        1
    shell:
        "touch {output}"

## This rule determines what is run in the directory building step
rule AGGREGATOR_build_directory_structure:
    input:
        "{path}operations/dir/main.built",
        "{path}operations/dir/operations.built",
        "{path}operations/dir/benchmark.built",
        "{path}operations/dir/metrics.built",
        "{path}operations/dir/preprocessing.built",
        "{path}operations/dir/saturation.built",
        "{path}operations/dir/footprints.built",
        "{path}operations/dir/peaks.built",
        "{path}operations/dir/figures.built"
    output:
        "{path}operations/dir/all.built"
    shell:
        "touch {output}"

## This rule determines what is done for peak calling
rule AGGREGATOR_peak_calling:
    input:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_peaks.narrowPeak",
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p001_peaks.narrowPeak",
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_p001_peaks.narrowPeak",
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p0001_peaks.narrowPeak",
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_p0001_peaks.narrowPeak",
        "{path}peaks/samplemerged/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
    output:
        "{path}operations/preprocessing/{sample}.rep{repnum}.ref{refgenome}.peak_calling.complete"
    priority:
        1
    shell:
        "touch {output}"

## This rule determines what is done for the preprocessing metrics
rule AGGREGATOR_preprocessing_metrics:
    input:
        "{path}metrics/genomecov/{sample}.rep{repnum}.ref{refgenome}.peak.globalnorm.genomecov.txt",
        "{path}metrics/genomecov/{sample}.rep{repnum}.ref{refgenome}.peak.localnorm.genomecov.txt",
        "{path}metrics/totalreads/{sample}.rep{repnum}.ref{refgenome}.totalreads.Rdata",
        "{path}operations/metrics/{sample}.rep{repnum}.ref{refgenome}.fragsizes.done",
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L1.myco.sam",
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L2.myco.sam",
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L3.myco.sam",
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L4.myco.sam"
    output:
        "{path}operations/preprocessing/{sample}.rep{repnum}.ref{refgenome}.preprocessing_metrics.complete"
    priority:
        1
    shell:
        "touch {output}"

########################################################################################################################################
#### FOOTPRINTING AGGREGATORS ##########################################################################################################
########################################################################################################################################

## This rule determines what is run for the footprinting analysis
rule AGGREGATOR_footprinting_analysis:
    input:
        "{path}operations/footprints/{sample}.rep{repnum}.ref{refgenome}.footprinting_raw_analysis.complete",
        "{path}operations/footprints/{sample}.rep{repnum}.ref{refgenome}.footprinting_raw_analysis.samplemergedpeaks.complete",
        "{path}operations/footprints/{sample}.rep{repnum}.ref{refgenome}.samplemergedpeaks.insertionmatrix.complete"
    output:
        "{path}operations/modules/{sample}.rep{repnum}.ref{refgenome}.footprinting_analysis.complete"
    shell:
        "touch {output}"

## This rule initiates the raw footprint analysis for all genes found in the config file
rule AGGREGATOR_footprinting_raw_analysis:
    input:
        expand("{{path}}operations/footprints/raw/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.rawFPanalysis.done", genename=config["geneNames"])
    output:
        "{path}operations/footprints/{sample}.rep{repnum}.ref{refgenome}.footprinting_raw_analysis.complete"
    shell:
        "touch {output}"

## Same as above but using the sample merged peak files
rule AGGREGATOR_footprinting_raw_analysis_sample_merged_peaks:
    input:
        expand("{{path}}operations/footprints/samplemergedpeaks/raw/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.rawFPanalysis.samplemergedpeaks.done", genename=config["geneNames"])
    output:
        "{path}operations/footprints/{sample}.rep{repnum}.ref{refgenome}.footprinting_raw_analysis.samplemergedpeaks.complete"
    shell:
        "touch {output}"

## Same as above but using the sample merged peak files
rule AGGREGATOR_footprinting_insertion_matrix_sample_merged_peaks:
    input:
        expand("{{path}}footprints/samplemergedpeaks/insertionmatrix/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.samplemergedpeaks.insertionmatrix.RData", genename=config["geneNames"])
    output:
        "{path}operations/footprints/{sample}.rep{repnum}.ref{refgenome}.samplemergedpeaks.insertionmatrix.complete"
    shell:
        "touch {output}"

########################################################################################################################################
#### FIGURES AGGREGATORS ###############################################################################################################
########################################################################################################################################

##
rule AGGREGATOR_generate_figures:
    input:
        "{path}figures/peakideogram/{sample}.rep{repnum}.ref{refgenome}.peakIdeogram.svg",
        "{path}operations/figures/{sample}.rep{repnum}.ref{refgenome}.motif_insertion_probability_graphs.complete",
        "{path}operations/figures/{sample}.rep{repnum}.ref{refgenome}.motif_aligned_heatmaps.complete"
    output:
        "{path}operations/modules/{sample}.rep{repnum}.ref{refgenome}.generate_figures.complete"
    shell:
        "touch {output}"

## 
rule AGGREGATOR_motif_insertion_probability_graph:
    input:
        expand("{{path}}figures/insertionprobability/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.insertionprobability.svg", genename=config["geneNames"])
    output:
        "{path}operations/figures/{sample}.rep{repnum}.ref{refgenome}.motif_insertion_probability_graphs.complete"
    shell:
        "touch {output}"

## 
rule AGGREGATOR_motif_aligned_heatmap:
    input:
        expand("{{path}}figures/motifalignedheatmap/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.motifalignedheatmap.svg", genename=config["geneNames"])
    output:
        "{path}operations/figures/{sample}.rep{repnum}.ref{refgenome}.motif_aligned_heatmaps.complete"
    shell:
        "touch {output}"

########################################################################################################################################
#### BUILD DIRECTORY STRUCTURE #########################################################################################################
########################################################################################################################################

rule DIR_main:
    output:
        "{path}operations/dir/main.built"
    shell:
        """
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}operations/dir
        mkdir -p -v {wildcards.path}benchmark
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}correlation
        mkdir -p -v {wildcards.path}bam
        mkdir -p -v {wildcards.path}bigwig
        mkdir -p -v {wildcards.path}figures
        touch {output}
        """

rule DIR_operations:
    output:
        "{path}operations/dir/operations.built"
    shell:
        """
        mkdir -p -v {wildcards.path}operations/modules
        mkdir -p -v {wildcards.path}operations/preprocessing
        mkdir -p -v {wildcards.path}operations/footprints
        mkdir -p -v {wildcards.path}operations/footprints/raw
        mkdir -p -v {wildcards.path}operations/footprints/samplemergedpeaks
        mkdir -p -v {wildcards.path}operations/footprints/samplemergedpeaks/raw {wildcards.path}operations/footprints/samplemergedpeaks/insertionmatrix
        mkdir -p -v {wildcards.path}operations/saturation
        mkdir -p -v {wildcards.path}operations/metrics
        mkdir -p -v {wildcards.path}operations/figures
        touch {output}
        """

rule DIR_benchmark:
    output:
        "{path}operations/dir/benchmark.built"
    shell:
        """
        mkdir -p -v {wildcards.path}benchmark/preprocessing
        mkdir -p -v {wildcards.path}benchmark/preprocessing/gunzip {wildcards.path}benchmark/preprocessing/fastp {wildcards.path}benchmark/preprocessing/mycoalign
        mkdir -p -v {wildcards.path}benchmark/preprocessing/align {wildcards.path}benchmark/preprocessing/coordsortsam {wildcards.path}benchmark/preprocessing/bamconversion
        mkdir -p -v {wildcards.path}benchmark/preprocessing/removemitochondrial {wildcards.path}benchmark/preprocessing/addRG {wildcards.path}benchmark/preprocessing/cleansam
        mkdir -p -v {wildcards.path}benchmark/preprocessing/mergelanes {wildcards.path}benchmark/preprocessing/purgeduplicates {wildcards.path}benchmark/preprocessing/mapqfilter
        mkdir -p -v {wildcards.path}benchmark/preprocessing/buildindex {wildcards.path}benchmark/preprocessing/bigwig {wildcards.path}benchmark/preprocessing/peaks
        mkdir -p -v {wildcards.path}benchmark/metrics
        mkdir -p -v {wildcards.path}benchmark/correlation
        mkdir -p -v {wildcards.path}benchmark/saturation
        mkdir -p -v {wildcards.path}benchmark/footprints
        mkdir -p -v {wildcards.path}benchmark/footprints/raw
        touch {output}
        """

rule DIR_metrics:
    output:
        "{path}operations/dir/metrics.built"
    shell:
        """
        mkdir -p -v {wildcards.path}metrics/saturation
        mkdir -p -v {wildcards.path}metrics/fastq
        mkdir -p -v {wildcards.path}metrics/myco
        mkdir -p -v {wildcards.path}metrics/align
        mkdir -p -v {wildcards.path}metrics/genomecov
        mkdir -p -v {wildcards.path}metrics/totalreads
        mkdir -p -v {wildcards.path}metrics/fragsize
        mkdir -p -v {wildcards.path}metrics/peakannotation
        mkdir -p -v {wildcards.path}metrics/duplication
        touch {output}
        """

rule DIR_preprocessing:
    output:
        "{path}operations/dir/preprocessing.built"
    shell:
        """
        mkdir -p -v {wildcards.path}preprocessing/2fastq
        mkdir -p -v {wildcards.path}preprocessing/3goodfastq
        mkdir -p -v {wildcards.path}preprocessing/4mycoalign
        mkdir -p -v {wildcards.path}preprocessing/5align
        mkdir -p -v {wildcards.path}preprocessing/6raw
        mkdir -p -v {wildcards.path}preprocessing/6raw/mitochondrial {wildcards.path}preprocessing/6raw/blacklist {wildcards.path}preprocessing/6raw/nonblacklist
        mkdir -p -v {wildcards.path}preprocessing/7rgsort
        mkdir -p -v {wildcards.path}preprocessing/8merged
        mkdir -p -v {wildcards.path}preprocessing/9dedup
        mkdir -p -v {wildcards.path}preprocessing/10unique
        mkdir -p -v {wildcards.path}preprocessing/11bigwig
        mkdir -p -v {wildcards.path}preprocessing/12saturation
        touch {output}
        """

rule DIR_saturation:
    output:
        "{path}operations/dir/saturation.built"
    shell:
        """
        mkdir -p -v {wildcards.path}preprocessing/12saturation/downsampled {wildcards.path}preprocessing/12saturation/downsampled/raw
        mkdir -p -v {wildcards.path}preprocessing/12saturation/downsampled/sorted {wildcards.path}preprocessing/12saturation/downsampled/deduplicated
        mkdir -p -v {wildcards.path}preprocessing/12saturation/duplication
        mkdir -p -v {wildcards.path}preprocessing/12saturation/peaks
        mkdir -p -v {wildcards.path}preprocessing/12saturation/footprints
        touch {output}
        """

rule DIR_footprints:
    output:
        "{path}operations/dir/footprints.built"
    shell:
        """
        mkdir -p -v {wildcards.path}footprints/raw
        mkdir -p -v {wildcards.path}footprints/aggregated
        mkdir -p -v {wildcards.path}footprints/figures
        mkdir -p -v {wildcards.path}footprints/samplemergedpeaks
        mkdir -p -v {wildcards.path}footprints/samplemergedpeaks/raw {wildcards.path}footprints/samplemergedpeaks/insertionmatrix
        touch {output}
        """

rule DIR_peaks:
    output:
        "{path}operations/dir/peaks.built"
    shell:
        """
        mkdir -p -v {wildcards.path}peaks/localnorm {wildcards.path}peaks/globalnorm {wildcards.path}peaks/samplemerged
        touch {output}
        """

rule DIR_figures:
    output:
        "{path}operations/dir/figures.built"
    shell:
        """
        mkdir -p -v {wildcards.path}figures/peakideogram
        mkdir -p -v {wildcards.path}figures/insertionprobability
        mkdir -p -v {wildcards.path}figures/motifalignedheatmap
        touch {output}
        """



########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################

## Gunzip the fastq files
rule STEP1_gunzip:
    input:
        a="{path}fastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R{read}.fastq.gz",
        b="{path}operations/dir/all.built"
    output:
        c="{path}preprocessing/2fastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R{read}.fastq"
    threads:
        1
    benchmark:
        '{path}benchmark/preprocessing/gunzip/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R{read}.gunzip.benchmark.txt'
    shell:
        "gunzip -c {input.a} > {output.c}"

## Fastq QC filtering
rule STEP2_fastp_filtering:
    input:
        a="{path}preprocessing/2fastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.fastq",
        b="{path}preprocessing/2fastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.fastq"
    output:
        c="{path}preprocessing/3goodfastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.good.fq",
        d="{path}preprocessing/3goodfastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.good.fq"
    benchmark:
        '{path}benchmark/preprocessing/fastp/{sample}.rep{repnum}.ref{refgenome}.{lane}.fastp.benchmark.txt'
    threads:
        6
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    conda:
    	"snakeResources/envs/fastp.yaml"
    shell:
        "fastp -i {input.a} -I {input.b} -o {output.c} -O {output.d} -w 6 -h {wildcards.path}metrics/fastq/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}.fastq.quality.html -j {wildcards.path}metrics/fastq/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}.fastq.quality.json"
    
## Check for mycoplasma contamination
rule STEP3_mycoalign:
    input:
        a="{path}preprocessing/3goodfastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.good.fq"
    output:
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L{lane}.myco.sam"
    benchmark:
        '{path}benchmark/preprocessing/mycoalign/{sample}.rep{repnum}.ref{refgenome}.{lane}.mycoalign.benchmark.txt'
    threads:
        6
    conda:
        "snakeResources/envs/bowtie2.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 2
    shell:
        "bowtie2 -q -p 12 -X2000 -x genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/myco/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}.myco.alignment.txt"
    
## Align reads to reference genome
rule STEP4_refgenome_align:
    input:
        a="{path}preprocessing/3goodfastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R2.good.fq"
    output:
        "{path}preprocessing/5align/{sample}.rep{repnum}.ref{refgenome}_L{lane}.sam"
    benchmark:
        '{path}benchmark/preprocessing/align/{sample}.rep{repnum}.ref{refgenome}.{lane}.align.benchmark.txt'
    threads:
        12
    conda:
        "snakeResources/envs/bowtie2.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "bowtie2 -q -p 12 -X2000 -x genomes/{wildcards.refgenome}/{wildcards.refgenome} -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/align/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}.alignment.txt"

## Coordinate sort the aligned reads. This is required for blacklist filtering
rule STEP5_coordsort_sam:
    input:
        "{path}preprocessing/5align/{sample}.rep{repnum}.ref{refgenome}_L{lane}.sam"
    output:
        "{path}preprocessing/5align/{sample}.rep{repnum}.ref{refgenome}_L{lane}.cs.sam"
    threads:
        1
    conda:
        "snakeResources/envs/samtools.yaml"
    benchmark:
        '{path}benchmark/preprocessing/coordsortsam/{sample}.rep{repnum}.ref{refgenome}.{lane}.coordsort.benchmark.txt'
    shell:
        "samtools sort {input} -o {output} -O sam"
    
## Remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
rule STEP6_blacklistfilter_bamconversion:
    input:
        "{path}preprocessing/5align/{sample}.rep{repnum}.ref{refgenome}_L{lane}.cs.sam"
    output:
        a="{path}preprocessing/6raw/blacklist/{sample}.rep{repnum}.ref{refgenome}_L{lane}.blacklist.bam",
        b="{path}preprocessing/6raw/nonblacklist/{sample}.rep{repnum}.ref{refgenome}_L{lane}.blrm.bam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        5
    benchmark:
        '{path}benchmark/preprocessing/bamconversion/{sample}.rep{repnum}.ref{refgenome}.{lane}.bamconvert.benchmark.txt'
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 4 {input}"
    
## Remove reads mapping to mitochondrial DNA
rule STEP7_chrM_contamination:
    input:
        "{path}preprocessing/6raw/nonblacklist/{sample}.rep{repnum}.ref{refgenome}_L{lane}.blrm.bam"
    output:
        a="{path}preprocessing/6raw/mitochondrial/{sample}.rep{repnum}.ref{refgenome}_L{lane}.mitochondrial.bam",
        b="{path}preprocessing/6raw/{sample}.rep{repnum}.ref{refgenome}_L{lane}.goodbam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        5
    benchmark:
        '{path}benchmark/preprocessing/removemitochondrial/{sample}.rep{repnum}.ref{refgenome}.{lane}.chrMfilter.benchmark.txt'
    shell:
        "samtools view -b -h -o {output.a} -L genomes/mtdna/mtdna.extents.bed -U {output.b} -@ 4 {input}"
    
## Add @RG tags to the reads and perform coordinate sorting
rule STEP8_addrgandcsbam:
    input:
        "{path}preprocessing/6raw/{sample}.rep{repnum}.ref{refgenome}_L{lane}.goodbam"
    output:
        "{path}preprocessing/7rgsort/{sample}.rep{repnum}.ref{refgenome}_L{lane}.rg.cs.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 1
    benchmark:
        '{path}benchmark/preprocessing/addRG/{sample}.rep{repnum}.ref{refgenome}.{lane}.addRGtag.benchmark.txt'
    shell:
        "picard AddOrReplaceReadGroups \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate \
        RGID=H5YHHBGX3.{wildcards.lane} \
        RGLB={wildcards.sample} \
        RGPL=ILLUMINA \
        RGPU=H5YHHBGX3.{wildcards.lane}.{wildcards.sample} \
        RGSM={wildcards.sample}"
    
## Clean the bam file
rule STEP9_cleansam:
    input:
        "{path}preprocessing/7rgsort/{sample}.rep{repnum}.ref{refgenome}_L{lane}.rg.cs.bam"
    output:
        "{path}preprocessing/7rgsort/{sample}.rep{repnum}.ref{refgenome}_L{lane}.clean.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/preprocessing/cleansam/{sample}.rep{repnum}.ref{refgenome}.{lane}.cleansam.benchmark.txt'
    shell:
        "picard CleanSam \
        I={input} \
        O={output}"
    
## Merge reads from different NextSeq lanes
rule STEP10_mergelanes:
    input:
        a="{path}preprocessing/7rgsort/{sample}.rep{repnum}.ref{refgenome}_L1.clean.bam",
        b="{path}preprocessing/7rgsort/{sample}.rep{repnum}.ref{refgenome}_L2.clean.bam",
        c="{path}preprocessing/7rgsort/{sample}.rep{repnum}.ref{refgenome}_L3.clean.bam",
        d="{path}preprocessing/7rgsort/{sample}.rep{repnum}.ref{refgenome}_L4.clean.bam"
    output:
        "{path}preprocessing/8merged/{sample}.rep{repnum}.ref{refgenome}.lanemerge.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 1
    benchmark:
        '{path}benchmark/preprocessing/mergelanes/{sample}.rep{repnum}.ref{refgenome}.mergelanes.benchmark.txt'
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        I={input.c} \
        I={input.d} \
        O={output} \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true"

## Purge PCR duplicate reads
rule STEP11_purgeduplicates:
    input:
        "{path}preprocessing/8merged/{sample}.rep{repnum}.ref{refgenome}.lanemerge.bam"
    output:
        "{path}preprocessing/9dedup/{sample}.rep{repnum}.ref{refgenome}.dp.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 3
    benchmark:
        '{path}benchmark/preprocessing/purgeduplicates/{sample}.rep{repnum}.ref{refgenome}.purgeduplicates.benchmark.txt'
    shell:
        "picard MarkDuplicates \
        I={input} \
        O={output} \
        M={wildcards.path}metrics/duplication/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}.duplication.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"
    
## Filter reads for only uniquely mapping
rule STEP12_mapqfilter:
    input:
        "{path}preprocessing/9dedup/{sample}.rep{repnum}.ref{refgenome}.dp.bam"
    output:
        "{path}preprocessing/10unique/{sample}.rep{repnum}.ref{refgenome}.u.bam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        1
    resources:
        run_time=lambda params, attempt: attempt * 4
    benchmark:
        '{path}benchmark/preprocessing/mapqfilter/{sample}.rep{repnum}.ref{refgenome}.mapqfilter.benchmark.txt'
    shell:
        "samtools view -h -q 2 -b {input} > {output}"
    
## Move bam to parent directory and rename
rule STEP13_move_bam:
    input:
        "{path}preprocessing/10unique/{sample}.rep{repnum}.ref{refgenome}.u.bam"
    output:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam"
    shell:
        """
        mv {wildcards.path}preprocessing/10unique/*rep{wildcards.repnum}*.bam {wildcards.path}bam/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}.bam
        touch {output}
        """

## Build the .bai index for the processed bam file
rule STEP14_build_bai_index:
    input:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam"
    output:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 3
    benchmark:
        '{path}benchmark/preprocessing/buildindex/{sample}.rep{repnum}.ref{refgenome}.buildindex.benchmark.txt'
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"
    
## Make a bigwig file from the bam file
rule STEP15_makebigwig_bamcov:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}bigwig/{sample}.rep{repnum}.ref{refgenome}.bw"
    conda:
        "snakeResources/envs/deeptools.yaml"
    threads:
        6
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    benchmark:
        '{path}benchmark/preprocessing/bigwig/{sample}.rep{repnum}.ref{refgenome}.makebigwig.benchmark.txt'
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 6 -v"
    
########################################################################################################################################
#### PEAK CALLING ######################################################################################################################
########################################################################################################################################

## Call peaks using global normalization. pvalue 0.01
rule STEP16_MACS2_peaks_global_normilization_p01:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/preprocessing/peaks/{sample}.rep{repnum}.ref{refgenome}.callpeaks.globalnorm.benchmark.txt'
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_globalnorm --outdir {wildcards.path}peaks/globalnorm --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
    
## Call peaks using local normalization. pvalue 0.01
rule STEP17_MACS2_peaks_local_normalization_p01:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/preprocessing/peaks/{sample}.rep{repnum}.ref{refgenome}.callpeaks.localnorm.benchmark.txt'
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_localnorm --outdir {wildcards.path}peaks/localnorm --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

## Call peaks using global normalization. pvalue 0.001
rule STEP16_MACS2_peaks_global_normilization_p001:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p001_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_globalnorm_p001 --outdir {wildcards.path}peaks/globalnorm --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.001"

## Call peaks using local normalization. pvalue 0.001
rule STEP17_MACS2_peaks_local_normalization_p001:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_p001_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_localnorm_p001 --outdir {wildcards.path}peaks/localnorm --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.001"

## Call peaks using global normalization. pvalue 0.0001
rule STEP16_MACS2_peaks_global_normilization_p0001:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p0001_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_globalnorm_p0001 --outdir {wildcards.path}peaks/globalnorm --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.0001"

## Call peaks using local normalization. pvalue 0.001
rule STEP17_MACS2_peaks_local_normalization_p0001:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_p0001_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_localnorm_p0001 --outdir {wildcards.path}peaks/localnorm --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.0001"

########################################################################################################################################
#### QC METRICS  #######################################################################################################################
########################################################################################################################################
    
## Calculate percent genome coverage from peaks with global normalization
rule METRICS_percent_peak_genome_coverage_globalnorm:
    input:
        a="{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/genomecov/{sample}.rep{repnum}.ref{refgenome}.peak.globalnorm.genomecov.txt"
    conda:
        "snakeResources/envs/bedops.yaml"
    threads:
        1  	
    benchmark:
        '{path}benchmark/metrics/{sample}.rep{repnum}.ref{refgenome}.genomecov.globalnorm.benchmark.txt'
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"
    
## Calculate percent genome coverage from peaks with local normalization
rule METRICS_percent_peak_genome_coverage_localnorm:
    input:
        a="{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/genomecov/{sample}.rep{repnum}.ref{refgenome}.peak.localnorm.genomecov.txt"
    conda:
        "snakeResources/envs/bedops.yaml"
    threads:
        1  
    benchmark:
        '{path}benchmark/metrics/{sample}.rep{repnum}.ref{refgenome}.genomecov.localnorm.benchmark.txt'
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"
    
## Generate the fragment size distribution graph
rule METRICS_fragment_size_distribution:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}operations/metrics/{sample}.rep{repnum}.ref{refgenome}.fragsizes.done"
    conda:
        "snakeResources/envs/fragSizeDistribution.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/metrics/{sample}.rep{repnum}.ref{refgenome}.fragsizes.benchmark.txt'
    script:
        "snakeResources/scripts/generateFragSizeDistribution.R"

## Annotate the peaks with global normalization
rule METRICS_annotate_peaks_global:
    input:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
    output:
        "{path}operations/metrics/peakannotation/{sample}.rep{repnum}.ref{refgenome}.globalpeak.annotations.done"
    conda:
        "snakeResources/envs/annotatePeaks.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/metrics/peakannotation/global/{sample}.rep{repnum}.ref{refgenome}.globalpeak.annotations.benchmark.txt'
    script:
        "snakeResources/scripts/annotatePeaks.R"

## Annotate the peaks with local normalization
rule METRICS_annotate_peaks_local:
    input:
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_peaks.narrowPeak",
    output:
        "{path}operations/metrics/peakannotation/{sample}.rep{repnum}.ref{refgenome}.localpeak.annotations.done"
    conda:
        "snakeResources/envs/annotatePeaks.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/metrics/peakannotation/local/{sample}.rep{repnum}.ref{refgenome}.localpeak.annotations.benchmark.txt'
    script:
        "snakeResources/scripts/annotatePeaks.R"

## Count the total number of reads in the sample
rule METRICS_sample_total_reads:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}metrics/totalreads/{sample}.rep{repnum}.ref{refgenome}.totalreads.Rdata"
    conda:
        "snakeResources/envs/countSampleReads.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/metrics/{sample}.rep{repnum}.ref{refgenome}.totalreads.benchmark.txt'
    script:
        "snakeResources/scripts/countTotalSampleReads.R"

########################################################################################################################################
#### SATURATION ANALYSIS ###############################################################################################################
########################################################################################################################################

## This rule determines what is run for the library saturation analysis
rule AGGREGATOR_saturation_analysis:
    input:
        "{path}operations/saturation/{sample}.rep{repnum}.ref{refgenome}.downsampling.done",
        "{path}metrics/saturation/{sample}.rep{repnum}.ref{refgenome}.downsampled_duplication_metrics.txt",
        "{path}metrics/saturation/{sample}.rep{repnum}.ref{refgenome}.downsampled_numpeaks.txt",
        "{path}operations/saturation/{sample}.rep{repnum}.ref{refgenome}.saturation_footprint_analysis.complete"
    output:
        "{path}operations/saturation/{sample}.rep{repnum}.ref{refgenome}.saturation_analysis.complete"
    shell:
    	"touch {output}"

## Downsample the processed but NOT duplicate purged .bam files
rule SATURATION_downsample_bam:
    input:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam"
    output:
        "{path}preprocessing/12saturation/downsampled/raw/{sample}.rep{repnum}.ref{refgenome}.{prob}.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 2
    benchmark:
        '{path}benchmark/saturation/{sample}.rep{repnum}.ref{refgenome}.{prob}.downsample.benchmark.txt'
    shell:
        "picard DownsampleSam \
        I={input} \
        O={output} \
        PROBABILITY=0.{wildcards.prob}"

## Coordinate sort the downsampled .bam files
rule SATURATION_sort_downsampled:
    input:
        "{path}preprocessing/12saturation/downsampled/raw/{sample}.rep{repnum}.ref{refgenome}.{prob}.bam"
    output:
        "{path}preprocessing/12saturation/downsampled/sorted/{sample}.rep{repnum}.ref{refgenome}.{prob}.sorted.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 3
    benchmark:
        '{path}benchmark/saturation/{sample}.rep{repnum}.ref{refgenome}.{prob}.sort.benchmark.txt'
    shell:
        "picard SortSam \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate"

## Purge duplicates from the downsampled .bam files
rule SATURATION_purge_duplicates:
    input:
        "{path}preprocessing/12saturation/downsampled/sorted/{sample}.rep{repnum}.ref{refgenome}.{prob}.sorted.bam"
    output:
        a="{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.{prob}.deduplicated.bam",
        b="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.{prob}.duplication-metrics.txt"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 3
    benchmark:
        '{path}benchmark/saturation/{sample}.rep{repnum}.ref{refgenome}.{prob}.purgeduplicates.benchmark.txt'
    shell:
        "picard MarkDuplicates \
        I={input} \
        O={output.a} \
        M={output.b} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"

## Generate .bai index for each downsampled .bam file
rule SATURATION_index_downsampled:
    input:
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.{prob}.deduplicated.bam"
    output:
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.{prob}.deduplicated.bam.bai"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 20000,
        run_time=lambda params, attempt: attempt * 1
    benchmark:
        '{path}benchmark/saturation/{sample}.rep{repnum}.ref{refgenome}.{prob}.index.benchmark.txt'
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"

## Aggregator rule for all the downsampling probabilities
rule SATURATION_downsample_aggregator:
    input:
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.9.deduplicated.bam.bai",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.8.deduplicated.bam.bai",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.7.deduplicated.bam.bai",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.6.deduplicated.bam.bai",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.5.deduplicated.bam.bai",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.4.deduplicated.bam.bai",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.3.deduplicated.bam.bai",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.2.deduplicated.bam.bai",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.1.deduplicated.bam.bai"    
    output:
        "{path}operations/saturation/{sample}.rep{repnum}.ref{refgenome}.downsampling.done"
    shell:
        "touch {output}"

## Determine the library complexity of the downsampled libraries and output to metrics
rule SATURATION_parse_duplication_metrics_downsampled:
    input:
        a="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.9.duplication-metrics.txt",
        b="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.8.duplication-metrics.txt",
        c="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.7.duplication-metrics.txt",
        d="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.6.duplication-metrics.txt",
        e="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.5.duplication-metrics.txt",
        f="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.4.duplication-metrics.txt",
        g="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.3.duplication-metrics.txt",
        h="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.2.duplication-metrics.txt",
        i="{path}preprocessing/12saturation/duplication/{sample}.rep{repnum}.ref{refgenome}.1.duplication-metrics.txt"
    output:
        "{path}metrics/saturation/{sample}.rep{repnum}.ref{refgenome}.downsampled_duplication_metrics.txt"
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

## Call peaks from downsampled libraries
rule SATURATION_peaks:
    input:
        a="{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.{prob}.deduplicated.bam",
        b="{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.{prob}.deduplicated.bam.bai"
    output:
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.{prob}_globalnorm_peaks.xls"
    conda:
        "snakeResources/envs/macs2.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.{prob}.downsampled.peak.benchmark.txt'
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}.{wildcards.prob}_globalnorm --outdir {wildcards.path}preprocessing/12saturation/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

## Count the number of peaks with global normalization from downsampled libraries
rule SATURATION_peak_calling_aggregator:
    input:
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.9_globalnorm_peaks.xls",
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.8_globalnorm_peaks.xls",
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.7_globalnorm_peaks.xls",
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.6_globalnorm_peaks.xls",
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.5_globalnorm_peaks.xls",
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.4_globalnorm_peaks.xls",
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.3_globalnorm_peaks.xls",
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.2_globalnorm_peaks.xls",
        "{path}preprocessing/12saturation/peaks/{sample}.rep{repnum}.ref{refgenome}.1_globalnorm_peaks.xls"
    output:
        "{path}metrics/saturation/{sample}.rep{repnum}.ref{refgenome}.downsampled_numpeaks.txt"
    shell:
        "wc -l < {input} >> {output}"

## Footprint analysis
rule SATURATION_footprint_raw_analysis:
    input:
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.{prob}.deduplicated.bam",
        "{path}preprocessing/12saturation/downsampled/deduplicated/{sample}.rep{repnum}.ref{refgenome}.{prob}.deduplicated.bam.bai",
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}operations/saturation/footprints/{sample}.rep{repnum}.ref{refgenome}.{gene}.{prob}.downsampled.rawFPanalysis.done"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/analyzeFPsaturation.R"

## An aggregator for the footprinting saturation analysis probabilities
rule AGGREGATOR_footprinting_saturation:
    input:
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.9.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"]),
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.8.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"]),
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.7.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"]),
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.6.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"]),
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.5.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"]),
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.4.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"]),
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.3.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"]),
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.2.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"]),
        expand("{{path}}operations/saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.1.downsampled.rawFPanalysis.done", genename=config["saturationGeneNames"])
    output:
        "{path}operations/saturation/{sample}.rep{repnum}.ref{refgenome}.saturation_footprint_analysis.complete"
    shell:
        "touch {output}"

########################################################################################################################################
#### FOOTPRINTING ######################################################################################################################
########################################################################################################################################

##
rule FOOTPRINTING_raw_analysis:
    input:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}operations/footprints/raw/{sample}.rep{repnum}.ref{refgenome}.{gene}.rawFPanalysis.done"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    benchmark:
        '{path}benchmark/footprints/raw/{sample}.rep{repnum}.ref{refgenome}.{gene}.rawFPanalysis.benchmark.txt'
    script:
        "snakeResources/scripts/analyzeFP.R"

##
rule FOOTPRINTING_raw_analysis_sample_merged_peaks:
    input:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}peaks/samplemerged/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}operations/footprints/samplemergedpeaks/raw/{sample}.rep{repnum}.ref{refgenome}.{gene}.rawFPanalysis.samplemergedpeaks.done"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    benchmark:
        '{path}benchmark/footprints/raw/{sample}.rep{repnum}.ref{refgenome}.{gene}.rawFPanalysis.samplemergedpeaks.benchmark.txt'
    script:
        "snakeResources/scripts/analyzeFPsampleMergedPeaks.R"

##
rule FOOTPRINTING_generate_insertion_matrix_sample_merged_peaks:
    input:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}peaks/samplemerged/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/samplemergedpeaks/insertionmatrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.samplemergedpeaks.insertionmatrix.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    benchmark:
        '{path}benchmark/footprints/raw/{sample}.rep{repnum}.ref{refgenome}.{gene}.generateinsertionmatrix.benchmark.txt'
    script:
        "snakeResources/scripts/generateInsertionMatrix.R"

########################################################################################################################################
#### GENERATE FIGURES ##################################################################################################################
########################################################################################################################################

##
rule FIGURES_generate_peak_idiogram:
    input:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p0001_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}figures/peakideogram/{sample}.rep{repnum}.ref{refgenome}.peakIdeogram.svg"
    conda:
        "snakeResources/envs/generatePeakIdeogram.yaml"
    script:
        "snakeResources/scripts/generatePeakIdeogram.R"

##
rule FIGURES_generate_motif_insertion_probability_graph_sample_merged_peaks:
    input:
        "{path}footprints/samplemergedpeaks/insertionmatrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.samplemergedpeaks.insertionmatrix.RData",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}figures/insertionprobability/{sample}.rep{repnum}.ref{refgenome}.{gene}.insertionprobability.svg"
    conda:
        "snakeResources/envs/generateMotifInsertionProbabilityGraph.yaml"
    script:
        "snakeResources/scripts/generateMotifInsertionProbabilityGraph.R"

##
rule FIGURES_generate_motif_aligned_heatmap_sample_merged_peaks:
    input:
        "{path}footprints/samplemergedpeaks/insertionmatrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.samplemergedpeaks.insertionmatrix.RData",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}figures/motifalignedheatmap/{sample}.rep{repnum}.ref{refgenome}.{gene}.motifalignedheatmap.svg"
    conda:
        "snakeResources/envs/generateMotifAlignedHeatmap.yaml"
    script:
        "snakeResources/scripts/generateMotifAlignedHeatmap.R"


########################################################################################################################################
#### MERGE SAMPLE PEAKS ################################################################################################################
########################################################################################################################################

## LNCAP
rule MERGE_sample_peaks_lncap:
    input:
        "data/pros/lncap/cr01/peaks/globalnorm/LNCaP-CR-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/pros/lncap/cr02/peaks/globalnorm/LNCaP-CR-02.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/pros/lncap/cr04/peaks/globalnorm/LNCaP-CR-04.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/pros/lncap/cr05/peaks/globalnorm/LNCaP-CR-05.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/pros/lncap/cr07/peaks/globalnorm/LNCaP-CR-07.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/pros/lncap/cr08/peaks/globalnorm/LNCaP-CR-08.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/pros/lncap/wt01/peaks/globalnorm/LNCaP-WT-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/pros/lncap/wt02/peaks/globalnorm/LNCaP-WT-02.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/pros/lncap/cr01/peaks/samplemerged/LNCaP-CR-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

## LNCAP
rule COPY_sample_merged_peaks_lncap:
    input:
        "data/pros/lncap/cr01/peaks/samplemerged/LNCaP-CR-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        "data/pros/lncap/cr02/peaks/samplemerged/LNCaP-CR-02.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr04/peaks/samplemerged/LNCaP-CR-04.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr05/peaks/samplemerged/LNCaP-CR-05.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr07/peaks/samplemerged/LNCaP-CR-07.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr08/peaks/samplemerged/LNCaP-CR-08.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/wt01/peaks/samplemerged/LNCaP-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/wt02/peaks/samplemerged/LNCaP-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} "data/pros/lncap/cr02/peaks/samplemerged/LNCaP-CR-02.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/cr04/peaks/samplemerged/LNCaP-CR-04.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/cr05/peaks/samplemerged/LNCaP-CR-05.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/cr07/peaks/samplemerged/LNCaP-CR-07.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/cr08/peaks/samplemerged/LNCaP-CR-08.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/wt01/peaks/samplemerged/LNCaP-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/wt02/peaks/samplemerged/LNCaP-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak"
        """

########################################################################################################################################
#### MERGE SAMPLE RUNS #################################################################################################################
########################################################################################################################################

#### ####
rule panc_merge_all_sample_runs:
    input:
        "data/panc/capan1/wt01/bam/CAPANI-WT-01-MERGED.rep1.refhg38.bam",
        "data/panc/capan1/wt02/bam/CAPANI-WT-02-MERGED.rep1.refhg38.bam",
        "data/panc/capan1/wt03/bam/CAPANI-WT-03-MERGED.rep1.refhg38.bam",
        "data/panc/hpafii/wt01/bam/HPAFII-WT-01-MERGED.rep1.refhg38.bam",
        "data/panc/hpafii/wt02/bam/HPAFII-WT-02-MERGED.rep1.refhg38.bam",
        "data/panc/hpafii/wt03/bam/HPAFII-WT-03-MERGED.rep1.refhg38.bam",
        "data/panc/kp4/wt01/bam/KP4-WT-01-MERGED.rep1.refhg38.bam",
        "data/panc/kp4/wt02/bam/KP4-WT-02-MERGED.rep1.refhg38.bam",
        "data/panc/kp4/wt03/bam/KP4-WT-03-MERGED.rep1.refhg38.bam",
        "data/panc/panc1/wt01/bam/PANC1-WT-01-MERGED.rep1.refhg38.bam",
        "data/panc/panc1/wt02/bam/PANC1-WT-02-MERGED.rep1.refhg38.bam",
        "data/panc/panc1/wt03/bam/PANC1-WT-03-MERGED.rep1.refhg38.bam",
        "data/panc/panc0403/wt01/bam/PANC0403-WT-01-MERGED.rep1.refhg38.bam",
        "data/panc/panc0403/wt02/bam/PANC0403-WT-02-MERGED.rep1.refhg38.bam",
        "data/panc/panc0403/wt03/bam/PANC0403-WT-03-MERGED.rep1.refhg38.bam",
        "data/panc/patu8ss89/wt01/bam/PATU8SS89-WT-01-MERGED.rep1.refhg38.bam",
        "data/panc/patu8ss89/wt02/bam/PATU8SS89-WT-02-MERGED.rep1.refhg38.bam",
        "data/panc/patu8ss89/wt03/bam/PATU8SS89-WT-03-MERGED.rep1.refhg38.bam",
        "data/panc/pk45h/wt01/bam/PK45H-WT-01-MERGED.rep1.refhg38.bam",
        "data/panc/pk45h/wt02/bam/PK45H-WT-02-MERGED.rep1.refhg38.bam",
        "data/panc/pk45h/wt03/bam/PK45H-WT-03-MERGED.rep1.refhg38.bam"

#### CAPANI ####
rule MERGE_sample_runs_capan1_wt01:
    input:
        a="data/panc/capan1/split/wt01r1/bam/CAPANI-WT-01-RUN1.rep1.refhg38.bam",
        b="data/panc/capan1/split/wt01r2/bam/CAPANI-WT-01-RUN2.rep1.refhg38.bam",
        c="data/panc/capan1/wt01/operations/dir/all.built"
    output:
        "data/panc/capan1/wt01/bam/CAPANI-WT-01-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_capan1_wt02:
    input:
        a="data/panc/capan1/split/wt02r1/bam/CAPANI-WT-02-RUN1.rep1.refhg38.bam",
        b="data/panc/capan1/split/wt02r2/bam/CAPANI-WT-02-RUN2.rep1.refhg38.bam",
        c="data/panc/capan1/wt02/operations/dir/all.built"
    output:
        "data/panc/capan1/wt02/bam/CAPANI-WT-02-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_capan1_wt03:
    input:
        a="data/panc/capan1/split/wt03r1/bam/CAPANI-WT-03-RUN1.rep1.refhg38.bam",
        b="data/panc/capan1/split/wt03r2/bam/CAPANI-WT-03-RUN2.rep1.refhg38.bam",
        c="data/panc/capan1/wt03/operations/dir/all.built"
    output:
        "data/panc/capan1/wt03/bam/CAPANI-WT-03-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

#### HPAFII ####
rule MERGE_sample_runs_hpafii_wt01:
    input:
        a="data/panc/hpafii/split/wt01r1/bam/HPAFII-WT-01-RUN1.rep1.refhg38.bam",
        b="data/panc/hpafii/split/wt01r2/bam/HPAFII-WT-01-RUN2.rep1.refhg38.bam",
        c="data/panc/hpafii/wt01/operations/dir/all.built"
    output:
        "data/panc/hpafii/wt01/bam/HPAFII-WT-01-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_hpafii_wt02:
    input:
        a="data/panc/hpafii/split/wt02r1/bam/HPAFII-WT-02-RUN1.rep1.refhg38.bam",
        b="data/panc/hpafii/split/wt02r2/bam/HPAFII-WT-02-RUN2.rep1.refhg38.bam",
         c="data/panc/hpafii/wt02/operations/dir/all.built"
    output:
        "data/panc/hpafii/wt02/bam/HPAFII-WT-02-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_hpafii_wt03:
    input:
        a="data/panc/hpafii/split/wt03r1/bam/HPAFII-WT-03-RUN1.rep1.refhg38.bam",
        b="data/panc/hpafii/split/wt03r2/bam/HPAFII-WT-03-RUN2.rep1.refhg38.bam",
         c="data/panc/hpafii/wt03/operations/dir/all.built"
    output:
        "data/panc/hpafii/wt03/bam/HPAFII-WT-03-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

#### KP4 ####
rule MERGE_sample_runs_kp4_wt01:
    input:
        a="data/panc/kp4/split/wt01r1/bam/KP4-WT-01-RUN1.rep1.refhg38.bam",
        b="data/panc/kp4/split/wt01r2/bam/KP4-WT-01-RUN2.rep1.refhg38.bam",
        c="data/panc/kp4/wt01/operations/dir/all.built"
    output:
        "data/panc/kp4/wt01/bam/KP4-WT-01-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_kp4_wt02:
    input:
        a="data/panc/kp4/split/wt02r1/bam/KP4-WT-02-RUN1.rep1.refhg38.bam",
        b="data/panc/kp4/split/wt02r2/bam/KP4-WT-02-RUN2.rep1.refhg38.bam",
        c="data/panc/kp4/wt02/operations/dir/all.built"
    output:
        "data/panc/kp4/wt02/bam/KP4-WT-02-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_kp4_wt03:
    input:
        a="data/panc/kp4/split/wt03r1/bam/KP4-WT-03-RUN1.rep1.refhg38.bam",
        b="data/panc/kp4/split/wt03r2/bam/KP4-WT-03-RUN2.rep1.refhg38.bam",
        c="data/panc/kp4/wt03/operations/dir/all.built"
    output:
        "data/panc/kp4/wt03/bam/KP4-WT-03-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

#### PANC1 ####
rule MERGE_sample_runs_panc1_wt01:
    input:
        a="data/panc/panc1/split/wt01r1/bam/PANC1-WT-01-RUN1.rep1.refhg38.bam",
        b="data/panc/panc1/split/wt01r2/bam/PANC1-WT-01-RUN2.rep1.refhg38.bam",
        c="data/panc/panc1/wt01/operations/dir/all.built"
    output:
        "data/panc/panc1/wt01/bam/PANC1-WT-01-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_panc1_wt02:
    input:
        a="data/panc/panc1/split/wt02r1/bam/PANC1-WT-02-RUN1.rep1.refhg38.bam",
        b="data/panc/panc1/split/wt02r2/bam/PANC1-WT-02-RUN2.rep1.refhg38.bam",
        c="data/panc/panc1/wt02/operations/dir/all.built"
    output:
        "data/panc/panc1/wt02/bam/PANC1-WT-02-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_panc1_wt03:
    input:
        a="data/panc/panc1/split/wt03r1/bam/PANC1-WT-03-RUN1.rep1.refhg38.bam",
        b="data/panc/panc1/split/wt03r2/bam/PANC1-WT-03-RUN2.rep1.refhg38.bam",
        c="data/panc/panc1/wt03/operations/dir/all.built"
    output:
        "data/panc/panc1/wt03/bam/PANC1-WT-03-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

#### PANC0403 ####
rule MERGE_sample_runs_panc0403_wt01:
    input:
        a="data/panc/panc0403/split/wt01r1/bam/PANC0403-WT-01-RUN1.rep1.refhg38.bam",
        b="data/panc/panc0403/split/wt01r2/bam/PANC0403-WT-01-RUN2.rep1.refhg38.bam",
        c="data/panc/panc0403/wt01/operations/dir/all.built"
    output:
        "data/panc/panc0403/wt01/bam/PANC0403-WT-01-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_panc0403_wt02:
    input:
        a="data/panc/panc0403/split/wt02r1/bam/PANC0403-WT-02-RUN1.rep1.refhg38.bam",
        b="data/panc/panc0403/split/wt02r2/bam/PANC0403-WT-02-RUN2.rep1.refhg38.bam",
        c="data/panc/panc0403/wt02/operations/dir/all.built"
    output:
        "data/panc/panc0403/wt02/bam/PANC0403-WT-02-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_panc0403_wt03:
    input:
        a="data/panc/panc0403/split/wt03r1/bam/PANC0403-WT-03-RUN1.rep1.refhg38.bam",
        b="data/panc/panc0403/split/wt03r2/bam/PANC0403-WT-03-RUN2.rep1.refhg38.bam",
        c="data/panc/panc0403/wt03/operations/dir/all.built"
    output:
        "data/panc/panc0403/wt03/bam/PANC0403-WT-03-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

#### PATU8SS89 ####
rule MERGE_sample_runs_patu8ss89_wt01:
    input:
        a="data/panc/patu8ss89/split/wt01r1/bam/PATU8SS89-WT-01-RUN1.rep1.refhg38.bam",
        b="data/panc/patu8ss89/split/wt01r2/bam/PATU8SS89-WT-01-RUN2.rep1.refhg38.bam",
        c="data/panc/patu8ss89/wt01/operations/dir/all.built"
    output:
        "data/panc/patu8ss89/wt01/bam/PATU8SS89-WT-01-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_patu8ss89_wt02:
    input:
        a="data/panc/patu8ss89/split/wt02r1/bam/PATU8SS89-WT-02-RUN1.rep1.refhg38.bam",
        b="data/panc/patu8ss89/split/wt02r2/bam/PATU8SS89-WT-02-RUN2.rep1.refhg38.bam",
        c="data/panc/patu8ss89/wt02/operations/dir/all.built"
    output:
        "data/panc/patu8ss89/wt02/bam/PATU8SS89-WT-02-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_patu8ss89_wt03:
    input:
        a="data/panc/patu8ss89/split/wt03r1/bam/PATU8SS89-WT-03-RUN1.rep1.refhg38.bam",
        b="data/panc/patu8ss89/split/wt03r2/bam/PATU8SS89-WT-03-RUN2.rep1.refhg38.bam",
        c="data/panc/patu8ss89/wt03/operations/dir/all.built"
    output:
        "data/panc/patu8ss89/wt03/bam/PATU8SS89-WT-03-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

#### PK45H ####
rule MERGE_sample_runs_pk45h_wt01:
    input:
        a="data/panc/pk45h/split/wt01r1/bam/PK45H-WT-01-RUN1.rep1.refhg38.bam",
        b="data/panc/pk45h/split/wt01r2/bam/PK45H-WT-01-RUN2.rep1.refhg38.bam",
        c="data/panc/pk45h/wt01/operations/dir/all.built"
    output:
        "data/panc/pk45h/wt01/bam/PK45H-WT-01-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_pk45h_wt02:
    input:
        a="data/panc/pk45h/split/wt02r1/bam/PK45H-WT-02-RUN1.rep1.refhg38.bam",
        b="data/panc/pk45h/split/wt02r2/bam/PK45H-WT-02-RUN2.rep1.refhg38.bam",
        c="data/panc/pk45h/wt02/operations/dir/all.built"
    output:
        "data/panc/pk45h/wt02/bam/PK45H-WT-02-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"

rule MERGE_sample_runs_pk45h_wt03:
    input:
        a="data/panc/pk45h/split/wt03r1/bam/PK45H-WT-03-RUN1.rep1.refhg38.bam",
        b="data/panc/pk45h/split/wt03r2/bam/PK45H-WT-03-RUN2.rep1.refhg38.bam",
        c="data/panc/pk45h/wt03/operations/dir/all.built"
    output:
        "data/panc/pk45h/wt03/bam/PK45H-WT-03-MERGED.rep1.refhg38.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 4
    shell:
        "picard MergeSamFiles \
        I={input.a} \
        I={input.b} \
        O={output} \
        ASSUME_SORTED=TRUE \
        MERGE_SEQUENCE_DICTIONARIES=TRUE \
        USE_THREADING=TRUE"
