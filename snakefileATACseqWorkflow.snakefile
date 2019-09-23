########################################################################################################################################
#### IMPORT MODULES AND CONFIG #########################################################################################################
########################################################################################################################################
configfile: "snakeResources/config/config.yaml"
include: "snakeResources/modules/spoolPreprocessing.snakefile"
include: "snakeResources/modules/spoolFootprinting.snakefile"
include: "snakeResources/modules/spoolFullAnalysis.snakefile"
include: "snakeResources/modules/spoolSampleCorrelation.snakefile"

########################################################################################################################################
#### DIRECTORY STRUCTURE ###############################################################################################################
########################################################################################################################################

rule DIR_main:
    output:
        "{path}operations/main_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}benchmark
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}correlation
        mkdir -p -v {wildcards.path}bam
        mkdir -p -v {wildcards.path}bigwig
        mkdir -p -v {wildcards.path}figures
        mkdir -p -v {wildcards.path}seqbias
        touch {output}
        """

rule DIR_benchmark:
    output:
        "{path}operations/benchmark_dir.built"
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
        "{path}operations/metrics_dir.built"
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
        "{path}operations/preprocessing_dir.built"
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
        "{path}operations/saturation_dir.built"
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
        "{path}operations/footprints_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}footprints/sample_specific
        mkdir -p -v {wildcards.path}footprints/sample_specific/sites
        mkdir -p -v {wildcards.path}footprints/sample_specific/insertion_matrix
        mkdir -p -v {wildcards.path}footprints/sample_specific/raw
        mkdir -p -v {wildcards.path}footprints/sample_specific/aggregated
        mkdir -p -v {wildcards.path}footprints/sample_merged
        mkdir -p -v {wildcards.path}footprints/sample_merged/sites
        mkdir -p -v {wildcards.path}footprints/sample_merged/insertion_matrix
        mkdir -p -v {wildcards.path}footprints/sample_merged/raw
        mkdir -p -v {wildcards.path}footprints/sample_merged/aggregated
        touch {output}
        """

rule DIR_peaks:
    output:
        "{path}operations/peaks_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}peaks/localnorm {wildcards.path}peaks/globalnorm {wildcards.path}peaks/sample_merged
        touch {output}
        """

rule DIR_figures:
    output:
        "{path}operations/figures_dir.built"
    shell:
        """
        mkdir -p -v {wildcards.path}figures/fragmentsizes
        mkdir -p -v {wildcards.path}figures/peakideogram
        mkdir -p -v {wildcards.path}figures/insertionprobability
        mkdir -p -v {wildcards.path}figures/motifalignedheatmap
        mkdir -p -v {wildcards.path}figures/seqbias
        touch {output}
        """

########################################################################################################################################
#### AGGREGATOR RULES ##################################################################################################################
########################################################################################################################################

## This rule determines what is run in the full analysis spooling option
rule AGGREGATOR_full_analysis:
    input:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.preprocessing.complete",
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.footprinting_analysis.complete",
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.seqbias_models.complete"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.full_analysis.finished"
    shell:
        "touch {output}"

## This rule determines what is run in the directory building step
rule AGGREGATOR_build_directory_structure:
    input:
        "{path}operations/main_dir.built",
        "{path}operations/benchmark_dir.built",
        "{path}operations/metrics_dir.built",
        "{path}operations/preprocessing_dir.built",
        "{path}operations/saturation_dir.built",
        "{path}operations/footprints_dir.built",
        "{path}operations/peaks_dir.built",
        "{path}operations/figures_dir.built"
    output:
        "{path}operations/all_dirs.built"
    shell:
        "touch {output}"

## This rule determines what is run in the preprocessing spooling option
rule AGGREGATOR_preprocessing:
    input:
        "{path}operations/all_dirs.built",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "{path}bigwig/{sample}.rep{repnum}.ref{refgenome}.bw",
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.peak_calling.complete",
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.preprocessing_metrics.complete",
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.generate_figures_preprocessing.complete"
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.saturation_analysis.complete"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.preprocessing_preclean.complete"
    priority:
        1
    shell:
        "touch {output}"

## Clean up all the intermediate files just before touching the final flag file
rule AGGREGATOR_preprocessing_clean_intermediate_data:
    input:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.preprocessing_preclean.complete"  
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.preprocessing.complete"
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

## This rule determines what is done for peak calling
rule AGGREGATOR_peak_calling:
    input:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_peaks.narrowPeak",
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p001_peaks.narrowPeak",
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_p001_peaks.narrowPeak",
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_p0001_peaks.narrowPeak",
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_p0001_peaks.narrowPeak",
        "{path}peaks/sample_merged/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.peak_calling.complete"
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
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L1.myco.sam",
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L2.myco.sam",
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L3.myco.sam",
        "{path}preprocessing/4mycoalign/{sample}.rep{repnum}.ref{refgenome}_L4.myco.sam"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.preprocessing_metrics.complete"
    priority:
        1
    shell:
        "touch {output}"

## This rule determines what is run for the footprinting analysis
rule AGGREGATOR_footprinting_analysis:
    input:
        "{path}footprints/sample_specific/aggregated/{sample}.rep{repnum}.ref{refgenome}.sample_specific_aggregated_footprint.RData",
        "{path}footprints/sample_merged/aggregated/{sample}.rep{repnum}.ref{refgenome}.sample_merged_aggregated_footprint.RData"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.footprinting_analysis.complete"
    shell:
        "touch {output}"

## This rule determines what is run for the footprinting analysis
rule AGGREGATOR_seqbias_models:
    input:
        "{path}seqbias/{sample}.rep{repnum}.ref{refgenome}.seqbias_sample_reads_against_sample_peaks.yml",
        "{path}seqbias/{sample}.rep{repnum}.ref{refgenome}.seqbias_sample_reads_against_hg38_chr1.yml"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.seqbias_models.complete"
    shell:
        "touch {output}"

## Aggregator for preprocessing figures
rule AGGREGATOR_generate_figures_preprocessing:
    input:
        "{path}figures/fragmentsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_4.svg",
        "{path}figures/peakideogram/{sample}.rep{repnum}.ref{refgenome}.peakIdeogram.svg"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.generate_figures_preprocessing.complete"
    shell:
        "touch {output}"

## Aggregator for footprinting figures
rule AGGREGATOR_generate_figures_footprinting:
    input:
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.motif_insertion_probability_graphs.complete",
        #"{path}operations/{sample}.rep{repnum}.ref{refgenome}.motif_aligned_heatmaps.complete"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.generate_figures_footprinting.complete"
    shell:
        "touch {output}"

########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################

## Gunzip the fastq files
rule STEP1_gunzip:
    input:
        a="{path}fastq/{sample}.rep{repnum}.ref{refgenome}_L{lane}_R{read}.fastq.gz",
        b="{path}operations/all_dirs.built"
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
        6
    conda:
        "snakeResources/envs/bowtie2.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 24
    shell:
        "bowtie2 -q -p 6 -X2000 -x genomes/{wildcards.refgenome}/{wildcards.refgenome} -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/align/{wildcards.sample}.rep{wildcards.repnum}.ref{wildcards.refgenome}_L{wildcards.lane}.alignment.txt"

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
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 12
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

## Annotate the peaks with global normalization
## Not currently working
rule METRICS_annotate_peaks_global:
    input:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
    output:
        "{path}operations/peakannotation/{sample}.rep{repnum}.ref{refgenome}.globalpeak.annotations.done"
    conda:
        "snakeResources/envs/annotatePeaks.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/metrics/peakannotation/global/{sample}.rep{repnum}.ref{refgenome}.globalpeak.annotations.benchmark.txt'
    script:
        "snakeResources/scripts/annotatePeaks.R"

## Annotate the peaks with local normalization
## not currently working
rule METRICS_annotate_peaks_local:
    input:
        "{path}peaks/localnorm/{sample}.rep{repnum}.ref{refgenome}_localnorm_peaks.narrowPeak",
    output:
        "{path}operations/peakannotation/{sample}.rep{repnum}.ref{refgenome}.localpeak.annotations.done"
    conda:
        "snakeResources/envs/annotatePeaks.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/metrics/peakannotation/local/{sample}.rep{repnum}.ref{refgenome}.localpeak.annotations.benchmark.txt'
    script:
        "snakeResources/scripts/annotatePeaks.R"

########################################################################################################################################
#### SATURATION ANALYSIS ###############################################################################################################
########################################################################################################################################

## This rule determines what is run for the library saturation analysis
rule AGGREGATOR_saturation_analysis:
    input:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.downsampling.done",
        "{path}metrics/saturation/{sample}.rep{repnum}.ref{refgenome}.downsampled_duplication_metrics.txt",
        "{path}metrics/saturation/{sample}.rep{repnum}.ref{refgenome}.downsampled_numpeaks.txt",
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.saturation_footprint_analysis.complete"
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.saturation_analysis.complete"
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
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.downsampling.done"
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
        "{path}saturation/footprints/{sample}.rep{repnum}.ref{refgenome}.{gene}.{prob}.downsampled_raw_footprint.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/analyzeRawFootprint.R"

rule EXPANDER_footprinting_saturation_raw_analysis:
    input:
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.9.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"]),
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.8.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"]),
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.7.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"]),
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.6.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"]),
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.5.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"]),
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.4.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"]),
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.3.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"]),
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.2.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"]),
        expand("{{path}}saturation/footprints/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.1.downsampled_raw_footprint.RData", genename=config["saturationGeneNames"])
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.saturation_footprint_analysis.complete"
    shell:
        "touch {output}"

########################################################################################################################################
#### FOOTPRINTING - SAMPLE SPECIFIC PEAKS ##############################################################################################
########################################################################################################################################

#### 
rule FOOTPRINTING_sample_specific_generate_binding_sites:
    input:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_specific/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_specific_binding_sites.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "snakeResources/scripts/generateBindingSites.R"

#### Generate the insertion matrix for sample specific peaks ####
rule FOOTPRINTING_sample_specific_generate_insertion_matrix:
    input:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "snakeResources/scripts/atacFunctions.R",
        "{path}footprints/sample_specific/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_specific_binding_sites.RData"
    output:
        "{path}footprints/sample_specific/insertion_matrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_specific_insertion_matrix.complete"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/generateInsertionMatrix.R"

#### Raw footprinting analysis, sample-specific peaks ####
rule FOOTPRINTING_sample_specific_raw_analysis:
    input:
        "{path}footprints/sample_specific/insertion_matrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_specific_insertion_matrix.complete",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_specific/raw/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_specific_raw_footprint.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/analyzeRawFootprint.R"

rule EXAPNDER_footprinting_analysis_sample_specific:
    input:
        expand("{{path}}footprints/sample_specific/raw/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.sample_specific_raw_footprint.RData", genename=config["geneNames"])
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.sample_specific_footprint_analysis.complete"
    shell:
        "touch {output}"

rule FOOTPRINTING_sample_specific_aggregate_footprint_data:
    input:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.sample_specific_footprint_analysis.complete",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_specific/aggregated/{sample}.rep{repnum}.ref{refgenome}.sample_specific_aggregated_footprint.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/aggregateFootprintData.R"

########################################################################################################################################
#### FOOTPRINTING - SAMPLE MERGED PEAKS ##############################################################################################
########################################################################################################################################

#### 
rule FOOTPRINTING_sample_merged_generate_binding_sites:
    input:
        "{path}peaks/sample_merged/{sample}.rep{repnum}.ref{refgenome}_sample_merged_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_merged/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_binding_sites.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 1
    script:
        "snakeResources/scripts/generateBindingSites.R"

#### Generate the insertion matrix for sample specific peaks ####
rule FOOTPRINTING_sample_merged_generate_insertion_matrix:
    input:
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai",
        "snakeResources/scripts/atacFunctions.R",
        "{path}footprints/sample_merged/sites/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_binding_sites.RData"
    output:
        "{path}footprints/sample_merged/insertion_matrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_insertion_matrix.complete"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/generateInsertionMatrix.R"

#### Raw footprinting analysis, sample-specific peaks ####
rule FOOTPRINTING_sample_merged_raw_analysis:
    input:
        "{path}footprints/sample_merged/insertion_matrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_insertion_matrix.complete",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_merged/raw/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_raw_footprint.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/analyzeRawFootprint.R"

rule EXAPNDER_footprinting_analysis_sample_merged:
    input:
        expand("{{path}}footprints/sample_merged/raw/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.sample_merged_raw_footprint.RData", genename=config["geneNames"])
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.sample_merged_footprint_analysis.complete"
    shell:
        "touch {output}"

rule FOOTPRINTING_sample_merged_aggregate_footprint_data:
    input:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.sample_merged_footprint_analysis.complete",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}footprints/sample_merged/aggregated/{sample}.rep{repnum}.ref{refgenome}.sample_merged_aggregated_footprint.RData"
    conda:
        "snakeResources/envs/footprintAnalysis.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 50000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/aggregateFootprintData.R"

########################################################################################################################################
#### SEQBIAS MODELS ####################################################################################################################
########################################################################################################################################

rule SEQBIAS_write_sample_peaks_to_BED:
    input:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.bed"
    conda:
        "snakeResources/envs/modelSeqbias.yaml"
    script:
        "snakeResources/scripts/convertPeaksToBEDinterval.R"

rule SEQBIAS_generate_seqbias_model_sample_reads_against_sample_peaks:
    input:
        "fasta/hg38/chr1.fa",
        "{path}peaks/globalnorm/{sample}.rep{repnum}.ref{refgenome}_globalnorm_peaks.bed",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}seqbias/{sample}.rep{repnum}.ref{refgenome}.seqbias_sample_reads_against_sample_peaks.yml",
        "{path}figures/seqbias/{sample}.rep{repnum}.ref{refgenome}.biasedPlot.svg",
        "{path}figures/seqbias/{sample}.rep{repnum}.ref{refgenome}.correctedPlot.svg"
    conda:
        "snakeResources/envs/modelSeqbias.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/generateSeqbiasModel.R"

rule SEQBIAS_generate_seqbias_model_sample_reads_against_hg38_chr1:
    input:
        "fasta/hg38/chr1.fa",
        "snakeResources/bed/hg38.chr1.bed",
        "{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}seqbias/{sample}.rep{repnum}.ref{refgenome}.seqbias_sample_reads_against_hg38_chr1.yml",
        "{path}figures/seqbias/{sample}.rep{repnum}.ref{refgenome}.biasedPlot.svg",
        "{path}figures/seqbias/{sample}.rep{repnum}.ref{refgenome}.correctedPlot.svg"
    conda:
        "snakeResources/envs/modelSeqbias.yaml"
    resources:
        mem_mb=lambda params, attempt: attempt * 30000,
        run_time=lambda params, attempt: attempt * 5
    script:
        "snakeResources/scripts/generateSeqbiasModel.R"

########################################################################################################################################
#### FIGURES ###########################################################################################################################
########################################################################################################################################

#### Generate the fragment size distribution graph
rule FIGURES_fragment_size_distribution:
    input:
        a="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam",
        b="{path}bam/{sample}.rep{repnum}.ref{refgenome}.bam.bai"
    output:
        "{path}figures/fragmentsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_1.svg",
        "{path}figures/fragmentsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_2.svg",
        "{path}figures/fragmentsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_3.svg",
        "{path}figures/fragmentsizes/{sample}.rep{repnum}.ref{refgenome}.fragment_sizes_4.svg"
    conda:
        "snakeResources/envs/fragSizeDistribution.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/metrics/{sample}.rep{repnum}.ref{refgenome}.fragsizes.benchmark.txt'
    script:
        "snakeResources/scripts/generateFragSizeDistribution.R"

#### Generate peak ideogram
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

#### Motif insertion probability graphs
rule FIGURES_generate_motif_insertion_probability_graph_sample_merged_peaks:
    input:
        "{path}footprints/sample_merged/insertion_matrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_insertion_matrix.RData",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}figures/insertionprobability/{sample}.rep{repnum}.ref{refgenome}.{gene}.insertionprobability.svg"
    conda:
        "snakeResources/envs/generateMotifInsertionProbabilityGraph.yaml"
    script:
        "snakeResources/scripts/generateMotifInsertionProbabilityGraph.R"

rule EXAPNDER_motif_insertion_probability_graph_sample_merged_peaks:
    input:
        expand("{{path}}figures/insertionprobability/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.insertionprobability.svg", genename=config["geneNames"])
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.motif_insertion_probability_graphs.complete"
    shell:
        "touch {output}"

#### Motif-aligned heatmaps
rule FIGURES_generate_motif_aligned_heatmap_sample_merged_peaks:
    input:
        "{path}footprints/sample_merged/insertion_matrix/{sample}.rep{repnum}.ref{refgenome}.{gene}.sample_merged_insertion_matrix.RData",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "{path}figures/motifalignedheatmap/{sample}.rep{repnum}.ref{refgenome}.{gene}.motifalignedheatmap.svg"
    conda:
        "snakeResources/envs/generateMotifAlignedHeatmap.yaml"
    script:
        "snakeResources/scripts/generateMotifAlignedHeatmap.R"

rule EXPANDER_motif_aligned_heatmap:
    input:
        expand("{{path}}figures/motifalignedheatmap/{{sample}}.rep{{repnum}}.ref{{refgenome}}.{genename}.motifalignedheatmap.svg", genename=config["geneNames"])
    output:
        "{path}operations/{sample}.rep{repnum}.ref{refgenome}.motif_aligned_heatmaps.complete"
    shell:
        "touch {output}"

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
        "data/pros/lncap/cr01/peaks/sample_merged/LNCaP-CR-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_lncap:
    input:
        "data/pros/lncap/cr01/peaks/sample_merged/LNCaP-CR-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        "data/pros/lncap/cr02/peaks/sample_merged/LNCaP-CR-02.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr04/peaks/sample_merged/LNCaP-CR-04.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr05/peaks/sample_merged/LNCaP-CR-05.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr07/peaks/sample_merged/LNCaP-CR-07.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/cr08/peaks/sample_merged/LNCaP-CR-08.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/wt01/peaks/sample_merged/LNCaP-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak",
        "data/pros/lncap/wt02/peaks/sample_merged/LNCaP-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} "data/pros/lncap/cr02/peaks/sample_merged/LNCaP-CR-02.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/cr04/peaks/sample_merged/LNCaP-CR-04.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/cr05/peaks/sample_merged/LNCaP-CR-05.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/cr07/peaks/sample_merged/LNCaP-CR-07.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/cr08/peaks/sample_merged/LNCaP-CR-08.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/wt01/peaks/sample_merged/LNCaP-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
        cp {input} "data/pros/lncap/wt02/peaks/sample_merged/LNCaP-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak"
        """

#### H508 #######################
rule MERGE_sample_peaks_h508:
    input:
        "data/coad/h508/wt01/peaks/globalnorm/H508-WT-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/h508/wt02/peaks/globalnorm/H508-WT-02.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/h508/wt03/peaks/globalnorm/H508-WT-03.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/h508/wt04/peaks/globalnorm/H508-WT-04.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/h508/wt05/peaks/globalnorm/H508-WT-05.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/h508/wt06/peaks/globalnorm/H508-WT-06.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/coad/h508/wt01/peaks/sample_merged/H508-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_h508:
    input:
        "data/coad/h508/wt01/peaks/sample_merged/H508-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/coad/h508/wt02/peaks/sample_merged/H508-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/coad/h508/wt03/peaks/sample_merged/H508-WT-03.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/coad/h508/wt04/peaks/sample_merged/H508-WT-04.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/coad/h508/wt05/peaks/sample_merged/H508-WT-05.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/coad/h508/wt06/peaks/sample_merged/H508-WT-06.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### hs675.t #######################
rule MERGE_sample_peaks_Hs675T:
    input:
        "data/coad/hs675t/wt01/peaks/globalnorm/Hs675T-WT-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/hs675t/wt02/peaks/globalnorm/Hs675T-WT-02.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/hs675t/wt03/peaks/globalnorm/Hs675T-WT-03.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/coad/hs675t/wt01/peaks/sample_merged/Hs675T-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_Hs675T:
    input:
        "data/coad/hs675t/wt01/peaks/sample_merged/Hs675T-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/coad/hs675t/wt02/peaks/sample_merged/Hs675T-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/coad/hs675t/wt03/peaks/sample_merged/Hs675T-WT-03.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        """

#### LS1034 #######################
rule MERGE_sample_peaks_ls1034:
    input:
        "data/coad/ls1034/wt01/peaks/globalnorm/LS1034-WT-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/ls1034/wt02/peaks/globalnorm/LS1034-WT-02.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/ls1034/wt03/peaks/globalnorm/LS1034-WT-03.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/ls1034/wt04/peaks/globalnorm/LS1034-WT-04.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/ls1034/wt05/peaks/globalnorm/LS1034-WT-05.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/ls1034/wt06/peaks/globalnorm/LS1034-WT-06.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/coad/ls1034/wt01/peaks/sample_merged/LS1034-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_ls1034:
    input:
        "data/coad/ls1034/wt01/peaks/sample_merged/LS1034-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/coad/ls1034/wt02/peaks/sample_merged/LS1034-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/coad/ls1034/wt03/peaks/sample_merged/LS1034-WT-03.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/coad/ls1034/wt04/peaks/sample_merged/LS1034-WT-04.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/coad/ls1034/wt05/peaks/sample_merged/LS1034-WT-05.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/coad/ls1034/wt06/peaks/sample_merged/LS1034-WT-06.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### MDST8 #######################
rule MERGE_sample_peaks_mdst8:
    input:
        "data/coad/mdst8/wt01/peaks/globalnorm/MDST8-WT-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/mdst8/wt02/peaks/globalnorm/MDST8-WT-02.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/coad/mdst8/wt01/peaks/sample_merged/MDST8-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_mdst8:
    input:
        "data/coad/mdst8/wt01/peaks/sample_merged/MDST8-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/coad/mdst8/wt02/peaks/sample_merged/MDST8-WT-02.rep1.refhg38_sample_merged_peaks.narrowPeak",
    shell:
        """
        cp {input} {output.b}
        """

#### SNU16 #######################
rule MERGE_sample_peaks_snu16:
    input:
        "data/coad/snu16/wt01/peaks/globalnorm/SNU16-WT-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/snu16/wt02/peaks/globalnorm/SNU16-WT-01.rep2.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/snu16/wt03/peaks/globalnorm/SNU16-WT-01.rep3.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/snu16/wt04/peaks/globalnorm/SNU16-WT-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/snu16/wt05/peaks/globalnorm/SNU16-WT-01.rep2.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/snu16/wt06/peaks/globalnorm/SNU16-WT-01.rep3.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/coad/snu16/wt01/peaks/sample_merged/SNU16-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_snu16:
    input:
        "data/coad/snu16/wt01/peaks/sample_merged/SNU16-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/coad/snu16/wt02/peaks/sample_merged/SNU16-WT-01.rep2.refhg38_sample_merged_peaks.narrowPeak",
        c="data/coad/snu16/wt03/peaks/sample_merged/SNU16-WT-01.rep3.refhg38_sample_merged_peaks.narrowPeak",
        d="data/coad/snu16/wt04/peaks/sample_merged/SNU16-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/coad/snu16/wt05/peaks/sample_merged/SNU16-WT-01.rep2.refhg38_sample_merged_peaks.narrowPeak",
        f="data/coad/snu16/wt06/peaks/sample_merged/SNU16-WT-01.rep3.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### SNU61 #######################
rule MERGE_sample_peaks_snu61:
    input:
        "data/coad/snu61/wt01/peaks/globalnorm/SNU61-WT-01.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/snu61/wt02/peaks/globalnorm/SNU61-WT-01.rep2.refhg38_globalnorm_peaks.narrowPeak",
        "data/coad/snu61/wt03/peaks/globalnorm/SNU61-WT-01.rep3.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/coad/snu61/wt01/peaks/sample_merged/SNU61-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_snu61:
    input:
        "data/coad/snu61/wt01/peaks/sample_merged/SNU61-WT-01.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/coad/snu61/wt02/peaks/sample_merged/SNU61-WT-01.rep2.refhg38_sample_merged_peaks.narrowPeak",
        c="data/coad/snu61/wt03/peaks/sample_merged/SNU61-WT-01.rep3.refhg38_sample_merged_peaks.narrowPeak",
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        """

#### CAPANI #######################
rule MERGE_sample_peaks_capan1:
    input:
        "data/panc/capan1/split/wt01r1/peaks/globalnorm/CAPANI-WT-01-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/capan1/split/wt02r1/peaks/globalnorm/CAPANI-WT-02-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/capan1/split/wt03r1/peaks/globalnorm/CAPANI-WT-03-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/capan1/split/wt01r2/peaks/globalnorm/CAPANI-WT-01-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/capan1/split/wt02r2/peaks/globalnorm/CAPANI-WT-02-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/capan1/split/wt03r2/peaks/globalnorm/CAPANI-WT-03-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/panc/capan1/split/wt01r1/peaks/sample_merged/CAPANI-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_capan1:
    input:
        "data/panc/capan1/split/wt01r1/peaks/sample_merged/CAPANI-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/panc/capan1/split/wt02r1/peaks/sample_merged/CAPANI-WT-02-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/panc/capan1/split/wt03r1/peaks/sample_merged/CAPANI-WT-03-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/panc/capan1/split/wt01r2/peaks/sample_merged/CAPANI-WT-01-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/panc/capan1/split/wt02r2/peaks/sample_merged/CAPANI-WT-02-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/panc/capan1/split/wt03r2/peaks/sample_merged/CAPANI-WT-03-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### HPAFII #######################
rule MERGE_sample_peaks_hpafii:
    input:
        "data/panc/hpafii/split/wt01r1/peaks/globalnorm/HPAFII-WT-01-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/hpafii/split/wt02r1/peaks/globalnorm/HPAFII-WT-02-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/hpafii/split/wt03r1/peaks/globalnorm/HPAFII-WT-03-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/hpafii/split/wt01r2/peaks/globalnorm/HPAFII-WT-01-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/hpafii/split/wt02r2/peaks/globalnorm/HPAFII-WT-02-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/hpafii/split/wt03r2/peaks/globalnorm/HPAFII-WT-03-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/panc/hpafii/split/wt01r1/peaks/sample_merged/HPAFII-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_hpafii:
    input:
        "data/panc/hpafii/split/wt01r1/peaks/sample_merged/HPAFII-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/panc/hpafii/split/wt02r1/peaks/sample_merged/HPAFII-WT-02-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/panc/hpafii/split/wt03r1/peaks/sample_merged/HPAFII-WT-03-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/panc/hpafii/split/wt01r2/peaks/sample_merged/HPAFII-WT-01-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/panc/hpafii/split/wt02r2/peaks/sample_merged/HPAFII-WT-02-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/panc/hpafii/split/wt03r2/peaks/sample_merged/HPAFII-WT-03-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### KP4 #######################
rule MERGE_sample_peaks_kp4:
    input:
        "data/panc/kp4/split/wt01r1/peaks/globalnorm/KP4-WT-01-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/kp4/split/wt02r1/peaks/globalnorm/KP4-WT-02-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/kp4/split/wt03r1/peaks/globalnorm/KP4-WT-03-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/kp4/split/wt01r2/peaks/globalnorm/KP4-WT-01-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/kp4/split/wt02r2/peaks/globalnorm/KP4-WT-02-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/kp4/split/wt03r2/peaks/globalnorm/KP4-WT-03-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/panc/kp4/split/wt01r1/peaks/sample_merged/KP4-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_kp4:
    input:
        "data/panc/kp4/split/wt01r1/peaks/sample_merged/KP4-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/panc/kp4/split/wt02r1/peaks/sample_merged/KP4-WT-02-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/panc/kp4/split/wt03r1/peaks/sample_merged/KP4-WT-03-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/panc/kp4/split/wt01r2/peaks/sample_merged/KP4-WT-01-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/panc/kp4/split/wt02r2/peaks/sample_merged/KP4-WT-02-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/panc/kp4/split/wt03r2/peaks/sample_merged/KP4-WT-03-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### PANC1 #######################
rule MERGE_sample_peaks_panc1:
    input:
        "data/panc/panc1/split/wt01r1/peaks/globalnorm/PANC1-WT-01-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc1/split/wt02r1/peaks/globalnorm/PANC1-WT-02-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc1/split/wt03r1/peaks/globalnorm/PANC1-WT-03-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc1/split/wt01r2/peaks/globalnorm/PANC1-WT-01-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc1/split/wt02r2/peaks/globalnorm/PANC1-WT-02-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc1/split/wt03r2/peaks/globalnorm/PANC1-WT-03-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/panc/panc1/split/wt01r1/peaks/sample_merged/PANC1-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_panc1:
    input:
        "data/panc/panc1/split/wt01r1/peaks/sample_merged/PANC1-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/panc/panc1/split/wt02r1/peaks/sample_merged/PANC1-WT-02-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/panc/panc1/split/wt03r1/peaks/sample_merged/PANC1-WT-03-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/panc/panc1/split/wt01r2/peaks/sample_merged/PANC1-WT-01-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/panc/panc1/split/wt02r2/peaks/sample_merged/PANC1-WT-02-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/panc/panc1/split/wt03r2/peaks/sample_merged/PANC1-WT-03-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### PANC0403 #######################
rule MERGE_sample_peaks_panc0403:
    input:
        "data/panc/panc0403/split/wt01r1/peaks/globalnorm/PANC0403-WT-01-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc0403/split/wt02r1/peaks/globalnorm/PANC0403-WT-02-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc0403/split/wt03r1/peaks/globalnorm/PANC0403-WT-03-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc0403/split/wt01r2/peaks/globalnorm/PANC0403-WT-01-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc0403/split/wt02r2/peaks/globalnorm/PANC0403-WT-02-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/panc0403/split/wt03r2/peaks/globalnorm/PANC0403-WT-03-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/panc/panc0403/split/wt01r1/peaks/sample_merged/PANC0403-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_panc0403:
    input:
        "data/panc/panc0403/split/wt01r1/peaks/sample_merged/PANC0403-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/panc/panc0403/split/wt02r1/peaks/sample_merged/PANC0403-WT-02-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/panc/panc0403/split/wt03r1/peaks/sample_merged/PANC0403-WT-03-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/panc/panc0403/split/wt01r2/peaks/sample_merged/PANC0403-WT-01-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/panc/panc0403/split/wt02r2/peaks/sample_merged/PANC0403-WT-02-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/panc/panc0403/split/wt03r2/peaks/sample_merged/PANC0403-WT-03-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### PATU8SS89 #######################
rule MERGE_sample_peaks_patu8ss89:
    input:
        "data/panc/patu8ss89/split/wt01r1/peaks/globalnorm/PATU8SS89-WT-01-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/patu8ss89/split/wt02r1/peaks/globalnorm/PATU8SS89-WT-02-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/patu8ss89/split/wt03r1/peaks/globalnorm/PATU8SS89-WT-03-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/patu8ss89/split/wt01r2/peaks/globalnorm/PATU8SS89-WT-01-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/patu8ss89/split/wt02r2/peaks/globalnorm/PATU8SS89-WT-02-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/patu8ss89/split/wt03r2/peaks/globalnorm/PATU8SS89-WT-03-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/panc/patu8ss89/split/wt01r1/peaks/sample_merged/PATU8SS89-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_patu8ss89:
    input:
        "data/panc/patu8ss89/split/wt01r1/peaks/sample_merged/PATU8SS89-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/panc/patu8ss89/split/wt02r1/peaks/sample_merged/PATU8SS89-WT-02-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/panc/patu8ss89/split/wt03r1/peaks/sample_merged/PATU8SS89-WT-03-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/panc/patu8ss89/split/wt01r2/peaks/sample_merged/PATU8SS89-WT-01-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/panc/patu8ss89/split/wt02r2/peaks/sample_merged/PATU8SS89-WT-02-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/panc/patu8ss89/split/wt03r2/peaks/sample_merged/PATU8SS89-WT-03-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

#### PK45H #######################
rule MERGE_sample_peaks_pk45h:
    input:
        "data/panc/pk45h/split/wt01r1/peaks/globalnorm/PK45H-WT-01-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/pk45h/split/wt02r1/peaks/globalnorm/PK45H-WT-02-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/pk45h/split/wt03r1/peaks/globalnorm/PK45H-WT-03-RUN1.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/pk45h/split/wt01r2/peaks/globalnorm/PK45H-WT-01-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/pk45h/split/wt02r2/peaks/globalnorm/PK45H-WT-02-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "data/panc/pk45h/split/wt03r2/peaks/globalnorm/PK45H-WT-03-RUN2.rep1.refhg38_globalnorm_peaks.narrowPeak",
        "snakeResources/scripts/atacFunctions.R"
    output:
        "data/panc/pk45h/split/wt01r1/peaks/sample_merged/PK45H-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    conda:
        "snakeResources/envs/mergeSamplePeaks.yaml"
    script:
        "snakeResources/scripts/mergeSamplePeaks.R"

rule COPY_sample_merged_peaks_pk45h:
    input:
        "data/panc/pk45h/split/wt01r1/peaks/sample_merged/PK45H-WT-01-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak"
    output:
        b="data/panc/pk45h/split/wt02r1/peaks/sample_merged/PK45H-WT-02-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        c="data/panc/pk45h/split/wt03r1/peaks/sample_merged/PK45H-WT-03-RUN1.rep1.refhg38_sample_merged_peaks.narrowPeak",
        d="data/panc/pk45h/split/wt01r2/peaks/sample_merged/PK45H-WT-01-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        e="data/panc/pk45h/split/wt02r2/peaks/sample_merged/PK45H-WT-02-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak",
        f="data/panc/pk45h/split/wt03r2/peaks/sample_merged/PK45H-WT-03-RUN2.rep1.refhg38_sample_merged_peaks.narrowPeak"
    shell:
        """
        cp {input} {output.b}
        cp {input} {output.c}
        cp {input} {output.d}
        cp {input} {output.e}
        cp {input} {output.f}
        """

########################################################################################################################################
#### MERGE SAMPLE RUNS #################################################################################################################
########################################################################################################################################

# #### ####
# rule panc_merge_all_sample_runs:
#     input:
#         "data/panc/capan1/wt01/bam/CAPANI-WT-01-MERGED.rep1.refhg38.bam",
#         "data/panc/capan1/wt02/bam/CAPANI-WT-02-MERGED.rep1.refhg38.bam",
#         "data/panc/capan1/wt03/bam/CAPANI-WT-03-MERGED.rep1.refhg38.bam",
#         "data/panc/hpafii/wt01/bam/HPAFII-WT-01-MERGED.rep1.refhg38.bam",
#         "data/panc/hpafii/wt02/bam/HPAFII-WT-02-MERGED.rep1.refhg38.bam",
#         "data/panc/hpafii/wt03/bam/HPAFII-WT-03-MERGED.rep1.refhg38.bam",
#         "data/panc/kp4/wt01/bam/KP4-WT-01-MERGED.rep1.refhg38.bam",
#         "data/panc/kp4/wt02/bam/KP4-WT-02-MERGED.rep1.refhg38.bam",
#         "data/panc/kp4/wt03/bam/KP4-WT-03-MERGED.rep1.refhg38.bam",
#         "data/panc/panc1/wt01/bam/PANC1-WT-01-MERGED.rep1.refhg38.bam",
#         "data/panc/panc1/wt02/bam/PANC1-WT-02-MERGED.rep1.refhg38.bam",
#         "data/panc/panc1/wt03/bam/PANC1-WT-03-MERGED.rep1.refhg38.bam",
#         "data/panc/panc0403/wt01/bam/PANC0403-WT-01-MERGED.rep1.refhg38.bam",
#         "data/panc/panc0403/wt02/bam/PANC0403-WT-02-MERGED.rep1.refhg38.bam",
#         "data/panc/panc0403/wt03/bam/PANC0403-WT-03-MERGED.rep1.refhg38.bam",
#         "data/panc/patu8ss89/wt01/bam/PATU8SS89-WT-01-MERGED.rep1.refhg38.bam",
#         "data/panc/patu8ss89/wt02/bam/PATU8SS89-WT-02-MERGED.rep1.refhg38.bam",
#         "data/panc/patu8ss89/wt03/bam/PATU8SS89-WT-03-MERGED.rep1.refhg38.bam",
#         "data/panc/pk45h/wt01/bam/PK45H-WT-01-MERGED.rep1.refhg38.bam",
#         "data/panc/pk45h/wt02/bam/PK45H-WT-02-MERGED.rep1.refhg38.bam",
#         "data/panc/pk45h/wt03/bam/PK45H-WT-03-MERGED.rep1.refhg38.bam"

# #### CAPANI ####
# rule MERGE_sample_runs_capan1_wt01:
#     input:
#         a="data/panc/capan1/split/wt01r1/bam/CAPANI-WT-01-RUN1.rep1.refhg38.bam",
#         b="data/panc/capan1/split/wt01r2/bam/CAPANI-WT-01-RUN2.rep1.refhg38.bam",
#         c="data/panc/capan1/wt01/operations/all_dirs.built"
#     output:
#         "data/panc/capan1/wt01/bam/CAPANI-WT-01-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_capan1_wt02:
#     input:
#         a="data/panc/capan1/split/wt02r1/bam/CAPANI-WT-02-RUN1.rep1.refhg38.bam",
#         b="data/panc/capan1/split/wt02r2/bam/CAPANI-WT-02-RUN2.rep1.refhg38.bam",
#         c="data/panc/capan1/wt02/operations/all_dirs.built"
#     output:
#         "data/panc/capan1/wt02/bam/CAPANI-WT-02-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_capan1_wt03:
#     input:
#         a="data/panc/capan1/split/wt03r1/bam/CAPANI-WT-03-RUN1.rep1.refhg38.bam",
#         b="data/panc/capan1/split/wt03r2/bam/CAPANI-WT-03-RUN2.rep1.refhg38.bam",
#         c="data/panc/capan1/wt03/operations/all_dirs.built"
#     output:
#         "data/panc/capan1/wt03/bam/CAPANI-WT-03-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# #### HPAFII ####
# rule MERGE_sample_runs_hpafii_wt01:
#     input:
#         a="data/panc/hpafii/split/wt01r1/bam/HPAFII-WT-01-RUN1.rep1.refhg38.bam",
#         b="data/panc/hpafii/split/wt01r2/bam/HPAFII-WT-01-RUN2.rep1.refhg38.bam",
#         c="data/panc/hpafii/wt01/operations/all_dirs.built"
#     output:
#         "data/panc/hpafii/wt01/bam/HPAFII-WT-01-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_hpafii_wt02:
#     input:
#         a="data/panc/hpafii/split/wt02r1/bam/HPAFII-WT-02-RUN1.rep1.refhg38.bam",
#         b="data/panc/hpafii/split/wt02r2/bam/HPAFII-WT-02-RUN2.rep1.refhg38.bam",
#          c="data/panc/hpafii/wt02/operations/all_dirs.built"
#     output:
#         "data/panc/hpafii/wt02/bam/HPAFII-WT-02-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_hpafii_wt03:
#     input:
#         a="data/panc/hpafii/split/wt03r1/bam/HPAFII-WT-03-RUN1.rep1.refhg38.bam",
#         b="data/panc/hpafii/split/wt03r2/bam/HPAFII-WT-03-RUN2.rep1.refhg38.bam",
#          c="data/panc/hpafii/wt03/operations/all_dirs.built"
#     output:
#         "data/panc/hpafii/wt03/bam/HPAFII-WT-03-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# #### KP4 ####
# rule MERGE_sample_runs_kp4_wt01:
#     input:
#         a="data/panc/kp4/split/wt01r1/bam/KP4-WT-01-RUN1.rep1.refhg38.bam",
#         b="data/panc/kp4/split/wt01r2/bam/KP4-WT-01-RUN2.rep1.refhg38.bam",
#         c="data/panc/kp4/wt01/operations/all_dirs.built"
#     output:
#         "data/panc/kp4/wt01/bam/KP4-WT-01-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_kp4_wt02:
#     input:
#         a="data/panc/kp4/split/wt02r1/bam/KP4-WT-02-RUN1.rep1.refhg38.bam",
#         b="data/panc/kp4/split/wt02r2/bam/KP4-WT-02-RUN2.rep1.refhg38.bam",
#         c="data/panc/kp4/wt02/operations/all_dirs.built"
#     output:
#         "data/panc/kp4/wt02/bam/KP4-WT-02-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_kp4_wt03:
#     input:
#         a="data/panc/kp4/split/wt03r1/bam/KP4-WT-03-RUN1.rep1.refhg38.bam",
#         b="data/panc/kp4/split/wt03r2/bam/KP4-WT-03-RUN2.rep1.refhg38.bam",
#         c="data/panc/kp4/wt03/operations/all_dirs.built"
#     output:
#         "data/panc/kp4/wt03/bam/KP4-WT-03-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# #### PANC1 ####
# rule MERGE_sample_runs_panc1_wt01:
#     input:
#         a="data/panc/panc1/split/wt01r1/bam/PANC1-WT-01-RUN1.rep1.refhg38.bam",
#         b="data/panc/panc1/split/wt01r2/bam/PANC1-WT-01-RUN2.rep1.refhg38.bam",
#         c="data/panc/panc1/wt01/operations/all_dirs.built"
#     output:
#         "data/panc/panc1/wt01/bam/PANC1-WT-01-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_panc1_wt02:
#     input:
#         a="data/panc/panc1/split/wt02r1/bam/PANC1-WT-02-RUN1.rep1.refhg38.bam",
#         b="data/panc/panc1/split/wt02r2/bam/PANC1-WT-02-RUN2.rep1.refhg38.bam",
#         c="data/panc/panc1/wt02/operations/all_dirs.built"
#     output:
#         "data/panc/panc1/wt02/bam/PANC1-WT-02-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_panc1_wt03:
#     input:
#         a="data/panc/panc1/split/wt03r1/bam/PANC1-WT-03-RUN1.rep1.refhg38.bam",
#         b="data/panc/panc1/split/wt03r2/bam/PANC1-WT-03-RUN2.rep1.refhg38.bam",
#         c="data/panc/panc1/wt03/operations/all_dirs.built"
#     output:
#         "data/panc/panc1/wt03/bam/PANC1-WT-03-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# #### PANC0403 ####
# rule MERGE_sample_runs_panc0403_wt01:
#     input:
#         a="data/panc/panc0403/split/wt01r1/bam/PANC0403-WT-01-RUN1.rep1.refhg38.bam",
#         b="data/panc/panc0403/split/wt01r2/bam/PANC0403-WT-01-RUN2.rep1.refhg38.bam",
#         c="data/panc/panc0403/wt01/operations/all_dirs.built"
#     output:
#         "data/panc/panc0403/wt01/bam/PANC0403-WT-01-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_panc0403_wt02:
#     input:
#         a="data/panc/panc0403/split/wt02r1/bam/PANC0403-WT-02-RUN1.rep1.refhg38.bam",
#         b="data/panc/panc0403/split/wt02r2/bam/PANC0403-WT-02-RUN2.rep1.refhg38.bam",
#         c="data/panc/panc0403/wt02/operations/all_dirs.built"
#     output:
#         "data/panc/panc0403/wt02/bam/PANC0403-WT-02-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_panc0403_wt03:
#     input:
#         a="data/panc/panc0403/split/wt03r1/bam/PANC0403-WT-03-RUN1.rep1.refhg38.bam",
#         b="data/panc/panc0403/split/wt03r2/bam/PANC0403-WT-03-RUN2.rep1.refhg38.bam",
#         c="data/panc/panc0403/wt03/operations/all_dirs.built"
#     output:
#         "data/panc/panc0403/wt03/bam/PANC0403-WT-03-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# #### PATU8SS89 ####
# rule MERGE_sample_runs_patu8ss89_wt01:
#     input:
#         a="data/panc/patu8ss89/split/wt01r1/bam/PATU8SS89-WT-01-RUN1.rep1.refhg38.bam",
#         b="data/panc/patu8ss89/split/wt01r2/bam/PATU8SS89-WT-01-RUN2.rep1.refhg38.bam",
#         c="data/panc/patu8ss89/wt01/operations/all_dirs.built"
#     output:
#         "data/panc/patu8ss89/wt01/bam/PATU8SS89-WT-01-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_patu8ss89_wt02:
#     input:
#         a="data/panc/patu8ss89/split/wt02r1/bam/PATU8SS89-WT-02-RUN1.rep1.refhg38.bam",
#         b="data/panc/patu8ss89/split/wt02r2/bam/PATU8SS89-WT-02-RUN2.rep1.refhg38.bam",
#         c="data/panc/patu8ss89/wt02/operations/all_dirs.built"
#     output:
#         "data/panc/patu8ss89/wt02/bam/PATU8SS89-WT-02-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_patu8ss89_wt03:
#     input:
#         a="data/panc/patu8ss89/split/wt03r1/bam/PATU8SS89-WT-03-RUN1.rep1.refhg38.bam",
#         b="data/panc/patu8ss89/split/wt03r2/bam/PATU8SS89-WT-03-RUN2.rep1.refhg38.bam",
#         c="data/panc/patu8ss89/wt03/operations/all_dirs.built"
#     output:
#         "data/panc/patu8ss89/wt03/bam/PATU8SS89-WT-03-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# #### PK45H ####
# rule MERGE_sample_runs_pk45h_wt01:
#     input:
#         a="data/panc/pk45h/split/wt01r1/bam/PK45H-WT-01-RUN1.rep1.refhg38.bam",
#         b="data/panc/pk45h/split/wt01r2/bam/PK45H-WT-01-RUN2.rep1.refhg38.bam",
#         c="data/panc/pk45h/wt01/operations/all_dirs.built"
#     output:
#         "data/panc/pk45h/wt01/bam/PK45H-WT-01-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_pk45h_wt02:
#     input:
#         a="data/panc/pk45h/split/wt02r1/bam/PK45H-WT-02-RUN1.rep1.refhg38.bam",
#         b="data/panc/pk45h/split/wt02r2/bam/PK45H-WT-02-RUN2.rep1.refhg38.bam",
#         c="data/panc/pk45h/wt02/operations/all_dirs.built"
#     output:
#         "data/panc/pk45h/wt02/bam/PK45H-WT-02-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"

# rule MERGE_sample_runs_pk45h_wt03:
#     input:
#         a="data/panc/pk45h/split/wt03r1/bam/PK45H-WT-03-RUN1.rep1.refhg38.bam",
#         b="data/panc/pk45h/split/wt03r2/bam/PK45H-WT-03-RUN2.rep1.refhg38.bam",
#         c="data/panc/pk45h/wt03/operations/all_dirs.built"
#     output:
#         "data/panc/pk45h/wt03/bam/PK45H-WT-03-MERGED.rep1.refhg38.bam"
#     conda:
#         "snakeResources/envs/picard.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda params, attempt: attempt * 30000,
#         run_time=lambda params, attempt: attempt * 4
#     shell:
#         "picard MergeSamFiles \
#         I={input.a} \
#         I={input.b} \
#         O={output} \
#         ASSUME_SORTED=TRUE \
#         MERGE_SEQUENCE_DICTIONARIES=TRUE \
#         USE_THREADING=TRUE"
