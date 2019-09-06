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
        "{path}operations/dir/all.built",
        "{path}operations/preprocessing/{sample}-REP{repnum}.preprocessing.complete",
        "{path}operations/footprints/{sample}-REP{repnum}.footprinting_raw_analysis.complete",
    output:
        "{path}operations/modules/{sample}-REP{repnum}.full_analysis.finished"
    shell:
        "touch {output}"

## This rule determines what is run in the directory building step
rule AGGREGATOR_builddirectories:
    input:
        "{path}operations/dir/main.built",
        "{path}operations/dir/operations.built",
        "{path}operations/dir/benchmark.built",
        "{path}operations/dir/metrics.built",
        "{path}operations/dir/preprocessing.built",
        "{path}operations/dir/footprints.built",
        "{path}operations/dir/peaks.built"
    output:
        "{path}operations/dir/all.built"
    shell:
        "touch {output}"

## This rule determines what is run in the preprocessing spooling option
rule AGGREGATOR_preprocessing:
    input:
        "{path}operations/dir/all.built",
        "{path}bam/{sample}-REP{repnum}.bai",
        "{path}bigwig/{sample}-REP{repnum}.bw",
        "{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak",
        "{path}peaks/localnorm/{sample}-REP{repnum}_localnorm_peaks.narrowPeak",
        "{path}metrics/genomecov/{sample}-REP{repnum}.peak.globalnorm.genomecov.txt",
        "{path}metrics/genomecov/{sample}-REP{repnum}.peak.localnorm.genomecov.txt",
        "{path}operations/metrics/{sample}-REP{repnum}.globalpeak.annotations.done",
        "{path}operations/metrics/{sample}-REP{repnum}.localpeak.annotations.done",
        "{path}metrics/totalreads/{sample}-REP{repnum}.totalreads.Rdata",
        "{path}operations/metrics/{sample}-REP{repnum}.fragsizes.done"
    output:
        "{path}operations/preprocessing/{sample}-REP{repnum}.preprocessing.complete"
    shell:
        "touch {output}"

## This rule initiates the raw footprint analysis for all genes found in the config file
rule AGGREGATOR_footprinting_raw_analysis:
    input:
        "{path}operations/preprocessing/{sample}-REP{repnum}.preprocessing.complete",
        expand("{{path}}operations/footprints/raw/{{sample}}-REP{{repnum}}.{genename}.rawFPanalysis.done", genename=config["geneNames"])
    output:
        "{path}operations/footprints/{sample}-REP{repnum}.footprinting_raw_analysis.complete"
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
        mkdir -p -v {wildcards.path}operations/saturation
        mkdir -p -v {wildcards.path}operations/metrics
        touch {output}
        """

rule DIR_benchmark:
    output:
        "{path}operations/dir/benchmark.built"
    shell:
        """
        mkdir -p -v {wildcards.path}benchmark/preprocessing
        mkdir -p -v {wildcards.path}benchmark/preprocessing/gunzip {wildcards.path}benchmark/preprocessing/fastp {wildcards.path}benchmark/preprocessing/mycoalign
        mkdir -p -v {wildcards.path}benchmark/preprocessing/hg38align {wildcards.path}benchmark/preprocessing/coordsortsam {wildcards.path}benchmark/preprocessing/bamconversion
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
        mkdir -p -v {wildcards.path}metrics/hg38
        mkdir -p -v {wildcards.path}metrics/genomecov
        mkdir -p -v {wildcards.path}metrics/duplication
        mkdir -p -v {wildcards.path}metrics/totalreads
        mkdir -p -v {wildcards.path}metrics/fragsize
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
        mkdir -p -v {wildcards.path}preprocessing/5hg38align
        mkdir -p -v {wildcards.path}preprocessing/6raw
        mkdir -p -v {wildcards.path}preprocessing/6raw/mitochondrial {wildcards.path}preprocessing/6raw/blacklist {wildcards.path}preprocessing/6raw/nonblacklist
        mkdir -p -v {wildcards.path}preprocessing/7rgsort
        mkdir -p -v {wildcards.path}preprocessing/8merged
        mkdir -p -v {wildcards.path}preprocessing/9dedup
        mkdir -p -v {wildcards.path}preprocessing/10unique
        mkdir -p -v {wildcards.path}preprocessing/11bigwig
        mkdir -p -v {wildcards.path}preprocessing/saturation
        mkdir -p -v {wildcards.path}preprocessing/saturation/complexity
        mkdir -p -v {wildcards.path}preprocessing/saturation/peaks 
        mkdir -p -v {wildcards.path}preprocessing/saturation/downsampled
        mkdir -p -v {wildcards.path}preprocessing/saturation/downsampled
        mkdir -p -v {wildcards.path}preprocessing/saturation/downsampled/raw {wildcards.path}preprocessing/saturation/downsampled/cs {wildcards.path}preprocessing/saturation/downsampled/md
        mkdir -p -v {wildcards.path}preprocessing/saturation/footprints
        mkdir -p -v {wildcards.path}preprocessing/saturation/footprints/raw {wildcards.path}preprocessing/saturation/footprints/parsed {wildcards.path}preprocessing/saturation/footprints/processed
        touch {output}
        """

rule DIR_footprints:
    output:
        "{path}operations/dir/footprints.built"
    shell:
        """
        mkdir -p -v {wildcards.path}footprints/data 
        mkdir -p -v {wildcards.path}footprints/data/raw {wildcards.path}footprints/data/parsed 
        mkdir -p -v {wildcards.path}footprints/data/processed {wildcards.path}footprints/data/aggregated
        mkdir -p -v {wildcards.path}footprints/data/temp
        mkdir -p -v {wildcards.path}footprints/graphs
        mkdir -p -v {wildcards.path}footprints/graphs/insprob {wildcards.path}footprints/graphs/heatmaps
        touch {output}
        """

rule DIR_peaks:
    output:
        "{path}operations/dir/peaks.built"
    shell:
        """
        mkdir -p -v {wildcards.path}peaks/localnorm {wildcards.path}peaks/globalnorm
        touch {output}
        """

########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################
## Gunzip the fastq files
rule STEP1_gunzip:
    input:
        a="{path}preprocessing/1gz/{sample}-REP{repnum}_L{lane}_R{read}.fastq.gz",
        b="{path}operations/dir/all.built"
    output:
        c="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R{read}.fastq"
    threads:
        1
    benchmark:
        '{path}benchmark/preprocessing/gunzip/{sample}-REP{repnum}_L{lane}_R{read}.gunzip.benchmark.txt'
    shell:
        "gunzip -k -c {input.a} > {output.c}"

## Fastq QC filtering
rule STEP2_fastp_filtering:
    input:
        a="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R1.fastq",
        b="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R2.fastq"
    output:
        c="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        d="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq"
    benchmark:
        '{path}benchmark/preprocessing/fastp/{sample}-REP{repnum}.{lane}.fastp.benchmark.txt'
    threads:
        2
    conda:
    	"snakeResources/envs/fastp.yaml"
    shell:
        "fastp -i {input.a} -I {input.b} -o {output.c} -O {output.d} -w 2 -h {wildcards.path}metrics/fastq/{wildcards.sample}-REP{wildcards.repnum}_L{wildcards.lane}.fastq.quality.html -j {wildcards.path}metrics/fastq/{wildcards.sample}-REP{wildcards.repnum}_L{wildcards.lane}.fastq.quality.json"
    
## Check for mycoplasma contamination
rule STEP3_mycoalign:
    input:
        a="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq"
    output:
        "{path}preprocessing/4mycoalign/{sample}-REP{repnum}_L{lane}.myco.sam"
    benchmark:
        '{path}benchmark/preprocessing/mycoalign/{sample}-REP{repnum}.{lane}.mycoalign.benchmark.txt'
    threads:
        10
    conda:
        "snakeResources/envs/bowtie2.yaml"
    resources:
        mem_mb=10000
    shell:
        "bowtie2 -q -p 10 -X2000 -x genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/myco/{wildcards.sample}-REP{wildcards.repnum}_L{wildcards.lane}.myco.alignment.txt"
    
## Align reads to human hg38 build
rule STEP4_hg38align:
    input:
        a="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq",
        c="{path}preprocessing/4mycoalign/{sample}-REP{repnum}_L{lane}.myco.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.sam"
    benchmark:
        '{path}benchmark/preprocessing/hg38align/{sample}-REP{repnum}.{lane}.hg38align.benchmark.txt'
    threads:
        10
    conda:
        "snakeResources/envs/bowtie2.yaml"
    resources:
        mem_mb=40000
    shell:
        "bowtie2 -q -p 10 -X2000 -x genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/hg38/{wildcards.sample}-REP{wildcards.repnum}_L{wildcards.lane}.hg38.alignment.txt"
    
## Coordinate sort the aligned reads. This is required for blacklist filtering
rule STEP5_coordsort_sam:
    input:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.cs.sam"
    threads:
        1
    conda:
        "snakeResources/envs/samtools.yaml"
    benchmark:
        '{path}benchmark/preprocessing/coordsortsam/{sample}-REP{repnum}.{lane}.coordsort.benchmark.txt'
    shell:
        "samtools sort {input} -o {output} -O sam"
    
## Remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
rule STEP6_blacklistfilter_bamconversion:
    input:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.cs.sam"
    output:
        a="{path}preprocessing/6raw/blacklist/{sample}-REP{repnum}_L{lane}.hg38blacklist.bam",
        b="{path}preprocessing/6raw/nonblacklist/{sample}-REP{repnum}_L{lane}.blrm.bam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        5
    benchmark:
        '{path}benchmark/preprocessing/bamconversion/{sample}-REP{repnum}.{lane}.bamconvert.benchmark.txt'
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 4 {input}"
    
## Remove reads mapping to mitochondrial DNA
rule STEP7_chrM_contamination:
    input:
        "{path}preprocessing/6raw/nonblacklist/{sample}-REP{repnum}_L{lane}.blrm.bam"
    output:
        a="{path}preprocessing/6raw/mitochondrial/{sample}-REP{repnum}_L{lane}.mitochondrial.bam",
        b="{path}preprocessing/6raw/{sample}-REP{repnum}_L{lane}.goodbam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        5
    benchmark:
        '{path}benchmark/preprocessing/removemitochondrial/{sample}-REP{repnum}.{lane}.chrMfilter.benchmark.txt'
    shell:
        "samtools view -b -h -o {output.a} -L genomes/hg38/hg38.blacklist.bed -U {output.b} -@ 4 {input}"
    
## Add @RG tags to the reads and perform coordinate sorting
rule STEP8_addrgandcsbam:
    input:
        "{path}preprocessing/6raw/{sample}-REP{repnum}_L{lane}.goodbam"
    output:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.rg.cs.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=10000
    benchmark:
        '{path}benchmark/preprocessing/addRG/{sample}-REP{repnum}.{lane}.addRGtag.benchmark.txt'
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
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.rg.cs.bam"
    output:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.clean.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=10000
    benchmark:
        '{path}benchmark/preprocessing/cleansam/{sample}-REP{repnum}.{lane}.cleansam.benchmark.txt'
    shell:
        "picard CleanSam \
        I={input} \
        O={output}"
    
## Merge reads from different NextSeq lanes
rule STEP10_mergelanes:
    input:
        a="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L1.clean.bam",
        b="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L2.clean.bam",
        c="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L3.clean.bam",
        d="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L4.clean.bam"
    output:
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        2
    resources:
        mem_mb=10000
    benchmark:
        '{path}benchmark/preprocessing/mergelanes/{sample}-REP{repnum}.mergelanes.benchmark.txt'
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

## After the lanes are merged, you can remove the intermediate data 
rule STEP10b_clean_intermediate_data:
    input:
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    output:
        "{path}operations/preprocessing/clean10b.{sample}.{repnum}.done"
    threads:
        1
    shell:
        """
        rm -f {wildcards.path}preprocessing/2fastq/*REP{wildcards.repnum}*.fastq
        rm -f {wildcards.path}preprocessing/3goodfastq/*REP{wildcards.repnum}*.fq
        rm -f {wildcards.path}preprocessing/4mycoalign/*REP{wildcards.repnum}*.sam
        rm -f {wildcards.path}preprocessing/5hg38align/*REP{wildcards.repnum}*.sam
        rm -f {wildcards.path}preprocessing/6raw/*REP{wildcards.repnum}*.goodbam
        rm -f {wildcards.path}preprocessing/6raw/blacklist/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/6raw/mitochondrial/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/6raw/nonblacklist/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/7rgsort/*REP{wildcards.repnum}*.bam
        #
        touch {output}
        """

## Purge PCR duplicate reads
rule STEP11_purgeduplicates:
    input:
        a="{path}operations/preprocessing/clean10b.{sample}.{repnum}.done",
        b="{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    output:
        "{path}preprocessing/9dedup/{sample}-REP{repnum}.dp.bam"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/preprocessing/purgeduplicates/{sample}-REP{repnum}.purgeduplicates.benchmark.txt'
    resources:
        mem_mb=25000
    shell:
        "picard MarkDuplicates \
        I={input.b} \
        O={output} \
        M={wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}.duplication.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"
    
## Filter reads for only uniquely mapping
rule STEP12_mapqfilter:
    input:
        "{path}preprocessing/9dedup/{sample}-REP{repnum}.dp.bam"
    output:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam"
    conda:
        "snakeResources/envs/samtools.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/preprocessing/mapqfilter/{sample}-REP{repnum}.mapqfilter.benchmark.txt'
    shell:
        "samtools view -h -q 2 -b {input} > {output}"
    
## Remove intermediate data to this point and copy the bam file out to the parent bam folder
rule STEP12b_clean_intermediate_data:
    input:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam"
    output:
        "{path}bam/{sample}-REP{repnum}.bam"
    shell:
        """
        rm -f {wildcards.path}preprocessing/8merged/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/9dedup/*REP{wildcards.repnum}*.bam
        mv {wildcards.path}preprocessing/10unique/*REP{wildcards.repnum}*.bam {wildcards.path}bam/{wildcards.sample}-REP{wildcards.repnum}.bam
        touch {output}
        """

## Build the .bai index for the processed bam file
rule STEP13_build_index:
    input:
        "{path}bam/{sample}-REP{repnum}.bam"
    output:
        "{path}bam/{sample}-REP{repnum}.bai"
    conda:
        "snakeResources/envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=10000
    benchmark:
        '{path}benchmark/preprocessing/buildindex/{sample}-REP{repnum}.buildindex.benchmark.txt'
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"
    
## Make a bigwig file from the bam file
rule STEP14_makebigwig_bamcov:
    input:
        a="{path}bam/{sample}-REP{repnum}.bam",
        b="{path}bam/{sample}-REP{repnum}.bai"
    output:
        "{path}bigwig/{sample}-REP{repnum}.bw"
    conda:
        "snakeResources/envs/deeptools.yaml"
    threads:
        10
    resources:
        mem_mb=10000
    benchmark:
        '{path}benchmark/preprocessing/bigwig/{sample}-REP{repnum}.makebigwig.benchmark.txt'
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 10 -v"
    
## Call peaks using global normalization
rule STEP15_MACS2_peaks_global_normilization:
    input:
        a="{path}bam/{sample}-REP{repnum}.bam",
        b="{path}bam/{sample}-REP{repnum}.bai"
    output:
        "{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    threads:
        4
    benchmark:
        '{path}benchmark/preprocessing/peaks/{sample}-REP{repnum}.callpeaks.globalnorm.benchmark.txt'
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}_globalnorm --outdir {wildcards.path}peaks/globalnorm --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"
    
## Call peaks using local normalization
rule STEP16_MACS2_peaks_local_normalization:
    input:
        a="{path}bam/{sample}-REP{repnum}.bam",
        b="{path}bam/{sample}-REP{repnum}.bai"
    output:
        "{path}peaks/localnorm/{sample}-REP{repnum}_localnorm_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    threads:
        4
    benchmark:
        '{path}benchmark/preprocessing/peaks/{sample}-REP{repnum}.callpeaks.localnorm.benchmark.txt'
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}_localnorm --outdir {wildcards.path}peaks/localnorm --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

########################################################################################################################################
#### QC METRICS  #######################################################################################################################
########################################################################################################################################
    
## Calculate percent genome coverage from peaks with global normalization
rule METRICS_percent_peak_genome_coverage_globalnorm:
    input:
        a="{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/genomecov/{sample}-REP{repnum}.peak.globalnorm.genomecov.txt"
    conda:
        "snakeResources/envs/bedops.yaml"
    threads:
        1  	
    benchmark:
        '{path}benchmark/metrics/{sample}-REP{repnum}.genomecov.globalnorm.benchmark.txt'
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"
    
## Calculate percent genome coverage from peaks with local normalization
rule METRICS_percent_peak_genome_coverage_localnorm:
    input:
        a="{path}peaks/localnorm/{sample}-REP{repnum}_localnorm_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/genomecov/{sample}-REP{repnum}.peak.localnorm.genomecov.txt"
    conda:
        "snakeResources/envs/bedops.yaml"
    threads:
        1  
    benchmark:
        '{path}benchmark/metrics/{sample}-REP{repnum}.genomecov.localnorm.benchmark.txt'
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"
    
## Generate the fragment size distribution graph
rule METRICS_fragment_size_distribution:
    input:
        a="{path}bam/{sample}-REP{repnum}.bam",
        b="{path}bam/{sample}-REP{repnum}.bai"
    output:
        "{path}operations/metrics/{sample}-REP{repnum}.fragsizes.done"
    conda:
        "snakeResources/envs/Rfragsizedistribution.yaml"
    threads:
        1
    resources:
        mem_mb=20000
    benchmark:
        '{path}benchmark/metrics/fragsize/{sample}-REP{repnum}.fragsizes.benchmark.txt'
    script:
        "snakeResources/scripts/generateFragSizeDistribution.R"

## Annotate the peaks with global normalization
rule METRICS_annotate_peaks_global:
    input:
        "{path}peaks/globalnorm/{sample}-REP{repnum}_globalnorm_peaks.narrowPeak",
    output:
        "{path}operations/metrics/{sample}-REP{repnum}.globalpeak.annotations.done"
    conda:
        "snakeResources/envs/Rannotatepeaks.yaml"
    threads:
        1
    resources:
        mem_mb=10000
    benchmark:
        '{path}benchmark/metrics/genomecov/{sample}-REP{repnum}.globalpeak.annotations.benchmark.txt'
    script:
        "snakeResources/scripts/annotatePeaks.R"

## Annotate the peaks with local normalization
rule METRICS_annotate_peaks_local:
    input:
        "{path}peaks/localnorm/{sample}-REP{repnum}_localnorm_peaks.narrowPeak",
    output:
        "{path}operations/metrics/{sample}-REP{repnum}.localpeak.annotations.done"
    conda:
        "snakeResources/envs/Rannotatepeaks.yaml"
    threads:
        1
    resources:
        mem_mb=10000
    benchmark:
        '{path}benchmark/metrics/genomecov/{sample}-REP{repnum}.localpeak.annotations.benchmark.txt'
    script:
        "snakeResources/scripts/annotatePeaks.R"

## Count the total number of reads in the sample
rule METRICS_sample_total_reads:
    input:
        a="{path}bam/{sample}-REP{repnum}.bam",
        b="{path}bam/{sample}-REP{repnum}.bai"
    output:
        "{path}metrics/totalreads/{sample}-REP{repnum}.totalreads.Rdata"
    conda:
        "snakeResources/envs/Rcountsamplereads.yaml"
    threads:
        1
    benchmark:
        '{path}benchmark/metrics/totalreads/{sample}-REP{repnum}.totalreads.benchmark.txt'
    script:
        "snakeResources/scripts/countTotalSampleReads.R"

########################################################################################################################################
#### SATURATION ANALYSIS ###############################################################################################################
########################################################################################################################################

## Analyze saturation of library in terms of library complexity, called peaks, and selected TF footprints
rule STEP21_saturation_analysis:
    input:
        "{path}operations/saturation/{sample}-REP{repnum}.downsample.done",
        "{path}operations/preprocessing/saturation/clean1.{sample}.{repnum}.done",
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_globalnorm_numpeaks.txt",
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_localnorm_numpeaks.txt",
        "{path}operations/saturation/footprints/{sample}-REP{repnum}.allgenes.footprint.downsampled.done"
    output:
        "{path}operations/saturation/{sample}-REP{repnum}.saturation_analysis.aggregator"
    shell:
    	"touch {output}"

## When all saturation analysis is done, remove the intermediate data, as it is quite large
rule SATURATION_clean_data_final:
    input:
        "{path}operations/saturation/{sample}-REP{repnum}.saturation_analysis.aggregator"
    output:
        "{path}operations/saturation/{sample}-REP{repnum}.saturation_analysis.done"
    shell:
        """
        rm -f {wildcards.path}preprocessing/saturation/downsampled/md/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/md/*REP{wildcards.repnum}*.bai
        rm -f {wildcards.path}preprocessing/saturation/peaks/*REP{wildcards.repnum}*
        touch {output}
        """

## This rule determines which genes will be analyzed for the footprinting saturation analysis
rule AGGREGATOR_saturation_footprints_genes:
    input:
        "{path}operations/saturation/footprints/{sample}-REP{repnum}.CTCF.footprint.downsampled.done"
    output:
        "{path}operations/saturation/footprints/{sample}-REP{repnum}.allgenes.footprint.downsampled.done"
    shell:
        """
        rm -f {wildcards.path}preprocessing/8merged/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/raw/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/cs/*REP{wildcards.repnum}*.bam
        touch {output}
        """

## Determines the downsampling levels of the libraries
rule AGGREGATOR_saturation_downsample:
    input:
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.9.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.8.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.7.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.6.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.5.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.4.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.3.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.2.md.bai",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.1.md.bai"      
    output:
        "{path}operations/saturation/{sample}-REP{repnum}.downsample.done"
    shell:
        "touch {output}"

## Downsample the processed but NOT duplicate purged .bam files
rule SATURATION_downsample_bam:
    input:
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    output:
        "{path}preprocessing/saturation/downsampled/raw/{sample}-REP{repnum}.{prob}.bam"
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.downsample.bam.benchmark.txt'
    shell:
        "picard DownsampleSam \
        I={input} \
        O={output} \
        PROBABILITY=0.{wildcards.prob}"

## Coordinate sort the downsampled .bam files
rule SATURATION_sort_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/raw/{sample}-REP{repnum}.{prob}.bam"
    output:
        "{path}preprocessing/saturation/downsampled/cs/{sample}-REP{repnum}.{prob}.cs.bam"
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.sort.downsampled.benchmark.txt'
    shell:
        "picard SortSam \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate"

## Purge duplicates from the downsampled .bam files
rule SATURATION_purge_duplicates_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/cs/{sample}-REP{repnum}.{prob}.cs.bam"
    output:
        a="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam",
        b="{path}metrics/saturation/{sample}-REP{repnum}.{prob}.duplication-metrics.txt"
    resources:
        purgeDuplicates=1,
        mem_mb=50000
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.purgeduplicates.downsampled.benchmark.txt'
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
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam"
    output:
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bai"
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{prob}.index.downsampled.benchmark.txt'
    shell:
        "picard BuildBamIndex \
        I={input} \
        O={output}"

## Determine the library complexity of the downsampled libraries and output to metrics
rule SATURATION_analyze_complexity_downsampled:
    input:
        a="{path}metrics/saturation/{sample}-REP{repnum}.9.duplication-metrics.txt",
        b="{path}metrics/saturation/{sample}-REP{repnum}.8.duplication-metrics.txt",
        c="{path}metrics/saturation/{sample}-REP{repnum}.7.duplication-metrics.txt",
        d="{path}metrics/saturation/{sample}-REP{repnum}.6.duplication-metrics.txt",
        e="{path}metrics/saturation/{sample}-REP{repnum}.5.duplication-metrics.txt",
        f="{path}metrics/saturation/{sample}-REP{repnum}.4.duplication-metrics.txt",
        g="{path}metrics/saturation/{sample}-REP{repnum}.3.duplication-metrics.txt",
        h="{path}metrics/saturation/{sample}-REP{repnum}.2.duplication-metrics.txt",
        i="{path}metrics/saturation/{sample}-REP{repnum}.1.duplication-metrics.txt",
        j="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.9.md.bai",
        k="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.8.md.bai",
        l="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.7.md.bai",
        m="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.6.md.bai",
        n="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.5.md.bai",
        o="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.4.md.bai",
        p="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.3.md.bai",
        q="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.2.md.bai",
        r="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.1.md.bai"
    output:
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_library_size.txt"
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

## Cleanup the uneeded intermediate files
rule SATURATION_clean_intermediate_data:
    input:
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_library_size.txt"
    output:
        "{path}operations/preprocessing/saturation/clean1.{sample}.{repnum}.done"
    shell:
        """
        rm -f {wildcards.path}preprocessing/8merged/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/raw/*REP{wildcards.repnum}*.bam
        rm -f {wildcards.path}preprocessing/saturation/downsampled/cs/*REP{wildcards.repnum}*.bam
        touch {output}
        """

## Call peaks with global normalization from downsampled libraries
rule SATURATION_peaks_globalnorm:
    input:
        a="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam",
        b="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bai"
    output:
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.{prob}_globalnorm_peaks.xls"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}.{wildcards.prob}_globalnorm --outdir {wildcards.path}preprocessing/saturation/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

## Call peaks with local normalization from downsampled libraries
rule SATURATION_peaks_localnorm:
    input:
        a="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam",
        b="{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bai"
    output:
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.{prob}_localnorm_peaks.xls"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}.{wildcards.prob}_localnorm --outdir {wildcards.path}preprocessing/saturation/peaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

## Count the number of peaks with global normalization from downsampled libraries and output to metrics
rule SATURATION_analyze_peak_saturation_globalnorm:
    input:
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.9_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.8_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.7_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.6_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.5_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.4_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.3_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.2_globalnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.1_globalnorm_peaks.xls"
    output:
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_globalnorm_numpeaks.txt"
    shell:
        "wc -l < {input} >> {output}"

## Count the number of peaks with local normalization from downsampled libraries and output to metrics
rule SATURATION_analyze_peak_saturation_localnorm:
    input:
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.9_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.8_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.7_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.6_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.5_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.4_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.3_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.2_localnorm_peaks.xls",
        "{path}preprocessing/saturation/peaks/{sample}-REP{repnum}.1_localnorm_peaks.xls"
    output:
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_localnorm_numpeaks.txt"
    shell:
        "wc -l < {input} >> {output}"

## An aggregator for the footprinting saturation analysis
rule AGGREGATOR_saturation_footprints:
    input:
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.9.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.8.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.7.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.6.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.5.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.4.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.3.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.2.done",
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.1.done"       
    output:
        "{path}operations/saturation/footprints/{sample}-REP{repnum}.{gene}.footprint.downsampled.done"
    shell:
        "touch {output}"

## Perform the footprint saturation analysis
rule SATURATION_analyze_raw_footprint_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bam",
        "{path}preprocessing/saturation/downsampled/md/{sample}-REP{repnum}.{prob}.md.bai",
        "snakeResources/sites/data/genes/{gene}.bindingSites.Rdata",
        "{path}metrics/saturation/{sample}-REP{repnum}.downsampled_library_size.txt"
    output:
        "{path}operations/saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.{prob}.done"
    benchmark:
        '{path}benchmark/saturation/{sample}-REP{repnum}.{gene}.rawfootprint.downsampled.{prob}.benchmark.txt'
    script:
        "snakeResources/scripts/saturation/snakeAnalyzeRawFootprintSaturation.R"

########################################################################################################################################
#### PEAK CALLING ######################################################################################################################
########################################################################################################################################
## Call peaks differentially across two samples
rule PEAKS_differential_peak_calling_2samples:
    input:
        a="{parentpath}{path1}bam/{sample1}.bam",
        b="{parentpath}{path2}bam/{sample2}.bam"
    output:
        "{parentpath}diffpeaks/ctrl-{sample1}.treat-{sample2}_globalnorm_peaks.narrowPeak"
    conda:
        "snakeResources/envs/macs2.yaml"
    threads:
        4
    shell:
        "macs2 callpeak -t {input.a} -c {input.b} -n ctrl-{wildcards.sample1}.treat-{wildcards.sample2}_globalnorm --outdir {wildcards.parentpath}diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

########################################################################################################################
#### FOOTPRINTING ######################################################################################################
########################################################################################################################
##
rule FOOTPRINTING_raw_analysis:
    input:
        "{path}bam/{sample}-REP{repnum}.bam",
        "{path}bam/{sample}-REP{repnum}.bai"
    output:
        "{path}operations/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.done"
    conda:
        "snakeResources/envs/RFootprint.yaml"
    threads:
        1
    resources:
        mem_mb=lambda params, attempt: attempt * 10000,
        run_time=lambda params, attempt: attempt * 4
    benchmark:
        '{path}benchmark/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.benchmark.txt'
    script:
        "snakeResources/scripts/analyzeFP.R"
