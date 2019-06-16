########################################################################################################################################
#### IMPORT MODULES AND CONFIG #########################################################################################################
########################################################################################################################################
include: "snakeModules/generateSites.snakefile"
include: "snakeModules/spoolPreprocessing.snakefile"
include: "snakeModules/spoolFootprinting.snakefile"
include: "snakeModules/rawFootprintGroups.snakefile"
########################################################################################################################################
#### PREPROCESSING RULES ###############################################################################################################
########################################################################################################################################

rule AGGREGATOR_preprocessing:
    input:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai",
        "{path}preprocessing/12bigwig/{sample}-REP{repnum}.bw",
        "{path}peaks/macs2/individual/{sample}-REP{repnum}_global_normalization_peaks.narrowPeak",
        "{path}peaks/macs2/individual/{sample}-REP{repnum}_local_normalization_peaks.narrowPeak",
        "{path}metrics/{sample}-REP{repnum}.peak.genome.coverage.txt",
        "{path}metrics/{sample}-REP{repnum}.u.fragsizes.svg",
        "{path}operations/{sample}-REP{repnum}.globalpeak.annotations.done.txt",
        "{path}operations/{sample}-REP{repnum}.localpeak.annotations.done.txt",
        "{path}metrics/{sample}-REP{repnum}.totalreads.Rdata",
        "{path}operations/{sample}-REP{repnum}.downsample.done.txt",
        "{path}metrics/{sample}-REP{repnum}.downsampled_library_size.txt",
        "{path}metrics/{sample}-REP{repnum}.downsampled_numpeaks.txt"
        #"{path}operations/{sample}-REP{repnum}.allgenes.footprint.downsampled.done.txt"
    output:
        "{path}operations/{sample}-REP{repnum}.preprocessing.complete.txt"
    shell:
        "touch {output}"

rule CLEAN_preprocessing:
    input:
        "{path}operations/{sample}-REP{repnum}.preprocessing.complete.txt"
    output:
        "{path}operations/{sample}-REP{repnum}.preprocessing.complete.clean.txt"
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
        rm -f {wildcards.path}preprocessing/8merged/*.bam
        rm -f {wildcards.path}preprocessing/9dedup/*.bam
        touch {output}
        """

rule AGGREGATOR_saturation:
    input:
        "{path}saturation/downsampled/{sample}-REP{repnum}.9.md.bai",
        "{path}saturation/downsampled/{sample}-REP{repnum}.8.md.bai",
        "{path}saturation/downsampled/{sample}-REP{repnum}.7.md.bai",
        "{path}saturation/downsampled/{sample}-REP{repnum}.6.md.bai",
        "{path}saturation/downsampled/{sample}-REP{repnum}.5.md.bai",
        "{path}saturation/downsampled/{sample}-REP{repnum}.4.md.bai",
        "{path}saturation/downsampled/{sample}-REP{repnum}.3.md.bai",
        "{path}saturation/downsampled/{sample}-REP{repnum}.2.md.bai",
        "{path}saturation/downsampled/{sample}-REP{repnum}.1.md.bai"      
    output:
        "{path}operations/{sample}-REP{repnum}.downsample.done.txt"
    shell:
        "touch {output}"

rule AGGREGATOR_saturation_footprints:
    input:
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.9.done",
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.8.done",
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.7.done",
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.6.done",
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.5.done",
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.4.done",
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.3.done",
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.2.done",
        "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.1.done",          
    output:
        "{path}operations/{sample}-REP{repnum}.{gene}.footprint.downsampled.done.txt"
    shell:
        "touch {output}"

rule AGGREGATOR_saturation_footprints_genes:
    input:
        "{path}operations/{sample}-REP{repnum}.CTCF.footprint.downsampled.done.txt",
        "{path}operations/{sample}-REP{repnum}.MNX1.footprint.downsampled.done.txt",
        "{path}operations/{sample}-REP{repnum}.CDX2.footprint.downsampled.done.txt"
    output:
        "{path}operations/{sample}-REP{repnum}.allgenes.footprint.downsampled.done.txt"
    shell:
        "touch {output}"

#### Preprocessing Rules #####

rule PREP_builddirstructure:
    # params: -p ignore error if existing, make parent dirs, -v verbose
    output:
        "{path}preprocessing/operations/dirtree.built"
    shell:
        """
        mkdir -p -v {wildcards.path}preprocessing
        mkdir -p -v {wildcards.path}preprocessing/2fastq {wildcards.path}preprocessing/3goodfastq {wildcards.path}preprocessing/4mycoalign {wildcards.path}preprocessing/5hg38align
        mkdir -p -v {wildcards.path}preprocessing/6rawbam {wildcards.path}preprocessing/7rgsort {wildcards.path}preprocessing/8merged {wildcards.path}preprocessing/9dedup
        mkdir -p -v {wildcards.path}preprocessing/10unique {wildcards.path}preprocessing/11repmerged {wildcards.path}preprocessing/12bigwig {wildcards.path}preprocessing/operations
        mkdir -p -v {wildcards.path}preprocessing/6rawbam/mitochondrial {wildcards.path}preprocessing/6rawbam/blacklist {wildcards.path}preprocessing/6rawbam/nonblacklist
        ##
        mkdir -p -v {wildcards.path}saturation
        mkdir -p -v {wildcards.path}saturation/footprints
        mkdir -p -v {wildcards.path}saturation/footprints/benchmark {wildcards.path}saturation/footprints/graphs {wildcards.path}saturation/footprints/operations
        mkdir -p -v {wildcards.path}saturation/footprints/data
        mkdir -p -v {wildcards.path}saturation/footprints/data/raw {wildcards.path}saturation/footprints/data/parsed 
        mkdir -p -v {wildcards.path}saturation/footprints/data/processed {wildcards.path}saturation/footprints/data/graphs
        mkdir -p -v {wildcards.path}saturation/complexity
        mkdir -p -v {wildcards.path}saturation/peaks
        mkdir -p -v {wildcards.path}saturation/downsampled
        ##
        mkdir -p -v {wildcards.path}footprints
        mkdir -p -v {wildcards.path}footprints/benchmark
        mkdir -p -v {wildcards.path}footprints/benchmark/parse {wildcards.path}footprints/benchmark/raw {wildcards.path}footprints/benchmark/processed
        mkdir -p -v {wildcards.path}footprints/data 
        mkdir -p -v {wildcards.path}footprints/data/parsed {wildcards.path}footprints/data/raw {wildcards.path}footprints/data/processed {wildcards.path}footprints/data/aggregated
        mkdir -p -v {wildcards.path}footprints/graphs
		mkdir -p -v {wildcards.path}footprints/graphs/insprob {wildcards.path}footprints/graphs/heatmaps
		mkdir -p -v {wildcards.path}footprints/operations
        mkdir -p -v {wildcards.path}footprints/operations/groups {wildcards.path}footprints/operations/parse {wildcards.path}footprints/operations/graphs
        mkdir -p -v {wildcards.path}footprints/operations/raw {wildcards.path}footprints/operations/processed {wildcards.path}footprints/operations/aggregated
        ##
        mkdir -p -v {wildcards.path}peaks
        mkdir -p -v {wildcards.path}peaks/genrich
        mkdir -p -v {wildcards.path}peaks/macs2
        mkdir -p -v {wildcards.path}peaks/macs2/individual {wildcards.path}peaks/macs2/merged
        ##
        mkdir -p -v {wildcards.path}operations
        mkdir -p -v {wildcards.path}metrics
        mkdir -p -v {wildcards.path}correlation
        ##
        touch {output}
        """

rule STEP1_gunzip:
    # params:
    # -k keep original files
    # -c write to standard output
    input:
        a="{path}preprocessing/1gz/{sample}-REP{repnum}_L{lane}_R{read}.fastq.gz",
        b="{path}preprocessing/operations/dirtree.built"
    output:
        c="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R{read}.fastq"
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
        a="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R1.fastq",
        b="{path}preprocessing/2fastq/{sample}-REP{repnum}_L{lane}_R2.fastq"
    output:
        c="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        d="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq"
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
        a="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq"
    output:
        "{path}preprocessing/4mycoalign/{sample}-REP{repnum}_L{lane}.myco.sam"
    shell:
        "bowtie2 -q -p 20 -X2000 -x genomes/myco/myco -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}_L{wildcards.lane}.myco.alignment.txt"

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
        a="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R1.good.fq",
        b="{path}preprocessing/3goodfastq/{sample}-REP{repnum}_L{lane}_R2.good.fq",
        c="{path}preprocessing/4mycoalign/{sample}-REP{repnum}_L{lane}.myco.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.sam"
    resources:
        hg38align=1
    shell:
        "bowtie2 -q -p 20 -X2000 -x genomes/hg38/hg38 -1 {input.a} -2 {input.b} -S {output} 2>{wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}_L{wildcards.lane}.hg38.alignment.txt"

rule STEP5_coordsort_sam:
    # coordinate sorting the sam files is required for blacklist filtering
    # params:
    # -o output file path
    # -O output file format
    input:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.sam"
    output:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.cs.sam"
    shell:
        "samtools sort {input} -o {output} -O sam"

rule STEP6_blacklistfilter_bamconversion:
    # remove aligned reads that map to hg38 blacklist regions as annotated by ENCODE
    input:
        "{path}preprocessing/5hg38align/{sample}-REP{repnum}_L{lane}.hg38.cs.sam"
    output:
        a="{path}preprocessing/6rawbam/blacklist/{sample}-REP{repnum}_L{lane}.hg38blacklist.bam",
        b="{path}preprocessing/6rawbam/nonblacklist/{sample}-REP{repnum}_L{lane}.blrm.bam"
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
        "{path}preprocessing/6rawbam/nonblacklist/{sample}-REP{repnum}_L{lane}.blrm.bam"
    output:
        a="{path}preprocessing/6rawbam/mitochondrial/{sample}-REP{repnum}_L{lane}.mitochondrial.bam",
        b="{path}preprocessing/6rawbam/{sample}-REP{repnum}_L{lane}.goodbam"
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
        "{path}preprocessing/6rawbam/{sample}-REP{repnum}_L{lane}.goodbam"
    output:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.rg.cs.bam"
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
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.rg.cs.bam"
    output:
        "{path}preprocessing/7rgsort/{sample}-REP{repnum}_L{lane}.clean.bam"
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
        a="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L1.clean.bam",
        b="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L2.clean.bam",
        c="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L3.clean.bam",
        d="{path}preprocessing/7rgsort/{sample}-REP{repnum}_L4.clean.bam"
    output:
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
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
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    output:
        "{path}preprocessing/9dedup/{sample}-REP{repnum}.dp.bam"
    shell:
        "java -jar programs/picard/picard.jar MarkDuplicates \
        I={input} \
        O={output} \
        M={wildcards.path}metrics/{wildcards.sample}-REP{wildcards.repnum}.duplication.txt \
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
        "{path}preprocessing/9dedup/{sample}-REP{repnum}.dp.bam"
    output:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam"
    shell:
        "samtools view -h -q 2 -b {input} > {output}"

rule STEP13_buildindex:
    # creates a bai index for the bam files
    # this is required for many downstream operations
    # the bai index allows other processes to access specific reads in the bam file without having to read through the entire bam contents to find them (its like a table of contents)
    # I specifies the input bam file
    # O specifies the output index file
    input:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam"
    output:
        "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    shell:
        "java -jar programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"

rule STEP14_makebigwig_bamcov:
    # params:
    # -b bam input
    # -o output file
    # -of output format
    # -bs binsize in bp
    # -p number of processors to use
    # -v verbose mode
    # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}preprocessing/12bigwig/{sample}-REP{repnum}.bw"
    resources:
        make_bigwig=1
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

rule STEP15_MACS2_peaks_global_normilization:
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
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}peaks/macs2/individual/{sample}-REP{repnum}_global_normalization_peaks.narrowPeak"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}_global_normalization --outdir {wildcards.path}peaks/macs2/individual --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule STEP16_MACS2_peaks_local_normalization:
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
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}peaks/macs2/individual/{sample}-REP{repnum}_local_normalization_peaks.narrowPeak"
    shell:
        "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}_local_normalization --outdir {wildcards.path}peaks/macs2/individual --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01"

rule STEP17_percent_peak_genome_coverage:
    # returns a fraction value of the basepairs of the genome covered by the merged peak file. multiple by 100 for percentages
    # parameters:
    # --echo output will be at least a three-column bed file
    # --bases-uniq the number of distinct bases from ref covered by overlap bed file
    # --delim change output delimeter from '|' to <delim>, e.g. '\t'
    input:
        a="{path}peaks/macs2/individual/{sample}-REP{repnum}_global_normalization_peaks.narrowPeak",
        b="genomes/hg38/hg38.extents.bed"
    output:
        "{path}metrics/{sample}-REP{repnum}.peak.genome.coverage.txt"
    shell:
        "bedmap --echo --bases-uniq --delim '\t' {input.b} {input.a} | awk 'BEGIN {{ genome_length = 0; masked_length = 0; }} {{ genome_length += ($3 - $2); masked_length += $4; }} END {{ print (masked_length / genome_length); }}' - > {output}"

rule STEP18_fragment_size_distribution:
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}metrics/{sample}-REP{repnum}.u.fragsizes.svg"
    script:
        "scripts/snakeFragSizeDist.R"

rule STEP19_annotate_peaks_global:
    input:
        "{path}peaks/macs2/individual/{sample}-REP{repnum}_global_normalization_peaks.narrowPeak"
    output:
        "{path}operations/{sample}-REP{repnum}.globalpeak.annotations.done.txt"
    script:
        "scripts/snakeAnnotatePeaks.R"

rule STEP19_annotate_peaks_local:
    input:
        "{path}peaks/macs2/individual/{sample}-REP{repnum}_local_normalization_peaks.narrowPeak"
    output:
        "{path}operations/{sample}-REP{repnum}.localpeak.annotations.done.txt"
    script:
        "scripts/snakeAnnotatePeaks.R"

rule STEP20_sample_total_reads:
    input:
        a="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam",
        b="{path}preprocessing/10unique/{sample}-REP{repnum}.u.bai"
    output:
        "{path}metrics/{sample}-REP{repnum}.totalreads.Rdata"
    script:
        "scripts/snakeCountSampleReads.R"

########################################################################################################################
#### FOOTPRINTING ######################################################################################################
########################################################################################################################

rule PANTF_raw_footprint_analysis:
    input:
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.{bamcopy}.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.{bamcopy}.bai",
        "sites/data/{gene}.bindingSites.Rdata",
        "{path}peaks/macs2/merged/{sample}-REP{repnum}-merged_global_normalization_peaks.narrowPeak",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.bamcopy.done"
    output:
        "{path}footprints/operations/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.bamcopy{bamcopy}.done"
    resources:
        analyzeRawFP=1
    benchmark:
        '{path}footprints/benchmark/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.bamcopy{bamcopy}.txt'
    script:
        "scripts/panTF/snakeAnalyzeRawFootprint.R"

rule PANTF_parse_and_generate_footprint_statistics:
    input:
        "{path}footprints/operations/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.bamcopy{bamcopy}.done",
        "{path}metrics/{sample}-REP{repnum}.totalreads.Rdata",
        "{path}peaks/macs2/merged/{sample}-REP{repnum}-merged_global_normalization_peaks.narrowPeak"
    output:
        "{path}footprints/operations/parse/{sample}-REP{repnum}.{gene}.parseFP.bamcopy{bamcopy}.done"
    resources:
        parseFootprint=1
    benchmark:
        '{path}footprints/benchmark/parse/{sample}-REP{repnum}.{gene}.bamcopy{bamcopy}.parseFP.txt'
    script:
        "scripts/panTF/snakeParseAndGenerateFootprintStats.R"

rule PANTF_process_footprint_analysis:
    input:
        "{path}footprints/operations/parse/{sample}-REP{repnum}.{gene}.parseFP.bamcopy{bamcopy}.done"
    output:
        "{path}footprints/operations/processed/{sample}-REP{repnum}.{gene}.processFP.bamcopy{bamcopy}.done"
    resources:
        processFootprint=1
    benchmark:
        '{path}footprints/benchmark/processed/{sample}-REP{repnum}.{gene}.bamcopy{bamcopy}.parseFP.txt'
    script:
        "scripts/panTF/snakeProcessFootprint.R"

rule PANTF_generate_tf_graphs:
    input:
        "{path}footprints/operations/processed/{sample}-REP{repnum}.{gene}.processFP.bamcopy{bamcopy}.done"
    output:
        "{path}footprints/operations/graphs/{sample}-REP{repnum}.{gene}.graphFP.bamcopy{bamcopy}.done"
    resources:
        graphFootprint=1
    script:
        "scripts/panTF/snakeGenerateFootprintGraphs.R"

rule PANTF_run_aggregator:
    input:
        "{path}footprints/data/processed/"
    output:
        "{path}footprints/operations/aggregated/{sample}-REP{repnum}.aggregated.done"
    script:
        "scripts/panTF/snakeAggregateProcessedFootprintData.R"

rule AGGREGATOR_copy_bam:
    input:
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.1.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.2.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.3.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.4.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.5.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.6.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.7.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.8.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.9.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.10.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.11.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.12.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.13.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.14.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.15.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.16.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.17.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.18.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.19.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.20.bam",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.1.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.2.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.3.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.4.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.5.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.6.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.7.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.8.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.9.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.10.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.11.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.12.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.13.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.14.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.15.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.16.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.17.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.18.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.19.bai",
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.20.bai"
    output:
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.bamcopy.done"
    shell:
        "touch {output}"

rule PANTF_copy_bam:
    # The TF analysis script runs in 20 simultaneous processes
    # Each process will need to access the bam file individually
    # To significantly speed this analysis up, temporarily make 20 copies of the bam file
    # And assign each individual process a unique file to access
    input:
        "{path}preprocessing/11repmerged/{sample}-REP{repnum}.bam"
    output:
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.{bamcopy}.bam"
    shell:
        "cp {input} {output}"

rule PANTF_copy_bai:
    # The TF analysis script runs in 20 simultaneous processes
    # Each process will need to access the bam file individually
    # To significantly speed this analysis up, temporarily make 20 copies of the bam file
    # And assign each individual process a unique file to access
    input:
        "{path}preprocessing/11repmerged/{sample}-REP{repnum}.bai"
    output:
        "{path}preprocessing/11repmerged/copy/{sample}-REP{repnum}.{bamcopy}.bai"
    shell:
        "cp {input} {output}"

rule PANTF_remove_bamcopy:
    input:
        "{path}footprints/operations/{sample}.rawTF.allgroups.done"
    output:
        "{path}footprints/operations/{sample}.rawTF.analysis.done"
    shell:
         """
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bam
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bai
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bamcopy.done
         touch {output}
         """