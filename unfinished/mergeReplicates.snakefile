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


rule STEP17_makebigwig_bamcov_merged_1replicate:
    # This rule will be used when only one replicate is present
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
        b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        c="{path}preprocessing/12bigwig/{mergedsample}-REP1of1.bw"
    output:
        "{path}preprocessing/12bigwig/{mergedsample}-repmerged.bw"
    resources:
        make_bigwig=1
    shell:
        "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

# rule STEP17_makebigwig_bamcov_merged_2replicates:
#     # This rule will be used when two replicates are present
#     # params:
#     # -b bam input
#     # -o output file
#     # -of output format
#     # -bs binsize in bp
#     # -p number of processors to use
#     # -v verbose mode
#     # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
#     input:
#         a="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
#         b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
#         c="{path}preprocessing/12bigwig/{mergedsample}-REP1of2.bw",
#         d="{path}preprocessing/12bigwig/{mergedsample}-REP2of2.bw"
#     output:
#         "{path}preprocessing/12bigwig/{mergedsample}-repmerged.bw"
#     resources:
#         make_bigwig=1
#     shell:
#         "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

# rule STEP17_makebigwig_bamcov_merged_3replicates:
#     # This rule will be used when three replicates are present
#     # params:
#     # -b bam input
#     # -o output file
#     # -of output format
#     # -bs binsize in bp
#     # -p number of processors to use
#     # -v verbose mode
#     # --normalizeUsing probably not useful for ATAC-seq normalization, need to find a good way (normalize to total library size)
#     input:
#         a="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
#         b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
#         c="{path}preprocessing/12bigwig/{mergedsample}-REP1of3.bw",
#         d="{path}preprocessing/12bigwig/{mergedsample}-REP2of3.bw",
#         e="{path}preprocessing/12bigwig/{mergedsample}-REP2of3.bw"
#     output:
#         "{path}preprocessing/12bigwig/{mergedsample}-repmerged.bw"
#     resources:
#         make_bigwig=1
#     shell:
#         "bamCoverage -b {input.a} -o {output} -of bigwig -bs 1 -p 20 -v"

rule STEP21_MACS2_peaks_merged_global_normilization_1replicate:
    # see above for notes applicable to MACS2 peak calling
    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP1of1.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP1of1.u.bai",
        c="{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt",
        d="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        e="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
        f="{path}peaks/macs2/individual/{mergedsample}-REP1of1_global_normalization_peaks.narrowPeak",
        g="{path}peaks/macs2/individual/{mergedsample}-REP1of1_local_normalization_peaks.narrowPeak"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
    shell:
        "macs2 callpeak -t {input.d} -n {wildcards.mergedsample}-merged_global_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# rule STEP21_MACS2_peaks_merged_global_normilization_2replicates:
#     # see above for notes applicable to MACS2 peak calling
#     input:
#         a="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bam",
#         b="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bai",
#         c="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bam",
#         d="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bai",
#         e="{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt",
#         f="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
#         g="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
#         h="{path}peaks/macs2/individual/{mergedsample}-REP1of2_global_normalization_peaks.narrowPeak",
#         i="{path}peaks/macs2/individual/{mergedsample}-REP1of2_local_normalization_peaks.narrowPeak",
#         j="{path}peaks/macs2/individual/{mergedsample}-REP2of2_global_normalization_peaks.narrowPeak",
#         k="{path}peaks/macs2/individual/{mergedsample}-REP2of2_local_normalization_peaks.narrowPeak"
#     output:
#         "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
#     shell:
#         "macs2 callpeak -t {input.f} -n {wildcards.mergedsample}-merged_global_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# rule STEP21_MACS2_peaks_merged_global_normilization_3replicates:
#     # see above for notes applicable to MACS2 peak calling
#     input:
#         a="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bam",
#         b="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bai",
#         c="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bam",
#         d="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bai",
#         e="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bam",
#         f="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bai",
#         g="{path}preprocessing/operations/{mergedsample}-preprocessing.done.txt",
#         h="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
#         i="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai",
#         j="{path}peaks/macs2/individual/{mergedsample}-REP1of3_global_normalization_peaks.narrowPeak",
#         k="{path}peaks/macs2/individual/{mergedsample}-REP1of3_local_normalization_peaks.narrowPeak",
#         l="{path}peaks/macs2/individual/{mergedsample}-REP2of3_global_normalization_peaks.narrowPeak",
#         m="{path}peaks/macs2/individual/{mergedsample}-REP2of3_local_normalization_peaks.narrowPeak",
#         n="{path}peaks/macs2/individual/{mergedsample}-REP3of3_global_normalization_peaks.narrowPeak",
#         o="{path}peaks/macs2/individual/{mergedsample}-REP3of3_local_normalization_peaks.narrowPeak"
#     output:
#         "{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak"
#     shell:
#         "macs2 callpeak -t {input.h} -n {wildcards.mergedsample}-merged_global_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

rule STEP22_MACS2_peaks_merged_local_normilization:
    # see above for notes applicable to MACS2 peak calling
    input:
        a="{path}peaks/macs2/merged/{mergedsample}-merged_global_normalization_peaks.narrowPeak",
        b="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam",
        c="{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai"
    output:
        "{path}peaks/macs2/merged/{mergedsample}-merged_local_normalization_peaks.narrowPeak"
    shell:
        "macs2 callpeak -t {input.b} -n {wildcards.mergedsample}-merged_local_normalization --outdir {wildcards.path}peaks/macs2/merged --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# rule FINISH_saturation_2rep:
#     input:
#         "{path}operations/{mergedsample}-REP1of2-downsample.done.txt",
#         "{path}operations/{mergedsample}-REP2of2-downsample.done.txt",
#     output:
#         "{path}operations/{mergedsample}-downsample.final.txt"
#     shell:
#         "touch {output}"

# rule FINISH_saturation_3rep:
#     input:
#         "{path}operations/{mergedsample}-REP1of3-downsample.done.txt",
#         "{path}operations/{mergedsample}-REP2of3-downsample.done.txt",
#         "{path}operations/{mergedsample}-REP3of3-downsample.done.txt",
#     output:
#         "{path}operations/{mergedsample}-downsample.final.txt"
#     shell:
#         "touch {output}"

# rule AGGREGATOR_fragsize_2reps:
#     input:
#         "{path}metrics/{mergedsample}-REP1of2.u.fragsizes.svg",
#         "{path}metrics/{mergedsample}-REP2of2.u.fragsizes.svg"
#     output:
#         "{path}operations/{mergedsample}.fragsizes.done.txt"
#     shell:
#         "touch {output}"

# rule AGGREGATOR_fragsize_3reps:
#     input:
#         "{path}metrics/{mergedsample}-REP1of3.u.fragsizes.svg",
#         "{path}metrics/{mergedsample}-REP2of3.u.fragsizes.svg",
#         "{path}metrics/{mergedsample}-REP3of3.u.fragsizes.svg"
#     output:
#         "{path}operations/{mergedsample}.fragsizes.done.txt"
#     shell:
#         "touch {output}"