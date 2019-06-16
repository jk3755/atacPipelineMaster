########################################################################################################################################
#### DOWNSAMPLE RULES ##################################################################################################################
########################################################################################################################################

#
rule AGGREGATOR_saturation:
    input:
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.9.md.bai",
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.8.md.bai",
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.7.md.bai",
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.6.md.bai",
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.5.md.bai",
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.4.md.bai",
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.3.md.bai",
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.2.md.bai",
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.1.md.bai"      
    output:
        "{path}operations/saturation/{sample}-REP{repnum}.downsample.done"
    shell:
        "touch {output}"

#
rule SATURATION_downsample_bam:
    input:
        "{path}preprocessing/8merged/{sample}-REP{repnum}.lanemerge.bam"
    output:
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.{prob}.bam"
    shell:
        "java -jar programs/picard/picard.jar DownsampleSam \
        I={input} \
        O={output} \
        PROBABILITY=0.{wildcards.prob}"

#
rule SATURATION_sort_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.{prob}.bam"
    output:
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.{prob}.cs.bam"
    shell:
        "java -jar programs/picard/picard.jar SortSam \
        I={input} \
        O={output} \
        SORT_ORDER=coordinate"

#
rule SATURATION_purge_duplicates_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.{prob}.cs.bam"
    output:
        a="{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.{prob}.md.bam",
        b="{path}metrics/saturation/{sample}-REP{repnum}.{prob}.duplication-metrics.txt"
    shell:
        "java -Xmx5g -jar programs/picard/picard.jar MarkDuplicates \
        I={input} \
        O={output.a} \
        M={output.b} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true"

#
rule SATURATION_index_downsampled:
    input:
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.{prob}.md.bam"
    output:
        "{path}preprocessing/saturation/downsampled/{sample}-REP{repnum}.{prob}.md.bai"
    shell:
        "java -jar programs/picard/picard.jar BuildBamIndex \
        I={input} \
        O={output}"

# # Clean up intermediate data to this point
# rule STEP12b_clean_intermediate_data:
#     input:
#         "{path}preprocessing/10unique/{sample}-REP{repnum}.u.bam"
#     output:
#         "{path}operations/preprocessing/clean12b.{sample}.{repnum}.done"
#     shell:
#         """
#         rm -f {wildcards.path}preprocessing/8merged/*REP{wildcards.repnum}*.bam
#         rm -f {wildcards.path}preprocessing/9dedup/*REP{wildcards.repnum}*.bam
#         touch {output}
#         """

########################################################################################################################################
#### LIBRARY COMPLEXITY ################################################################################################################
########################################################################################################################################

# #
# rule STEP25_analyze_complexity_downsampled:
#     input:
#         a="{path}metrics/{sample}-REP{repnum}.9.duplication-metrics.txt",
#         b="{path}metrics/{sample}-REP{repnum}.8.duplication-metrics.txt",
#         c="{path}metrics/{sample}-REP{repnum}.7.duplication-metrics.txt",
#         d="{path}metrics/{sample}-REP{repnum}.6.duplication-metrics.txt",
#         e="{path}metrics/{sample}-REP{repnum}.5.duplication-metrics.txt",
#         f="{path}metrics/{sample}-REP{repnum}.4.duplication-metrics.txt",
#         g="{path}metrics/{sample}-REP{repnum}.3.duplication-metrics.txt",
#         h="{path}metrics/{sample}-REP{repnum}.2.duplication-metrics.txt",
#         i="{path}metrics/{sample}-REP{repnum}.1.duplication-metrics.txt",
#         j="{path}saturation/downsampled/{sample}-REP{repnum}.9.md.bai",
#         k="{path}saturation/downsampled/{sample}-REP{repnum}.8.md.bai",
#         l="{path}saturation/downsampled/{sample}-REP{repnum}.7.md.bai",
#         m="{path}saturation/downsampled/{sample}-REP{repnum}.6.md.bai",
#         n="{path}saturation/downsampled/{sample}-REP{repnum}.5.md.bai",
#         o="{path}saturation/downsampled/{sample}-REP{repnum}.4.md.bai",
#         p="{path}saturation/downsampled/{sample}-REP{repnum}.3.md.bai",
#         q="{path}saturation/downsampled/{sample}-REP{repnum}.2.md.bai",
#         r="{path}saturation/downsampled/{sample}-REP{repnum}.1.md.bai"
#     output:
#         "{path}metrics/{sample}-REP{repnum}.downsampled_library_size.txt"
#     shell:
#         """
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.a} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.b} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.c} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.d} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.e} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.f} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.g} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.h} >> {output}
#         awk '/ESTIMATED_LIBRARY_SIZE/ {{ getline; print $10; }}' {input.i} >> {output}
#         """

# ########################################################################################################################################
# #### PEAKS #############################################################################################################################
# ########################################################################################################################################
    
# rule STEP26_MACS2_peaks_downsampled:
#     input:
#         a="{path}saturation/downsampled/{sample}-REP{repnum}.{prob}.md.bam",
#         b="{path}saturation/downsampled/{sample}-REP{repnum}.{prob}.md.bai"
#     output:
#         "{path}saturation/peaks/{sample}-REP{repnum}.{prob}_global_normalization_peaks.xls"
#     shell:
#         "macs2 callpeak -t {input.a} -n {wildcards.sample}-REP{wildcards.repnum}.{wildcards.prob}_global_normalization --outdir {wildcards.path}saturation/peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01"

# rule STEP27_analyze_peak_saturation_downsampled:
#     input:
#         "{path}saturation/peaks/{sample}-REP{repnum}.9_global_normalization_peaks.xls",
#         "{path}saturation/peaks/{sample}-REP{repnum}.8_global_normalization_peaks.xls",
#         "{path}saturation/peaks/{sample}-REP{repnum}.7_global_normalization_peaks.xls",
#         "{path}saturation/peaks/{sample}-REP{repnum}.6_global_normalization_peaks.xls",
#         "{path}saturation/peaks/{sample}-REP{repnum}.5_global_normalization_peaks.xls",
#         "{path}saturation/peaks/{sample}-REP{repnum}.4_global_normalization_peaks.xls",
#         "{path}saturation/peaks/{sample}-REP{repnum}.3_global_normalization_peaks.xls",
#         "{path}saturation/peaks/{sample}-REP{repnum}.2_global_normalization_peaks.xls",
#         "{path}saturation/peaks/{sample}-REP{repnum}.1_global_normalization_peaks.xls"
#     output:
#         "{path}metrics/{sample}-REP{repnum}.downsampled_numpeaks.txt"
#     shell:
#         "wc -l < {input} >> {output}"



# ########################################################################################################################################
# #### FOOTPRINTS ########################################################################################################################
# ########################################################################################################################################

# rule STEP28_analyze_raw_footprint_downsampled:
#     input:
#         "{path}saturation/downsampled/{sample}-REP{repnum}.{prob}.bam",
#         "{path}saturation/downsampled/{sample}-REP{repnum}.{prob}.bai",
#         "sites/data/{gene}.bindingSites.Rdata",
#         "{path}operations/{sample}-REP{repnum}.downsample.done.txt"
#     output:
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.{prob}.done"
#     script:
#         "scripts/panTF/snakeAnalyzeRawFootprint.R"

# rule AGGREGATOR_saturation_footprints:
#     input:
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.9.done",
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.8.done",
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.7.done",
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.6.done",
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.5.done",
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.4.done",
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.3.done",
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.2.done",
#         "{path}saturation/footprints/raw/{sample}-REP{repnum}.{gene}.rawFPanalysis.downsampled.1.done",          
#     output:
#         "{path}operations/{sample}-REP{repnum}.{gene}.footprint.downsampled.done.txt"
#     shell:
#         "touch {output}"


# rule AGGREGATOR_saturation_footprints_genes:
#     input:
#         "{path}operations/{sample}-REP{repnum}.CTCF.footprint.downsampled.done.txt",
#         "{path}operations/{sample}-REP{repnum}.MNX1.footprint.downsampled.done.txt",
#         "{path}operations/{sample}-REP{repnum}.CDX2.footprint.downsampled.done.txt"
#     output:
#         "{path}operations/{sample}-REP{repnum}.allgenes.footprint.downsampled.done.txt"
#     shell:
#         "touch {output}"