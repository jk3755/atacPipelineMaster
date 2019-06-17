########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################
# This analysis relies on the deepTools package: https://deeptools.readthedocs.io/en/develop/index.html
#
# parameters:
# -b input bam files
# -o output file name
# -bs set the bin size used for comparison, default is 10000 bp
# -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
# -p set the number of computing processors to use
# -v verbose mode

########################################################################################################################################
#### SAMPLE CORRELATION ANALYSIS RULES #################################################################################################
########################################################################################################################################

rule CORRELATION_spearman_2samples:
    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bam",
        c="{path}preprocessing/10unique/{mergedsample}-REP1of2.u.bai",
        d="{path}preprocessing/10unique/{mergedsample}-REP2of2.u.bai"
    output:
        "{path}correlation/{mergedsample}.spearman.corrTest"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} -o {output} -bs 10000 -p 20 -v"

rule STEP23_sample_correlation_spearman_3replicates:

    input:
        a="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bam",
        b="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bam",
        c="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bam",
        d="{path}preprocessing/10unique/{mergedsample}-REP1of3.u.bai",
        e="{path}preprocessing/10unique/{mergedsample}-REP2of3.u.bai",
        f="{path}preprocessing/10unique/{mergedsample}-REP3of3.u.bai"
    output:
        "{path}correlation/{mergedsample}.spearman.corrTest"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} {input.c} -o {output} -bs 10000 -p 20 -v"

rule STEP24_makecorrheatmap:
    input:
        "{path}correlation/{sample}.spearman.corrTest"
    output:
        "{path}correlation/{sample}.spearman.heatmap.svg"
    shell:
        "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"

rule AGGREGATOR_correlation:
    input:
        "{path}correlation/{mergedsample}.spearman.heatmap.svg"
    output:
        "{path}operations/{mergedsample}-correlation.done.txt"
    shell:
        "touch {output}"

rule run_xsample_corr_h508_snu61_ls1034:
    input:
        "xsample_analysis/correlation/H508-wt-01.LS1034-wt-01.SNU61-wt-01.spearman.heatmap.svg"

rule run_xsample_corr_replicates_h508_snu61_ls1034:
    input:
        "xsample_analysis/correlation/H508-wt-01.LS1034-wt-01.SNU61-wt-01.spearman.heatmap.svg"