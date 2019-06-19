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
#### SPOOL #############################################################################################################################
########################################################################################################################################
rule correlation_mdst8_wt01:
    input:
        "mdst8/wt01/correlation/MDST8-WT-01-REP1.MDST8-WT-01-REP2.spearman.heatmap.svg"

########################################################################################################################################
#### CORRELATION RULES #################################################################################################################
########################################################################################################################################
rule CORRELATION_spearman_2samples:
    # parameters:
    # -b input bam files
    # -o output file name
    # -bs set the bin size used for comparison, default is 10000 bp
    # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
    # -p set the number of computing processors to use
    # -v verbose mode
    input:
        a="{path}preprocessing/10unique/{sample1}.u.bam",
        b="{path}preprocessing/10unique/{sample2}.u.bam"
    output:
        "{path}correlation/{sample1}.{sample2}.spearman.corrTest"
    benchmark:
        "{path}benchmark/correlation/{sample1}.{sample2}.spearman.corrTest.benchmark.txt"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} -o {output} -bs 10000 -p 20 -v"

rule CORRELATION_make_heatmap_2samples:
    input:
        "{path}correlation/{sample1}.{sample2}.spearman.corrTest"
    output:
        "{path}correlation/{sample1}.{sample2}.spearman.heatmap.svg"
    shell:
        "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"

rule CORRELATION_spearman_3samples:
    # parameters:
    # -b input bam files
    # -o output file name
    # -bs set the bin size used for comparison, default is 10000 bp
    # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
    # -p set the number of computing processors to use
    # -v verbose mode
    input:
        a="{path}preprocessing/10unique/{sample1}.u.bam",
        b="{path}preprocessing/10unique/{sample2}.u.bam",
        c="{path}preprocessing/10unique/{sample3}.u.bam",
    output:
        "{path}correlation/{sample1}.{sample2}.{sample2}.spearman.corrTest"
    benchmark:
        "{path}benchmark/correlation/{sample1}.{sample2}.{sample2}.spearman.corrTest.benchmark.txt"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} {input.c} -o {output} -bs 10000 -p 20 -v"

rule CORRELATION_make_heatmap_3samples:
    input:
        "{path}correlation/{sample1}.{sample2}.{sample3}.spearman.corrTest"
    output:
        "{path}correlation/{sample1}.{sample2}.{sample3}.spearman.heatmap.svg"
    shell:
        "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"
