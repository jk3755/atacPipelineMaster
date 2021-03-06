########################################################################################################################################
#### SPOOL #############################################################################################################################
########################################################################################################################################
rule correlation_mdst8_wt01:
    input:
        "mdst8/wt01/correlation/MDST8-WT-01-REP1.MDST8-WT-01-REP2.spearman.heatmap.svg"

rule correlation_lncap_ex01:
    input:
        "lncap/ex01/correlation/LNCaP-WT-01-REP1.LNCaP-WT-02-REP1.LNCaP-CR-01-REP1.LNCaP-CR-02-REP1.LNCaP-CR-04-REP1.LNCaP-CR-05-REP1.LNCaP-CR-07-REP1.LNCaP-CR-08-REP1.spearman.heatmap.svg"

##
rule CORRELATION_spearman_8samples:
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
        d="{path}preprocessing/10unique/{sample4}.u.bam",
        e="{path}preprocessing/10unique/{sample5}.u.bam",
        f="{path}preprocessing/10unique/{sample6}.u.bam",
        g="{path}preprocessing/10unique/{sample7}.u.bam",
        h="{path}preprocessing/10unique/{sample8}.u.bam"
    output:
        "{path}correlation/{sample1}.{sample2}.{sample3}.{sample4}.{sample5}.{sample6}.{sample7}.{sample8}.spearman.corrTest"
    benchmark:
        "{path}benchmark/correlation/{sample1}.{sample2}.{sample3}.{sample4}.{sample5}.{sample6}.{sample7}.{sample8}.spearman.corrTest.benchmark.txt"
    shell:
        "multiBamSummary bins -b {input.a} {input.b} {input.c} {input.d} {input.e} {input.f} {input.g} {input.h} -o {output} -bs 10000 -p 20 -v"

rule CORRELATION_make_heatmap_8samples:
    input:
        "{path}correlation/{sample1}.{sample2}.{sample3}.{sample4}.{sample5}.{sample6}.{sample7}.{sample8}.spearman.corrTest"
    output:
        "{path}correlation/{sample1}.{sample2}.{sample3}.{sample4}.{sample5}.{sample6}.{sample7}.{sample8}.spearman.heatmap.svg"
    shell:
        "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"
