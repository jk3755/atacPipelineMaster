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
