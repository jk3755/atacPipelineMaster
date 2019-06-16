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