rule AGGREGATOR_copy_bam:
    input:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.1.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.2.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.3.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.4.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.5.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.6.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.7.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.8.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.9.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.10.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.11.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.12.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.13.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.14.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.15.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.16.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.17.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.18.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.19.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.20.bam",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.1.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.2.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.3.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.4.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.5.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.6.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.7.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.8.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.9.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.10.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.11.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.12.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.13.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.14.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.15.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.16.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.17.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.18.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.19.bai",
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.20.bai"
    output:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.bamcopy.done"
    shell:
        "touch {output}"

rule PANTF_copy_bam:
    # The TF analysis script runs in 20 simultaneous processes
    # Each process will need to access the bam file individually
    # To significantly speed this analysis up, temporarily make 20 copies of the bam file
    # And assign each individual process a unique file to access
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bam"
    output:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.{bamcopy}.bam"
    shell:
        "cp {input} {output}"

rule PANTF_copy_bai:
    # The TF analysis script runs in 20 simultaneous processes
    # Each process will need to access the bam file individually
    # To significantly speed this analysis up, temporarily make 20 copies of the bam file
    # And assign each individual process a unique file to access
    input:
        "{path}preprocessing/11repmerged/{mergedsample}-repmerged.bai"
    output:
        "{path}preprocessing/11repmerged/copy/{mergedsample}-repmerged.{bamcopy}.bai"
    shell:
        "cp {input} {output}"

rule PANTF_remove_bamcopy:
    input:
        "{path}footprints/operations/{mergedsample}.rawTF.allgroups.done"
    output:
        "{path}footprints/operations/{mergedsample}.rawTF.analysis.done"
    shell:
         """
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bam
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bai
         rm -f {wildcards.path}preprocessing/11repmerged/copy/*.bamcopy.done
         touch {output}
         """