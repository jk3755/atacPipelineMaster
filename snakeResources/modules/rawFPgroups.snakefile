rule rawFPanalysis_group1:
    input:
        '{path}operations/footprints/raw/{sample}.CTCF.rawFPanalysis.bamcopy1.done'
    output:
        '{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group1.done'
    shell:
        'touch {output}'

rule rawFPanalysis_group2:
    input:
        '{path}operations/footprints/raw/{sample}.MUSC.rawFPanalysis.bamcopy1.done'
    output:
        '{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group2.done'
    shell:
        'touch {output}'
