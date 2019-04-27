########################################################################################################################################
################################ GENERAL INFO ##########################################################################################
########################################################################################################################################
# This script is provided for modularization purposes to the ATACseq snakemake workflow
# It should be included in the main workflow using 'include'
# 
# The flow goes: copy bam -> raw FP analysis -> parse FP -> process FP

########################################################################################################################################
#### COPY BAM FILES ####################################################################################################################
########################################################################################################################################

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

########################################################################################################################################
#### GROUP AGGREGATOR ##################################################################################################################
########################################################################################################################################

rule PANTF_TFgroup_aggregator:
    input:
        '{path}footprints/operations/{mergedsample}.parseTF.group1.done',
        '{path}footprints/operations/{mergedsample}.parseTF.group2.done',
        '{path}footprints/operations/{mergedsample}.parseTF.group3.done',
        '{path}footprints/operations/{mergedsample}.parseTF.group4.done',
        '{path}footprints/operations/{mergedsample}.parseTF.group5.done',
        '{path}footprints/operations/{mergedsample}.parseTF.group6.done'
        # '{path}footprints/operations/{mergedsample}.parseTF.group7.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group8.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group9.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group10.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group11.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group12.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group13.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group14.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group15.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group16.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group17.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group18.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group19.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group20.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group21.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group22.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group23.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group24.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group25.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group26.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group27.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group28.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group29.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group30.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group31.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group32.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group33.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group34.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group35.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group36.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group37.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group38.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group39.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group40.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group41.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group42.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group43.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group44.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group45.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group46.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group47.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group48.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group49.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group50.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group51.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group52.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group53.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group54.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group55.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group56.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group57.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group58.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group59.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group60.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group61.done',
        # '{path}footprints/operations/{mergedsample}.parseTF.group62.done'
    output:
        "{path}footprints/operations/{mergedsample}.parseFP.allgroups.done"
    shell:
        "touch {output}"

########################################################################################################################################
#### PARSE FP ANALYSIS #################################################################################################################
########################################################################################################################################

rule panTF_group1:
    input:
        '{path}footprints/operations/processed/{mergedsample}.TFAP2A.{type}.bamcopy1.done', 
        '{path}footprints/operations/processed/{mergedsample}.NFIL3.{type}.bamcopy2.done', 
        '{path}footprints/operations/processed/{mergedsample}.HLF.{type}.bamcopy3.done', 
        '{path}footprints/operations/processed/{mergedsample}.NHLH1.{type}.bamcopy4.done', 
        '{path}footprints/operations/processed/{mergedsample}.MAX.{type}.bamcopy5.done', 
        '{path}footprints/operations/processed/{mergedsample}.USF1.{type}.bamcopy6.done', 
        '{path}footprints/operations/processed/{mergedsample}.CEBPA.{type}.bamcopy7.done', 
        '{path}footprints/operations/processed/{mergedsample}.EBF1.{type}.bamcopy8.done', 
        '{path}footprints/operations/processed/{mergedsample}.CEBPB.{type}.bamcopy9.done', 
        '{path}footprints/operations/processed/{mergedsample}.FOS.{type}.bamcopy10.done', 
        '{path}footprints/operations/processed/{mergedsample}.FOSL1.{type}.bamcopy11.done', 
        '{path}footprints/operations/processed/{mergedsample}.FOSL2.{type}.bamcopy12.done', 
        '{path}footprints/operations/processed/{mergedsample}.JUN.{type}.bamcopy13.done', 
        '{path}footprints/operations/processed/{mergedsample}.JUNB.{type}.bamcopy14.done', 
        '{path}footprints/operations/processed/{mergedsample}.JUND.{type}.bamcopy15.done', 
        '{path}footprints/operations/processed/{mergedsample}.MAFF.{type}.bamcopy16.done', 
        '{path}footprints/operations/processed/{mergedsample}.MAFK.{type}.bamcopy17.done', 
        '{path}footprints/operations/processed/{mergedsample}.TFAP2C.{type}.bamcopy18.done', 
        '{path}footprints/operations/processed/{mergedsample}.USF2.{type}.bamcopy19.done', 
        '{path}footprints/operations/processed/{mergedsample}.SREBF1.{type}.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.{type}.group1.done'
    shell:
        'touch {output}'