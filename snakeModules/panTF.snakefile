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

rule processFP_COAD_MRs:
    input:
        '{path}footprints/operations/processed/{mergedsample}.TCF7.processFP.bamcopy1.done',
        '{path}footprints/operations/processed/{mergedsample}.MNX1.processFP.bamcopy2.done',
        '{path}footprints/operations/processed/{mergedsample}.POU5F1B.processFP.bamcopy3.done',
        '{path}footprints/operations/processed/{mergedsample}.ESRRA.processFP.bamcopy4.done',
        '{path}footprints/operations/processed/{mergedsample}.CDX2.processFP.bamcopy5.done',
        '{path}footprints/operations/processed/{mergedsample}.HNF4A.processFP.bamcopy6.done',
        '{path}footprints/operations/processed/{mergedsample}.GMEB2.processFP.bamcopy7.done',
        '{path}footprints/operations/processed/{mergedsample}.HOXA3.processFP.bamcopy8.done',
        '{path}footprints/operations/processed/{mergedsample}.OVOL1.processFP.bamcopy9.done',
        '{path}footprints/operations/processed/{mergedsample}.ASCL2.processFP.bamcopy10.done',
        '{path}footprints/operations/processed/{mergedsample}.ZSWIM1.processFP.bamcopy11.done',
        '{path}footprints/operations/processed/{mergedsample}.CBFA2T2.processFP.bamcopy12.done',
        '{path}footprints/operations/processed/{mergedsample}.PAX6.processFP.bamcopy13.done',
        '{path}footprints/operations/processed/{mergedsample}.ADNP.processFP.bamcopy14.done',
        '{path}footprints/operations/processed/{mergedsample}.TAF4.processFP.bamcopy15.done',
        '{path}footprints/operations/processed/{mergedsample}.ZMYND8.processFP.bamcopy16.done',
        '{path}footprints/operations/processed/{mergedsample}.ZNF696.processFP.bamcopy17.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.processFP.COADMR.done'
    shell:
        'touch {output}'
rule processFP_group1:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFAP2A.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFIL3.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HLF.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NHLH1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAX.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.USF1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CEBPA.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EBF1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CEBPB.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOS.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOSL1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOSL2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.JUN.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.JUNB.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.JUND.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAFF.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAFK.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFAP2C.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.USF2.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SREBF1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group1.done'
    shell:
        'touch {output}'
rule processFP_group2:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.SREBF2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AHR.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFAP4.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARNT.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATF6.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BACH1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BACH2.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CREB1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATF2.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCF3.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.XBP1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARID5B.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYOD1.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFE2.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYCN.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFE2L1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TEF.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATF3.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BATF.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCF12.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group2.done'
    shell:
        'touch {output}'
rule processFP_group3:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYC.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MXI1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHLHE40.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARNTL.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATF4.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATF7.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BATF3.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHLHA15.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHLHE41.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHLHE22.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHLHE23.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CEBPD.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CEBPE.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CEBPG.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CLOCK.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CREB3.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CREB3L1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DBP.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FIGLA.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HES5.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group3.done'
    shell:
        'touch {output}'
rule processFP_group4:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.HES7.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HEY1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HEY2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ID4.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.JDP2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAFG.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MESP1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MGA.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MLX.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MLXIPL.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MNT.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MSC.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYF6.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NEUROD2.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NEUROG2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NRL.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OLIG1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OLIG2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OLIG3.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCF4.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group4.done'
    shell:
        'touch {output}'
rule processFP_group5:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFAP2B.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFE3.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFEB.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFEC.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFAP2D.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARID3A.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARNT2.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATF1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATF5.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CREM.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DDIT3.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EPAS1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOSB.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HAND1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HES1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HIF1A.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HMGA1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HMGA2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAFA.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAFB.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group5.done'
    shell:
        'touch {output}'
rule processFP_group6:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAF.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MITF.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYOG.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NEUROD1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFE2L2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PTF1A.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TAL1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TWIST1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AIRE.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ALX1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ALX3.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ALX4.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ANDR.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AP2A.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AP2B.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AP2C.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AP2D.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARI3A.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARI5B.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARX.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group6.done'
    shell:
        'touch {output}'
rule processFP_group7:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ASCL2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATF6A.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ATOH1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BARH1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BARH2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BARX1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BARX2.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BC11A.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BCL6B.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BCL6.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHA15.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHE22.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHE23.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHE40.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHE41.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BMAL1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BPTF.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BRAC.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BRCA1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BSH.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group7.done'
    shell:
        'touch {output}'
rule processFP_group8:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.CDC5L.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CDX1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CDX2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CEBPZ.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CENPB.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.COE1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.COT1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.COT2.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CPEB1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CR3L1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CR3L2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CREB5.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CRX.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CTCFL.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CTCF.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CUX1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CUX2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CXXC1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DLX1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DLX2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group8.done'
    shell:
        'touch {output}'
rule processFP_group9:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.DLX3.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DLX4.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DLX5.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DLX6.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DMBX1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DPRX.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DRGX.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DUXA.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E2F1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E2F2.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E2F3.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E2F4.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E2F5.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E2F6.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E2F7.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E2F8.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.E4F1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EGR1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EGR2.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EGR3.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group9.done'
    shell:
        'touch {output}'
rule processFP_group10:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.EGR4.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EHF.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ELF1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ELF2.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ELF3.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ELF5.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ELK1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ELK3.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ELK4.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EMX1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EMX2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EOMES.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ERF.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ERG.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ERR1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ERR2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ERR3.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ESR1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ESR2.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ESX1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group10.done'
    shell:
        'touch {output}'
rule processFP_group11:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETS1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETS2.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETV1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETV2.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETV3.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETV4.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETV5.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETV6.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETV7.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EVI1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EVX1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EVX2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FEV.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FLI1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXA1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXA2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXA3.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXB1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXC1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXC2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group11.done'
    shell:
        'touch {output}'
rule processFP_group12:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXD1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXD2.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXD3.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXF1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXF2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXG1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXH1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXI1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXJ2.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXJ3.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXK1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXL1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXM1.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXO1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXO3.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXO4.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXO6.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXP2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXP3.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXQ1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group12.done'
    shell:
        'touch {output}'
rule processFP_group13:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.FUBP1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GABP1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GABPA.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GATA1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GATA2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GATA3.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GATA4.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GATA5.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GATA6.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GBX1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GBX2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GCM1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GCM2.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GCR.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GFI1B.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GFI1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GLI1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GLI2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GLI3.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GLIS1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group13.done'
    shell:
        'touch {output}'
rule processFP_group14:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.GLIS2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GLIS3.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GMEB2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GRHL1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GSC2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GSC.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GSX1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GSX2.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HBP1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HEN1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HESX1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HIC1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HIC2.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HINFP.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HLTF.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HMBX1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HME1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HME2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HMX1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HMX2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group14.done'
    shell:
        'touch {output}'
rule processFP_group15:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.HMX3.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNF1A.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNF1B.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNF4A.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNF4G.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNF6.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOMEZ.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HSF1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HSF2.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HSF4.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HSFY1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HTF4.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXA10.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXA11.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXA13.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXA1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXA2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXA5.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXA7.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXA9.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group15.done'
    shell:
        'touch {output}'
rule processFP_group16:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXB13.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXB1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXB2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXB3.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXB6.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXB7.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXB8.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXC10.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXC11.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXC12.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXC13.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXC6.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXC8.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXD10.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXD11.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXD12.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXD13.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXD3.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXD4.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXD8.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group16.done'
    shell:
        'touch {output}'
rule processFP_group17:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.HXD9.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IKZF1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.INSM1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF3.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF4.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF5.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF7.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF8.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF9.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRX2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRX3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ISL1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ISL2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ISX.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ITF2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KAISO.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF12.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF13.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group17.done'
    shell:
        'touch {output}'
rule processFP_group18:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF14.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF15.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF16.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF3.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF4.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF6.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF8.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LBX2.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LEF1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LHX2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LHX3.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LHX4.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LHX6.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LHX8.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LHX9.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LMX1A.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LMX1B.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAZ.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MBD2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group18.done'
    shell:
        'touch {output}'
rule processFP_group19:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.MCR.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MECP2.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEF2A.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEF2B.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEF2C.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEF2D.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEIS1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEIS2.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEIS3.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEOX1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEOX2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MGAP.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MIXL1.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MLXPL.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MNX1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MSX1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MSX2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MTF1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MUSC.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYBA.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group19.done'
    shell:
        'touch {output}'
rule processFP_group20:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYBB.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYB.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MZF1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NANOG.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NDF1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NDF2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NF2L1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NF2L2.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFAC1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFAC2.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFAC3.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFAC4.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFAT5.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFIA.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFIC.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFKB1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFKB2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFYA.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFYB.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFYC.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group20.done'
    shell:
        'touch {output}'
rule processFP_group21:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.NGN2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX21.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX22.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX23.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX25.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX28.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX31.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX32.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX61.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NKX62.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NOBOX.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NOTO.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR0B1.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR1D1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR1H2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR1H4.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR1I2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR1I3.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR2C1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR2C2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group21.done'
    shell:
        'touch {output}'
rule processFP_group22:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR2E1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR2E3.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR2F6.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR4A1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR4A2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR4A3.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR5A2.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR6A1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NRF1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ONEC2.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ONEC3.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OTX1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OTX2.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OVOL1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.P53.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.P5F1B.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.P63.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.P73.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group22.done'
    shell:
        'touch {output}'
rule processFP_group23:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX3.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX4.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX5.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX6.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX7.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX8.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PBX1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PBX2.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PBX3.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PDX1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PEBB.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PHX2A.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PHX2B.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PIT1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PITX1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PITX2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PITX3.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PKNX1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PKNX2.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PLAG1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group23.done'
    shell:
        'touch {output}'
rule processFP_group24:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.PLAL1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO2F1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO2F2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO2F3.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO3F1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO3F2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO3F3.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO3F4.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO4F1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO4F2.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO4F3.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO5F1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO6F1.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PO6F2.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PPARA.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PPARD.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PPARG.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRD14.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRDM1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRDM4.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group24.done'
    shell:
        'touch {output}'
rule processFP_group25:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRGR.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PROP1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PROX1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRRX1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRRX2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PURA.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RARA.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RARB.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RARG.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAX2.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RELB.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.REL.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.REST.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFX1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFX2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFX3.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFX4.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFX5.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RHXF1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RORA.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group25.done'
    shell:
        'touch {output}'
rule processFP_group26:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.RORG.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RREB1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RUNX1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RUNX2.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RUNX3.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RXRA.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RXRB.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RXRG.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RX.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SCRT1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SCRT2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SHOX2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SHOX.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMAD1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMAD2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMAD3.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMAD4.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMRC1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SNAI1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SNAI2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group26.done'
    shell:
        'touch {output}'
rule processFP_group27:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX10.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX11.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX13.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX15.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX17.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX18.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX21.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX2.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX3.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX4.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX5.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX7.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX8.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX9.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SP1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SP2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SP3.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SP4.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SPDEF.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group27.done'
    shell:
        'touch {output}'
rule processFP_group28:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.SPI1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SPIB.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SPIC.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SPZ1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SRBP1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SRBP2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SRF.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SRY.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STA5A.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STA5B.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STAT1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STAT2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STAT3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STAT4.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STAT6.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STF1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SUH.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBP.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBR1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX15.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group28.done'
    shell:
        'touch {output}'
rule processFP_group29:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX19.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX20.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX21.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX3.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX4.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBX5.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCF7.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TEAD1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TEAD3.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TEAD4.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TF2LX.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TF65.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TF7L1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TF7L2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFCP2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFDP1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFE2.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TGIF1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group29.done'
    shell:
        'touch {output}'
rule processFP_group30:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.TGIF2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.THAP1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.THA.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.THB.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TLX1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TWST1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TYY1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TYY2.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UBIP1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UNC4.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.VAX1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.VAX2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.VDR.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.VENTX.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.VSX1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.VSX2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.WT1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.YBOX1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBED1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBT18.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group30.done'
    shell:
        'touch {output}'
rule processFP_group31:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBT49.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBT7A.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBT7B.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB4.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB6.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZEB1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZEP1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZEP2.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZFHX3.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZFX.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZIC1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZIC2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZIC3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZIC4.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZKSC1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZKSC3.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN143.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN148.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN219.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN232.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group31.done'
    shell:
        'touch {output}'
rule processFP_group32:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN282.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN333.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN350.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN384.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN410.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN423.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN524.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN589.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN639.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN652.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN713.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN740.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZN784.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZSC16.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZSCA4.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ABCF2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.A1CF.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ACO1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ADARB1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AFF4.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group32.done'
    shell:
        'touch {output}'
rule processFP_group33:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.AGGF1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AKR1A1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ANXA1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ANXA11.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.APEX2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARFGAP1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ASCC1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ASPSCR1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AVEN.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BAD.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GPANK1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BAX.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BCL11A.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BOLL.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CELF4.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CELF5.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CELF6.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.C19orf25.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.C19orf40.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EXO5.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group33.done'
    shell:
        'touch {output}'
rule processFP_group34:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.LINC00471.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.C9orf156.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CANX.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CAT.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CBFA2T2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CBFB.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CBX7.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF830.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CCDC25.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CD59.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CDK2AP1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AGAP2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CFL2.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXN3.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CKMT1B.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CLK1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CNOT6.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NELFB.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CPSF4.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CSNK2B.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group34.done'
    shell:
        'touch {output}'
rule processFP_group35:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.CSTF2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CYB5R1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CYCS.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DAB2.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DAZAP1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ASAP3.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DDX20.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DDX4.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DDX43.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DDX53.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DGCR8.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DHX36.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DIABLO.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DIS3.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DNMT3A.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DTL.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DUS3L.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DUSP22.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DUSP26.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ECSIT.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group35.done'
    shell:
        'touch {output}'
rule processFP_group36:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.EDN1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EEF1D.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EIF5A2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ENO1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ESRRA.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ETFB.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EWSR1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EXOSC3.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.METTL21B.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FAM127B.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FEZ1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FEZF2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FGF19.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FHL2.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FIP1L1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SRRM3.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXP4.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GADD45A.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GIT2.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GLYCTK.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group36.done'
    shell:
        'touch {output}'
rule processFP_group37:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.GOT1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GPAM.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GPD1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GRHPR.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GTF2B.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GTF2H3.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GTF3C2.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GTF3C5.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GTPBP1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GTPBP6.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.H1FX.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.H2AFY.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.H2AFZ.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HCFC2.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HCLS1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HDAC8.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HHAT.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HHEX.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UBE2K.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HIRIP3.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group37.done'
    shell:
        'touch {output}'
rule processFP_group38:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.HIST1H2BN.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HIST2H2AB.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HIST2H2BE.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HLCS.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HMG20A.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNRNPA0.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNRNPA1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNRNPC.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNRNPH3.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HNRNPLL.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB13.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB9.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HP1BP3.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HSPA1L.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HSPA5.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HTATIP2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ID2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IL24.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ING3.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group38.done'
    shell:
        'touch {output}'
rule processFP_group39:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRF6.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IVD.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KDM5A.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KDM5D.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KCNIP1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KIAA0907.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KIF22.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LARP1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LARP4.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LAS1L.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CERS4.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UBXN1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CBX3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LRRFIP1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LSM6.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LUZP1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LUZP2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAGEA8.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAGED4B.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAGEF1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group39.done'
    shell:
        'touch {output}'
rule processFP_group40:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAGOH.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAP4K2.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MAPK1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MBTPS2.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MCTP2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MDM2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GLTPD1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBM42.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.WDR83.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MORN1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MRPL1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MRPL2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MRPS25.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MSI1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MSI2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MSRA.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MSRB3.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MTHFD1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MXD4.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYEF2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group40.done'
    shell:
        'touch {output}'
rule processFP_group41:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYLK.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NANOS1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NAP1L1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NCALD.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NCBP2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFATC3.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFATC4.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFIB.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFIX.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NME1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NMI.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NMRAL1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NNT.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NOC2L.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GAR1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NONO.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR2F1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NUCB1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NUP107.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NUP133.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group41.done'
    shell:
        'touch {output}'
rule processFP_group42:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.NXPH3.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ODC1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OTUD4.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.P4HB.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAXIP1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PCK2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PDCD11.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PDE6H.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PDLIM5.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PGAM2.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PHLDA2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PHOX2A.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PHTF1.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PICK1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PIK3C3.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PIR.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PKM.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PKNOX2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PLAGL1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PLG.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group42.done'
    shell:
        'touch {output}'
rule processFP_group43:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.POLE3.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POLI.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU3F2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU4F3.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PPP1R10.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PPP2R3B.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PPP5C.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PQBP1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRDX5.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRKRIR.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRNP.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PSMA6.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PSMC2.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PTCD1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PTPMT1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PURG.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.R3HDM2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAB14.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAB18.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAB2A.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group43.done'
    shell:
        'touch {output}'
rule processFP_group44:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAB7A.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAN.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAX.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBBP5.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBBP9.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBM17.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBM22.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBM3.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ESRP1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ESRP2.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBM7.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBM8A.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBFOX2.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBMS1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFC2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFC3.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFXANK.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RIOK2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MEX3C.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RNASEH2C.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group44.done'
    shell:
        'touch {output}'
rule processFP_group45:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.RNF138.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RPL35.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RPL6.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RPP25.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RPS10.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RPS4X.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RPS6KA5.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RUFY3.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RUVBL1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SCAND2P.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PDS5A.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SCMH1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SEMA4A.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SF1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SF3B1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SFT2D1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SLC18A1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMAP2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMCR7L.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMPX.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group45.done'
    shell:
        'touch {output}'
rule processFP_group46:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMUG1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SNAPC4.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SNAPC5.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SND1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SNRNP70.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SNRPB2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOCS4.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOD1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX14.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SPAG7.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SPATS2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SPR.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SRBD1.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SRP9.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SSBP3.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SSX2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SSX3.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STAU2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STUB1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SUCLG1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group46.done'
    shell:
        'touch {output}'
rule processFP_group47:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.TAF1A.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TAF9.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TAGLN2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBPL1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCEAL2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCEAL6.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFAM.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TGIF2LX.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.THAP5.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.THRA.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MED30.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TIA1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TIMELESS.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TIMM44.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TIMM8A.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TMSB4XP8.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TOB2.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TP73.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TPI1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TPPP.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group47.done'
    shell:
        'touch {output}'
rule processFP_group48:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.TRIM21.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TRIM69.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TRIP10.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TRMT1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TROVE2.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TSC22D4.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TSN.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TSNAX.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TULP1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.U2AF1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UBB.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UBE2V1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UGP2.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UQCRB.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.USP39.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UTP18.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.VAMP3.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EZR.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.VPS4B.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NELFA.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group48.done'
    shell:
        'touch {output}'
rule processFP_group49:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.WISP2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.XG.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.XRCC1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.YEATS4.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.YWHAE.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.YWHAZ.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB12.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB25.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB43.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB46.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZC3H7A.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZCCHC14.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZCCHC17.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZDHHC15.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZDHHC5.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZFP3.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZHX3.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZMAT2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZMAT4.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF124.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group49.done'
    shell:
        'touch {output}'
rule processFP_group50:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF131.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF160.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZKSCAN8.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZSCAN9.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF205.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF207.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB18.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF250.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF26.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF3.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF304.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RNF114.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZSCAN31.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF326.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF385A.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF503.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF510.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF655.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF671.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF695.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group50.done'
    shell:
        'touch {output}'
rule processFP_group51:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF706.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF71.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF720.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF76.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF766.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZRSR2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZSWIM1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.Myf.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.Pax6.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RORA_1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RORA_2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.YY1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TP53.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RELA.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF354C.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MIZF.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AP1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DUX4.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXP1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU2F2.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group51.done'
    shell:
        'touch {output}'
rule processFP_group52:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCF7L2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TP63.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB33.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF263.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.AR.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF5.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.T.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EN1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF143.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR3C1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ESRRB.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA5.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DMRT3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LBX1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU6F1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BARHL2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ELF4.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EN2.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA13.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC11.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group52.done'
    shell:
        'touch {output}'
rule processFP_group53:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ONECUT1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU4F2.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB7B.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB7C.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RHOXF1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UNCX.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR3C2.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SP8.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.YY2.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB7A.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF410.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF740.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ONECUT2.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ONECUT3.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYBL1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MYBL2.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PAX9.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PKNOX1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU1F1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU2F1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group53.done'
    shell:
        'touch {output}'
rule processFP_group54:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU3F1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU3F3.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU3F4.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU4F1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU5F1B.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU6F2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD12.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BSX.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HMBOX1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA10.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB2.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC10.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC12.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC13.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD11.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD13.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFATC2.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ASCL1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group54.done'
    shell:
        'touch {output}'
rule processFP_group55:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXK2.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GRHL2.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF9.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NR2F2.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU5F1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RBPJ.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SIX1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SIX2.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TEAD2.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF24.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF384.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF282.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZSCAN4.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RORB.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RORC.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCF7L1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HINFP1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF238.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF306.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF524.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group55.done'
    shell:
        'touch {output}'
rule processFP_group56:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF75A.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF784.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HSFY2.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.NFATC1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU2F3.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.POU5F1P1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHLHB2.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BHLHB3.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CART1.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB5.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD8.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRX5.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PHOX2B.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAXL1.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ESRRG.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.THRB.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.Trp53.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.Trp73.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB49.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group56.done'
    shell:
        'touch {output}'
rule processFP_group57:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF232.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF435.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF713.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ARID5A.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BARHL1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BBX.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.BCL3.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CHD1.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CHD2.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.CREB3L2.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DBX2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.DMC1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EBF3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EP300.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.EZH2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXJ1.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.FOXN1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GMEB1.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GTF2F1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.GTF2I.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group57.done'
    shell:
        'touch {output}'
rule processFP_group58:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.GZF1.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HCFC1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HDX.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HIVEP1.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HLX.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA11.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA3.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA4.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA6.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA7.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXA9.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB1.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB4.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB6.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB7.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXB8.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC4.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC5.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC6.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC8.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group58.done'
    shell:
        'touch {output}'
rule processFP_group59:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXC9.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD10.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD4.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.HOXD9.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IKZF2.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRX4.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.IRX6.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.KLF7.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LHX1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.LHX5.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MECOM.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.MTA3.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OSR1.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OSR2.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.OTP.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PATZ1.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PGR.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PML.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.PRDM14.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group59.done'
    shell:
        'touch {output}'
rule processFP_group60:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.RAD21.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RCOR1.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RFX7.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.RHOXF2.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SIN3A.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SIX3.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SIX4.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SIX5.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SIX6.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMARCC1.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMARCC2.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SMC3.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX12.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX30.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SOX6.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.SP100.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STAT5A.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.STAT5B.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TAF1.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TBL1XR1.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group60.done'
    shell:
        'touch {output}'
rule processFP_group61:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.TCF21.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFAP2E.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TFCP2L1.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.TLX2.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.UBP1.processFP.bamcopy5.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.WRNIP1.processFP.bamcopy6.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.YBX1.processFP.bamcopy7.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB14.processFP.bamcopy8.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB16.processFP.bamcopy9.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZBTB3.processFP.bamcopy10.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZKSCAN1.processFP.bamcopy11.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZKSCAN3.processFP.bamcopy12.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF148.processFP.bamcopy13.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF219.processFP.bamcopy14.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF274.processFP.bamcopy15.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF281.processFP.bamcopy16.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF333.processFP.bamcopy17.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF350.processFP.bamcopy18.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF35.processFP.bamcopy19.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF423.processFP.bamcopy20.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group61.done'
    shell:
        'touch {output}'
rule processFP_group62:
    input:
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF652.processFP.bamcopy1.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF691.processFP.bamcopy2.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF711.processFP.bamcopy3.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.ZNF8.processFP.bamcopy4.done', 
        '{path}footprints/operations/peaks/processed/{mergedsample}.Sox4.processFP.bamcopy5.done'
    output:
        '{path}footprints/operations/peaks/groups/{mergedsample}.processFP.group62.done'
    shell:
        'touch {output}'


