########################################################################################################################################
################################ GENERAL INFO ##########################################################################################
########################################################################################################################################
# This script is provided for modularization purposes to the ATACseq snakemake workflow
# It should be included in the main workflow using 'include'

########################################################################################################################################
#### PAN TF FOOTPRINTING ANALYSIS ######################################################################################################
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

rule PANTF_TFgroup_aggregator:
    input:
        '{path}footprints/operations/{mergedsample}.rawTF.group1.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group2.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group3.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group4.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group5.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group6.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group7.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group8.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group9.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group10.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group11.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group12.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group13.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group14.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group15.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group16.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group17.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group18.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group19.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group20.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group21.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group22.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group23.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group24.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group25.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group26.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group27.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group28.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group29.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group30.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group31.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group32.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group33.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group34.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group35.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group36.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group37.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group38.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group39.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group40.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group41.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group42.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group43.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group44.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group45.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group46.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group47.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group48.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group49.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group50.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group51.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group52.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group53.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group54.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group55.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group56.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group57.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group58.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group59.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group60.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group61.done',
        '{path}footprints/operations/{mergedsample}.rawTF.group62.done'
    output:
        "{path}footprints/operations/{mergedsample}.rawTF.allgroups.done"
    shell:
        "touch {output}"

rule run_group1:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group1.done'
rule run_group2:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group2.done'
rule run_group3:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group3.done'
rule run_group4:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group4.done'
rule run_group5:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group5.done'
rule run_group6:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group6.done'
rule run_group7:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group7.done'
rule run_group8:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group8.done'
rule run_group9:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group9.done'
rule run_group10:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group10.done'
rule run_group11:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group11.done'
rule run_group12:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group12.done'
rule run_group13:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group13.done'
rule run_group14:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group14.done'
rule run_group15:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group15.done'
rule run_group16:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group16.done'
rule run_group17:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group17.done'
rule run_group18:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group18.done'
rule run_group19:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group19.done'
rule run_group20:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group20.done'
rule run_group21:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group21.done'
rule run_group22:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group22.done'
rule run_group23:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group23.done'
rule run_group24:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group24.done'
rule run_group25:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group25.done'
rule run_group26:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group26.done'
rule run_group27:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group27.done'
rule run_group28:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group28.done'
rule run_group29:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group29.done'
rule run_group30:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group30.done'
rule run_group31:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group31.done'
rule run_group32:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group32.done'
rule run_group33:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group33.done'
rule run_group34:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group34.done'
rule run_group35:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group35.done'
rule run_group36:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group36.done'
rule run_group37:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group37.done'
rule run_group38:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group38.done'
rule run_group39:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group39.done'
rule run_group40:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group40.done'
rule run_group41:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group41.done'
rule run_group42:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group42.done'
rule run_group43:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group43.done'
rule run_group44:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group44.done'
rule run_group45:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group45.done'
rule run_group46:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group46.done'
rule run_group47:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group47.done'
rule run_group48:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group48.done'
rule run_group49:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group49.done'
rule run_group50:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group50.done'
rule run_group51:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group51.done'
rule run_group52:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group52.done'
rule run_group53:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group53.done'
rule run_group54:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group54.done'
rule run_group55:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group55.done'
rule run_group56:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group56.done'
rule run_group57:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group57.done'
rule run_group58:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group58.done'
rule run_group59:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group59.done'
rule run_group60:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group60.done'
rule run_group61:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group61.done'
rule run_group62:
	input:
		'snu61/wt01/footprints/operations/SNU61-WT-01.rawTF.group62.done'

rule rawTF_group1:
	input:
		'{path}footprints/operations/{mergedsample}.TFAP2A.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.NFIL3.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.HLF.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.NHLH1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.MAX.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.USF1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.CEBPA.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.EBF1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.CEBPB.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.FOS.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.FOSL1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.FOSL2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.JUN.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.JUNB.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.JUND.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.MAFF.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.MAFK.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.TFAP2C.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.USF2.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.SREBF1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group1.done'
	shell:
		'touch {output}'
rule rawTF_group2:
	input:
		'{path}footprints/operations/{mergedsample}.SREBF2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.AHR.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.TFAP4.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ARNT.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ATF6.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.BACH1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.BACH2.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.CREB1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ATF2.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.TCF3.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.XBP1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.ARID5B.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.MYOD1.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.NFE2.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.MYCN.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.NFE2L1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.TEF.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ATF3.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.BATF.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.TCF12.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group2.done'
	shell:
		'touch {output}'
rule rawTF_group3:
	input:
		'{path}footprints/operations/{mergedsample}.MYC.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.MXI1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.BHLHE40.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ARNTL.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ATF4.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ATF7.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.BATF3.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.BHLHA15.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.BHLHE41.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.BHLHE22.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.BHLHE23.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.CEBPD.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.CEBPE.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.CEBPG.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.CLOCK.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.CREB3.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.CREB3L1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.DBP.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.FIGLA.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.HES5.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group3.done'
	shell:
		'touch {output}'
rule rawTF_group4:
	input:
		'{path}footprints/operations/{mergedsample}.HES7.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.HEY1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.HEY2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ID4.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.JDP2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.MAFG.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.MESP1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.MGA.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.MLX.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.MLXIPL.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.MNT.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.MSC.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.MYF6.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.NEUROD2.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.NEUROG2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.NRL.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.OLIG1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.OLIG2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.OLIG3.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.TCF4.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group4.done'
	shell:
		'touch {output}'
rule rawTF_group5:
	input:
		'{path}footprints/operations/{mergedsample}.TFAP2B.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.TFE3.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.TFEB.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.TFEC.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.TFAP2D.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ARID3A.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ARNT2.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ATF1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ATF5.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.CREM.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.DDIT3.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.EPAS1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.FOSB.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.HAND1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.HES1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.HIF1A.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HMGA1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.HMGA2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.MAFA.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.MAFB.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group5.done'
	shell:
		'touch {output}'
rule rawTF_group6:
	input:
		'{path}footprints/operations/{mergedsample}.MAF.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.MITF.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.MYOG.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.NEUROD1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.NFE2L2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.PTF1A.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.TAL1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.TWIST1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.AIRE.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ALX1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ALX3.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.ALX4.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ANDR.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.AP2A.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.AP2B.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.AP2C.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.AP2D.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ARI3A.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ARI5B.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ARX.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group6.done'
	shell:
		'touch {output}'
rule rawTF_group7:
	input:
		'{path}footprints/operations/{mergedsample}.ASCL2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ATF6A.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ATOH1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.BARH1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.BARH2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.BARX1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.BARX2.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.BC11A.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.BCL6B.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.BCL6.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.BHA15.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.BHE22.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.BHE23.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.BHE40.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.BHE41.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.BMAL1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.BPTF.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.BRAC.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.BRCA1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.BSH.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group7.done'
	shell:
		'touch {output}'
rule rawTF_group8:
	input:
		'{path}footprints/operations/{mergedsample}.CDC5L.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.CDX1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.CDX2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.CEBPZ.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.CENPB.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.COE1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.COT1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.COT2.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.CPEB1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.CR3L1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.CR3L2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.CREB5.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.CRX.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.CTCFL.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.CTCF.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.CUX1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.CUX2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.CXXC1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.DLX1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.DLX2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group8.done'
	shell:
		'touch {output}'
rule rawTF_group9:
	input:
		'{path}footprints/operations/{mergedsample}.DLX3.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.DLX4.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.DLX5.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.DLX6.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.DMBX1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.DPRX.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.DRGX.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.DUXA.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.E2F1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.E2F2.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.E2F3.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.E2F4.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.E2F5.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.E2F6.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.E2F7.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.E2F8.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.E4F1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.EGR1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.EGR2.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.EGR3.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group9.done'
	shell:
		'touch {output}'
rule rawTF_group10:
	input:
		'{path}footprints/operations/{mergedsample}.EGR4.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.EHF.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ELF1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ELF2.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ELF3.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ELF5.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ELK1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ELK3.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ELK4.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.EMX1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.EMX2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.EOMES.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ERF.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.ERG.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.ERR1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.ERR2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.ERR3.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ESR1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ESR2.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ESX1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group10.done'
	shell:
		'touch {output}'
rule rawTF_group11:
	input:
		'{path}footprints/operations/{mergedsample}.ETS1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ETS2.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ETV1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ETV2.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ETV3.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ETV4.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ETV5.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ETV6.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ETV7.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.EVI1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.EVX1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.EVX2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.FEV.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.FLI1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.FOXA1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.FOXA2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.FOXA3.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.FOXB1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.FOXC1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.FOXC2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group11.done'
	shell:
		'touch {output}'
rule rawTF_group12:
	input:
		'{path}footprints/operations/{mergedsample}.FOXD1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.FOXD2.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.FOXD3.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.FOXF1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.FOXF2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.FOXG1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.FOXH1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.FOXI1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.FOXJ2.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.FOXJ3.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.FOXK1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.FOXL1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.FOXM1.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.FOXO1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.FOXO3.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.FOXO4.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.FOXO6.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.FOXP2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.FOXP3.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.FOXQ1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group12.done'
	shell:
		'touch {output}'
rule rawTF_group13:
	input:
		'{path}footprints/operations/{mergedsample}.FUBP1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.GABP1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.GABPA.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.GATA1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.GATA2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.GATA3.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.GATA4.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.GATA5.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.GATA6.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.GBX1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.GBX2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.GCM1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.GCM2.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.GCR.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.GFI1B.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.GFI1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.GLI1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.GLI2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.GLI3.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.GLIS1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group13.done'
	shell:
		'touch {output}'
rule rawTF_group14:
	input:
		'{path}footprints/operations/{mergedsample}.GLIS2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.GLIS3.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.GMEB2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.GRHL1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.GSC2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.GSC.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.GSX1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.GSX2.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.HBP1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.HEN1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.HESX1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.HIC1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.HIC2.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.HINFP.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.HLTF.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.HMBX1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HME1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.HME2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.HMX1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.HMX2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group14.done'
	shell:
		'touch {output}'
rule rawTF_group15:
	input:
		'{path}footprints/operations/{mergedsample}.HMX3.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.HNF1A.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.HNF1B.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.HNF4A.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.HNF4G.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.HNF6.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.HOMEZ.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.HSF1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.HSF2.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.HSF4.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.HSFY1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.HTF4.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.HXA10.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.HXA11.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.HXA13.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.HXA1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HXA2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.HXA5.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.HXA7.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.HXA9.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group15.done'
	shell:
		'touch {output}'
rule rawTF_group16:
	input:
		'{path}footprints/operations/{mergedsample}.HXB13.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.HXB1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.HXB2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.HXB3.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.HXB6.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.HXB7.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.HXB8.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.HXC10.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.HXC11.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.HXC12.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.HXC13.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.HXC6.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.HXC8.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.HXD10.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.HXD11.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.HXD12.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HXD13.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.HXD3.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.HXD4.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.HXD8.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group16.done'
	shell:
		'touch {output}'
rule rawTF_group17:
	input:
		'{path}footprints/operations/{mergedsample}.HXD9.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.IKZF1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.INSM1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.IRF1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.IRF2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.IRF3.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.IRF4.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.IRF5.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.IRF7.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.IRF8.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.IRF9.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.IRX2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.IRX3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.ISL1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.ISL2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.ISX.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.ITF2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.KAISO.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.KLF12.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.KLF13.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group17.done'
	shell:
		'touch {output}'
rule rawTF_group18:
	input:
		'{path}footprints/operations/{mergedsample}.KLF14.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.KLF15.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.KLF16.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.KLF1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.KLF3.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.KLF4.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.KLF6.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.KLF8.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.LBX2.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.LEF1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.LHX2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.LHX3.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.LHX4.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.LHX6.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.LHX8.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.LHX9.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.LMX1A.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.LMX1B.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.MAZ.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.MBD2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group18.done'
	shell:
		'touch {output}'
rule rawTF_group19:
	input:
		'{path}footprints/operations/{mergedsample}.MCR.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.MECP2.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.MEF2A.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.MEF2B.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.MEF2C.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.MEF2D.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.MEIS1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.MEIS2.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.MEIS3.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.MEOX1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.MEOX2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.MGAP.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.MIXL1.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.MLXPL.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.MNX1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.MSX1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.MSX2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.MTF1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.MUSC.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.MYBA.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group19.done'
	shell:
		'touch {output}'
rule rawTF_group20:
	input:
		'{path}footprints/operations/{mergedsample}.MYBB.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.MYB.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.MZF1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.NANOG.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.NDF1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.NDF2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.NF2L1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.NF2L2.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.NFAC1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.NFAC2.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.NFAC3.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.NFAC4.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.NFAT5.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.NFIA.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.NFIC.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.NFKB1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.NFKB2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.NFYA.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.NFYB.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.NFYC.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group20.done'
	shell:
		'touch {output}'
rule rawTF_group21:
	input:
		'{path}footprints/operations/{mergedsample}.NGN2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.NKX21.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.NKX22.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.NKX23.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.NKX25.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.NKX28.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.NKX31.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.NKX32.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.NKX61.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.NKX62.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.NOBOX.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.NOTO.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.NR0B1.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.NR1D1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.NR1H2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.NR1H4.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.NR1I2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.NR1I3.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.NR2C1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.NR2C2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group21.done'
	shell:
		'touch {output}'
rule rawTF_group22:
	input:
		'{path}footprints/operations/{mergedsample}.NR2E1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.NR2E3.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.NR2F6.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.NR4A1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.NR4A2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.NR4A3.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.NR5A2.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.NR6A1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.NRF1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ONEC2.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ONEC3.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.OTX1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.OTX2.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.OVOL1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.P53.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.P5F1B.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.P63.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.P73.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.PAX1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.PAX2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group22.done'
	shell:
		'touch {output}'
rule rawTF_group23:
	input:
		'{path}footprints/operations/{mergedsample}.PAX3.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.PAX4.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.PAX5.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.PAX6.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.PAX7.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.PAX8.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.PBX1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.PBX2.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.PBX3.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.PDX1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.PEBB.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.PHX2A.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.PHX2B.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.PIT1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.PITX1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.PITX2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.PITX3.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.PKNX1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.PKNX2.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.PLAG1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group23.done'
	shell:
		'touch {output}'
rule rawTF_group24:
	input:
		'{path}footprints/operations/{mergedsample}.PLAL1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.PO2F1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.PO2F2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.PO2F3.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.PO3F1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.PO3F2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.PO3F3.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.PO3F4.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.PO4F1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.PO4F2.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.PO4F3.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.PO5F1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.PO6F1.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.PO6F2.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.PPARA.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.PPARD.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.PPARG.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.PRD14.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.PRDM1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.PRDM4.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group24.done'
	shell:
		'touch {output}'
rule rawTF_group25:
	input:
		'{path}footprints/operations/{mergedsample}.PRGR.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.PROP1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.PROX1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.PRRX1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.PRRX2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.PURA.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.RARA.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.RARB.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.RARG.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.RAX2.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.RELB.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.REL.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.REST.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.RFX1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.RFX2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.RFX3.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.RFX4.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.RFX5.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.RHXF1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.RORA.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group25.done'
	shell:
		'touch {output}'
rule rawTF_group26:
	input:
		'{path}footprints/operations/{mergedsample}.RORG.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.RREB1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.RUNX1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.RUNX2.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.RUNX3.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.RXRA.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.RXRB.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.RXRG.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.RX.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.SCRT1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.SCRT2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.SHOX2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.SHOX.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.SMAD1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.SMAD2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.SMAD3.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.SMAD4.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.SMRC1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.SNAI1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.SNAI2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group26.done'
	shell:
		'touch {output}'
rule rawTF_group27:
	input:
		'{path}footprints/operations/{mergedsample}.SOX10.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.SOX11.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.SOX13.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.SOX15.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.SOX17.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.SOX18.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.SOX1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.SOX21.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.SOX2.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.SOX3.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.SOX4.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.SOX5.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.SOX7.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.SOX8.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.SOX9.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.SP1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.SP2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.SP3.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.SP4.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.SPDEF.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group27.done'
	shell:
		'touch {output}'
rule rawTF_group28:
	input:
		'{path}footprints/operations/{mergedsample}.SPI1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.SPIB.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.SPIC.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.SPZ1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.SRBP1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.SRBP2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.SRF.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.SRY.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.STA5A.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.STA5B.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.STAT1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.STAT2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.STAT3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.STAT4.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.STAT6.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.STF1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.SUH.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.TBP.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.TBR1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.TBX15.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group28.done'
	shell:
		'touch {output}'
rule rawTF_group29:
	input:
		'{path}footprints/operations/{mergedsample}.TBX19.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.TBX1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.TBX20.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.TBX21.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.TBX2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.TBX3.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.TBX4.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.TBX5.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.TCF7.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.TEAD1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.TEAD3.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.TEAD4.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.TF2LX.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.TF65.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.TF7L1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.TF7L2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.TFCP2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.TFDP1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.TFE2.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.TGIF1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group29.done'
	shell:
		'touch {output}'
rule rawTF_group30:
	input:
		'{path}footprints/operations/{mergedsample}.TGIF2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.THAP1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.THA.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.THB.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.TLX1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.TWST1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.TYY1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.TYY2.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.UBIP1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.UNC4.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.VAX1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.VAX2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.VDR.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.VENTX.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.VSX1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.VSX2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.WT1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.YBOX1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ZBED1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ZBT18.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group30.done'
	shell:
		'touch {output}'
rule rawTF_group31:
	input:
		'{path}footprints/operations/{mergedsample}.ZBT49.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ZBT7A.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ZBT7B.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB4.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB6.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ZEB1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ZEP1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ZEP2.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ZFHX3.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ZFX.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ZIC1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.ZIC2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ZIC3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.ZIC4.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.ZKSC1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.ZKSC3.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.ZN143.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ZN148.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ZN219.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ZN232.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group31.done'
	shell:
		'touch {output}'
rule rawTF_group32:
	input:
		'{path}footprints/operations/{mergedsample}.ZN282.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ZN333.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ZN350.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ZN384.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ZN410.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ZN423.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ZN524.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ZN589.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ZN639.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ZN652.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ZN713.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.ZN740.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ZN784.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.ZSC16.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.ZSCA4.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.ABCF2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.A1CF.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ACO1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ADARB1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.AFF4.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group32.done'
	shell:
		'touch {output}'
rule rawTF_group33:
	input:
		'{path}footprints/operations/{mergedsample}.AGGF1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.AKR1A1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ANXA1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ANXA11.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.APEX2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ARFGAP1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ASCC1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ASPSCR1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.AVEN.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.BAD.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.GPANK1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.BAX.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.BCL11A.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.BOLL.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.CELF4.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.CELF5.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.CELF6.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.C19orf25.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.C19orf40.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.EXO5.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group33.done'
	shell:
		'touch {output}'
rule rawTF_group34:
	input:
		'{path}footprints/operations/{mergedsample}.LINC00471.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.C9orf156.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.CANX.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.CAT.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.CBFA2T2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.CBFB.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.CBX7.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ZNF830.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.CCDC25.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.CD59.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.CDK2AP1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.AGAP2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.CFL2.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.FOXN3.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.CKMT1B.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.CLK1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.CNOT6.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.NELFB.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.CPSF4.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.CSNK2B.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group34.done'
	shell:
		'touch {output}'
rule rawTF_group35:
	input:
		'{path}footprints/operations/{mergedsample}.CSTF2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.CYB5R1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.CYCS.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.DAB2.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.DAZAP1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ASAP3.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.DDX20.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.DDX4.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.DDX43.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.DDX53.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.DGCR8.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.DHX36.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.DIABLO.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.DIS3.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.DNMT3A.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.DTL.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.DUS3L.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.DUSP22.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.DUSP26.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ECSIT.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group35.done'
	shell:
		'touch {output}'
rule rawTF_group36:
	input:
		'{path}footprints/operations/{mergedsample}.EDN1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.EEF1D.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.EIF5A2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ENO1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ESRRA.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ETFB.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.EWSR1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.EXOSC3.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.METTL21B.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.FAM127B.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.FEZ1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.FEZF2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.FGF19.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.FHL2.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.FIP1L1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.SRRM3.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.FOXP4.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.GADD45A.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.GIT2.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.GLYCTK.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group36.done'
	shell:
		'touch {output}'
rule rawTF_group37:
	input:
		'{path}footprints/operations/{mergedsample}.GOT1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.GPAM.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.GPD1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.GRHPR.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.GTF2B.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.GTF2H3.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.GTF3C2.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.GTF3C5.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.GTPBP1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.GTPBP6.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.H1FX.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.H2AFY.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.H2AFZ.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.HCFC2.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.HCLS1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.HDAC8.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HHAT.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.HHEX.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.UBE2K.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.HIRIP3.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group37.done'
	shell:
		'touch {output}'
rule rawTF_group38:
	input:
		'{path}footprints/operations/{mergedsample}.HIST1H2BN.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.HIST2H2AB.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.HIST2H2BE.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.HLCS.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.HMG20A.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.HNRNPA0.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.HNRNPA1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.HNRNPC.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.HNRNPH3.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.HNRNPLL.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.HOXB13.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.HOXB9.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.HOXD3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.HP1BP3.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.HSPA1L.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.HSPA5.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HTATIP2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ID2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.IL24.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ING3.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group38.done'
	shell:
		'touch {output}'
rule rawTF_group39:
	input:
		'{path}footprints/operations/{mergedsample}.IRF6.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.IVD.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.KDM5A.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.KDM5D.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.KCNIP1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.KIAA0907.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.KIF22.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.LARP1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.LARP4.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.LAS1L.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.CERS4.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.UBXN1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.CBX3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.LRRFIP1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.LSM6.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.LUZP1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.LUZP2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.MAGEA8.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.MAGED4B.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.MAGEF1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group39.done'
	shell:
		'touch {output}'
rule rawTF_group40:
	input:
		'{path}footprints/operations/{mergedsample}.MAGOH.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.MAP4K2.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.MAPK1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.MBTPS2.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.MCTP2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.MDM2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.GLTPD1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.RBM42.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.WDR83.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.MORN1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.MRPL1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.MRPL2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.MRPS25.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.MSI1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.MSI2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.MSRA.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.MSRB3.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.MTHFD1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.MXD4.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.MYEF2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group40.done'
	shell:
		'touch {output}'
rule rawTF_group41:
	input:
		'{path}footprints/operations/{mergedsample}.MYLK.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.NANOS1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.NAP1L1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.NCALD.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.NCBP2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.NFATC3.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.NFATC4.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.NFIB.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.NFIX.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.NME1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.NMI.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.NMRAL1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.NNT.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.NOC2L.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.GAR1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.NONO.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.NR2F1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.NUCB1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.NUP107.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.NUP133.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group41.done'
	shell:
		'touch {output}'
rule rawTF_group42:
	input:
		'{path}footprints/operations/{mergedsample}.NXPH3.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ODC1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.OTUD4.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.P4HB.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.PAXIP1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.PCK2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.PDCD11.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.PDE6H.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.PDLIM5.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.PGAM2.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.PHLDA2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.PHOX2A.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.PHTF1.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.PICK1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.PIK3C3.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.PIR.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.PKM.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.PKNOX2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.PLAGL1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.PLG.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group42.done'
	shell:
		'touch {output}'
rule rawTF_group43:
	input:
		'{path}footprints/operations/{mergedsample}.POLE3.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.POLI.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.POU3F2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.POU4F3.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.PPP1R10.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.PPP2R3B.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.PPP5C.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.PQBP1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.PRDX5.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.PRKRIR.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.PRNP.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.PSMA6.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.PSMC2.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.PTCD1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.PTPMT1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.PURG.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.R3HDM2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.RAB14.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.RAB18.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.RAB2A.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group43.done'
	shell:
		'touch {output}'
rule rawTF_group44:
	input:
		'{path}footprints/operations/{mergedsample}.RAB7A.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.RAN.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.RAX.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.RBBP5.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.RBBP9.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.RBM17.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.RBM22.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.RBM3.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ESRP1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ESRP2.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.RBM7.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.RBM8A.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.RBFOX2.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.RBMS1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.RFC2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.RFC3.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.RFXANK.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.RIOK2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.MEX3C.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.RNASEH2C.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group44.done'
	shell:
		'touch {output}'
rule rawTF_group45:
	input:
		'{path}footprints/operations/{mergedsample}.RNF138.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.RPL35.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.RPL6.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.RPP25.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.RPS10.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.RPS4X.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.RPS6KA5.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.RUFY3.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.RUVBL1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.SCAND2P.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.PDS5A.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.SCMH1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.SEMA4A.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.SF1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.SF3B1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.SFT2D1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.SLC18A1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.SMAP2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.SMCR7L.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.SMPX.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group45.done'
	shell:
		'touch {output}'
rule rawTF_group46:
	input:
		'{path}footprints/operations/{mergedsample}.SMUG1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.SNAPC4.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.SNAPC5.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.SND1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.SNRNP70.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.SNRPB2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.SOCS4.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.SOD1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.SOX14.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.SPAG7.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.SPATS2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.SPR.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.SRBD1.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.SRP9.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.SSBP3.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.SSX2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.SSX3.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.STAU2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.STUB1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.SUCLG1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group46.done'
	shell:
		'touch {output}'
rule rawTF_group47:
	input:
		'{path}footprints/operations/{mergedsample}.TAF1A.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.TAF9.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.TAGLN2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.TBPL1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.TCEAL2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.TCEAL6.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.TFAM.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.TGIF2LX.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.THAP5.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.THRA.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.MED30.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.TIA1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.TIMELESS.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.TIMM44.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.TIMM8A.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.TMSB4XP8.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.TOB2.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.TP73.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.TPI1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.TPPP.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group47.done'
	shell:
		'touch {output}'
rule rawTF_group48:
	input:
		'{path}footprints/operations/{mergedsample}.TRIM21.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.TRIM69.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.TRIP10.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.TRMT1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.TROVE2.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.TSC22D4.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.TSN.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.TSNAX.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.TULP1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.U2AF1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.UBB.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.UBE2V1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.UGP2.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.UQCRB.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.USP39.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.UTP18.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.VAMP3.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.EZR.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.VPS4B.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.NELFA.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group48.done'
	shell:
		'touch {output}'
rule rawTF_group49:
	input:
		'{path}footprints/operations/{mergedsample}.WISP2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.XG.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.XRCC1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.YEATS4.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.YWHAE.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.YWHAZ.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB12.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB25.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB43.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB46.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ZC3H7A.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.ZCCHC14.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ZCCHC17.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.ZDHHC15.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.ZDHHC5.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.ZFP3.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.ZHX3.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ZMAT2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ZMAT4.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ZNF124.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group49.done'
	shell:
		'touch {output}'
rule rawTF_group50:
	input:
		'{path}footprints/operations/{mergedsample}.ZNF131.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ZNF160.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ZKSCAN8.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ZSCAN9.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ZNF205.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ZNF207.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB18.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ZNF250.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ZNF26.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ZNF3.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ZNF304.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.RNF114.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ZSCAN31.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.ZNF326.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.ZNF385A.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.ZNF503.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.ZNF510.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ZNF655.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ZNF671.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ZNF695.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group50.done'
	shell:
		'touch {output}'
rule rawTF_group51:
	input:
		'{path}footprints/operations/{mergedsample}.ZNF706.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ZNF71.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ZNF720.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ZNF76.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.ZNF766.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.ZRSR2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.ZSWIM1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.Myf.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.Pax6.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.RORA_1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.RORA_2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.YY1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.TP53.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.RELA.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.ZNF354C.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.MIZF.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.AP1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.DUX4.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.FOXP1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.POU2F2.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group51.done'
	shell:
		'touch {output}'
rule rawTF_group52:
	input:
		'{path}footprints/operations/{mergedsample}.TCF7L2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.TP63.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB33.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ZNF263.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.AR.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.KLF5.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.T.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.EN1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ZNF143.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.NR3C1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ESRRB.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.HOXA5.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.DMRT3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.LBX1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.POU6F1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.BARHL2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.ELF4.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.EN2.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.HOXA13.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.HOXC11.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group52.done'
	shell:
		'touch {output}'
rule rawTF_group53:
	input:
		'{path}footprints/operations/{mergedsample}.ONECUT1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.POU4F2.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB7B.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB7C.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.RHOXF1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.UNCX.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.NR3C2.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.SP8.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.YY2.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB7A.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ZNF410.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.ZNF740.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ONECUT2.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.ONECUT3.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.MYBL1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.MYBL2.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.PAX9.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.PKNOX1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.POU1F1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.POU2F1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group53.done'
	shell:
		'touch {output}'
rule rawTF_group54:
	input:
		'{path}footprints/operations/{mergedsample}.POU3F1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.POU3F3.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.POU3F4.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.POU4F1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.POU5F1B.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.POU6F2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.HOXD12.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.BSX.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.HMBOX1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.HOXA10.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.HOXA2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.HOXB2.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.HOXB3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.HOXC10.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.HOXC12.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.HOXC13.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HOXD11.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.HOXD13.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.NFATC2.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ASCL1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group54.done'
	shell:
		'touch {output}'
rule rawTF_group55:
	input:
		'{path}footprints/operations/{mergedsample}.FOXK2.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.GRHL2.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.KLF9.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.NR2F2.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.POU5F1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.RBPJ.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.SIX1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.SIX2.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.TEAD2.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ZNF24.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ZNF384.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.ZNF282.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ZSCAN4.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.RORB.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.RORC.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.TCF7L1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HINFP1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ZNF238.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ZNF306.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ZNF524.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group55.done'
	shell:
		'touch {output}'
rule rawTF_group56:
	input:
		'{path}footprints/operations/{mergedsample}.ZNF75A.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ZNF784.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.HSFY2.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.NFATC1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.POU2F3.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.POU5F1P1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.BHLHB2.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.BHLHB3.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.CART1.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.HOXA1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.HOXB5.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.HOXD8.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.IRX5.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.PHOX2B.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.RAXL1.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.ESRRG.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.THRB.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.Trp53.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.Trp73.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB49.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group56.done'
	shell:
		'touch {output}'
rule rawTF_group57:
	input:
		'{path}footprints/operations/{mergedsample}.ZNF232.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ZNF435.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ZNF713.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ARID5A.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.BARHL1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.BBX.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.BCL3.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.CHD1.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.CHD2.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.CREB3L2.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.DBX2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.DMC1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.EBF3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.EP300.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.EZH2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.FOXJ1.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.FOXN1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.GMEB1.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.GTF2F1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.GTF2I.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group57.done'
	shell:
		'touch {output}'
rule rawTF_group58:
	input:
		'{path}footprints/operations/{mergedsample}.GZF1.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.HCFC1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.HDX.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.HIVEP1.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.HLX.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.HOXA11.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.HOXA3.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.HOXA4.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.HOXA6.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.HOXA7.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.HOXA9.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.HOXB1.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.HOXB4.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.HOXB6.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.HOXB7.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.HOXB8.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.HOXC4.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.HOXC5.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.HOXC6.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.HOXC8.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group58.done'
	shell:
		'touch {output}'
rule rawTF_group59:
	input:
		'{path}footprints/operations/{mergedsample}.HOXC9.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.HOXD10.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.HOXD1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.HOXD4.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.HOXD9.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.IKZF2.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.IRX4.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.IRX6.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.KLF7.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.LHX1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.LHX5.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.MECOM.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.MTA3.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.OSR1.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.OSR2.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.OTP.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.PATZ1.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.PGR.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.PML.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.PRDM14.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group59.done'
	shell:
		'touch {output}'
rule rawTF_group60:
	input:
		'{path}footprints/operations/{mergedsample}.RAD21.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.RCOR1.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.RFX7.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.RHOXF2.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.SIN3A.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.SIX3.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.SIX4.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.SIX5.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.SIX6.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.SMARCC1.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.SMARCC2.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.SMC3.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.SOX12.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.SOX30.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.SOX6.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.SP100.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.STAT5A.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.STAT5B.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.TAF1.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.TBL1XR1.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group60.done'
	shell:
		'touch {output}'
rule rawTF_group61:
	input:
		'{path}footprints/operations/{mergedsample}.TCF21.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.TFAP2E.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.TFCP2L1.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.TLX2.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.UBP1.rawFPanalysis.bamcopy5.done', 
		'{path}footprints/operations/{mergedsample}.WRNIP1.rawFPanalysis.bamcopy6.done', 
		'{path}footprints/operations/{mergedsample}.YBX1.rawFPanalysis.bamcopy7.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB14.rawFPanalysis.bamcopy8.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB16.rawFPanalysis.bamcopy9.done', 
		'{path}footprints/operations/{mergedsample}.ZBTB3.rawFPanalysis.bamcopy10.done', 
		'{path}footprints/operations/{mergedsample}.ZKSCAN1.rawFPanalysis.bamcopy11.done', 
		'{path}footprints/operations/{mergedsample}.ZKSCAN3.rawFPanalysis.bamcopy12.done', 
		'{path}footprints/operations/{mergedsample}.ZNF148.rawFPanalysis.bamcopy13.done', 
		'{path}footprints/operations/{mergedsample}.ZNF219.rawFPanalysis.bamcopy14.done', 
		'{path}footprints/operations/{mergedsample}.ZNF274.rawFPanalysis.bamcopy15.done', 
		'{path}footprints/operations/{mergedsample}.ZNF281.rawFPanalysis.bamcopy16.done', 
		'{path}footprints/operations/{mergedsample}.ZNF333.rawFPanalysis.bamcopy17.done', 
		'{path}footprints/operations/{mergedsample}.ZNF350.rawFPanalysis.bamcopy18.done', 
		'{path}footprints/operations/{mergedsample}.ZNF35.rawFPanalysis.bamcopy19.done', 
		'{path}footprints/operations/{mergedsample}.ZNF423.rawFPanalysis.bamcopy20.done' 
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group61.done'
	shell:
		'touch {output}'
rule rawTF_group62:
	input:
		'{path}footprints/operations/{mergedsample}.ZNF652.rawFPanalysis.bamcopy1.done', 
		'{path}footprints/operations/{mergedsample}.ZNF691.rawFPanalysis.bamcopy2.done', 
		'{path}footprints/operations/{mergedsample}.ZNF711.rawFPanalysis.bamcopy3.done', 
		'{path}footprints/operations/{mergedsample}.ZNF8.rawFPanalysis.bamcopy4.done', 
		'{path}footprints/operations/{mergedsample}.Sox4.rawFPanalysis.bamcopy5.done'
	output:
		'{path}footprints/operations/{mergedsample}.rawTF.group62.done'
	shell:
		'touch {output}'