########################################################################################################################################
################################ GENERATE SITES DATABASE ###############################################################################
########################################################################################################################################

# Build the directory structure for the sites database
rule SITES_build_dir_structure:
    # params: -p ignore error if existing, make parent dirs, -v verbose
    output:
        "snakeResources/sites/operations/dirtree.built"
    shell:
        """
        mkdir -p -v snakeResources/sites/operations snakeResources/sites/scripts snakeResources/sites/data snakeResources/sites/genes
        ##
        touch {output}
        """

# Generate the motif data Rdata file
rule generate_motifData:
    input:
        "snakeResources/sites/operations/dirtree.built"
    output:
        "snakeResources/sites/data/motifData.Rdata"
    script:
        "snakeResources/sites/scripts/generateMotifData.R"

# Generate the gene names file
rule generate_geneNames:
    input:
        "snakeResources/sites/data/motifData.Rdata"
    output:
        "snakeResources/sites/data/geneNames.txt"
    script:
        "snakeResources/sites/scripts/generateNames.R"

# Use the PWM information to scan the genome for matches and store the data
rule scanPWM:
    input:
        "snakeResources/sites/data/motifData.Rdata"
    output:
        "snakeResources/sites/operations/genes/{gene}.PWMscan.done"
    script:
        'snakeResources/sites/scripts/scanPWM/snakeScanPWM.R'

####
rule PWMscan_group_aggregator:
    input:
        'snakeResources/sites/operations/groups/PWMscan.group1.done',
        'snakeResources/sites/operations/groups/PWMscan.group2.done',
        'snakeResources/sites/operations/groups/PWMscan.group3.done',
        'snakeResources/sites/operations/groups/PWMscan.group4.done',
        'snakeResources/sites/operations/groups/PWMscan.group5.done',
        'snakeResources/sites/operations/groups/PWMscan.group6.done',
        'snakeResources/sites/operations/groups/PWMscan.group7.done',
        'snakeResources/sites/operations/groups/PWMscan.group8.done',
        'snakeResources/sites/operations/groups/PWMscan.group9.done',
        'snakeResources/sites/operations/groups/PWMscan.group10.done',
        'snakeResources/sites/operations/groups/PWMscan.group11.done',
        'snakeResources/sites/operations/groups/PWMscan.group12.done',
        'snakeResources/sites/operations/groups/PWMscan.group13.done',
        'snakeResources/sites/operations/groups/PWMscan.group14.done',
        'snakeResources/sites/operations/groups/PWMscan.group15.done',
        'snakeResources/sites/operations/groups/PWMscan.group16.done',
        'snakeResources/sites/operations/groups/PWMscan.group17.done',
        'snakeResources/sites/operations/groups/PWMscan.group18.done',
        'snakeResources/sites/operations/groups/PWMscan.group19.done',
        'snakeResources/sites/operations/groups/PWMscan.group20.done',
        'snakeResources/sites/operations/groups/PWMscan.group21.done',
        'snakeResources/sites/operations/groups/PWMscan.group22.done',
        'snakeResources/sites/operations/groups/PWMscan.group23.done',
        'snakeResources/sites/operations/groups/PWMscan.group24.done',
        'snakeResources/sites/operations/groups/PWMscan.group25.done',
        'snakeResources/sites/operations/groups/PWMscan.group26.done',
        'snakeResources/sites/operations/groups/PWMscan.group27.done',
        'snakeResources/sites/operations/groups/PWMscan.group28.done',
        'snakeResources/sites/operations/groups/PWMscan.group29.done',
        'snakeResources/sites/operations/groups/PWMscan.group30.done',
        'snakeResources/sites/operations/groups/PWMscan.group31.done',
        'snakeResources/sites/operations/groups/PWMscan.group32.done',
        'snakeResources/sites/operations/groups/PWMscan.group33.done',
        'snakeResources/sites/operations/groups/PWMscan.group34.done',
        'snakeResources/sites/operations/groups/PWMscan.group35.done',
        'snakeResources/sites/operations/groups/PWMscan.group36.done',
        'snakeResources/sites/operations/groups/PWMscan.group37.done',
        'snakeResources/sites/operations/groups/PWMscan.group38.done',
        'snakeResources/sites/operations/groups/PWMscan.group39.done',
        'snakeResources/sites/operations/groups/PWMscan.group40.done',
        'snakeResources/sites/operations/groups/PWMscan.group41.done',
        'snakeResources/sites/operations/groups/PWMscan.group42.done',
        'snakeResources/sites/operations/groups/PWMscan.group43.done',
        'snakeResources/sites/operations/groups/PWMscan.group44.done',
        'snakeResources/sites/operations/groups/PWMscan.group45.done',
        'snakeResources/sites/operations/groups/PWMscan.group46.done',
        'snakeResources/sites/operations/groups/PWMscan.group47.done',
        'snakeResources/sites/operations/groups/PWMscan.group48.done',
        'snakeResources/sites/operations/groups/PWMscan.group49.done',
        'snakeResources/sites/operations/groups/PWMscan.group50.done',
        'snakeResources/sites/operations/groups/PWMscan.group51.done',
        'snakeResources/sites/operations/groups/PWMscan.group52.done',
        'snakeResources/sites/operations/groups/PWMscan.group53.done',
        'snakeResources/sites/operations/groups/PWMscan.group54.done',
        'snakeResources/sites/operations/groups/PWMscan.group55.done',
        'snakeResources/sites/operations/groups/PWMscan.group56.done',
        'snakeResources/sites/operations/groups/PWMscan.group57.done',
        'snakeResources/sites/operations/groups/PWMscan.group58.done',
        'snakeResources/sites/operations/groups/PWMscan.group59.done',
        'snakeResources/sites/operations/groups/PWMscan.group60.done',
        'snakeResources/sites/operations/groups/PWMscan.group61.done',
        'snakeResources/sites/operations/groups/PWMscan.group62.done'
    output:
        "snakeResources/sites/operations/PWMscan.allgroups.done"
    shell:
        "touch {output}"

####
rule PWMscan_group1:
	input:
		'snakeResources/sites/operations/genes/TFAP2A.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFIL3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HLF.PWMscan.done', 
		'snakeResources/sites/operations/genes/NHLH1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAX.PWMscan.done', 
		'snakeResources/sites/operations/genes/USF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CEBPA.PWMscan.done', 
		'snakeResources/sites/operations/genes/EBF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CEBPB.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOS.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOSL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOSL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/JUN.PWMscan.done', 
		'snakeResources/sites/operations/genes/JUNB.PWMscan.done', 
		'snakeResources/sites/operations/genes/JUND.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAFF.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAFK.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFAP2C.PWMscan.done', 
		'snakeResources/sites/operations/genes/USF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SREBF1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group1.done'
	shell:
		'touch {output}'

rule PWMscan_group2:
	input:
		'snakeResources/sites/operations/genes/SREBF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/AHR.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFAP4.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARNT.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATF6.PWMscan.done', 
		'snakeResources/sites/operations/genes/BACH1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BACH2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CREB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TCF3.PWMscan.done', 
		'snakeResources/sites/operations/genes/XBP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARID5B.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYOD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFE2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYCN.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFE2L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TEF.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATF3.PWMscan.done', 
		'snakeResources/sites/operations/genes/BATF.PWMscan.done', 
		'snakeResources/sites/operations/genes/TCF12.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group2.done'
	shell:
		'touch {output}'

rule PWMscan_group3:
	input:
		'snakeResources/sites/operations/genes/MYC.PWMscan.done', 
		'snakeResources/sites/operations/genes/MXI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHLHE40.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARNTL.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATF4.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATF7.PWMscan.done', 
		'snakeResources/sites/operations/genes/BATF3.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHLHA15.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHLHE41.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHLHE22.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHLHE23.PWMscan.done', 
		'snakeResources/sites/operations/genes/CEBPD.PWMscan.done', 
		'snakeResources/sites/operations/genes/CEBPE.PWMscan.done', 
		'snakeResources/sites/operations/genes/CEBPG.PWMscan.done', 
		'snakeResources/sites/operations/genes/CLOCK.PWMscan.done', 
		'snakeResources/sites/operations/genes/CREB3.PWMscan.done', 
		'snakeResources/sites/operations/genes/CREB3L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/DBP.PWMscan.done', 
		'snakeResources/sites/operations/genes/FIGLA.PWMscan.done', 
		'snakeResources/sites/operations/genes/HES5.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group3.done'
	shell:
		'touch {output}'

rule PWMscan_group4:
	input:
		'snakeResources/sites/operations/genes/HES7.PWMscan.done', 
		'snakeResources/sites/operations/genes/HEY1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HEY2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ID4.PWMscan.done', 
		'snakeResources/sites/operations/genes/JDP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAFG.PWMscan.done', 
		'snakeResources/sites/operations/genes/MESP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MGA.PWMscan.done', 
		'snakeResources/sites/operations/genes/MLX.PWMscan.done', 
		'snakeResources/sites/operations/genes/MLXIPL.PWMscan.done', 
		'snakeResources/sites/operations/genes/MNT.PWMscan.done', 
		'snakeResources/sites/operations/genes/MSC.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYF6.PWMscan.done', 
		'snakeResources/sites/operations/genes/NEUROD2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NEUROG2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NRL.PWMscan.done', 
		'snakeResources/sites/operations/genes/OLIG1.PWMscan.done', 
		'snakeResources/sites/operations/genes/OLIG2.PWMscan.done', 
		'snakeResources/sites/operations/genes/OLIG3.PWMscan.done', 
		'snakeResources/sites/operations/genes/TCF4.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group4.done'
	shell:
		'touch {output}'

rule PWMscan_group5:
	input:
		'snakeResources/sites/operations/genes/TFAP2B.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFE3.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFEB.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFEC.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFAP2D.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARID3A.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARNT2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATF5.PWMscan.done', 
		'snakeResources/sites/operations/genes/CREM.PWMscan.done', 
		'snakeResources/sites/operations/genes/DDIT3.PWMscan.done', 
		'snakeResources/sites/operations/genes/EPAS1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOSB.PWMscan.done', 
		'snakeResources/sites/operations/genes/HAND1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HES1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HIF1A.PWMscan.done', 
		'snakeResources/sites/operations/genes/HMGA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HMGA2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAFA.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAFB.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group5.done'
	shell:
		'touch {output}'

rule PWMscan_group6:
	input:
		'snakeResources/sites/operations/genes/MAF.PWMscan.done', 
		'snakeResources/sites/operations/genes/MITF.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYOG.PWMscan.done', 
		'snakeResources/sites/operations/genes/NEUROD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFE2L2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PTF1A.PWMscan.done', 
		'snakeResources/sites/operations/genes/TAL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TWIST1.PWMscan.done', 
		'snakeResources/sites/operations/genes/AIRE.PWMscan.done', 
		'snakeResources/sites/operations/genes/ALX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ALX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ALX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/ANDR.PWMscan.done', 
		'snakeResources/sites/operations/genes/AP2A.PWMscan.done', 
		'snakeResources/sites/operations/genes/AP2B.PWMscan.done', 
		'snakeResources/sites/operations/genes/AP2C.PWMscan.done', 
		'snakeResources/sites/operations/genes/AP2D.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARI3A.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARI5B.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARX.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group6.done'
	shell:
		'touch {output}'

rule PWMscan_group7:
	input:
		'snakeResources/sites/operations/genes/ASCL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATF6A.PWMscan.done', 
		'snakeResources/sites/operations/genes/ATOH1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BARH1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BARH2.PWMscan.done', 
		'snakeResources/sites/operations/genes/BARX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BARX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/BC11A.PWMscan.done', 
		'snakeResources/sites/operations/genes/BCL6B.PWMscan.done', 
		'snakeResources/sites/operations/genes/BCL6.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHA15.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHE22.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHE23.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHE40.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHE41.PWMscan.done', 
		'snakeResources/sites/operations/genes/BMAL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BPTF.PWMscan.done', 
		'snakeResources/sites/operations/genes/BRAC.PWMscan.done', 
		'snakeResources/sites/operations/genes/BRCA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BSH.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group7.done'
	shell:
		'touch {output}'

rule PWMscan_group8:
	input:
		'snakeResources/sites/operations/genes/CDC5L.PWMscan.done', 
		'snakeResources/sites/operations/genes/CDX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CDX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CEBPZ.PWMscan.done', 
		'snakeResources/sites/operations/genes/CENPB.PWMscan.done', 
		'snakeResources/sites/operations/genes/COE1.PWMscan.done', 
		'snakeResources/sites/operations/genes/COT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/COT2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CPEB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CR3L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CR3L2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CREB5.PWMscan.done', 
		'snakeResources/sites/operations/genes/CRX.PWMscan.done', 
		'snakeResources/sites/operations/genes/CTCFL.PWMscan.done', 
		'snakeResources/sites/operations/genes/CTCF.PWMscan.done', 
		'snakeResources/sites/operations/genes/CUX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CUX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CXXC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/DLX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/DLX2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group8.done'
	shell:
		'touch {output}'

rule PWMscan_group9:
	input:
		'snakeResources/sites/operations/genes/DLX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/DLX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/DLX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/DLX6.PWMscan.done', 
		'snakeResources/sites/operations/genes/DMBX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/DPRX.PWMscan.done', 
		'snakeResources/sites/operations/genes/DRGX.PWMscan.done', 
		'snakeResources/sites/operations/genes/DUXA.PWMscan.done', 
		'snakeResources/sites/operations/genes/E2F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/E2F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/E2F3.PWMscan.done', 
		'snakeResources/sites/operations/genes/E2F4.PWMscan.done', 
		'snakeResources/sites/operations/genes/E2F5.PWMscan.done', 
		'snakeResources/sites/operations/genes/E2F6.PWMscan.done', 
		'snakeResources/sites/operations/genes/E2F7.PWMscan.done', 
		'snakeResources/sites/operations/genes/E2F8.PWMscan.done', 
		'snakeResources/sites/operations/genes/E4F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/EGR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/EGR2.PWMscan.done', 
		'snakeResources/sites/operations/genes/EGR3.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group9.done'
	shell:
		'touch {output}'

rule PWMscan_group10:
	input:
		'snakeResources/sites/operations/genes/EGR4.PWMscan.done', 
		'snakeResources/sites/operations/genes/EHF.PWMscan.done', 
		'snakeResources/sites/operations/genes/ELF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ELF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ELF3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ELF5.PWMscan.done', 
		'snakeResources/sites/operations/genes/ELK1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ELK3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ELK4.PWMscan.done', 
		'snakeResources/sites/operations/genes/EMX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/EMX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/EOMES.PWMscan.done', 
		'snakeResources/sites/operations/genes/ERF.PWMscan.done', 
		'snakeResources/sites/operations/genes/ERG.PWMscan.done', 
		'snakeResources/sites/operations/genes/ERR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ERR2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ERR3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ESR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ESR2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ESX1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group10.done'
	shell:
		'touch {output}'

rule PWMscan_group11:
	input:
		'snakeResources/sites/operations/genes/ETS1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETS2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETV1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETV2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETV3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETV4.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETV5.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETV6.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETV7.PWMscan.done', 
		'snakeResources/sites/operations/genes/EVI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/EVX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/EVX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FEV.PWMscan.done', 
		'snakeResources/sites/operations/genes/FLI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXA2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXA3.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXC2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group11.done'
	shell:
		'touch {output}'

rule PWMscan_group12:
	input:
		'snakeResources/sites/operations/genes/FOXD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXD2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXD3.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXG1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXH1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXJ2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXJ3.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXK1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXM1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXO1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXO3.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXO4.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXO6.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXP3.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXQ1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group12.done'
	shell:
		'touch {output}'

rule PWMscan_group13:
	input:
		'snakeResources/sites/operations/genes/FUBP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GABP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GABPA.PWMscan.done', 
		'snakeResources/sites/operations/genes/GATA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GATA2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GATA3.PWMscan.done', 
		'snakeResources/sites/operations/genes/GATA4.PWMscan.done', 
		'snakeResources/sites/operations/genes/GATA5.PWMscan.done', 
		'snakeResources/sites/operations/genes/GATA6.PWMscan.done', 
		'snakeResources/sites/operations/genes/GBX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GBX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GCM1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GCM2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GCR.PWMscan.done', 
		'snakeResources/sites/operations/genes/GFI1B.PWMscan.done', 
		'snakeResources/sites/operations/genes/GFI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GLI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GLI2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GLI3.PWMscan.done', 
		'snakeResources/sites/operations/genes/GLIS1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group13.done'
	shell:
		'touch {output}'

rule PWMscan_group14:
	input:
		'snakeResources/sites/operations/genes/GLIS2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GLIS3.PWMscan.done', 
		'snakeResources/sites/operations/genes/GMEB2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GRHL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GSC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GSC.PWMscan.done', 
		'snakeResources/sites/operations/genes/GSX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GSX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HBP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HEN1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HESX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HIC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HIC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HINFP.PWMscan.done', 
		'snakeResources/sites/operations/genes/HLTF.PWMscan.done', 
		'snakeResources/sites/operations/genes/HMBX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HME1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HME2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HMX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HMX2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group14.done'
	shell:
		'touch {output}'

rule PWMscan_group15:
	input:
		'snakeResources/sites/operations/genes/HMX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNF1A.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNF1B.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNF4A.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNF4G.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNF6.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOMEZ.PWMscan.done', 
		'snakeResources/sites/operations/genes/HSF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HSF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HSF4.PWMscan.done', 
		'snakeResources/sites/operations/genes/HSFY1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HTF4.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXA10.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXA11.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXA13.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXA2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXA5.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXA7.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXA9.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group15.done'
	shell:
		'touch {output}'

rule PWMscan_group16:
	input:
		'snakeResources/sites/operations/genes/HXB13.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXB2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXB3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXB6.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXB7.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXB8.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXC10.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXC11.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXC12.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXC13.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXC6.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXC8.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXD10.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXD11.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXD12.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXD13.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXD3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXD4.PWMscan.done', 
		'snakeResources/sites/operations/genes/HXD8.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group16.done'
	shell:
		'touch {output}'

rule PWMscan_group17:
	input:
		'snakeResources/sites/operations/genes/HXD9.PWMscan.done', 
		'snakeResources/sites/operations/genes/IKZF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/INSM1.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRF3.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRF4.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRF5.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRF7.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRF8.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRF9.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ISL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ISL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ISX.PWMscan.done', 
		'snakeResources/sites/operations/genes/ITF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/KAISO.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF12.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF13.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group17.done'
	shell:
		'touch {output}'

rule PWMscan_group18:
	input:
		'snakeResources/sites/operations/genes/KLF14.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF15.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF16.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF3.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF4.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF6.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF8.PWMscan.done', 
		'snakeResources/sites/operations/genes/LBX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/LEF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/LHX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/LHX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/LHX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/LHX6.PWMscan.done', 
		'snakeResources/sites/operations/genes/LHX8.PWMscan.done', 
		'snakeResources/sites/operations/genes/LHX9.PWMscan.done', 
		'snakeResources/sites/operations/genes/LMX1A.PWMscan.done', 
		'snakeResources/sites/operations/genes/LMX1B.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAZ.PWMscan.done', 
		'snakeResources/sites/operations/genes/MBD2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group18.done'
	shell:
		'touch {output}'

rule PWMscan_group19:
	input:
		'snakeResources/sites/operations/genes/MCR.PWMscan.done', 
		'snakeResources/sites/operations/genes/MECP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEF2A.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEF2B.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEF2C.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEF2D.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEIS1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEIS2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEIS3.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEOX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEOX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MGAP.PWMscan.done', 
		'snakeResources/sites/operations/genes/MIXL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MLXPL.PWMscan.done', 
		'snakeResources/sites/operations/genes/MNX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MSX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MSX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MTF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MUSC.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYBA.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group19.done'
	shell:
		'touch {output}'

rule PWMscan_group20:
	input:
		'snakeResources/sites/operations/genes/MYBB.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYB.PWMscan.done', 
		'snakeResources/sites/operations/genes/MZF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NANOG.PWMscan.done', 
		'snakeResources/sites/operations/genes/NDF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NDF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NF2L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NF2L2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFAC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFAC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFAC3.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFAC4.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFAT5.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFIA.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFIC.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFKB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFKB2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFYA.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFYB.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFYC.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group20.done'
	shell:
		'touch {output}'

rule PWMscan_group21:
	input:
		'snakeResources/sites/operations/genes/NGN2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX21.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX22.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX23.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX25.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX28.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX31.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX32.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX61.PWMscan.done', 
		'snakeResources/sites/operations/genes/NKX62.PWMscan.done', 
		'snakeResources/sites/operations/genes/NOBOX.PWMscan.done', 
		'snakeResources/sites/operations/genes/NOTO.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR0B1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR1D1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR1H2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR1H4.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR1I2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR1I3.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR2C1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR2C2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group21.done'
	shell:
		'touch {output}'

rule PWMscan_group22:
	input:
		'snakeResources/sites/operations/genes/NR2E1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR2E3.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR2F6.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR4A1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR4A2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR4A3.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR5A2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR6A1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NRF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ONEC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ONEC3.PWMscan.done', 
		'snakeResources/sites/operations/genes/OTX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/OTX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/OVOL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/P53.PWMscan.done', 
		'snakeResources/sites/operations/genes/P5F1B.PWMscan.done', 
		'snakeResources/sites/operations/genes/P63.PWMscan.done', 
		'snakeResources/sites/operations/genes/P73.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAX2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group22.done'
	shell:
		'touch {output}'

rule PWMscan_group23:
	input:
		'snakeResources/sites/operations/genes/PAX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAX6.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAX7.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAX8.PWMscan.done', 
		'snakeResources/sites/operations/genes/PBX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PBX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PBX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/PDX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PEBB.PWMscan.done', 
		'snakeResources/sites/operations/genes/PHX2A.PWMscan.done', 
		'snakeResources/sites/operations/genes/PHX2B.PWMscan.done', 
		'snakeResources/sites/operations/genes/PIT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PITX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PITX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PITX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/PKNX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PKNX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PLAG1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group23.done'
	shell:
		'touch {output}'

rule PWMscan_group24:
	input:
		'snakeResources/sites/operations/genes/PLAL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO2F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO2F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO2F3.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO3F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO3F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO3F3.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO3F4.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO4F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO4F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO4F3.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO5F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO6F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PO6F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PPARA.PWMscan.done', 
		'snakeResources/sites/operations/genes/PPARD.PWMscan.done', 
		'snakeResources/sites/operations/genes/PPARG.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRD14.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRDM1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRDM4.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group24.done'
	shell:
		'touch {output}'

rule PWMscan_group25:
	input:
		'snakeResources/sites/operations/genes/PRGR.PWMscan.done', 
		'snakeResources/sites/operations/genes/PROP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PROX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRRX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRRX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PURA.PWMscan.done', 
		'snakeResources/sites/operations/genes/RARA.PWMscan.done', 
		'snakeResources/sites/operations/genes/RARB.PWMscan.done', 
		'snakeResources/sites/operations/genes/RARG.PWMscan.done', 
		'snakeResources/sites/operations/genes/RAX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/RELB.PWMscan.done', 
		'snakeResources/sites/operations/genes/REL.PWMscan.done', 
		'snakeResources/sites/operations/genes/REST.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/RHXF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RORA.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group25.done'
	shell:
		'touch {output}'

rule PWMscan_group26:
	input:
		'snakeResources/sites/operations/genes/RORG.PWMscan.done', 
		'snakeResources/sites/operations/genes/RREB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RUNX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RUNX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/RUNX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/RXRA.PWMscan.done', 
		'snakeResources/sites/operations/genes/RXRB.PWMscan.done', 
		'snakeResources/sites/operations/genes/RXRG.PWMscan.done', 
		'snakeResources/sites/operations/genes/RX.PWMscan.done', 
		'snakeResources/sites/operations/genes/SCRT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SCRT2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SHOX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SHOX.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMAD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMAD2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMAD3.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMAD4.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMRC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SNAI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SNAI2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group26.done'
	shell:
		'touch {output}'

rule PWMscan_group27:
	input:
		'snakeResources/sites/operations/genes/SOX10.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX11.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX13.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX15.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX17.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX18.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX21.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX7.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX8.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX9.PWMscan.done', 
		'snakeResources/sites/operations/genes/SP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SP3.PWMscan.done', 
		'snakeResources/sites/operations/genes/SP4.PWMscan.done', 
		'snakeResources/sites/operations/genes/SPDEF.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group27.done'
	shell:
		'touch {output}'

rule PWMscan_group28:
	input:
		'snakeResources/sites/operations/genes/SPI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SPIB.PWMscan.done', 
		'snakeResources/sites/operations/genes/SPIC.PWMscan.done', 
		'snakeResources/sites/operations/genes/SPZ1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SRBP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SRBP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SRF.PWMscan.done', 
		'snakeResources/sites/operations/genes/SRY.PWMscan.done', 
		'snakeResources/sites/operations/genes/STA5A.PWMscan.done', 
		'snakeResources/sites/operations/genes/STA5B.PWMscan.done', 
		'snakeResources/sites/operations/genes/STAT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/STAT2.PWMscan.done', 
		'snakeResources/sites/operations/genes/STAT3.PWMscan.done', 
		'snakeResources/sites/operations/genes/STAT4.PWMscan.done', 
		'snakeResources/sites/operations/genes/STAT6.PWMscan.done', 
		'snakeResources/sites/operations/genes/STF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SUH.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBP.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBX15.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group28.done'
	shell:
		'touch {output}'

rule PWMscan_group29:
	input:
		'snakeResources/sites/operations/genes/TBX19.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBX20.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBX21.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/TCF7.PWMscan.done', 
		'snakeResources/sites/operations/genes/TEAD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TEAD3.PWMscan.done', 
		'snakeResources/sites/operations/genes/TEAD4.PWMscan.done', 
		'snakeResources/sites/operations/genes/TF2LX.PWMscan.done', 
		'snakeResources/sites/operations/genes/TF65.PWMscan.done', 
		'snakeResources/sites/operations/genes/TF7L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TF7L2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFCP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFDP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFE2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TGIF1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group29.done'
	shell:
		'touch {output}'

rule PWMscan_group30:
	input:
		'snakeResources/sites/operations/genes/TGIF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/THAP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/THA.PWMscan.done', 
		'snakeResources/sites/operations/genes/THB.PWMscan.done', 
		'snakeResources/sites/operations/genes/TLX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TWST1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TYY1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TYY2.PWMscan.done', 
		'snakeResources/sites/operations/genes/UBIP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/UNC4.PWMscan.done', 
		'snakeResources/sites/operations/genes/VAX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/VAX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/VDR.PWMscan.done', 
		'snakeResources/sites/operations/genes/VENTX.PWMscan.done', 
		'snakeResources/sites/operations/genes/VSX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/VSX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/WT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/YBOX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBED1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBT18.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group30.done'
	shell:
		'touch {output}'

rule PWMscan_group31:
	input:
		'snakeResources/sites/operations/genes/ZBT49.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBT7A.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBT7B.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB4.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB6.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZEB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZEP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZEP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZFHX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZFX.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZIC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZIC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZIC3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZIC4.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZKSC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZKSC3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN143.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN148.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN219.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN232.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group31.done'
	shell:
		'touch {output}'

rule PWMscan_group32:
	input:
		'snakeResources/sites/operations/genes/ZN282.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN333.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN350.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN384.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN410.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN423.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN524.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN589.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN639.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN652.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN713.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN740.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZN784.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZSC16.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZSCA4.PWMscan.done', 
		'snakeResources/sites/operations/genes/ABCF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/A1CF.PWMscan.done', 
		'snakeResources/sites/operations/genes/ACO1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ADARB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/AFF4.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group32.done'
	shell:
		'touch {output}'

rule PWMscan_group33:
	input:
		'snakeResources/sites/operations/genes/AGGF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/AKR1A1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ANXA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ANXA11.PWMscan.done', 
		'snakeResources/sites/operations/genes/APEX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARFGAP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ASCC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ASPSCR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/AVEN.PWMscan.done', 
		'snakeResources/sites/operations/genes/BAD.PWMscan.done', 
		'snakeResources/sites/operations/genes/GPANK1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BAX.PWMscan.done', 
		'snakeResources/sites/operations/genes/BCL11A.PWMscan.done', 
		'snakeResources/sites/operations/genes/BOLL.PWMscan.done', 
		'snakeResources/sites/operations/genes/CELF4.PWMscan.done', 
		'snakeResources/sites/operations/genes/CELF5.PWMscan.done', 
		'snakeResources/sites/operations/genes/CELF6.PWMscan.done', 
		'snakeResources/sites/operations/genes/C19orf25.PWMscan.done', 
		'snakeResources/sites/operations/genes/C19orf40.PWMscan.done', 
		'snakeResources/sites/operations/genes/EXO5.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group33.done'
	shell:
		'touch {output}'

rule PWMscan_group34:
	input:
		'snakeResources/sites/operations/genes/LINC00471.PWMscan.done', 
		'snakeResources/sites/operations/genes/C9orf156.PWMscan.done', 
		'snakeResources/sites/operations/genes/CANX.PWMscan.done', 
		'snakeResources/sites/operations/genes/CAT.PWMscan.done', 
		'snakeResources/sites/operations/genes/CBFA2T2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CBFB.PWMscan.done', 
		'snakeResources/sites/operations/genes/CBX7.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF830.PWMscan.done', 
		'snakeResources/sites/operations/genes/CCDC25.PWMscan.done', 
		'snakeResources/sites/operations/genes/CD59.PWMscan.done', 
		'snakeResources/sites/operations/genes/CDK2AP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/AGAP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CFL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXN3.PWMscan.done', 
		'snakeResources/sites/operations/genes/CKMT1B.PWMscan.done', 
		'snakeResources/sites/operations/genes/CLK1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CNOT6.PWMscan.done', 
		'snakeResources/sites/operations/genes/NELFB.PWMscan.done', 
		'snakeResources/sites/operations/genes/CPSF4.PWMscan.done', 
		'snakeResources/sites/operations/genes/CSNK2B.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group34.done'
	shell:
		'touch {output}'

rule PWMscan_group35:
	input:
		'snakeResources/sites/operations/genes/CSTF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CYB5R1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CYCS.PWMscan.done', 
		'snakeResources/sites/operations/genes/DAB2.PWMscan.done', 
		'snakeResources/sites/operations/genes/DAZAP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ASAP3.PWMscan.done', 
		'snakeResources/sites/operations/genes/DDX20.PWMscan.done', 
		'snakeResources/sites/operations/genes/DDX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/DDX43.PWMscan.done', 
		'snakeResources/sites/operations/genes/DDX53.PWMscan.done', 
		'snakeResources/sites/operations/genes/DGCR8.PWMscan.done', 
		'snakeResources/sites/operations/genes/DHX36.PWMscan.done', 
		'snakeResources/sites/operations/genes/DIABLO.PWMscan.done', 
		'snakeResources/sites/operations/genes/DIS3.PWMscan.done', 
		'snakeResources/sites/operations/genes/DNMT3A.PWMscan.done', 
		'snakeResources/sites/operations/genes/DTL.PWMscan.done', 
		'snakeResources/sites/operations/genes/DUS3L.PWMscan.done', 
		'snakeResources/sites/operations/genes/DUSP22.PWMscan.done', 
		'snakeResources/sites/operations/genes/DUSP26.PWMscan.done', 
		'snakeResources/sites/operations/genes/ECSIT.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group35.done'
	shell:
		'touch {output}'

rule PWMscan_group36:
	input:
		'snakeResources/sites/operations/genes/EDN1.PWMscan.done', 
		'snakeResources/sites/operations/genes/EEF1D.PWMscan.done', 
		'snakeResources/sites/operations/genes/EIF5A2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ENO1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ESRRA.PWMscan.done', 
		'snakeResources/sites/operations/genes/ETFB.PWMscan.done', 
		'snakeResources/sites/operations/genes/EWSR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/EXOSC3.PWMscan.done', 
		'snakeResources/sites/operations/genes/METTL21B.PWMscan.done', 
		'snakeResources/sites/operations/genes/FAM127B.PWMscan.done', 
		'snakeResources/sites/operations/genes/FEZ1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FEZF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FGF19.PWMscan.done', 
		'snakeResources/sites/operations/genes/FHL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FIP1L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SRRM3.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXP4.PWMscan.done', 
		'snakeResources/sites/operations/genes/GADD45A.PWMscan.done', 
		'snakeResources/sites/operations/genes/GIT2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GLYCTK.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group36.done'
	shell:
		'touch {output}'

rule PWMscan_group37:
	input:
		'snakeResources/sites/operations/genes/GOT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GPAM.PWMscan.done', 
		'snakeResources/sites/operations/genes/GPD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GRHPR.PWMscan.done', 
		'snakeResources/sites/operations/genes/GTF2B.PWMscan.done', 
		'snakeResources/sites/operations/genes/GTF2H3.PWMscan.done', 
		'snakeResources/sites/operations/genes/GTF3C2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GTF3C5.PWMscan.done', 
		'snakeResources/sites/operations/genes/GTPBP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GTPBP6.PWMscan.done', 
		'snakeResources/sites/operations/genes/H1FX.PWMscan.done', 
		'snakeResources/sites/operations/genes/H2AFY.PWMscan.done', 
		'snakeResources/sites/operations/genes/H2AFZ.PWMscan.done', 
		'snakeResources/sites/operations/genes/HCFC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HCLS1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HDAC8.PWMscan.done', 
		'snakeResources/sites/operations/genes/HHAT.PWMscan.done', 
		'snakeResources/sites/operations/genes/HHEX.PWMscan.done', 
		'snakeResources/sites/operations/genes/UBE2K.PWMscan.done', 
		'snakeResources/sites/operations/genes/HIRIP3.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group37.done'
	shell:
		'touch {output}'

rule PWMscan_group38:
	input:
		'snakeResources/sites/operations/genes/HIST1H2BN.PWMscan.done', 
		'snakeResources/sites/operations/genes/HIST2H2AB.PWMscan.done', 
		'snakeResources/sites/operations/genes/HIST2H2BE.PWMscan.done', 
		'snakeResources/sites/operations/genes/HLCS.PWMscan.done', 
		'snakeResources/sites/operations/genes/HMG20A.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNRNPA0.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNRNPA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNRNPC.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNRNPH3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HNRNPLL.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB13.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB9.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HP1BP3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HSPA1L.PWMscan.done', 
		'snakeResources/sites/operations/genes/HSPA5.PWMscan.done', 
		'snakeResources/sites/operations/genes/HTATIP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ID2.PWMscan.done', 
		'snakeResources/sites/operations/genes/IL24.PWMscan.done', 
		'snakeResources/sites/operations/genes/ING3.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group38.done'
	shell:
		'touch {output}'

rule PWMscan_group39:
	input:
		'snakeResources/sites/operations/genes/IRF6.PWMscan.done', 
		'snakeResources/sites/operations/genes/IVD.PWMscan.done', 
		'snakeResources/sites/operations/genes/KDM5A.PWMscan.done', 
		'snakeResources/sites/operations/genes/KDM5D.PWMscan.done', 
		'snakeResources/sites/operations/genes/KCNIP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/KIAA0907.PWMscan.done', 
		'snakeResources/sites/operations/genes/KIF22.PWMscan.done', 
		'snakeResources/sites/operations/genes/LARP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/LARP4.PWMscan.done', 
		'snakeResources/sites/operations/genes/LAS1L.PWMscan.done', 
		'snakeResources/sites/operations/genes/CERS4.PWMscan.done', 
		'snakeResources/sites/operations/genes/UBXN1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CBX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/LRRFIP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/LSM6.PWMscan.done', 
		'snakeResources/sites/operations/genes/LUZP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/LUZP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAGEA8.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAGED4B.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAGEF1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group39.done'
	shell:
		'touch {output}'

rule PWMscan_group40:
	input:
		'snakeResources/sites/operations/genes/MAGOH.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAP4K2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MAPK1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MBTPS2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MCTP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MDM2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GLTPD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBM42.PWMscan.done', 
		'snakeResources/sites/operations/genes/WDR83.PWMscan.done', 
		'snakeResources/sites/operations/genes/MORN1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MRPL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MRPL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MRPS25.PWMscan.done', 
		'snakeResources/sites/operations/genes/MSI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MSI2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MSRA.PWMscan.done', 
		'snakeResources/sites/operations/genes/MSRB3.PWMscan.done', 
		'snakeResources/sites/operations/genes/MTHFD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MXD4.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYEF2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group40.done'
	shell:
		'touch {output}'

rule PWMscan_group41:
	input:
		'snakeResources/sites/operations/genes/MYLK.PWMscan.done', 
		'snakeResources/sites/operations/genes/NANOS1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NAP1L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NCALD.PWMscan.done', 
		'snakeResources/sites/operations/genes/NCBP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFATC3.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFATC4.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFIB.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFIX.PWMscan.done', 
		'snakeResources/sites/operations/genes/NME1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NMI.PWMscan.done', 
		'snakeResources/sites/operations/genes/NMRAL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NNT.PWMscan.done', 
		'snakeResources/sites/operations/genes/NOC2L.PWMscan.done', 
		'snakeResources/sites/operations/genes/GAR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NONO.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR2F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NUCB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/NUP107.PWMscan.done', 
		'snakeResources/sites/operations/genes/NUP133.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group41.done'
	shell:
		'touch {output}'

rule PWMscan_group42:
	input:
		'snakeResources/sites/operations/genes/NXPH3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ODC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/OTUD4.PWMscan.done', 
		'snakeResources/sites/operations/genes/P4HB.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAXIP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PCK2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PDCD11.PWMscan.done', 
		'snakeResources/sites/operations/genes/PDE6H.PWMscan.done', 
		'snakeResources/sites/operations/genes/PDLIM5.PWMscan.done', 
		'snakeResources/sites/operations/genes/PGAM2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PHLDA2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PHOX2A.PWMscan.done', 
		'snakeResources/sites/operations/genes/PHTF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PICK1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PIK3C3.PWMscan.done', 
		'snakeResources/sites/operations/genes/PIR.PWMscan.done', 
		'snakeResources/sites/operations/genes/PKM.PWMscan.done', 
		'snakeResources/sites/operations/genes/PKNOX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PLAGL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PLG.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group42.done'
	shell:
		'touch {output}'

rule PWMscan_group43:
	input:
		'snakeResources/sites/operations/genes/POLE3.PWMscan.done', 
		'snakeResources/sites/operations/genes/POLI.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU3F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU4F3.PWMscan.done', 
		'snakeResources/sites/operations/genes/PPP1R10.PWMscan.done', 
		'snakeResources/sites/operations/genes/PPP2R3B.PWMscan.done', 
		'snakeResources/sites/operations/genes/PPP5C.PWMscan.done', 
		'snakeResources/sites/operations/genes/PQBP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRDX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRKRIR.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRNP.PWMscan.done', 
		'snakeResources/sites/operations/genes/PSMA6.PWMscan.done', 
		'snakeResources/sites/operations/genes/PSMC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PTCD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PTPMT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PURG.PWMscan.done', 
		'snakeResources/sites/operations/genes/R3HDM2.PWMscan.done', 
		'snakeResources/sites/operations/genes/RAB14.PWMscan.done', 
		'snakeResources/sites/operations/genes/RAB18.PWMscan.done', 
		'snakeResources/sites/operations/genes/RAB2A.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group43.done'
	shell:
		'touch {output}'

rule PWMscan_group44:
	input:
		'snakeResources/sites/operations/genes/RAB7A.PWMscan.done', 
		'snakeResources/sites/operations/genes/RAN.PWMscan.done', 
		'snakeResources/sites/operations/genes/RAX.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBBP5.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBBP9.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBM17.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBM22.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBM3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ESRP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ESRP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBM7.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBM8A.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBFOX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBMS1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFC3.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFXANK.PWMscan.done', 
		'snakeResources/sites/operations/genes/RIOK2.PWMscan.done', 
		'snakeResources/sites/operations/genes/MEX3C.PWMscan.done', 
		'snakeResources/sites/operations/genes/RNASEH2C.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group44.done'
	shell:
		'touch {output}'

rule PWMscan_group45:
	input:
		'snakeResources/sites/operations/genes/RNF138.PWMscan.done', 
		'snakeResources/sites/operations/genes/RPL35.PWMscan.done', 
		'snakeResources/sites/operations/genes/RPL6.PWMscan.done', 
		'snakeResources/sites/operations/genes/RPP25.PWMscan.done', 
		'snakeResources/sites/operations/genes/RPS10.PWMscan.done', 
		'snakeResources/sites/operations/genes/RPS4X.PWMscan.done', 
		'snakeResources/sites/operations/genes/RPS6KA5.PWMscan.done', 
		'snakeResources/sites/operations/genes/RUFY3.PWMscan.done', 
		'snakeResources/sites/operations/genes/RUVBL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SCAND2P.PWMscan.done', 
		'snakeResources/sites/operations/genes/PDS5A.PWMscan.done', 
		'snakeResources/sites/operations/genes/SCMH1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SEMA4A.PWMscan.done', 
		'snakeResources/sites/operations/genes/SF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SF3B1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SFT2D1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SLC18A1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMAP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMCR7L.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMPX.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group45.done'
	shell:
		'touch {output}'

rule PWMscan_group46:
	input:
		'snakeResources/sites/operations/genes/SMUG1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SNAPC4.PWMscan.done', 
		'snakeResources/sites/operations/genes/SNAPC5.PWMscan.done', 
		'snakeResources/sites/operations/genes/SND1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SNRNP70.PWMscan.done', 
		'snakeResources/sites/operations/genes/SNRPB2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOCS4.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX14.PWMscan.done', 
		'snakeResources/sites/operations/genes/SPAG7.PWMscan.done', 
		'snakeResources/sites/operations/genes/SPATS2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SPR.PWMscan.done', 
		'snakeResources/sites/operations/genes/SRBD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SRP9.PWMscan.done', 
		'snakeResources/sites/operations/genes/SSBP3.PWMscan.done', 
		'snakeResources/sites/operations/genes/SSX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SSX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/STAU2.PWMscan.done', 
		'snakeResources/sites/operations/genes/STUB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SUCLG1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group46.done'
	shell:
		'touch {output}'

rule PWMscan_group47:
	input:
		'snakeResources/sites/operations/genes/TAF1A.PWMscan.done', 
		'snakeResources/sites/operations/genes/TAF9.PWMscan.done', 
		'snakeResources/sites/operations/genes/TAGLN2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBPL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TCEAL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TCEAL6.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFAM.PWMscan.done', 
		'snakeResources/sites/operations/genes/TGIF2LX.PWMscan.done', 
		'snakeResources/sites/operations/genes/THAP5.PWMscan.done', 
		'snakeResources/sites/operations/genes/THRA.PWMscan.done', 
		'snakeResources/sites/operations/genes/MED30.PWMscan.done', 
		'snakeResources/sites/operations/genes/TIA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TIMELESS.PWMscan.done', 
		'snakeResources/sites/operations/genes/TIMM44.PWMscan.done', 
		'snakeResources/sites/operations/genes/TIMM8A.PWMscan.done', 
		'snakeResources/sites/operations/genes/TMSB4XP8.PWMscan.done', 
		'snakeResources/sites/operations/genes/TOB2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TP73.PWMscan.done', 
		'snakeResources/sites/operations/genes/TPI1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TPPP.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group47.done'
	shell:
		'touch {output}'

rule PWMscan_group48:
	input:
		'snakeResources/sites/operations/genes/TRIM21.PWMscan.done', 
		'snakeResources/sites/operations/genes/TRIM69.PWMscan.done', 
		'snakeResources/sites/operations/genes/TRIP10.PWMscan.done', 
		'snakeResources/sites/operations/genes/TRMT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TROVE2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TSC22D4.PWMscan.done', 
		'snakeResources/sites/operations/genes/TSN.PWMscan.done', 
		'snakeResources/sites/operations/genes/TSNAX.PWMscan.done', 
		'snakeResources/sites/operations/genes/TULP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/U2AF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/UBB.PWMscan.done', 
		'snakeResources/sites/operations/genes/UBE2V1.PWMscan.done', 
		'snakeResources/sites/operations/genes/UGP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/UQCRB.PWMscan.done', 
		'snakeResources/sites/operations/genes/USP39.PWMscan.done', 
		'snakeResources/sites/operations/genes/UTP18.PWMscan.done', 
		'snakeResources/sites/operations/genes/VAMP3.PWMscan.done', 
		'snakeResources/sites/operations/genes/EZR.PWMscan.done', 
		'snakeResources/sites/operations/genes/VPS4B.PWMscan.done', 
		'snakeResources/sites/operations/genes/NELFA.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group48.done'
	shell:
		'touch {output}'

rule PWMscan_group49:
	input:
		'snakeResources/sites/operations/genes/WISP2.PWMscan.done', 
		'snakeResources/sites/operations/genes/XG.PWMscan.done', 
		'snakeResources/sites/operations/genes/XRCC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/YEATS4.PWMscan.done', 
		'snakeResources/sites/operations/genes/YWHAE.PWMscan.done', 
		'snakeResources/sites/operations/genes/YWHAZ.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB12.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB25.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB43.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB46.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZC3H7A.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZCCHC14.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZCCHC17.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZDHHC15.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZDHHC5.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZFP3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZHX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZMAT2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZMAT4.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF124.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group49.done'
	shell:
		'touch {output}'

rule PWMscan_group50:
	input:
		'snakeResources/sites/operations/genes/ZNF131.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF160.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZKSCAN8.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZSCAN9.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF205.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF207.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB18.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF250.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF26.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF304.PWMscan.done', 
		'snakeResources/sites/operations/genes/RNF114.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZSCAN31.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF326.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF385A.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF503.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF510.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF655.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF671.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF695.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group50.done'
	shell:
		'touch {output}'

rule PWMscan_group51:
	input:
		'snakeResources/sites/operations/genes/ZNF706.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF71.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF720.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF76.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF766.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZRSR2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZSWIM1.PWMscan.done', 
		'snakeResources/sites/operations/genes/Myf.PWMscan.done', 
		'snakeResources/sites/operations/genes/Pax6.PWMscan.done', 
		'snakeResources/sites/operations/genes/RORA_1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RORA_2.PWMscan.done', 
		'snakeResources/sites/operations/genes/YY1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TP53.PWMscan.done', 
		'snakeResources/sites/operations/genes/RELA.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF354C.PWMscan.done', 
		'snakeResources/sites/operations/genes/MIZF.PWMscan.done', 
		'snakeResources/sites/operations/genes/AP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/DUX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU2F2.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group51.done'
	shell:
		'touch {output}'

rule PWMscan_group52:
	input:
		'snakeResources/sites/operations/genes/TCF7L2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TP63.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB33.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF263.PWMscan.done', 
		'snakeResources/sites/operations/genes/AR.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF5.PWMscan.done', 
		'snakeResources/sites/operations/genes/T.PWMscan.done', 
		'snakeResources/sites/operations/genes/EN1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF143.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR3C1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ESRRB.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA5.PWMscan.done', 
		'snakeResources/sites/operations/genes/DMRT3.PWMscan.done', 
		'snakeResources/sites/operations/genes/LBX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU6F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BARHL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ELF4.PWMscan.done', 
		'snakeResources/sites/operations/genes/EN2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA13.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXC11.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group52.done'
	shell:
		'touch {output}'

rule PWMscan_group53:
	input:
		'snakeResources/sites/operations/genes/ONECUT1.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU4F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB7B.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB7C.PWMscan.done', 
		'snakeResources/sites/operations/genes/RHOXF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/UNCX.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR3C2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SP8.PWMscan.done', 
		'snakeResources/sites/operations/genes/YY2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB7A.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF410.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF740.PWMscan.done', 
		'snakeResources/sites/operations/genes/ONECUT2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ONECUT3.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYBL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/MYBL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/PAX9.PWMscan.done', 
		'snakeResources/sites/operations/genes/PKNOX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU1F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU2F1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group53.done'
	shell:
		'touch {output}'

rule PWMscan_group54:
	input:
		'snakeResources/sites/operations/genes/POU3F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU3F3.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU3F4.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU4F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU5F1B.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU6F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD12.PWMscan.done', 
		'snakeResources/sites/operations/genes/BSX.PWMscan.done', 
		'snakeResources/sites/operations/genes/HMBOX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA10.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB2.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXC10.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXC12.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXC13.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD11.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD13.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFATC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ASCL1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group54.done'
	shell:
		'touch {output}'

rule PWMscan_group55:
	input:
		'snakeResources/sites/operations/genes/FOXK2.PWMscan.done', 
		'snakeResources/sites/operations/genes/GRHL2.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF9.PWMscan.done', 
		'snakeResources/sites/operations/genes/NR2F2.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU5F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RBPJ.PWMscan.done', 
		'snakeResources/sites/operations/genes/SIX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SIX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/TEAD2.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF24.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF384.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF282.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZSCAN4.PWMscan.done', 
		'snakeResources/sites/operations/genes/RORB.PWMscan.done', 
		'snakeResources/sites/operations/genes/RORC.PWMscan.done', 
		'snakeResources/sites/operations/genes/TCF7L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HINFP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF238.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF306.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF524.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group55.done'
	shell:
		'touch {output}'

rule PWMscan_group56:
	input:
		'snakeResources/sites/operations/genes/ZNF75A.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF784.PWMscan.done', 
		'snakeResources/sites/operations/genes/HSFY2.PWMscan.done', 
		'snakeResources/sites/operations/genes/NFATC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU2F3.PWMscan.done', 
		'snakeResources/sites/operations/genes/POU5F1P1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHLHB2.PWMscan.done', 
		'snakeResources/sites/operations/genes/BHLHB3.PWMscan.done', 
		'snakeResources/sites/operations/genes/CART1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB5.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD8.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/PHOX2B.PWMscan.done', 
		'snakeResources/sites/operations/genes/RAXL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ESRRG.PWMscan.done', 
		'snakeResources/sites/operations/genes/THRB.PWMscan.done', 
		'snakeResources/sites/operations/genes/Trp53.PWMscan.done', 
		'snakeResources/sites/operations/genes/Trp73.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB49.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group56.done'
	shell:
		'touch {output}'

rule PWMscan_group57:
	input:
		'snakeResources/sites/operations/genes/ZNF232.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF435.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF713.PWMscan.done', 
		'snakeResources/sites/operations/genes/ARID5A.PWMscan.done', 
		'snakeResources/sites/operations/genes/BARHL1.PWMscan.done', 
		'snakeResources/sites/operations/genes/BBX.PWMscan.done', 
		'snakeResources/sites/operations/genes/BCL3.PWMscan.done', 
		'snakeResources/sites/operations/genes/CHD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/CHD2.PWMscan.done', 
		'snakeResources/sites/operations/genes/CREB3L2.PWMscan.done', 
		'snakeResources/sites/operations/genes/DBX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/DMC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/EBF3.PWMscan.done', 
		'snakeResources/sites/operations/genes/EP300.PWMscan.done', 
		'snakeResources/sites/operations/genes/EZH2.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXJ1.PWMscan.done', 
		'snakeResources/sites/operations/genes/FOXN1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GMEB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GTF2F1.PWMscan.done', 
		'snakeResources/sites/operations/genes/GTF2I.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group57.done'
	shell:
		'touch {output}'

rule PWMscan_group58:
	input:
		'snakeResources/sites/operations/genes/GZF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HCFC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HDX.PWMscan.done', 
		'snakeResources/sites/operations/genes/HIVEP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HLX.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA11.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA3.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA4.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA6.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA7.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXA9.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB4.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB6.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB7.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXB8.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXC4.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXC5.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXC6.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXC8.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group58.done'
	shell:
		'touch {output}'

rule PWMscan_group59:
	input:
		'snakeResources/sites/operations/genes/HOXC9.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD10.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD1.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD4.PWMscan.done', 
		'snakeResources/sites/operations/genes/HOXD9.PWMscan.done', 
		'snakeResources/sites/operations/genes/IKZF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/IRX6.PWMscan.done', 
		'snakeResources/sites/operations/genes/KLF7.PWMscan.done', 
		'snakeResources/sites/operations/genes/LHX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/LHX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/MECOM.PWMscan.done', 
		'snakeResources/sites/operations/genes/MTA3.PWMscan.done', 
		'snakeResources/sites/operations/genes/OSR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/OSR2.PWMscan.done', 
		'snakeResources/sites/operations/genes/OTP.PWMscan.done', 
		'snakeResources/sites/operations/genes/PATZ1.PWMscan.done', 
		'snakeResources/sites/operations/genes/PGR.PWMscan.done', 
		'snakeResources/sites/operations/genes/PML.PWMscan.done', 
		'snakeResources/sites/operations/genes/PRDM14.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group59.done'
	shell:
		'touch {output}'

rule PWMscan_group60:
	input:
		'snakeResources/sites/operations/genes/RAD21.PWMscan.done', 
		'snakeResources/sites/operations/genes/RCOR1.PWMscan.done', 
		'snakeResources/sites/operations/genes/RFX7.PWMscan.done', 
		'snakeResources/sites/operations/genes/RHOXF2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SIN3A.PWMscan.done', 
		'snakeResources/sites/operations/genes/SIX3.PWMscan.done', 
		'snakeResources/sites/operations/genes/SIX4.PWMscan.done', 
		'snakeResources/sites/operations/genes/SIX5.PWMscan.done', 
		'snakeResources/sites/operations/genes/SIX6.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMARCC1.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMARCC2.PWMscan.done', 
		'snakeResources/sites/operations/genes/SMC3.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX12.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX30.PWMscan.done', 
		'snakeResources/sites/operations/genes/SOX6.PWMscan.done', 
		'snakeResources/sites/operations/genes/SP100.PWMscan.done', 
		'snakeResources/sites/operations/genes/STAT5A.PWMscan.done', 
		'snakeResources/sites/operations/genes/STAT5B.PWMscan.done', 
		'snakeResources/sites/operations/genes/TAF1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TBL1XR1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group60.done'
	shell:
		'touch {output}'

rule PWMscan_group61:
	input:
		'snakeResources/sites/operations/genes/TCF21.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFAP2E.PWMscan.done', 
		'snakeResources/sites/operations/genes/TFCP2L1.PWMscan.done', 
		'snakeResources/sites/operations/genes/TLX2.PWMscan.done', 
		'snakeResources/sites/operations/genes/UBP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/WRNIP1.PWMscan.done', 
		'snakeResources/sites/operations/genes/YBX1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB14.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB16.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZBTB3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZKSCAN1.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZKSCAN3.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF148.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF219.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF274.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF281.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF333.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF350.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF35.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF423.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group61.done'
	shell:
		'touch {output}'

rule PWMscan_group62:
	input:
		'snakeResources/sites/operations/genes/ZNF652.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF691.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF711.PWMscan.done', 
		'snakeResources/sites/operations/genes/ZNF8.PWMscan.done', 
		'snakeResources/sites/operations/genes/Sox4.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group62.done'
	shell:
		'touch {output}'
