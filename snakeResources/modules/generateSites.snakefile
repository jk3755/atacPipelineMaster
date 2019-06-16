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
		'sites/operations/scans/TFAP2A.PWMscan.done', 
		'sites/operations/scans/NFIL3.PWMscan.done', 
		'sites/operations/scans/HLF.PWMscan.done', 
		'sites/operations/scans/NHLH1.PWMscan.done', 
		'sites/operations/scans/MAX.PWMscan.done', 
		'sites/operations/scans/USF1.PWMscan.done', 
		'sites/operations/scans/CEBPA.PWMscan.done', 
		'sites/operations/scans/EBF1.PWMscan.done', 
		'sites/operations/scans/CEBPB.PWMscan.done', 
		'sites/operations/scans/FOS.PWMscan.done', 
		'sites/operations/scans/FOSL1.PWMscan.done', 
		'sites/operations/scans/FOSL2.PWMscan.done', 
		'sites/operations/scans/JUN.PWMscan.done', 
		'sites/operations/scans/JUNB.PWMscan.done', 
		'sites/operations/scans/JUND.PWMscan.done', 
		'sites/operations/scans/MAFF.PWMscan.done', 
		'sites/operations/scans/MAFK.PWMscan.done', 
		'sites/operations/scans/TFAP2C.PWMscan.done', 
		'sites/operations/scans/USF2.PWMscan.done', 
		'sites/operations/scans/SREBF1.PWMscan.done'
	output:
		'snakeResources/sites/operations/groups/PWMscan.group1.done'
	shell:
		'touch {output}'

