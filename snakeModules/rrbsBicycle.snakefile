####################################################################################################
##### METHYLATION ANALYSIS WITH BICYCLE PIPELINE ###################################################
####################################################################################################
# hg38 reference sequence in .fa (fasta) format can be found here:
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/



rule generate_motifData:
    output:
        "sites/motifData.Rdata"
    script:
        "scripts/scanPWM/generateMotifData.R"

rule generate_geneNames:
    output:
        "sites/geneNames.txt"
    script:
        "scripts/scanPWM/generateNames.R"
