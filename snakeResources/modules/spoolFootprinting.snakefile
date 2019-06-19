########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################

########################################################################################################################################
#### TEST ##############################################################################################################################
########################################################################################################################################
# rule run_raw_footprint_test:
#     input:
#         "test/operations/preprocessing/test-REP1.preprocessing.complete",
#         expand("test/operations/footprints/groups/raw/test-REP1.rawFPanalysis.group{param}.done", param=config["group"])

########################################################################################################################################
#### MDST8 #############################################################################################################################
########################################################################################################################################
rule rawFP_mdst8_wt01:
    input:
        "mdst8/wt01/operations/preprocessing/MDST8-WT-01-REP1.preprocessing.complete",
        "mdst8/wt01/operations/preprocessing/MDST8-WT-01-REP2.preprocessing.complete",
        expand("mdst8/wt01/operations/footprints/raw/MDST8-WT-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("mdst8/wt01/operations/footprints/raw/MDST8-WT-01-REP2.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("mdst8/wt01/operations/footprints/temp/MDST8-WT-01-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
        expand("mdst8/wt01/operations/footprints/temp/MDST8-WT-01-REP2.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"])