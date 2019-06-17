########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################
# Spool the footprinting as below. Preprocessing will be completed first, if not already done
# for i in {1..55}; do snakemake -j 20 run_raw_footprint_test --resources hg38align=1 --config group=$i; done

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
        expand("mdst8/wt01/operations/footprints/groups/raw/MDST8-WT-01-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("mdst8/wt01/operations/footprints/groups/raw/MDST8-WT-01-REP2.rawFPanalysis.group{param}.done", param=config["group"])

rule testFP:
    input:
        "mdst8/wt01/operations/preprocessing/MDST8-WT-01-REP1.preprocessing.complete",
        expand("mdst8/wt01/operations/footprints/raw/MDST8-WT-01-REP1.{genename}.rawFPanalysis.bamcopy1.done", genename=config["geneNames"]),



