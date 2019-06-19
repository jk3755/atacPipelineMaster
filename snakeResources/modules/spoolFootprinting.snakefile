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

########################################################################################################################################
#### LNCAP #############################################################################################################################
########################################################################################################################################
rule rawFP_lncap_ex01:
    input:
        "lncap/ex01/operations/preprocessing/LNCaP-WT-01-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-WT-02-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-01-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-02-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-04-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-05-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-07-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-08-REP1.preprocessing.complete",
        expand("lncap/ex01/operations/footprints/raw/LNCaP-WT-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/operations/footprints/raw/LNCaP-WT-02-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/operations/footprints/raw/LNCaP-CR-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/operations/footprints/raw/LNCaP-CR-02-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/operations/footprints/raw/LNCaP-CR-04-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/operations/footprints/raw/LNCaP-CR-05-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/operations/footprints/raw/LNCaP-CR-07-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/operations/footprints/raw/LNCaP-CR-08-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        # expand("lncap/ex01/operations/footprints/temp/LNCaP-WT-01-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
        # expand("lncap/ex01/operations/footprints/temp/LNCaP-WT-02-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
        # expand("lncap/ex01/operations/footprints/temp/LNCaP-CR-01-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
        # expand("lncap/ex01/operations/footprints/temp/LNCaP-CR-02-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
        # expand("lncap/ex01/operations/footprints/temp/LNCaP-CR-04-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
        # expand("lncap/ex01/operations/footprints/temp/LNCaP-CR-05-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
        # expand("lncap/ex01/operations/footprints/temp/LNCaP-CR-07-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
        # expand("lncap/ex01/operations/footprints/temp/LNCaP-CR-08-REP1.{genename}.rawFPanalysis.large.done", genename=config["geneNamesLarge"]),
