########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################

########################################################################################################################################
#### MANUAL ############################################################################################################################
########################################################################################################################################
rule manual_raw_footprints:
    input:
        "lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr01/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr02/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr04/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr05/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr07/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr08/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        #
        "lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.SOX2.rawFPanalysis.done",
        "lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-01-REP1.SOX2.rawFPanalysis.done",
        "lncap/ex01/cr01/operations/footprints/raw/LNCaP-WT-01-REP1.SOX2.rawFPanalysis.done",
        "lncap/ex01/cr02/operations/footprints/raw/LNCaP-WT-01-REP1.SOX2.rawFPanalysis.done",
        "lncap/ex01/cr04/operations/footprints/raw/LNCaP-WT-01-REP1.SOX2.rawFPanalysis.done",
        "lncap/ex01/cr05/operations/footprints/raw/LNCaP-WT-01-REP1.SOX2.rawFPanalysis.done",
        "lncap/ex01/cr07/operations/footprints/raw/LNCaP-WT-01-REP1.SOX2.rawFPanalysis.done",
        "lncap/ex01/cr08/operations/footprints/raw/LNCaP-WT-01-REP1.SOX2.rawFPanalysis.done",
        #
        "lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.FOXM1.rawFPanalysis.done",
        "lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-01-REP1.FOXM1.rawFPanalysis.done",
        "lncap/ex01/cr01/operations/footprints/raw/LNCaP-WT-01-REP1.FOXM1.rawFPanalysis.done",
        "lncap/ex01/cr02/operations/footprints/raw/LNCaP-WT-01-REP1.FOXM1.rawFPanalysis.done",
        "lncap/ex01/cr04/operations/footprints/raw/LNCaP-WT-01-REP1.FOXM1.rawFPanalysis.done",
        "lncap/ex01/cr05/operations/footprints/raw/LNCaP-WT-01-REP1.FOXM1.rawFPanalysis.done",
        "lncap/ex01/cr07/operations/footprints/raw/LNCaP-WT-01-REP1.FOXM1.rawFPanalysis.done",
        "lncap/ex01/cr08/operations/footprints/raw/LNCaP-WT-01-REP1.FOXM1.rawFPanalysis.done",
        #
        "lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.EZH2.rawFPanalysis.done",
        "lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-01-REP1.EZH2.rawFPanalysis.done",
        "lncap/ex01/cr01/operations/footprints/raw/LNCaP-WT-01-REP1.EZH2.rawFPanalysis.done",
        "lncap/ex01/cr02/operations/footprints/raw/LNCaP-WT-01-REP1.EZH2.rawFPanalysis.done",
        "lncap/ex01/cr04/operations/footprints/raw/LNCaP-WT-01-REP1.EZH2.rawFPanalysis.done",
        "lncap/ex01/cr05/operations/footprints/raw/LNCaP-WT-01-REP1.EZH2.rawFPanalysis.done",
        "lncap/ex01/cr07/operations/footprints/raw/LNCaP-WT-01-REP1.EZH2.rawFPanalysis.done",
        "lncap/ex01/cr08/operations/footprints/raw/LNCaP-WT-01-REP1.EZH2.rawFPanalysis.done",
        #
        "lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.MYCN.rawFPanalysis.done",
        "lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-01-REP1.MYCN.rawFPanalysis.done",
        "lncap/ex01/cr01/operations/footprints/raw/LNCaP-WT-01-REP1.MYCN.rawFPanalysis.done",
        "lncap/ex01/cr02/operations/footprints/raw/LNCaP-WT-01-REP1.MYCN.rawFPanalysis.done",
        "lncap/ex01/cr04/operations/footprints/raw/LNCaP-WT-01-REP1.MYCN.rawFPanalysis.done",
        "lncap/ex01/cr05/operations/footprints/raw/LNCaP-WT-01-REP1.MYCN.rawFPanalysis.done",
        "lncap/ex01/cr07/operations/footprints/raw/LNCaP-WT-01-REP1.MYCN.rawFPanalysis.done",
        "lncap/ex01/cr08/operations/footprints/raw/LNCaP-WT-01-REP1.MYCN.rawFPanalysis.done",
        #
        "lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.MYC.rawFPanalysis.done",
        "lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-01-REP1.MYC.rawFPanalysis.done",
        "lncap/ex01/cr01/operations/footprints/raw/LNCaP-WT-01-REP1.MYC.rawFPanalysis.done",
        "lncap/ex01/cr02/operations/footprints/raw/LNCaP-WT-01-REP1.MYC.rawFPanalysis.done",
        "lncap/ex01/cr04/operations/footprints/raw/LNCaP-WT-01-REP1.MYC.rawFPanalysis.done",
        "lncap/ex01/cr05/operations/footprints/raw/LNCaP-WT-01-REP1.MYC.rawFPanalysis.done",
        "lncap/ex01/cr07/operations/footprints/raw/LNCaP-WT-01-REP1.MYC.rawFPanalysis.done",
        "lncap/ex01/cr08/operations/footprints/raw/LNCaP-WT-01-REP1.MYC.rawFPanalysis.done"

rule manual_process_raw_footprints:
    input:
        "lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr01/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr02/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr04/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr05/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr07/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done",
        "lncap/ex01/cr08/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done"

########################################################################################################################################
#### TEST ##############################################################################################################################
########################################################################################################################################
# rule run_raw_footprint_test:
#     input:
#         expand("test/operations/footprints/raw/test-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
#         expand("test/operations/footprints/raw/test-REP2.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
#         expand("test/operations/footprints/raw/test-REP3.{genename}.rawFPanalysis.done", genename=config["geneNames"])

rule PREP_footprinting:
    input:
        "{path}preprocessing/footprint_dirtree.built",
        "{path}operations/preprocessing/{sample}-REP{repnum}.preprocessing.complete"
    output:
        "{path}operations/footprints/{sample}-REP{repnum}.footprint.prep.complete"
    shell:
        "touch {output}"

rule run_raw_footprint_test:
    input:
        "data/test/operations/footprints/test-REP1.footprinting_raw_analysis.complete",
        "data/test/operations/footprints/test-REP2.footprinting_raw_analysis.complete",
        "data/test/operations/footprints/test-REP3.footprinting_raw_analysis.complete"

########################################################################################################################################
#### LNCAP #############################################################################################################################
########################################################################################################################################
rule rawFP_lncap:
    input:
        expand("data/pros/lncap/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("data/pros/lncap/wt02/operations/footprints/raw/LNCaP-WT-02-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("data/pros/lncap/cr01/operations/footprints/raw/LNCaP-CR-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("data/pros/lncap/cr02/operations/footprints/raw/LNCaP-CR-02-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("data/pros/lncap/cr04/operations/footprints/raw/LNCaP-CR-04-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("data/pros/lncap/cr05/operations/footprints/raw/LNCaP-CR-05-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("data/pros/lncap/cr07/operations/footprints/raw/LNCaP-CR-07-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("data/pros/lncap/cr08/operations/footprints/raw/LNCaP-CR-08-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"])
        
