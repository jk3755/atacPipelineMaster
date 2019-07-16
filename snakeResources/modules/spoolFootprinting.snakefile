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
        "lncap/ex01/cr08/operations/footprints/raw/LNCaP-WT-01-REP1.AR.rawFPanalysis.done"

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
rule run_raw_footprint_test:
    input:
        expand("test/operations/footprints/raw/test-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("test/operations/footprints/raw/test-REP2.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("test/operations/footprints/raw/test-REP3.{genename}.rawFPanalysis.done", genename=config["geneNames"])

########################################################################################################################################
#### LNCAP #############################################################################################################################
########################################################################################################################################
rule rawFP_lncap_ex01:
    input:
        expand("lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-02-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr01/operations/footprints/raw/LNCaP-CR-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr02/operations/footprints/raw/LNCaP-CR-02-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr04/operations/footprints/raw/LNCaP-CR-04-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr05/operations/footprints/raw/LNCaP-CR-05-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr07/operations/footprints/raw/LNCaP-CR-07-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr08/operations/footprints/raw/LNCaP-CR-08-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"])
        
