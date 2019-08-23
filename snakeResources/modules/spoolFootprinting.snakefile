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
