########################################################################################################################################
#### TEST DATA #########################################################################################################################
########################################################################################################################################
rule full_analysis_test:
    input:
        "data/test/rep1/operations/modules/test-REP1.full_analysis.finished"
        #"data/test/rep2/operations/modules/test-REP2.full_analysis.finished",
        #"data/test/rep3/operations/modules/test-REP3.full_analysis.finished"

rule full_lncap:
    input:
        "data/pros/lncap/cr01/operations/modules/LNCaP-CR-01-REP1.full_analysis.finished",
        "data/pros/lncap/cr02/operations/modules/LNCaP-CR-02-REP1.full_analysis.finished",
        "data/pros/lncap/cr04/operations/modules/LNCaP-CR-04-REP1.full_analysis.finished",
        "data/pros/lncap/cr05/operations/modules/LNCaP-CR-05-REP1.full_analysis.finished",
        "data/pros/lncap/cr07/operations/modules/LNCaP-CR-07-REP1.full_analysis.finished",
        "data/pros/lncap/cr08/operations/modules/LNCaP-CR-08-REP1.full_analysis.finished",
        "data/pros/lncap/wt01/operations/modules/LNCaP-WT-01-REP1.full_analysis.finished",
        "data/pros/lncap/wt02/operations/modules/LNCaP-WT-02-REP1.full_analysis.finished"
