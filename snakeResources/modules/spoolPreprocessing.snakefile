########################################################################################################################################
#### SPOOL PREPROCESSING ###############################################################################################################
########################################################################################################################################

rule preprocessing_h508_wt01:
    input:
        "h508/wt01/operations/H508-WT-01-REP1.preprocessing.complete.clean.txt",
        "h508/wt01/operations/H508-WT-01-REP2.preprocessing.complete.clean.txt",
        "h508/wt01/operations/H508-WT-01-REP3.preprocessing.complete.clean.txt"

rule preprocessing_lncap_group1:
    input:
        "lncap/wt02/operations/LNCaP-WT-02-pipeline.complete.txt",
        "lncap/cr01/operations/LNCaP-CR-01-pipeline.complete.txt",
        "lncap/cr04/operations/LNCaP-CR-04-pipeline.complete.txt",
        "lncap/cr07/operations/LNCaP-CR-07-pipeline.complete.txt"

rule preprocessing_lncap_group2:
    input:
        "lncap/wt01/operations/LNCaP-WT-01-pipeline.complete.txt",
        "lncap/cr02/operations/LNCaP-CR-02-pipeline.complete.txt",
        "lncap/cr05/operations/LNCaP-CR-05-pipeline.complete.txt",
         "lncap/cr08/operations/LNCaP-CR-08-pipeline.complete.txt"