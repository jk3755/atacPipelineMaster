########################################################################################################################################
#### TEST DATA #########################################################################################################################
########################################################################################################################################
rule preprocessing_test:
    input:
        "test/operations/test-REP1.preprocessing.complete.txt",
        "test/operations/test-REP2.preprocessing.complete.txt",
        "test/operations/test-REP3.preprocessing.complete.txt"

########################################################################################################################################
#### H508 ##############################################################################################################################
########################################################################################################################################
rule preprocessing_h508_wt01:
    input:
        "h508/wt01/operations/H508-WT-01-REP1.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP2.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP3.preprocessing.complete.txt"

########################################################################################################################################
#### LNCaP #############################################################################################################################
########################################################################################################################################
rule preprocessing_lncap_all:
    input:
        "lncap/wt01/operations/LNCaP-WT-01-REP1.preprocessing.complete.txt",
        "lncap/wt02/operations/LNCaP-WT-02-REP1.preprocessing.complete.txt",
        "lncap/cr01/operations/LNCaP-CR-01-REP1.preprocessing.complete.txt",
        "lncap/cr02/operations/LNCaP-CR-02-REP1.preprocessing.complete.txt",
        "lncap/cr04/operations/LNCaP-CR-04-REP1.preprocessing.complete.txt",
        "lncap/cr05/operations/LNCaP-CR-05-REP1.preprocessing.complete.txt",
        "lncap/cr07/operations/LNCaP-CR-07-REP1.preprocessing.complete.txt",
        "lncap/cr08/operations/LNCaP-CR-08-REP1.preprocessing.complete.txt"