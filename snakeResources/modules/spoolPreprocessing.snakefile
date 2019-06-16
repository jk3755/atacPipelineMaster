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
        "lncap/ex01/wt01/operations/LNCaP-WT-01-REP1.preprocessing.complete.txt",
        "lncap/ex01/wt02/operations/LNCaP-WT-02-REP1.preprocessing.complete.txt",
        "lncap/ex01/cr01/operations/LNCaP-CR-01-REP1.preprocessing.complete.txt",
        "lncap/ex01/cr02/operations/LNCaP-CR-02-REP1.preprocessing.complete.txt",
        "lncap/ex01/cr04/operations/LNCaP-CR-04-REP1.preprocessing.complete.txt",
        "lncap/ex01/cr05/operations/LNCaP-CR-05-REP1.preprocessing.complete.txt",
        "lncap/ex01/cr07/operations/LNCaP-CR-07-REP1.preprocessing.complete.txt",
        "lncap/ex01/cr08/operations/LNCaP-CR-08-REP1.preprocessing.complete.txt"