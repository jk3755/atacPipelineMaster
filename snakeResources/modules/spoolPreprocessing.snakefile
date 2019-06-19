########################################################################################################################################
#### TEST DATA #########################################################################################################################
########################################################################################################################################
rule preprocessing_test:
    input:
        "test/operations/preprocessing/test-REP1.preprocessing.complete",
        "test/operations/preprocessing/test-REP2.preprocessing.complete",
        "test/operations/preprocessing/test-REP3.preprocessing.complete"

########################################################################################################################################
#### LNCaP #############################################################################################################################
########################################################################################################################################
rule preprocessing_lncap_ex01:
    input:
        "lncap/ex01/wt01/operations/preprocessing/LNCaP-WT-01-REP1.preprocessing.complete",
        "lncap/ex01/wt02/operations/preprocessing/LNCaP-WT-02-REP1.preprocessing.complete",
        "lncap/ex01/cr01/operations/preprocessing/LNCaP-CR-01-REP1.preprocessing.complete",
        "lncap/ex01/cr02/operations/preprocessing/LNCaP-CR-02-REP1.preprocessing.complete",
        "lncap/ex01/cr04/operations/preprocessing/LNCaP-CR-04-REP1.preprocessing.complete",
        "lncap/ex01/cr05/operations/preprocessing/LNCaP-CR-05-REP1.preprocessing.complete",
        "lncap/ex01/cr07/operations/preprocessing/LNCaP-CR-07-REP1.preprocessing.complete",
        "lncap/ex01/cr08/operations/preprocessing/LNCaP-CR-08-REP1.preprocessing.complete"

########################################################################################################################################
#### SNU-16 ############################################################################################################################
########################################################################################################################################
rule preprocessing_snu16_wt01:
    input:
        "snu16/wt01/operations/preprocessing/SNU16-WT-01-REP1A.preprocessing.complete",
        "snu16/wt01/operations/preprocessing/SNU16-WT-01-REP2A.preprocessing.complete",
        "snu16/wt01/operations/preprocessing/SNU16-WT-01-REP3A.preprocessing.complete",
        "snu16/wt01/operations/preprocessing/SNU16-WT-01-REP1F.preprocessing.complete",
        "snu16/wt01/operations/preprocessing/SNU16-WT-01-REP2F.preprocessing.complete",
        "snu16/wt01/operations/preprocessing/SNU16-WT-01-REP3F.preprocessing.complete"

########################################################################################################################################
#### MDST8 #############################################################################################################################
########################################################################################################################################
rule preprocessing_mdst8_wt01:
    input:
        "mdst8/wt01/operations/preprocessing/MDST8-WT-01-REP1.preprocessing.complete",
        "mdst8/wt01/operations/preprocessing/MDST8-WT-01-REP2.preprocessing.complete"

########################################################################################################################################
#### Hs675.T ###########################################################################################################################
########################################################################################################################################
rule preprocessing_hs675t_wt01:
    input:
        "hs675t/wt01/operations/preprocessing/Hs675T-WT-01-REP1.preprocessing.complete",
        "hs675t/wt01/operations/preprocessing/Hs675T-WT-01-REP2.preprocessing.complete",
        "hs675t/wt01/operations/preprocessing/Hs675T-WT-01-REP3.preprocessing.complete"


