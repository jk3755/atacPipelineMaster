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

########################################################################################################################################
#### H508 ##############################################################################################################################
########################################################################################################################################

rule preprocessing_h508_wt01:
    input:
        "h508/wt01/operations/H508-WT-01-REP1.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP2.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP3.preprocessing.complete.txt"

########################################################################################################################################
#### Hs675.T ###########################################################################################################################
########################################################################################################################################

rule preprocessing_hs675t:
    input:
        "hs675t/r1/operations/Hs675T-WT-01-REP1.preprocessing.complete.txt",
        "hs675t/r2/operations/Hs675T-WT-01-REP2.preprocessing.complete.txt",
        "hs675t/r3/operations/Hs675T-WT-01-REP3.preprocessing.complete.txt"

########################################################################################################################################
#### LS1034 ############################################################################################################################
########################################################################################################################################

rule preprocessing_h508_wt01:
    input:
        "h508/wt01/operations/H508-WT-01-REP1.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP2.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP3.preprocessing.complete.txt"

########################################################################################################################################
#### MDST8 #############################################################################################################################
########################################################################################################################################

rule preprocessing_h508_wt01:
    input:
        "h508/wt01/operations/H508-WT-01-REP1.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP2.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP3.preprocessing.complete.txt"

########################################################################################################################################
#### SNU-61 ############################################################################################################################
########################################################################################################################################

rule preprocessing_h508_wt01:
    input:
        "h508/wt01/operations/H508-WT-01-REP1.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP2.preprocessing.complete.txt",
        "h508/wt01/operations/H508-WT-01-REP3.preprocessing.complete.txt"

########################################################################################################################################
#### SNU-16 ############################################################################################################################
########################################################################################################################################

rule preprocessing_snu16:
    input:
        "snu16/ad01/operations/SNU16-WT-01-AD-REP1.preprocessing.complete.txt",
        "snu16/ad02/operations/SNU16-WT-01-AD-REP2.preprocessing.complete.txt",
        "snu16/ad03/operations/SNU16-WT-01-AD-REP3.preprocessing.complete.txt",
        "snu16/fl01/operations/SNU16-WT-01-FL-REP1.preprocessing.complete.txt",
        "snu16/fl02/operations/SNU16-WT-01-FL-REP2.preprocessing.complete.txt",
        "snu16/fl03/operations/SNU16-WT-01-FL-REP3.preprocessing.complete.txt"
