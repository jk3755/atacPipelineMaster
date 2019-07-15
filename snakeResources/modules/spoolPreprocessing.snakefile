########################################################################################################################################
#### TEST DATA #########################################################################################################################
########################################################################################################################################
rule preprocessing_test:
    input:
        "test/operations/preprocessing/test-REP1.preprocessing.complete",
        "test/operations/preprocessing/test-REP2.preprocessing.complete",
        "test/operations/preprocessing/test-REP3.preprocessing.complete"

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

########################################################################################################################################
#### H508 ##############################################################################################################################
########################################################################################################################################
rule preprocessing_h508_wt01:
    input:
        "h508/wt01/operations/preprocessing/H508-WT-01-REP1.preprocessing.complete",
        "h508/wt01/operations/preprocessing/H508-WT-01-REP2.preprocessing.complete",
        "h508/wt01/operations/preprocessing/H508-WT-01-REP3.preprocessing.complete"

rule preprocessing_h508_wt02:
    input:
        "h508/wt02/operations/preprocessing/H508-WT-02-REP1A.preprocessing.complete",
        "h508/wt02/operations/preprocessing/H508-WT-02-REP2A.preprocessing.complete",
        "h508/wt02/operations/preprocessing/H508-WT-02-REP1F.preprocessing.complete"

########################################################################################################################################
#### LS1034 ############################################################################################################################
########################################################################################################################################
rule preprocessing_ls1034_wt01:
    input:
        "ls1034/wt01/operations/preprocessing/LS1034-WT-01-REP1.preprocessing.complete",
        "ls1034/wt01/operations/preprocessing/LS1034-WT-01-REP2.preprocessing.complete",
        "ls1034/wt01/operations/preprocessing/LS1034-WT-01-REP3.preprocessing.complete"

rule preprocessing_ls1034_wt02:
    input:
        "ls1034/wt02/operations/preprocessing/LS1034-WT-02-REP1.preprocessing.complete",
        "ls1034/wt02/operations/preprocessing/LS1034-WT-02-REP2.preprocessing.complete",
        "ls1034/wt02/operations/preprocessing/LS1034-WT-02-REP3.preprocessing.complete"

########################################################################################################################################
#### SNU61 #############################################################################################################################
########################################################################################################################################
rule preprocessing_snu61_wt01:
    input:
        "snu61/wt01/operations/preprocessing/SNU61-WT-01-REP1.preprocessing.complete",
        "snu61/wt01/operations/preprocessing/SNU61-WT-01-REP2.preprocessing.complete",
        "snu61/wt01/operations/preprocessing/SNU61-WT-01-REP3.preprocessing.complete"

########################################################################################################################################
#### LNCaP #############################################################################################################################
########################################################################################################################################
rule preprocessing_lncap_ex01:
    input:
        "lncap/ex01/operations/preprocessing/LNCaP-WT-01-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-WT-02-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-01-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-02-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-04-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-05-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-07-REP1.preprocessing.complete",
        "lncap/ex01/operations/preprocessing/LNCaP-CR-08-REP1.preprocessing.complete"

########################################################################################################################################
#### PANC1 #############################################################################################################################
########################################################################################################################################
rule preprocessing_capani_wt01:
    input:
        "panc/capani/operations/preprocessing/CAPANI-WT-01-RUN1-REP1.preprocessing.complete",
        "panc/capani/operations/preprocessing/CAPANI-WT-01-RUN1-REP2.preprocessing.complete",
        "panc/capani/operations/preprocessing/CAPANI-WT-01-RUN1-REP3.preprocessing.complete",
        "panc/capani/operations/preprocessing/CAPANI-WT-01-RUN2-REP1.preprocessing.complete",
        "panc/capani/operations/preprocessing/CAPANI-WT-01-RUN2-REP2.preprocessing.complete",
        "panc/capani/operations/preprocessing/CAPANI-WT-01-RUN2-REP3.preprocessing.complete"