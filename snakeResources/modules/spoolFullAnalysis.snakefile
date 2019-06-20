########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################

########################################################################################################################################
#### TEST DATA #########################################################################################################################
########################################################################################################################################
rule full_analysis_test:
    input:
        "test/operations/modules/test-REP1.full_analysis.finished",
        "test/operations/modules/test-REP2.full_analysis.finished",
        "test/operations/modules/test-REP3.full_analysis.finished"

########################################################################################################################################
#### H508 ##############################################################################################################################
########################################################################################################################################
rule full_analysis_h508_wt01:
    input:
        "h508/wt01/operations/modules/H508-WT-01-REP1.full_analysis.finished",
        "h508/wt01/operations/modules/H508-WT-01-REP2.full_analysis.finished",
        "h508/wt01/operations/modules/H508-WT-01-REP3.full_analysis.finished"

rule full_analysis_h508_wt02:
    input:
        "h508/wt02/operations/modules/H508-WT-02-REP1A.full_analysis.finished",
        "h508/wt02/operations/modules/H508-WT-02-REP2A.full_analysis.finished",
        "h508/wt02/operations/modules/H508-WT-02-REP1F.full_analysis.finished"

########################################################################################################################################
#### SNU-16 ############################################################################################################################
########################################################################################################################################
rule full_analysis_snu16_wt01:
    input:
        "snu16/wt01/operations/modules/SNU16-WT-01-REP1A.full_analysis.finished",
        "snu16/wt01/operations/modules/SNU16-WT-01-REP2A.full_analysis.finished",
        "snu16/wt01/operations/modules/SNU16-WT-01-REP3A.full_analysis.finished",
        "snu16/wt01/operations/modules/SNU16-WT-01-REP1F.full_analysis.finished",
        "snu16/wt01/operations/modules/SNU16-WT-01-REP2F.full_analysis.finished",
        "snu16/wt01/operations/modules/SNU16-WT-01-REP3F.full_analysis.finished"

########################################################################################################################################
#### MDST8 #############################################################################################################################
########################################################################################################################################
rule full_analysis_mdst8_wt01:
    input:
        "mdst8/wt01/operations/modules/MDST8-WT-01-REP1.full_analysis.finished",
        "mdst8/wt01/operations/modules/MDST8-WT-01-REP2.full_analysis.finished"

########################################################################################################################################
#### Hs675.T ###########################################################################################################################
########################################################################################################################################
rule full_analysis_hs675t_wt01:
    input:
        "hs675t/wt01/operations/modules/Hs675T-WT-01-REP1.full_analysis.finished",
        "hs675t/wt01/operations/modules/Hs675T-WT-01-REP2.full_analysis.finished",
        "hs675t/wt01/operations/modules/Hs675T-WT-01-REP3.full_analysis.finished"

########################################################################################################################################
#### LS1034 ############################################################################################################################
########################################################################################################################################
rule full_analysis_ls1034_wt01:
    input:
        "ls1034/wt01/operations/modules/LS1034-WT-01-REP1.full_analysis.finished",
        "ls1034/wt01/operations/modules/LS1034-WT-01-REP2.full_analysis.finished",
        "ls1034/wt01/operations/modules/LS1034-WT-01-REP3.full_analysis.finished"

rule full_analysis_ls1034_wt02:
    input:
        "ls1034/wt02/operations/modules/LS1034-WT-02-REP1.full_analysis.finished",
        "ls1034/wt02/operations/modules/LS1034-WT-02-REP2.full_analysis.finished",
        "ls1034/wt02/operations/modules/LS1034-WT-02-REP3.full_analysis.finished"

########################################################################################################################################
#### SNU-61 #############################################################################################################################
########################################################################################################################################
rule full_analysis_snu61_wt01:
    input:
        "snu61/wt01/operations/modules/SNU61-WT-01-REP1.full_analysis.finished",
        "snu61/wt01/operations/modules/SNU61-WT-01-REP2.full_analysis.finished",
        "snu61/wt01/operations/modules/SNU61-WT-01-REP3.full_analysis.finished"

########################################################################################################################################
#### LNCaP #############################################################################################################################
########################################################################################################################################
rule full_analysis_lncap_ex01:
    input:
        "lncap/ex01/operations/modules/LNCaP-WT-01-REP1.full_analysis.finished",
        "lncap/ex01/operations/modules/LNCaP-WT-02-REP1.full_analysis.finished",
        "lncap/ex01/operations/modules/LNCaP-CR-01-REP1.full_analysis.finished",
        "lncap/ex01/operations/modules/LNCaP-CR-02-REP1.full_analysis.finished",
        "lncap/ex01/operations/modules/LNCaP-CR-04-REP1.full_analysis.finished",
        "lncap/ex01/operations/modules/LNCaP-CR-05-REP1.full_analysis.finished",
        "lncap/ex01/operations/modules/LNCaP-CR-07-REP1.full_analysis.finished",
        "lncap/ex01/operations/modules/LNCaP-CR-08-REP1.full_analysis.finished"
