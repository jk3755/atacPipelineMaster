########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################
# Spool the footprinting as below. Preprocessing will be completed first, if not already done
# for i in {1..55}; do snakemake -j 20 run_raw_footprint_test --resources hg38align=1 --config group=$i; done

########################################################################################################################################
#### TEST ##############################################################################################################################
########################################################################################################################################
rule run_raw_footprint_test:
    input:
        "test/operations/preprocessing/test-REP1.preprocessing.complete",
        expand("test/operations/footprints/groups/raw/test-REP1.rawFPanalysis.group{param}.done", param=config["group"])

########################################################################################################################################
#### LNCaP #############################################################################################################################
########################################################################################################################################
rule run_raw_footprint_lncap_ex01:
    input:
        "lncap/ex01/wt01/operations/preprocessing/LNCaP-WT-01-REP1.preprocessing.complete",
        "lncap/ex01/wt02/operations/preprocessing/LNCaP-WT-02-REP1.preprocessing.complete",
        "lncap/ex01/cr01/operations/preprocessing/LNCaP-CR-01-REP1.preprocessing.complete",
        "lncap/ex01/cr02/operations/preprocessing/LNCaP-CR-02-REP1.preprocessing.complete",
        "lncap/ex01/cr04/operations/preprocessing/LNCaP-CR-04-REP1.preprocessing.complete",
        "lncap/ex01/cr05/operations/preprocessing/LNCaP-CR-05-REP1.preprocessing.complete",
        "lncap/ex01/cr07/operations/preprocessing/LNCaP-CR-07-REP1.preprocessing.complete",
        "lncap/ex01/cr08/operations/preprocessing/LNCaP-CR-08-REP1.preprocessing.complete",
        ##
        expand("lncap/ex01/wt01/operations/footprints/groups/raw/LNCaP-WT-01-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/ex01/wt02/operations/footprints/groups/raw/LNCaP-WT-02-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/ex01/cr01/operations/footprints/groups/raw/LNCaP-CR-01-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/ex01/cr02/operations/footprints/groups/raw/LNCaP-CR-02-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/ex01/cr04/operations/footprints/groups/raw/LNCaP-CR-04-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/ex01/cr05/operations/footprints/groups/raw/LNCaP-CR-05-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/ex01/cr07/operations/footprints/groups/raw/LNCaP-CR-07-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/ex01/cr08/operations/footprints/groups/raw/LNCaP-CR-08-REP1.rawFPanalysis.group{param}.done", param=config["group"])

########################################################################################################################################
#### MDST8 #############################################################################################################################
########################################################################################################################################
rule run_raw_footprint_mdst8_wt01:
    input:
        "mdst8/wt01/operations/preprocessing/MDST8-WT-01-REP1.preprocessing.complete",
        "mdst8/wt01/operations/preprocessing/MDST8-WT-01-REP2.preprocessing.complete",
        ##
        expand("mdst8/wt01/operations/footprints/groups/raw/MDST8-WT-01-REP1.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("mdst8/wt01/operations/footprints/groups/raw/MDST8-WT-01-REP2.rawFPanalysis.group{param}.done", param=config["group"])
