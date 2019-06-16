########################################################################################################################################
#### SPOOL FOOTPRINTING ################################################################################################################
########################################################################################################################################

rule run_pantf_lncap_group1:
    input:
        expand("lncap/wt02/footprints/operations/groups/LNCaP-WT-02.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr01/footprints/operations/groups/LNCaP-CR-01.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr04/footprints/operations/groups/LNCaP-CR-04.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr07/footprints/operations/groups/LNCaP-CR-07.rawFPanalysis.group{param}.done", param=config["group"])

rule run_pantf_lncap_group2:
    input:
        expand("lncap/wt01/footprints/operations/groups/LNCaP-WT-01.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr02/footprints/operations/groups/LNCaP-CR-02.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr05/footprints/operations/groups/LNCaP-CR-05.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr08/footprints/operations/groups/LNCaP-CR-08.rawFPanalysis.group{param}.done", param=config["group"])

rule run_pantf_lncap_group1:
    input:
        expand("lncap/wt02/footprints/operations/groups/LNCaP-WT-02.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr01/footprints/operations/groups/LNCaP-CR-01.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr04/footprints/operations/groups/LNCaP-CR-04.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr07/footprints/operations/groups/LNCaP-CR-07.rawFPanalysis.group{param}.done", param=config["group"])

rule parse_pantf_lncap_group1:
    input:
        expand("lncap/wt02/footprints/operations/groups/LNCaP-WT-02.parseFP.group{param}.done", param=config["group"]),
        expand("lncap/cr01/footprints/operations/groups/LNCaP-CR-01.parseFP.group{param}.done", param=config["group"]),
        expand("lncap/cr04/footprints/operations/groups/LNCaP-CR-04.parseFP.group{param}.done", param=config["group"]),
        expand("lncap/cr07/footprints/operations/groups/LNCaP-CR-07.parseFP.group{param}.done", param=config["group"])

rule run_pantf_lncap_group2:
    input:
        expand("lncap/wt01/footprints/operations/groups/LNCaP-WT-01.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr02/footprints/operations/groups/LNCaP-CR-02.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr05/footprints/operations/groups/LNCaP-CR-05.rawFPanalysis.group{param}.done", param=config["group"]),
        expand("lncap/cr08/footprints/operations/groups/LNCaP-CR-08.rawFPanalysis.group{param}.done", param=config["group"])

rule parse_pantf_lncap_group2:
    input:
        expand("lncap/wt01/footprints/operations/groups/LNCaP-WT-01.parseFP.group{param}.done", param=config["group"]),
        expand("lncap/cr02/footprints/operations/groups/LNCaP-CR-02.parseFP.group{param}.done", param=config["group"]),
        expand("lncap/cr05/footprints/operations/groups/LNCaP-CR-05.parseFP.group{param}.done", param=config["group"]),
        expand("lncap/cr08/footprints/operations/groups/LNCaP-CR-08.parseFP.group{param}.done", param=config["group"])