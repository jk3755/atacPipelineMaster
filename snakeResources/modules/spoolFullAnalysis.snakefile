## run preprocessing and footprinting analysis for all samples to date
rule all_lncap:
    input:
        "data/pros/lncap/wt01/operations/modules/LNCaP-WT-01.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/wt02/operations/modules/LNCaP-WT-02.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr01/operations/modules/LNCaP-CR-01.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr02/operations/modules/LNCaP-CR-02.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr04/operations/modules/LNCaP-CR-04.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr05/operations/modules/LNCaP-CR-05.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr07/operations/modules/LNCaP-CR-07.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr08/operations/modules/LNCaP-CR-08.rep1.refhg38.full_analysis.finished"

rule full_ls1034:
    input:
        "data/coad/ls1034/wt01/operations/modules/LS1034-WT-01.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt02/operations/modules/LS1034-WT-02.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt03/operations/modules/LS1034-WT-03.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt04/operations/modules/LS1034-WT-04.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt05/operations/modules/LS1034-WT-05.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt06/operations/modules/LS1034-WT-06.rep1.refhg38.full_analysis.finished"
