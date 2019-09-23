## run preprocessing and footprinting analysis for all samples to date
rule full_lncap_refhg38:
    input:
        "data/pros/lncap/wt01/operations/LNCaP-WT-01.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/wt02/operations/LNCaP-WT-02.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr01/operations/LNCaP-CR-01.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr02/operations/LNCaP-CR-02.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr04/operations/LNCaP-CR-04.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr05/operations/LNCaP-CR-05.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr07/operations/LNCaP-CR-07.rep1.refhg38.full_analysis.finished",
        "data/pros/lncap/cr08/operations/LNCaP-CR-08.rep1.refhg38.full_analysis.finished"

rule full_ls1034:
    input:
        "data/coad/ls1034/wt01/operations/LS1034-WT-01.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt02/operations/LS1034-WT-02.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt03/operations/LS1034-WT-03.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt04/operations/LS1034-WT-04.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt05/operations/LS1034-WT-05.rep1.refhg38.full_analysis.finished",
        "data/coad/ls1034/wt06/operations/LS1034-WT-06.rep1.refhg38.full_analysis.finished"

rule test:
    input:
        "data/coad/ls1034/wt01/operations/LS1034-WT-01.rep1.refhg38.full_analysis.finished"
