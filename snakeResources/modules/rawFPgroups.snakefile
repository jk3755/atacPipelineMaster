rule rawFPanalysis_group1:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MUSC.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SPZ1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ONEC2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.MLXPL.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZSC16.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN350.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.BRAC.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF435.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.HEN1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN639.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF350.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN524.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF282.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX9.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB49.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.DMC1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF713.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF8.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.STAT5B.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB16.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group1.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group2:
    input:
        '{path}footprints/operations/raw/{mergedsample}.TBX19.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PTF1A.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.AP2D.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.MCR.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHE41.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.TF2LX.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.TF7L1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN423.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN232.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN143.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBT7A.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXA7.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN589.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX6.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.Pax6.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXC6.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXC8.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.KCNIP1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.AP2A.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF232.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group2.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group3:
    input:
        '{path}footprints/operations/raw/{mergedsample}.HME1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO5F1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.MGAP.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.P5F1B.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC8.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.CR3L2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.RAD21.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA7.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBT7B.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMC3.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZKSCAN3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF306.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.HSFY1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB7C.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.PURA.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF410.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.RFX7.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.HLX.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC4.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.BC11A.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group3.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group4:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZEP2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR3C2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZSCAN4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBT18.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB4.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.GMEB1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HCFC1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.P53.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF7.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF691.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.OTP.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC5.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF423.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.SIX3.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.Trp73.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.RHOXF2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.BARH2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO3F4.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRD14.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.SIN3A.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group4.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group5:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NKX61.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PSMC2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SIX6.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.COE1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ANDR.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.GZF1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXD1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.EVI1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.SIX5.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB14.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXA9.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCF21.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF148.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.BCL3.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.SP100.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.SIX4.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC9.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO4F3.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN713.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group5.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group6:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZN282.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFCP2L1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO3F3.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO2F3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXB2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.CHD2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.LHX1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.LHX5.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.BBX.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.HDX.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.BMAL1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO2F2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.CREB3L2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.PKNX1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX30.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.BSH.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.ERR1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO3F1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.OSR2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYBA.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group6.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group7:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NDF1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO4F1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO6F1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.GRHL2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RORA_2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATF6A.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.CYCS.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.GTF2F1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN410.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXC10.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.AP2C.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA9.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZEP1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRX4.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.HIVEP1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF281.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU4F1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN384.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHE23.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.OSR1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group7.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group8:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZKSC3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFAP2E.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC6.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRX6.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZKSC1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZSCA4.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMARCC2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.RORG.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXC11.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.PKNX2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF274.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBL1XR1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXJ1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR1H2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.EP300.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.TYY1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF9.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRDM14.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group8.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group9:
    input:
        '{path}footprints/operations/raw/{mergedsample}.FOXN1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARI5B.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXC13.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.SP8.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXB6.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.TLX2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.MECOM.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN148.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFAP2D.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.CR3L1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXB3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETV7.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.KAISO.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBX3.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX17.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.AIRE.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN740.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB7B.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.PIT1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.NGN2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group9.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group10:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GLIS1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.MIZF.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF143.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZIC4.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.PATZ1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXD3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX8.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO4F2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHE22.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRDM4.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBT49.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.CENPB.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.PPARD.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF524.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR1D1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.SIX2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF75A.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.PROX1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.SCRT1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group10.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group11:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MAZ.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.HME2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATF4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOMEZ.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ASCL1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX12.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCF7.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.TAF1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ELF4.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFYC.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.NF2L2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.SRBP2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.TYY2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.E2F8.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZFX.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.BCL6.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.PLAL1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.CTCF.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.GLIS2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.SCRT2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group11.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group12:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZBTB6.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB33.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF15.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.STAT2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TF7L2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.HMGA2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.RREB1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF14.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ERF.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAFA.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.HINFP1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.MTF1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.BARHL1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB7A.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR3C1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.RHXF1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.SP4.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXA10.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARID5A.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX11.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group12.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group13:
    input:
        '{path}footprints/operations/raw/{mergedsample}.TP63.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.Trp53.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ELF3.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF9.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.T.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.SIX1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.THRB.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HES1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.E2F7.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.CREB5.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.PLAG1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF4.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR6A1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF8.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN784.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMRC1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.SP2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATF7.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group13.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group14:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NKX31.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU1F1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.DBX2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ASCL2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.EPAS1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNF1A.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXA11.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR2C1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.COT2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.E2F5.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.CTCFL.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.LHX8.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.REST.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.PKNOX1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYOG.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.TP53.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.ONECUT3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF5.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.HSF4.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF652.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group14.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group15:
    input:
        '{path}footprints/operations/raw/{mergedsample}.SOX7.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.TEAD2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN219.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HES7.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF219.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF16.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HINFP.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.AR.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFCP2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRX5.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.PML.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.STAT4.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.BACH2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF13.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU4F2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFKB2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX18.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.CREM.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEF2B.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARID5B.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group15.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group16:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NFAC1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.RFX2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.EBF3.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETV5.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFE2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.EGR2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.NDF2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.GCM2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.E2F3.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.BATF3.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.CREB3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.VENTX.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETV3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRX2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.NKX32.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETV2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.CDC5L.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFAT5.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.RFX5.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.RFX1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group16.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group17:
    input:
        '{path}footprints/operations/raw/{mergedsample}.BHA15.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.Sox4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR2F6.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN652.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.E4F1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR2C2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETV1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ESR1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF35.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.WT1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.HLTF.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.GCM1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.NKX22.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.EGR3.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCF7L1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETV6.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.GABP1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.DMBX1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.CLOCK.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group17.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group18:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZKSCAN1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX21.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.CUX2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU3F1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.Myf.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNF6.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.YBX1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF238.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.PBX1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNF4G.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAFF.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBX2.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.DMRT3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.MLXIPL.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.CYB5R1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.SPDEF.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.SRF.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYBL1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.ELK3.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.TLX1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group18.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group19:
    input:
        '{path}footprints/operations/raw/{mergedsample}.DDIT3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.BARH1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.NRF1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HIF1A.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF784.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.GLI2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.SRBP1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXA2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.NME1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHE40.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.CHD1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.EOMES.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZSCAN31.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHLHB2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYBB.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXB8.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.SPAG7.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.YY2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB8.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.P73.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group19.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group20:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RORB.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.CEBPE.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATF5.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXA13.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.YBOX1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.RXRG.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNF1B.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.THB.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBX20.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAGOH.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.ONEC3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXD8.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.XBP1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF385A.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.GABPA.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR0B1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXC12.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ESRRB.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO6F2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.HMBOX1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group20.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group21:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ELF5.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX8.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMAD1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PBX2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFYB.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.RORA_1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNF4A.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.RORC.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.NKX21.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHLHA15.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.AP2B.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.HAND1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.ELF1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEF2C.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.HLF.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.PBX3.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.NKX23.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.HSFY2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.E2F4.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group21.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group22:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GLI1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.HTF4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.TGIF2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFEC.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBP.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.ONECUT2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.GSC2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HIST1H2BN.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA10.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXA2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.THAP1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU2F1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXA1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.MBD2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.BCL6B.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.E2F6.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.CART1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.TF65.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.LHX3.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.GLIS3.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group22.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group23:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NFE2L2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.EGR4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXB1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ONECUT1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.RELA.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.GBX1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.NHLH1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR1H4.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.UBE2V1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.E2F2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXB7.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXA5.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB7.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.DUX4.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBR1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.HIC1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.TRIM69.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.SP3.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.DUXA.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group23.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group24:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MSC.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.GLI3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.RORA.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXD12.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SPIC.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.PHX2B.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.RELB.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.PHX2A.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.UNC4.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATF1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.ELK4.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARX.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.EGR1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX7.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.STAU2.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.HESX1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.STAT6.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ESRRG.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.STAT1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.EMX2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group24.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group25:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NEUROD1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZHX3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.RX.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCF7L2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.STAT5A.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.P63.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXD13.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HMX2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.DRGX.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.HIST2H2BE.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRDM1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZIC1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.NKX62.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.DPRX.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.CEBPZ.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.WRNIP1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.LHX6.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF24.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOSL2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXD10.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group25.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group26:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NEUROD2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXD10.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOSB.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXD11.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.REL.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.PROP1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.EVX2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ELK1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ERG.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZIC3.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.THA.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXP4.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXA3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXK2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.LRRFIP1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZDHHC15.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.CUX1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFYA.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARNTL.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.ALX1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group26.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group27:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZNF740.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR4A3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFE2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.TWST1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.STAT3.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.SUH.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBX21.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO3F2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.HSF2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ALX4.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR2F2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATF6.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZFHX3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU3F3.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.HMBX1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.GBX2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.CEBPD.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATF2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.NKX28.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.GFI1B.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group27.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group28:
    input:
        '{path}footprints/operations/raw/{mergedsample}.DLX1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFKB1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HMX3.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC13.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHLHE40.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXQ1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.E2F1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATOH1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXH1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETS2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEOX1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.GSC.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ASAP3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU6F2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.MLX.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.SREBF2.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.NOTO.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB3.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.TEAD4.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFAC4.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group28.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group29:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GMEB2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.BARX2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.RFX3.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU2F2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RXRB.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.EZH2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.LBX2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.MNT.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFAP2B.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXF2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.YWHAE.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHLHB3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF263.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.STA5B.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF333.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZIC2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.RUVBL1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.HMX1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.GSX2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group29.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group30:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RUNX1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.DNMT3A.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA11.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.BPTF.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.COT1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.PHOX2B.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.TAGLN2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHLHE41.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ID4.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFATC1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXB13.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRX3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXA1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXC2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBM22.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF5.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXD8.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.CEBPA.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXD13.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXF1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group30.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group31:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ERR3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.STF1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXP1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ERR2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.EMX1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.JUND.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX6.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.EN2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.STA5A.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.HBP1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.MRPS25.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.PEBB.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHLHE22.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.CERS4.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXD11.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.XG.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.H2AFZ.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR5A2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.HEY1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.VSX2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group31.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group32:
    input:
        '{path}footprints/operations/raw/{mergedsample}.FIGLA.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.NCALD.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SREBF1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.TEF.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU6F1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.GSX1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.MTA3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.EN1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.PITX2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.BARHL2.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEF2D.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYC.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF8.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.THAP5.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.TWIST1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.JUNB.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB13.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.NEUROG2.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAFK.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.BATF.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group32.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group33:
    input:
        '{path}footprints/operations/raw/{mergedsample}.HOXB2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.RCOR1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.OVOL1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.MRPL1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZN333.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.ATF3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.IKZF2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.PCK2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.HSPA5.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.GATA1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA13.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.BHLHE23.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMARCC1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZC3H7A.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.MITF.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEF2A.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.MNX1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.EWSR1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYBL2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group33.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group34:
    input:
        '{path}footprints/operations/raw/{mergedsample}.CEBPB.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOS.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF711.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.H2AFY.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX4.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC12.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.ISL2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.APEX2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.NOBOX.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU5F1B.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.PO2F1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.PITX3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.OTX1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFAC2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNRNPLL.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU2F3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU5F1P1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.SEMA4A.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR2E1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.RUNX2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group34.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group35:
    input:
        '{path}footprints/operations/raw/{mergedsample}.FLI1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMUG1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SRP9.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HES5.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.DUSP26.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAX.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HEY2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF7.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.EIF5A2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU3F4.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.C9orf156.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAFB.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA5.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.OLIG2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFAP4.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.ALX3.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBPL1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMPX.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXC1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.AGGF1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group35.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group36:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GLYCTK.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SND1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PPARG.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC10.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TIMM8A.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.NNT.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.RPS4X.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX9.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNRNPA1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBPJ.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.SNRPB2.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.BACH1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARNT.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.DDX20.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.SP1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB4.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.MESP1.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.GTF2I.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.VAX2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF6.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group36.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group37:
    input:
        '{path}footprints/operations/raw/{mergedsample}.PTCD1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXB1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ODC1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.JUN.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.PURG.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.ELF2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.MSRA.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HIC2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.GATA5.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRKRIR.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.NMI.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.OTX2.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.HCLS1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFE3.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.LARP4.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.FEZF2.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.TEAD3.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.LBX1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.LUZP1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group37.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group38:
    input:
        '{path}footprints/operations/raw/{mergedsample}.HOXD12.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SNAPC5.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.NANOG.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXJ3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TP73.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.XRCC1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEIS3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.UBE2K.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.CNOT6.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.SSBP3.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFAC3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.OLIG3.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.IL24.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.MXI1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.TAL1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMAD4.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXC11.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.GFI1.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.GOT1.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group38.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group39:
    input:
        '{path}footprints/operations/raw/{mergedsample}.DHX36.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.RAN.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.TRIP10.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.LMX1A.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNRNPA0.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXJ2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.RUNX3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.FIP1L1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF695.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.CXXC1.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.VAX1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZFP3.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZEB1.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.FHL2.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.PDX1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.ESRP1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.FEV.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF766.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX2.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.JDP2.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group39.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group40:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MIXL1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAF.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXO1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYCN.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.FUBP1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.ITF2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.UQCRB.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.LARP1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBX15.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBX4.rawFPanalysis.bamcopy10.done', 
        '{path}footprints/operations/raw/{mergedsample}.ANXA1.rawFPanalysis.bamcopy11.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBX1.rawFPanalysis.bamcopy12.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF384.rawFPanalysis.bamcopy13.done', 
        '{path}footprints/operations/raw/{mergedsample}.EDN1.rawFPanalysis.bamcopy14.done', 
        '{path}footprints/operations/raw/{mergedsample}.SNAI1.rawFPanalysis.bamcopy15.done', 
        '{path}footprints/operations/raw/{mergedsample}.CREB1.rawFPanalysis.bamcopy16.done', 
        '{path}footprints/operations/raw/{mergedsample}.BSX.rawFPanalysis.bamcopy17.done', 
        '{path}footprints/operations/raw/{mergedsample}.GPD1.rawFPanalysis.bamcopy18.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRDX5.rawFPanalysis.bamcopy19.done', 
        '{path}footprints/operations/raw/{mergedsample}.RAB14.rawFPanalysis.bamcopy20.done'
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group40.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group41:
    input:
        '{path}footprints/operations/raw/{mergedsample}.TFDP1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX14.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.DLX4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR2E3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.FAM127B.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.LINC00471.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEOX2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFAP2A.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCF12.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB18.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group41.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group42:
    input:
        '{path}footprints/operations/raw/{mergedsample}.PDE6H.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.CBX3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBED1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXD1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRNP.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFEB.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARNT2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.AHR.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.RAB2A.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group42.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group43:
    input:
        '{path}footprints/operations/raw/{mergedsample}.FOXO3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.CRX.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZMYND8.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.SSX2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.KDM5A.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.VSX1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.LUZP2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXD4.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.RFXANK.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.PGAM2.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group43.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group44:
    input:
        '{path}footprints/operations/raw/{mergedsample}.FOSL1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETV4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXD4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.LMX1B.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TSC22D4.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXD9.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HXD9.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR1I3.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.DIABLO.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ESRRA.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group44.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group45:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MSX1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.GRHL1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAGEF1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB5.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF124.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.PDCD11.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.TEAD1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYB.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.CFL2.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.CD59.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group45.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group46:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GATA4.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.AP1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PDS5A.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NELFB.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXN3.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCF3.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.EZR.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF71.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXO6.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXK1.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group46.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group47:
    input:
        '{path}footprints/operations/raw/{mergedsample}.FOXD2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.OLIG1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXI1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NANOS1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ACO1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFAP2C.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.BOLL.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFIA.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.DTL.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.POLI.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group47.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group48:
    input:
        '{path}footprints/operations/raw/{mergedsample}.FOXP2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF12.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.INSM1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFIB.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.GRHPR.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.BRCA1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBMS1.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.HIRIP3.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.UNCX.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.RAXL1.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group48.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group49:
    input:
        '{path}footprints/operations/raw/{mergedsample}.SHOX2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.LHX9.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.RAX2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.DBP.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF250.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.MRPL2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.ADNP.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ADARB1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.TAF4.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.TMSB4XP8.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group49.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group50:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MEIS2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR1I2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZSWIM1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.CREB3L1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.MGA.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.TAF1A.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.USF2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.OTUD4.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.WDR83.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.KDM5D.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group50.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group51:
    input:
        '{path}footprints/operations/raw/{mergedsample}.TULP1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PQBP1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.RPP25.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXP3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.DLX3.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXG1.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.TIMELESS.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCEAL2.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEX3C.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.EEF1D.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group51.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group52:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NUCB1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.TSN.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.MXD4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ISX.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARFGAP1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.TBX5.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.DAB2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.ISL1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.TPI1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB46.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group52.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group53:
    input:
        '{path}footprints/operations/raw/{mergedsample}.PGR.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.GCR.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRGR.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.MSX2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.HTATIP2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.ESR2.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.MBTPS2.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXO4.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZKSCAN8.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNRNPC.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group53.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group54:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NKX25.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.HHEX.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.NUP107.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAFG.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SUCLG1.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOCS4.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXD3.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.SHOX.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.USF1.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF76.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group54.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group55:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZNF326.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SRBD1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.DGCR8.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.SNAPC4.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.MDM2.rawFPanalysis.bamcopy5.done', 
        '{path}footprints/operations/raw/{mergedsample}.TPPP.rawFPanalysis.bamcopy6.done', 
        '{path}footprints/operations/raw/{mergedsample}.NELFA.rawFPanalysis.bamcopy7.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXL1.rawFPanalysis.bamcopy8.done', 
        '{path}footprints/operations/raw/{mergedsample}.GTPBP6.rawFPanalysis.bamcopy9.done', 
        '{path}footprints/operations/raw/{mergedsample}.H1FX.rawFPanalysis.bamcopy10.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group55.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group56:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MAGEA8.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.HP1BP3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU3F2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.EBF1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.GTF3C5.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group56.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group57:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NONO.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.CDX2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYOD1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PPP2R3B.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SNRNP70.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group57.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group58:
    input:
        '{path}footprints/operations/raw/{mergedsample}.TOB2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB12.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF6.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ETFB.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.DLX2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group58.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group59:
    input:
        '{path}footprints/operations/raw/{mergedsample}.SPIB.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF304.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.VAMP3.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF205.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU5F1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group59.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group60:
    input:
        '{path}footprints/operations/raw/{mergedsample}.EHF.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.RAX.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SRY.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.MEIS1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TIA1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group60.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group61:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RXRA.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMAD3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ESRP2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NRL.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBM42.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group61.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group62:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZNF3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.CPEB1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.CDX1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NF2L1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TGIF2LX.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group62.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group63:
    input:
        '{path}footprints/operations/raw/{mergedsample}.PKNOX2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.LEF1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PPARA.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PPP1R10.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SLC18A1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group63.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group64:
    input:
        '{path}footprints/operations/raw/{mergedsample}.DIS3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.MAP4K2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.USP39.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.GTF2B.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZCCHC14.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group64.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group65:
    input:
        '{path}footprints/operations/raw/{mergedsample}.FOXD3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.EXO5.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.MSI2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.TSNAX.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF696.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group65.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group66:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NFATC2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCF4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.TIMM44.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.PIK3C3.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group66.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group67:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GTF3C2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX15.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PKM.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.RIOK2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF720.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group67.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group68:
    input:
        '{path}footprints/operations/raw/{mergedsample}.DLX5.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PIR.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ECSIT.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.IVD.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.NCBP2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group68.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group69:
    input:
        '{path}footprints/operations/raw/{mergedsample}.SF3B1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.KLF4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF160.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF354C.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMCR7L.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group69.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group70:
    input:
        '{path}footprints/operations/raw/{mergedsample}.HCFC2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.TAF9.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HHAT.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.RARB.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RHOXF1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group70.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group71:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ARI3A.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ARID3A.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR4A2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF830.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.TRMT1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group71.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group72:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RFX4.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SRRM3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFIC.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ASPSCR1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ABCF2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group72.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group73:
    input:
        '{path}footprints/operations/raw/{mergedsample}.YEATS4.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF26.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SNAI2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HSF1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA4.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group73.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group74:
    input:
        '{path}footprints/operations/raw/{mergedsample}.AGAP2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB9.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PDLIM5.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ASCC1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SFT2D1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group74.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group75:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RARA.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.BAX.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRRX1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBFOX2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.GIT2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group75.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group76:
    input:
        '{path}footprints/operations/raw/{mergedsample}.CBFA2T2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF706.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.RAB7A.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.RPS10.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF655.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group76.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group77:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NUP133.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.CBX7.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAXIP1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NMRAL1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.MED30.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group77.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group78:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MYEF2.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.FEZ1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.PAX5.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBM8A.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.MSI1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group78.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group79:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ANXA11.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.RARG.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB25.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.THRA.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RNASEH2C.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group79.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group80:
    input:
        '{path}footprints/operations/raw/{mergedsample}.CAT.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.RFC3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.RPS6KA5.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX5.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFE2L1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group80.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group81:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NXPH3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.TCEAL6.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.YY1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PSMA6.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFIL3.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group81.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group82:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MZF1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.VDR.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.LHX4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PITX1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ID2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group82.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group83:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ETS1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.LSM6.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.TGIF1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMAD2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.UTP18.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group83.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group84:
    input:
        '{path}footprints/operations/raw/{mergedsample}.PLAGL1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.POLE3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF503.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PHOX2A.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.CLK1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group84.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group85:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RAB18.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.HNRNPH3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ING3.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ENO1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.U2AF1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group85.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group86:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RPL6.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.IKZF1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.P4HB.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.CCDC25.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.NOC2L.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group86.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group87:
    input:
        '{path}footprints/operations/raw/{mergedsample}.HDAC8.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.DDX43.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SCAND2P.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.CSNK2B.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.FOXM1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group87.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group88:
    input:
        '{path}footprints/operations/raw/{mergedsample}.SCMH1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.POU4F3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.UGP2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.GATA6.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.C19orf40.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group88.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group89:
    input:
        '{path}footprints/operations/raw/{mergedsample}.HIST2H2AB.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SMAP2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.TRIM21.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.BCL11A.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.CSTF2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group89.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group90:
    input:
        '{path}footprints/operations/raw/{mergedsample}.NR2F1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.CELF5.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.NAP1L1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.UBIP1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.UBP1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group90.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group91:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MYF6.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.DLX6.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SPI1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.MTHFD1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBBP5.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group91.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group92:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MAGED4B.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PHLDA2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.TFAM.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.LAS1L.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.MYLK.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group92.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group93:
    input:
        '{path}footprints/operations/raw/{mergedsample}.AVEN.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZMAT4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.EVX1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFATC4.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.PICK1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group93.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group94:
    input:
        '{path}footprints/operations/raw/{mergedsample}.CPSF4.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.VPS4B.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.R3HDM2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.FGF19.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.AFF4.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group94.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group95:
    input:
        '{path}footprints/operations/raw/{mergedsample}.DUSP22.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFATC3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.BARX1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HSPA1L.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.MSRB3.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group95.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group96:
    input:
        '{path}footprints/operations/raw/{mergedsample}.CEBPG.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZCCHC17.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.RNF114.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBM3.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.GATA2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group96.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group97:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GTPBP1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PLG.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.MECP2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.CELF6.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.GAR1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group97.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group98:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GTF2H3.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.MORN1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.NFIX.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.YWHAZ.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBM17.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group98.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group99:
    input:
        '{path}footprints/operations/raw/{mergedsample}.CKMT1B.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.GATA3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX10.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PPP5C.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOD1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group99.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group100:
    input:
        '{path}footprints/operations/raw/{mergedsample}.PTPMT1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.HMGA1.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.CANX.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.SOX13.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.NR4A1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group100.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group101:
    input:
        '{path}footprints/operations/raw/{mergedsample}.ZDHHC5.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZMAT2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.WISP2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.SF1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF510.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group101.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group102:
    input:
        '{path}footprints/operations/raw/{mergedsample}.GPAM.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF131.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HMG20A.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF671.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RBBP9.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group102.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group103:
    input:
        '{path}footprints/operations/raw/{mergedsample}.UBXN1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.C19orf25.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.STUB1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.PHTF1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.RUFY3.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group103.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group104:
    input:
        '{path}footprints/operations/raw/{mergedsample}.METTL21B.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.CELF4.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.A1CF.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.GPANK1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.KIF22.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group104.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group105:
    input:
        '{path}footprints/operations/raw/{mergedsample}.HLCS.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.RFC2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.IRF1.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.CDK2AP1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.GADD45A.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group105.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group106:
    input:
        '{path}footprints/operations/raw/{mergedsample}.DUS3L.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.LHX2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZNF207.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.KIAA0907.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SPATS2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group106.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group107:
    input:
        '{path}footprints/operations/raw/{mergedsample}.MAPK1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.EXOSC3.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.DDX53.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZBTB43.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.GLTPD1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group107.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group108:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RNF138.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.UBB.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.BAD.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.DAZAP1.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.SSX3.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group108.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group109:
    input:
        '{path}footprints/operations/raw/{mergedsample}.AKR1A1.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.SPR.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.DDX4.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.RPL35.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.ESX1.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group109.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group110:
    input:
        '{path}footprints/operations/raw/{mergedsample}.CBFB.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.PRRX2.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.TROVE2.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZRSR2.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.MCTP2.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group110.done'
    shell:
        'touch {output}'
rule rawFPanalysis_group111:
    input:
        '{path}footprints/operations/raw/{mergedsample}.RBM7.rawFPanalysis.bamcopy1.done', 
        '{path}footprints/operations/raw/{mergedsample}.ZSCAN9.rawFPanalysis.bamcopy2.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXA6.rawFPanalysis.bamcopy3.done', 
        '{path}footprints/operations/raw/{mergedsample}.HOXB6.rawFPanalysis.bamcopy4.done', 
        '{path}footprints/operations/raw/{mergedsample}.NA.rawFPanalysis.bamcopy5.done', 
    output:
        '{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group111.done'
    shell:
        'touch {output}'
