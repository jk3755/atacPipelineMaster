rule parseFP_group1:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MUSC.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SPZ1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ONEC2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.MLXPL.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZSC16.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN350.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.BRAC.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF435.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.HEN1.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN639.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF350.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN524.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF282.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX9.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB49.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.DMC1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF713.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF8.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.STAT5B.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB16.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group1.done'
	shell:
		'touch {output}'
rule parseFP_group2:
	input:
		'{path}footprints/operations/parse/{mergedsample}.TBX19.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PTF1A.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.AP2D.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.MCR.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHE41.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.TF2LX.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.TF7L1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN423.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN232.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN143.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBT7A.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXA7.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN589.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX6.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.Pax6.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXC6.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXC8.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.KCNIP1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.AP2A.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF232.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group2.done'
	shell:
		'touch {output}'
rule parseFP_group3:
	input:
		'{path}footprints/operations/parse/{mergedsample}.HME1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO5F1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.MGAP.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.P5F1B.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC8.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.CR3L2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.RAD21.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA7.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBT7B.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMC3.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZKSCAN3.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF306.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.HSFY1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB7C.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.PURA.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF410.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.RFX7.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.HLX.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC4.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.BC11A.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group3.done'
	shell:
		'touch {output}'
rule parseFP_group4:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZEP2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR3C2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZSCAN4.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBT18.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB4.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.GMEB1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HCFC1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.P53.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF7.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF691.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.OTP.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC5.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF423.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.SIX3.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.Trp73.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.RHOXF2.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.BARH2.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO3F4.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRD14.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.SIN3A.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group4.done'
	shell:
		'touch {output}'
rule parseFP_group5:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NKX61.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PSMC2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SIX6.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.COE1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ANDR.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.GZF1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXD1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.EVI1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.SIX5.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB14.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXA9.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB3.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCF21.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF148.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.BCL3.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.SP100.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.SIX4.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC9.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO4F3.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN713.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group5.done'
	shell:
		'touch {output}'
rule parseFP_group6:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZN282.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFCP2L1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO3F3.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO2F3.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXB2.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.CHD2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.LHX1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.LHX5.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.BBX.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.HDX.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.BMAL1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO2F2.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.CREB3L2.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.PKNX1.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX30.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.BSH.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.ERR1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO3F1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.OSR2.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYBA.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group6.done'
	shell:
		'touch {output}'
rule parseFP_group7:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NDF1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO4F1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO6F1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.GRHL2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RORA_2.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATF6A.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.CYCS.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.GTF2F1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN410.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXC10.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.AP2C.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA9.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZEP1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRX4.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.HIVEP1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF281.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU4F1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN384.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHE23.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.OSR1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group7.done'
	shell:
		'touch {output}'
rule parseFP_group8:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZKSC3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFAP2E.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC6.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRX6.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZKSC1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZSCA4.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMARCC2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.RORG.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXC11.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.PKNX2.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF274.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBL1XR1.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXJ1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR1H2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.EP300.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA3.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.TYY1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF9.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRDM14.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group8.done'
	shell:
		'touch {output}'
rule parseFP_group9:
	input:
		'{path}footprints/operations/parse/{mergedsample}.FOXN1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARI5B.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXC13.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.SP8.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXB6.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.TLX2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.MECOM.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN148.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFAP2D.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.CR3L1.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXB3.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETV7.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.KAISO.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBX3.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX17.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.AIRE.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN740.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB7B.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.PIT1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.NGN2.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group9.done'
	shell:
		'touch {output}'
rule parseFP_group10:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GLIS1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.MIZF.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF143.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF3.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZIC4.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.PATZ1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXD3.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX8.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO4F2.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHE22.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRDM4.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBT49.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.CENPB.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.PPARD.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF524.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR1D1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.SIX2.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF75A.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.PROX1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.SCRT1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group10.done'
	shell:
		'touch {output}'
rule parseFP_group11:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MAZ.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.HME2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATF4.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOMEZ.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ASCL1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX12.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCF7.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.TAF1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ELF4.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFYC.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.NF2L2.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.SRBP2.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.TYY2.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.E2F8.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZFX.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.BCL6.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.PLAL1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.CTCF.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.GLIS2.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.SCRT2.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group11.done'
	shell:
		'touch {output}'
rule parseFP_group12:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZBTB6.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB33.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF15.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.STAT2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TF7L2.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.HMGA2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.RREB1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF14.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ERF.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAFA.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.HINFP1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.MTF1.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.BARHL1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB7A.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR3C1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.RHXF1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.SP4.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXA10.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARID5A.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX11.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group12.done'
	shell:
		'touch {output}'
rule parseFP_group13:
	input:
		'{path}footprints/operations/parse/{mergedsample}.TP63.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.Trp53.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ELF3.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF9.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.T.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.SIX1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.THRB.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HES1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.E2F7.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.CREB5.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.PLAG1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF4.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR6A1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX1.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF8.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF2.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN784.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMRC1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.SP2.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATF7.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group13.done'
	shell:
		'touch {output}'
rule parseFP_group14:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NKX31.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU1F1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.DBX2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ASCL2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.EPAS1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNF1A.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXA11.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR2C1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.COT2.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.E2F5.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.CTCFL.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.LHX8.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.REST.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.PKNOX1.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYOG.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.TP53.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.ONECUT3.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF5.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.HSF4.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF652.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group14.done'
	shell:
		'touch {output}'
rule parseFP_group15:
	input:
		'{path}footprints/operations/parse/{mergedsample}.SOX7.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.TEAD2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN219.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HES7.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF219.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF16.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HINFP.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.AR.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFCP2.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRX5.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.PML.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.STAT4.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.BACH2.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF13.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU4F2.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFKB2.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX18.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.CREM.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEF2B.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARID5B.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group15.done'
	shell:
		'touch {output}'
rule parseFP_group16:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NFAC1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.RFX2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.EBF3.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETV5.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFE2.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.EGR2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.NDF2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.GCM2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.E2F3.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.BATF3.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.CREB3.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.VENTX.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETV3.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRX2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.NKX32.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETV2.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.CDC5L.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFAT5.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.RFX5.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.RFX1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group16.done'
	shell:
		'touch {output}'
rule parseFP_group17:
	input:
		'{path}footprints/operations/parse/{mergedsample}.BHA15.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX4.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.Sox4.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR2F6.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN652.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.E4F1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR2C2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETV1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ESR1.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF35.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.WT1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.HLTF.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.GCM1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.NKX22.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.EGR3.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCF7L1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETV6.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.GABP1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.DMBX1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.CLOCK.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group17.done'
	shell:
		'touch {output}'
rule parseFP_group18:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZKSCAN1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX21.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.CUX2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU3F1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.Myf.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNF6.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.YBX1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF238.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.PBX1.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNF4G.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAFF.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBX2.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.DMRT3.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.MLXIPL.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.CYB5R1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.SPDEF.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.SRF.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYBL1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.ELK3.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.TLX1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group18.done'
	shell:
		'touch {output}'
rule parseFP_group19:
	input:
		'{path}footprints/operations/parse/{mergedsample}.DDIT3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.BARH1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.NRF1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HIF1A.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF784.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.GLI2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.SRBP1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXA2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.NME1.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHE40.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.CHD1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.EOMES.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZSCAN31.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHLHB2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYBB.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXB8.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.SPAG7.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.YY2.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB8.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.P73.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group19.done'
	shell:
		'touch {output}'
rule parseFP_group20:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RORB.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.CEBPE.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATF5.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXA13.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.YBOX1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.RXRG.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNF1B.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.THB.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBX20.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAGOH.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.ONEC3.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXD8.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.XBP1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF385A.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.GABPA.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR0B1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXC12.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ESRRB.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO6F2.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.HMBOX1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group20.done'
	shell:
		'touch {output}'
rule parseFP_group21:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ELF5.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX8.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMAD1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PBX2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFYB.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.RORA_1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNF4A.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.RORC.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.NKX21.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHLHA15.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.AP2B.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.HAND1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.ELF1.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEF2C.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.HLF.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.PBX3.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.NKX23.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.HSFY2.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.E2F4.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group21.done'
	shell:
		'touch {output}'
rule parseFP_group22:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GLI1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.HTF4.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.TGIF2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFEC.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBP.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.ONECUT2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.GSC2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HIST1H2BN.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA10.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXA2.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.THAP1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU2F1.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXA1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.MBD2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.BCL6B.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.E2F6.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.CART1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.TF65.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.LHX3.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.GLIS3.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group22.done'
	shell:
		'touch {output}'
rule parseFP_group23:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NFE2L2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.EGR4.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXB1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ONECUT1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.RELA.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.GBX1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.NHLH1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR1H4.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.UBE2V1.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.E2F2.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXB7.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXA5.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB7.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.DUX4.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBR1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.HIC1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.TRIM69.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.SP3.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.DUXA.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group23.done'
	shell:
		'touch {output}'
rule parseFP_group24:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MSC.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.GLI3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.RORA.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXD12.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SPIC.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.PHX2B.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.RELB.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.PHX2A.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.UNC4.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATF1.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.ELK4.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARX.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.EGR1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX7.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.STAU2.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.HESX1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.STAT6.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ESRRG.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.STAT1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.EMX2.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group24.done'
	shell:
		'touch {output}'
rule parseFP_group25:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NEUROD1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZHX3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.RX.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCF7L2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.STAT5A.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.P63.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXD13.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HMX2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.DRGX.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.HIST2H2BE.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRDM1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZIC1.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.NKX62.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.DPRX.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.CEBPZ.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.WRNIP1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.LHX6.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF24.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOSL2.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXD10.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group25.done'
	shell:
		'touch {output}'
rule parseFP_group26:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NEUROD2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXD10.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOSB.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXD11.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.REL.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.PROP1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.EVX2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ELK1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ERG.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZIC3.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.THA.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXP4.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXA3.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXK2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.LRRFIP1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZDHHC15.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.CUX1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFYA.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARNTL.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.ALX1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group26.done'
	shell:
		'touch {output}'
rule parseFP_group27:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZNF740.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR4A3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFE2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.TWST1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.STAT3.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.SUH.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBX21.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO3F2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.HSF2.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ALX4.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR2F2.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATF6.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZFHX3.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU3F3.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.HMBX1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.GBX2.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.CEBPD.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATF2.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.NKX28.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.GFI1B.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group27.done'
	shell:
		'touch {output}'
rule parseFP_group28:
	input:
		'{path}footprints/operations/parse/{mergedsample}.DLX1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFKB1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HMX3.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC13.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHLHE40.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXQ1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.E2F1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATOH1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXH1.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETS2.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEOX1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.GSC.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ASAP3.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU6F2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.MLX.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.SREBF2.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.NOTO.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB3.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.TEAD4.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFAC4.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group28.done'
	shell:
		'touch {output}'
rule parseFP_group29:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GMEB2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.BARX2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.RFX3.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU2F2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RXRB.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.EZH2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.LBX2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.MNT.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFAP2B.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXF2.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.YWHAE.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHLHB3.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA2.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF263.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.STA5B.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF333.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZIC2.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.RUVBL1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.HMX1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.GSX2.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group29.done'
	shell:
		'touch {output}'
rule parseFP_group30:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RUNX1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.DNMT3A.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA11.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.BPTF.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.COT1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.PHOX2B.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.TAGLN2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHLHE41.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ID4.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFATC1.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXB13.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRX3.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXA1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXC2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBM22.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF5.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXD8.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.CEBPA.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXD13.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXF1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group30.done'
	shell:
		'touch {output}'
rule parseFP_group31:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ERR3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.STF1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXP1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ERR2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.EMX1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.JUND.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX6.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.EN2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.STA5A.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.HBP1.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.MRPS25.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.PEBB.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHLHE22.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.CERS4.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXD11.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.XG.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.H2AFZ.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR5A2.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.HEY1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.VSX2.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group31.done'
	shell:
		'touch {output}'
rule parseFP_group32:
	input:
		'{path}footprints/operations/parse/{mergedsample}.FIGLA.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.NCALD.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SREBF1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.TEF.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU6F1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.GSX1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.MTA3.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.EN1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.PITX2.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.BARHL2.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEF2D.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYC.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF8.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.THAP5.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.TWIST1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.JUNB.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB13.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.NEUROG2.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAFK.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.BATF.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group32.done'
	shell:
		'touch {output}'
rule parseFP_group33:
	input:
		'{path}footprints/operations/parse/{mergedsample}.HOXB2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.RCOR1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.OVOL1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.MRPL1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZN333.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.ATF3.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.IKZF2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.PCK2.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.HSPA5.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.GATA1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA13.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.BHLHE23.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMARCC1.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZC3H7A.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.MITF.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEF2A.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.MNX1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.EWSR1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYBL2.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group33.done'
	shell:
		'touch {output}'
rule parseFP_group34:
	input:
		'{path}footprints/operations/parse/{mergedsample}.CEBPB.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOS.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF711.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.H2AFY.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX4.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC12.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.ISL2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.APEX2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.NOBOX.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU5F1B.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.PO2F1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.PITX3.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.OTX1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFAC2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNRNPLL.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU2F3.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU5F1P1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.SEMA4A.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR2E1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.RUNX2.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group34.done'
	shell:
		'touch {output}'
rule parseFP_group35:
	input:
		'{path}footprints/operations/parse/{mergedsample}.FLI1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMUG1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SRP9.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HES5.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.DUSP26.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAX.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HEY2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF7.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.EIF5A2.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU3F4.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.C9orf156.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAFB.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA5.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.OLIG2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFAP4.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.ALX3.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBPL1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMPX.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXC1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.AGGF1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group35.done'
	shell:
		'touch {output}'
rule parseFP_group36:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GLYCTK.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SND1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PPARG.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC10.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TIMM8A.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.NNT.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.RPS4X.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX9.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNRNPA1.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBPJ.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.SNRPB2.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.BACH1.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARNT.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.DDX20.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.SP1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB4.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.MESP1.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.GTF2I.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.VAX2.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF6.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group36.done'
	shell:
		'touch {output}'
rule parseFP_group37:
	input:
		'{path}footprints/operations/parse/{mergedsample}.PTCD1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXB1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ODC1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.JUN.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.PURG.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.ELF2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.MSRA.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HIC2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.GATA5.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRKRIR.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.NMI.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX3.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.OTX2.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.HCLS1.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFE3.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.LARP4.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.FEZF2.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.TEAD3.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.LBX1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.LUZP1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group37.done'
	shell:
		'touch {output}'
rule parseFP_group38:
	input:
		'{path}footprints/operations/parse/{mergedsample}.HOXD12.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SNAPC5.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.NANOG.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXJ3.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TP73.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.XRCC1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEIS3.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.UBE2K.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.CNOT6.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.SSBP3.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFAC3.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.OLIG3.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.IL24.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.MXI1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.TAL1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMAD4.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXC11.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.GFI1.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.GOT1.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group38.done'
	shell:
		'touch {output}'
rule parseFP_group39:
	input:
		'{path}footprints/operations/parse/{mergedsample}.DHX36.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.RAN.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.TRIP10.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.LMX1A.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNRNPA0.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXJ2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.RUNX3.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.FIP1L1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF695.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.CXXC1.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.VAX1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZFP3.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZEB1.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.FHL2.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.PDX1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.ESRP1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.FEV.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF766.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX2.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.JDP2.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group39.done'
	shell:
		'touch {output}'
rule parseFP_group40:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MIXL1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAF.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXO1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYCN.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.FUBP1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.ITF2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.UQCRB.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.LARP1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBX15.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBX4.parseFP.bamcopy10.done', 
		'{path}footprints/operations/parse/{mergedsample}.ANXA1.parseFP.bamcopy11.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBX1.parseFP.bamcopy12.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF384.parseFP.bamcopy13.done', 
		'{path}footprints/operations/parse/{mergedsample}.EDN1.parseFP.bamcopy14.done', 
		'{path}footprints/operations/parse/{mergedsample}.SNAI1.parseFP.bamcopy15.done', 
		'{path}footprints/operations/parse/{mergedsample}.CREB1.parseFP.bamcopy16.done', 
		'{path}footprints/operations/parse/{mergedsample}.BSX.parseFP.bamcopy17.done', 
		'{path}footprints/operations/parse/{mergedsample}.GPD1.parseFP.bamcopy18.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRDX5.parseFP.bamcopy19.done', 
		'{path}footprints/operations/parse/{mergedsample}.RAB14.parseFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group40.done'
	shell:
		'touch {output}'
rule parseFP_group41:
	input:
		'{path}footprints/operations/parse/{mergedsample}.TFDP1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX14.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.DLX4.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR2E3.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.FAM127B.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.LINC00471.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEOX2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFAP2A.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCF12.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB18.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group41.done'
	shell:
		'touch {output}'
rule parseFP_group42:
	input:
		'{path}footprints/operations/parse/{mergedsample}.PDE6H.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.CBX3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBED1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXD1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRNP.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX3.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFEB.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARNT2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.AHR.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.RAB2A.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group42.done'
	shell:
		'touch {output}'
rule parseFP_group43:
	input:
		'{path}footprints/operations/parse/{mergedsample}.FOXO3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.CRX.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZMYND8.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.SSX2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.KDM5A.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.VSX1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.LUZP2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXD4.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.RFXANK.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.PGAM2.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group43.done'
	shell:
		'touch {output}'
rule parseFP_group44:
	input:
		'{path}footprints/operations/parse/{mergedsample}.FOSL1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETV4.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXD4.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.LMX1B.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TSC22D4.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXD9.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HXD9.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR1I3.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.DIABLO.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ESRRA.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group44.done'
	shell:
		'touch {output}'
rule parseFP_group45:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MSX1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.GRHL1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAGEF1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB5.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF124.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.PDCD11.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.TEAD1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYB.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.CFL2.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.CD59.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group45.done'
	shell:
		'touch {output}'
rule parseFP_group46:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GATA4.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.AP1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PDS5A.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NELFB.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXN3.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCF3.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.EZR.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF71.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXO6.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXK1.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group46.done'
	shell:
		'touch {output}'
rule parseFP_group47:
	input:
		'{path}footprints/operations/parse/{mergedsample}.FOXD2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.OLIG1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXI1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NANOS1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ACO1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFAP2C.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.BOLL.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFIA.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.DTL.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.POLI.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group47.done'
	shell:
		'touch {output}'
rule parseFP_group48:
	input:
		'{path}footprints/operations/parse/{mergedsample}.FOXP2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF12.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.INSM1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFIB.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.GRHPR.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.BRCA1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBMS1.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.HIRIP3.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.UNCX.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.RAXL1.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group48.done'
	shell:
		'touch {output}'
rule parseFP_group49:
	input:
		'{path}footprints/operations/parse/{mergedsample}.SHOX2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.LHX9.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.RAX2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.DBP.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF250.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.MRPL2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.ADNP.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ADARB1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.TAF4.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.TMSB4XP8.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group49.done'
	shell:
		'touch {output}'
rule parseFP_group50:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MEIS2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR1I2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZSWIM1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.CREB3L1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.MGA.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.TAF1A.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.USF2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.OTUD4.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.WDR83.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.KDM5D.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group50.done'
	shell:
		'touch {output}'
rule parseFP_group51:
	input:
		'{path}footprints/operations/parse/{mergedsample}.TULP1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PQBP1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.RPP25.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXP3.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.DLX3.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXG1.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.TIMELESS.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCEAL2.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEX3C.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.EEF1D.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group51.done'
	shell:
		'touch {output}'
rule parseFP_group52:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NUCB1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.TSN.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.MXD4.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ISX.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARFGAP1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.TBX5.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.DAB2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.ISL1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.TPI1.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB46.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group52.done'
	shell:
		'touch {output}'
rule parseFP_group53:
	input:
		'{path}footprints/operations/parse/{mergedsample}.PGR.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.GCR.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRGR.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.MSX2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.HTATIP2.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.ESR2.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.MBTPS2.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXO4.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZKSCAN8.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNRNPC.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group53.done'
	shell:
		'touch {output}'
rule parseFP_group54:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NKX25.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.HHEX.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.NUP107.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAFG.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SUCLG1.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOCS4.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXD3.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.SHOX.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.USF1.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF76.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group54.done'
	shell:
		'touch {output}'
rule parseFP_group55:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZNF326.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SRBD1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.DGCR8.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.SNAPC4.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.MDM2.parseFP.bamcopy5.done', 
		'{path}footprints/operations/parse/{mergedsample}.TPPP.parseFP.bamcopy6.done', 
		'{path}footprints/operations/parse/{mergedsample}.NELFA.parseFP.bamcopy7.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXL1.parseFP.bamcopy8.done', 
		'{path}footprints/operations/parse/{mergedsample}.GTPBP6.parseFP.bamcopy9.done', 
		'{path}footprints/operations/parse/{mergedsample}.H1FX.parseFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group55.done'
	shell:
		'touch {output}'
rule parseFP_group56:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MAGEA8.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.HP1BP3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU3F2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.EBF1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.GTF3C5.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group56.done'
	shell:
		'touch {output}'
rule parseFP_group57:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NONO.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.CDX2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYOD1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PPP2R3B.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SNRNP70.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group57.done'
	shell:
		'touch {output}'
rule parseFP_group58:
	input:
		'{path}footprints/operations/parse/{mergedsample}.TOB2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB12.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF6.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ETFB.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.DLX2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group58.done'
	shell:
		'touch {output}'
rule parseFP_group59:
	input:
		'{path}footprints/operations/parse/{mergedsample}.SPIB.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF304.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.VAMP3.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF205.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU5F1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group59.done'
	shell:
		'touch {output}'
rule parseFP_group60:
	input:
		'{path}footprints/operations/parse/{mergedsample}.EHF.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.RAX.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SRY.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.MEIS1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TIA1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group60.done'
	shell:
		'touch {output}'
rule parseFP_group61:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RXRA.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMAD3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ESRP2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NRL.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBM42.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group61.done'
	shell:
		'touch {output}'
rule parseFP_group62:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZNF3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.CPEB1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.CDX1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NF2L1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TGIF2LX.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group62.done'
	shell:
		'touch {output}'
rule parseFP_group63:
	input:
		'{path}footprints/operations/parse/{mergedsample}.PKNOX2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.LEF1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PPARA.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PPP1R10.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SLC18A1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group63.done'
	shell:
		'touch {output}'
rule parseFP_group64:
	input:
		'{path}footprints/operations/parse/{mergedsample}.DIS3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.MAP4K2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.USP39.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.GTF2B.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZCCHC14.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group64.done'
	shell:
		'touch {output}'
rule parseFP_group65:
	input:
		'{path}footprints/operations/parse/{mergedsample}.FOXD3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.EXO5.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.MSI2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.TSNAX.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF696.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group65.done'
	shell:
		'touch {output}'
rule parseFP_group66:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NFATC2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCF4.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.TIMM44.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF3.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.PIK3C3.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group66.done'
	shell:
		'touch {output}'
rule parseFP_group67:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GTF3C2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX15.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PKM.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.RIOK2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF720.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group67.done'
	shell:
		'touch {output}'
rule parseFP_group68:
	input:
		'{path}footprints/operations/parse/{mergedsample}.DLX5.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PIR.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ECSIT.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.IVD.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.NCBP2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group68.done'
	shell:
		'touch {output}'
rule parseFP_group69:
	input:
		'{path}footprints/operations/parse/{mergedsample}.SF3B1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.KLF4.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF160.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF354C.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMCR7L.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group69.done'
	shell:
		'touch {output}'
rule parseFP_group70:
	input:
		'{path}footprints/operations/parse/{mergedsample}.HCFC2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.TAF9.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HHAT.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.RARB.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RHOXF1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group70.done'
	shell:
		'touch {output}'
rule parseFP_group71:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ARI3A.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ARID3A.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR4A2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF830.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.TRMT1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group71.done'
	shell:
		'touch {output}'
rule parseFP_group72:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RFX4.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SRRM3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFIC.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ASPSCR1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ABCF2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group72.done'
	shell:
		'touch {output}'
rule parseFP_group73:
	input:
		'{path}footprints/operations/parse/{mergedsample}.YEATS4.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF26.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SNAI2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HSF1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA4.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group73.done'
	shell:
		'touch {output}'
rule parseFP_group74:
	input:
		'{path}footprints/operations/parse/{mergedsample}.AGAP2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB9.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PDLIM5.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ASCC1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SFT2D1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group74.done'
	shell:
		'touch {output}'
rule parseFP_group75:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RARA.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.BAX.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRRX1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBFOX2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.GIT2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group75.done'
	shell:
		'touch {output}'
rule parseFP_group76:
	input:
		'{path}footprints/operations/parse/{mergedsample}.CBFA2T2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF706.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.RAB7A.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.RPS10.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF655.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group76.done'
	shell:
		'touch {output}'
rule parseFP_group77:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NUP133.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.CBX7.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAXIP1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NMRAL1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.MED30.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group77.done'
	shell:
		'touch {output}'
rule parseFP_group78:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MYEF2.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.FEZ1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.PAX5.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBM8A.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.MSI1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group78.done'
	shell:
		'touch {output}'
rule parseFP_group79:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ANXA11.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.RARG.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB25.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.THRA.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RNASEH2C.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group79.done'
	shell:
		'touch {output}'
rule parseFP_group80:
	input:
		'{path}footprints/operations/parse/{mergedsample}.CAT.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.RFC3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.RPS6KA5.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX5.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFE2L1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group80.done'
	shell:
		'touch {output}'
rule parseFP_group81:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NXPH3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.TCEAL6.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.YY1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PSMA6.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFIL3.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group81.done'
	shell:
		'touch {output}'
rule parseFP_group82:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MZF1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.VDR.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.LHX4.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PITX1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ID2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group82.done'
	shell:
		'touch {output}'
rule parseFP_group83:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ETS1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.LSM6.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.TGIF1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMAD2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.UTP18.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group83.done'
	shell:
		'touch {output}'
rule parseFP_group84:
	input:
		'{path}footprints/operations/parse/{mergedsample}.PLAGL1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.POLE3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF503.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PHOX2A.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.CLK1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group84.done'
	shell:
		'touch {output}'
rule parseFP_group85:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RAB18.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.HNRNPH3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ING3.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ENO1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.U2AF1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group85.done'
	shell:
		'touch {output}'
rule parseFP_group86:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RPL6.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.IKZF1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.P4HB.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.CCDC25.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.NOC2L.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group86.done'
	shell:
		'touch {output}'
rule parseFP_group87:
	input:
		'{path}footprints/operations/parse/{mergedsample}.HDAC8.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.DDX43.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SCAND2P.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.CSNK2B.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.FOXM1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group87.done'
	shell:
		'touch {output}'
rule parseFP_group88:
	input:
		'{path}footprints/operations/parse/{mergedsample}.SCMH1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.POU4F3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.UGP2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.GATA6.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.C19orf40.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group88.done'
	shell:
		'touch {output}'
rule parseFP_group89:
	input:
		'{path}footprints/operations/parse/{mergedsample}.HIST2H2AB.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SMAP2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.TRIM21.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.BCL11A.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.CSTF2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group89.done'
	shell:
		'touch {output}'
rule parseFP_group90:
	input:
		'{path}footprints/operations/parse/{mergedsample}.NR2F1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.CELF5.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.NAP1L1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.UBIP1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.UBP1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group90.done'
	shell:
		'touch {output}'
rule parseFP_group91:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MYF6.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.DLX6.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SPI1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.MTHFD1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBBP5.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group91.done'
	shell:
		'touch {output}'
rule parseFP_group92:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MAGED4B.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PHLDA2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.TFAM.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.LAS1L.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.MYLK.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group92.done'
	shell:
		'touch {output}'
rule parseFP_group93:
	input:
		'{path}footprints/operations/parse/{mergedsample}.AVEN.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZMAT4.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.EVX1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFATC4.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.PICK1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group93.done'
	shell:
		'touch {output}'
rule parseFP_group94:
	input:
		'{path}footprints/operations/parse/{mergedsample}.CPSF4.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.VPS4B.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.R3HDM2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.FGF19.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.AFF4.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group94.done'
	shell:
		'touch {output}'
rule parseFP_group95:
	input:
		'{path}footprints/operations/parse/{mergedsample}.DUSP22.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFATC3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.BARX1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HSPA1L.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.MSRB3.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group95.done'
	shell:
		'touch {output}'
rule parseFP_group96:
	input:
		'{path}footprints/operations/parse/{mergedsample}.CEBPG.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZCCHC17.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.RNF114.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBM3.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.GATA2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group96.done'
	shell:
		'touch {output}'
rule parseFP_group97:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GTPBP1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PLG.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.MECP2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.CELF6.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.GAR1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group97.done'
	shell:
		'touch {output}'
rule parseFP_group98:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GTF2H3.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.MORN1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.NFIX.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.YWHAZ.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBM17.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group98.done'
	shell:
		'touch {output}'
rule parseFP_group99:
	input:
		'{path}footprints/operations/parse/{mergedsample}.CKMT1B.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.GATA3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX10.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PPP5C.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOD1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group99.done'
	shell:
		'touch {output}'
rule parseFP_group100:
	input:
		'{path}footprints/operations/parse/{mergedsample}.PTPMT1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.HMGA1.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.CANX.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.SOX13.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.NR4A1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group100.done'
	shell:
		'touch {output}'
rule parseFP_group101:
	input:
		'{path}footprints/operations/parse/{mergedsample}.ZDHHC5.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZMAT2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.WISP2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.SF1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF510.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group101.done'
	shell:
		'touch {output}'
rule parseFP_group102:
	input:
		'{path}footprints/operations/parse/{mergedsample}.GPAM.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF131.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HMG20A.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF671.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RBBP9.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group102.done'
	shell:
		'touch {output}'
rule parseFP_group103:
	input:
		'{path}footprints/operations/parse/{mergedsample}.UBXN1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.C19orf25.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.STUB1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.PHTF1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.RUFY3.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group103.done'
	shell:
		'touch {output}'
rule parseFP_group104:
	input:
		'{path}footprints/operations/parse/{mergedsample}.METTL21B.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.CELF4.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.A1CF.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.GPANK1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.KIF22.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group104.done'
	shell:
		'touch {output}'
rule parseFP_group105:
	input:
		'{path}footprints/operations/parse/{mergedsample}.HLCS.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.RFC2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.IRF1.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.CDK2AP1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.GADD45A.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group105.done'
	shell:
		'touch {output}'
rule parseFP_group106:
	input:
		'{path}footprints/operations/parse/{mergedsample}.DUS3L.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.LHX2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZNF207.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.KIAA0907.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SPATS2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group106.done'
	shell:
		'touch {output}'
rule parseFP_group107:
	input:
		'{path}footprints/operations/parse/{mergedsample}.MAPK1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.EXOSC3.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.DDX53.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZBTB43.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.GLTPD1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group107.done'
	shell:
		'touch {output}'
rule parseFP_group108:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RNF138.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.UBB.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.BAD.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.DAZAP1.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.SSX3.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group108.done'
	shell:
		'touch {output}'
rule parseFP_group109:
	input:
		'{path}footprints/operations/parse/{mergedsample}.AKR1A1.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.SPR.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.DDX4.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.RPL35.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.ESX1.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group109.done'
	shell:
		'touch {output}'
rule parseFP_group110:
	input:
		'{path}footprints/operations/parse/{mergedsample}.CBFB.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.PRRX2.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.TROVE2.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZRSR2.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.MCTP2.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group110.done'
	shell:
		'touch {output}'
rule parseFP_group111:
	input:
		'{path}footprints/operations/parse/{mergedsample}.RBM7.parseFP.bamcopy1.done', 
		'{path}footprints/operations/parse/{mergedsample}.ZSCAN9.parseFP.bamcopy2.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXA6.parseFP.bamcopy3.done', 
		'{path}footprints/operations/parse/{mergedsample}.HOXB6.parseFP.bamcopy4.done', 
		'{path}footprints/operations/parse/{mergedsample}.NA.parseFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.parseFP.group111.done'
	shell:
		'touch {output}'
