rule graphFP_group1:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MUSC.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SPZ1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ONEC2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MLXPL.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZSC16.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN350.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BRAC.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF435.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HEN1.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN639.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF350.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN524.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF282.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX9.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB49.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DMC1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF713.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF8.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STAT5B.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB16.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group1.done'
	shell:
		'touch {output}'
rule graphFP_group2:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.TBX19.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PTF1A.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AP2D.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MCR.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHE41.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TF2LX.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TF7L1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN423.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN232.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN143.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBT7A.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXA7.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN589.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX6.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.Pax6.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXC6.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXC8.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KCNIP1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AP2A.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF232.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group2.done'
	shell:
		'touch {output}'
rule graphFP_group3:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.HME1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO5F1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MGAP.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.P5F1B.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC8.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CR3L2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RAD21.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA7.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBT7B.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMC3.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZKSCAN3.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF306.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HSFY1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB7C.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PURA.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF410.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RFX7.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HLX.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC4.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BC11A.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group3.done'
	shell:
		'touch {output}'
rule graphFP_group4:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZEP2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR3C2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZSCAN4.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBT18.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB4.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GMEB1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HCFC1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.P53.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF7.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF691.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OTP.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC5.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF423.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SIX3.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.Trp73.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RHOXF2.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BARH2.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO3F4.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRD14.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SIN3A.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group4.done'
	shell:
		'touch {output}'
rule graphFP_group5:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NKX61.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PSMC2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SIX6.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.COE1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ANDR.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GZF1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXD1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EVI1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SIX5.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB14.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXA9.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB3.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCF21.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF148.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BCL3.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SP100.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SIX4.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC9.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO4F3.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN713.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group5.done'
	shell:
		'touch {output}'
rule graphFP_group6:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZN282.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFCP2L1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO3F3.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO2F3.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXB2.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CHD2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LHX1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LHX5.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BBX.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HDX.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BMAL1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO2F2.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CREB3L2.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PKNX1.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX30.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BSH.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ERR1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO3F1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OSR2.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYBA.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group6.done'
	shell:
		'touch {output}'
rule graphFP_group7:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NDF1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO4F1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO6F1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GRHL2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RORA_2.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATF6A.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CYCS.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GTF2F1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN410.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXC10.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AP2C.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA9.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZEP1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRX4.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HIVEP1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF281.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU4F1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN384.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHE23.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OSR1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group7.done'
	shell:
		'touch {output}'
rule graphFP_group8:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZKSC3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFAP2E.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC6.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRX6.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZKSC1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZSCA4.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMARCC2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RORG.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXC11.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PKNX2.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF274.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBL1XR1.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXJ1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR1H2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EP300.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA3.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TYY1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF9.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRDM14.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group8.done'
	shell:
		'touch {output}'
rule graphFP_group9:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.FOXN1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARI5B.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXC13.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SP8.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXB6.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TLX2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MECOM.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN148.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFAP2D.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CR3L1.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXB3.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETV7.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KAISO.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBX3.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX17.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AIRE.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN740.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB7B.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PIT1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NGN2.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group9.done'
	shell:
		'touch {output}'
rule graphFP_group10:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GLIS1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MIZF.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF143.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF3.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZIC4.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PATZ1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXD3.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX8.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO4F2.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHE22.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRDM4.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBT49.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CENPB.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PPARD.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF524.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR1D1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SIX2.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF75A.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PROX1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SCRT1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group10.done'
	shell:
		'touch {output}'
rule graphFP_group11:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MAZ.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HME2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATF4.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOMEZ.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ASCL1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX12.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCF7.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TAF1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ELF4.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFYC.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NF2L2.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SRBP2.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TYY2.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E2F8.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZFX.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BCL6.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PLAL1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CTCF.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GLIS2.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SCRT2.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group11.done'
	shell:
		'touch {output}'
rule graphFP_group12:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB6.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB33.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF15.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STAT2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TF7L2.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HMGA2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RREB1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF14.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ERF.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAFA.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HINFP1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MTF1.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BARHL1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB7A.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR3C1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RHXF1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SP4.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXA10.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARID5A.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX11.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group12.done'
	shell:
		'touch {output}'
rule graphFP_group13:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.TP63.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.Trp53.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ELF3.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF9.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.T.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SIX1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.THRB.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HES1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E2F7.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CREB5.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PLAG1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF4.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR6A1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX1.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF8.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF2.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN784.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMRC1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SP2.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATF7.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group13.done'
	shell:
		'touch {output}'
rule graphFP_group14:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NKX31.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU1F1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DBX2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ASCL2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EPAS1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNF1A.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXA11.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR2C1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.COT2.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E2F5.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CTCFL.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LHX8.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.REST.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PKNOX1.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYOG.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TP53.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ONECUT3.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF5.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HSF4.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF652.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group14.done'
	shell:
		'touch {output}'
rule graphFP_group15:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.SOX7.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TEAD2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN219.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HES7.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF219.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF16.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HINFP.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AR.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFCP2.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRX5.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PML.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STAT4.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BACH2.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF13.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU4F2.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFKB2.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX18.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CREM.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEF2B.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARID5B.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group15.done'
	shell:
		'touch {output}'
rule graphFP_group16:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NFAC1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RFX2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EBF3.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETV5.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFE2.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EGR2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NDF2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GCM2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E2F3.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BATF3.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CREB3.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.VENTX.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETV3.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRX2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NKX32.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETV2.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CDC5L.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFAT5.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RFX5.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RFX1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group16.done'
	shell:
		'touch {output}'
rule graphFP_group17:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.BHA15.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX4.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.Sox4.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR2F6.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN652.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E4F1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR2C2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETV1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ESR1.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF35.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.WT1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HLTF.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GCM1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NKX22.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EGR3.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCF7L1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETV6.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GABP1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DMBX1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CLOCK.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group17.done'
	shell:
		'touch {output}'
rule graphFP_group18:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZKSCAN1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX21.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CUX2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU3F1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.Myf.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNF6.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.YBX1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF238.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PBX1.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNF4G.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAFF.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBX2.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DMRT3.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MLXIPL.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CYB5R1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SPDEF.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SRF.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYBL1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ELK3.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TLX1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group18.done'
	shell:
		'touch {output}'
rule graphFP_group19:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.DDIT3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BARH1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NRF1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HIF1A.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF784.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GLI2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SRBP1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXA2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NME1.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHE40.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CHD1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EOMES.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZSCAN31.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHLHB2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYBB.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXB8.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SPAG7.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.YY2.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB8.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.P73.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group19.done'
	shell:
		'touch {output}'
rule graphFP_group20:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RORB.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CEBPE.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATF5.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXA13.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.YBOX1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RXRG.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNF1B.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.THB.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBX20.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAGOH.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ONEC3.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXD8.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.XBP1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF385A.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GABPA.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR0B1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXC12.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ESRRB.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO6F2.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HMBOX1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group20.done'
	shell:
		'touch {output}'
rule graphFP_group21:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ELF5.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX8.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMAD1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PBX2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFYB.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RORA_1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNF4A.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RORC.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NKX21.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHLHA15.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AP2B.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HAND1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ELF1.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEF2C.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HLF.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PBX3.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NKX23.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HSFY2.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E2F4.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group21.done'
	shell:
		'touch {output}'
rule graphFP_group22:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GLI1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HTF4.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TGIF2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFEC.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBP.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ONECUT2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GSC2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HIST1H2BN.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA10.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXA2.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.THAP1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU2F1.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXA1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MBD2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BCL6B.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E2F6.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CART1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TF65.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LHX3.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GLIS3.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group22.done'
	shell:
		'touch {output}'
rule graphFP_group23:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NFE2L2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EGR4.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXB1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ONECUT1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RELA.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GBX1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NHLH1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR1H4.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UBE2V1.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E2F2.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXB7.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXA5.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB7.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DUX4.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBR1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HIC1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TRIM69.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SP3.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DUXA.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group23.done'
	shell:
		'touch {output}'
rule graphFP_group24:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MSC.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GLI3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RORA.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXD12.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SPIC.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PHX2B.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RELB.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PHX2A.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UNC4.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATF1.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ELK4.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARX.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EGR1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX7.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STAU2.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HESX1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STAT6.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ESRRG.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STAT1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EMX2.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group24.done'
	shell:
		'touch {output}'
rule graphFP_group25:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NEUROD1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZHX3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RX.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCF7L2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STAT5A.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.P63.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXD13.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HMX2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DRGX.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HIST2H2BE.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRDM1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZIC1.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NKX62.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DPRX.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CEBPZ.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.WRNIP1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LHX6.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF24.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOSL2.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXD10.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group25.done'
	shell:
		'touch {output}'
rule graphFP_group26:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NEUROD2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXD10.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOSB.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXD11.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.REL.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PROP1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EVX2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ELK1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ERG.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZIC3.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.THA.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXP4.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXA3.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXK2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LRRFIP1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZDHHC15.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CUX1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFYA.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARNTL.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ALX1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group26.done'
	shell:
		'touch {output}'
rule graphFP_group27:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZNF740.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR4A3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFE2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TWST1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STAT3.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SUH.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBX21.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO3F2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HSF2.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ALX4.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR2F2.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATF6.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZFHX3.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU3F3.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HMBX1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GBX2.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CEBPD.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATF2.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NKX28.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GFI1B.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group27.done'
	shell:
		'touch {output}'
rule graphFP_group28:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.DLX1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFKB1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HMX3.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC13.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHLHE40.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXQ1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.E2F1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATOH1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXH1.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETS2.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEOX1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GSC.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ASAP3.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU6F2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MLX.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SREBF2.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NOTO.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB3.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TEAD4.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFAC4.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group28.done'
	shell:
		'touch {output}'
rule graphFP_group29:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GMEB2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BARX2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RFX3.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU2F2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RXRB.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EZH2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LBX2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MNT.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFAP2B.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXF2.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.YWHAE.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHLHB3.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA2.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF263.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STA5B.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF333.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZIC2.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RUVBL1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HMX1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GSX2.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group29.done'
	shell:
		'touch {output}'
rule graphFP_group30:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RUNX1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DNMT3A.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA11.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BPTF.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.COT1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PHOX2B.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TAGLN2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHLHE41.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ID4.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFATC1.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXB13.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRX3.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXA1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXC2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBM22.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF5.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXD8.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CEBPA.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXD13.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXF1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group30.done'
	shell:
		'touch {output}'
rule graphFP_group31:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ERR3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STF1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXP1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ERR2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EMX1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.JUND.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX6.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EN2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STA5A.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HBP1.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MRPS25.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PEBB.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHLHE22.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CERS4.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXD11.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.XG.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.H2AFZ.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR5A2.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HEY1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.VSX2.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group31.done'
	shell:
		'touch {output}'
rule graphFP_group32:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.FIGLA.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NCALD.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SREBF1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TEF.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU6F1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GSX1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MTA3.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EN1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PITX2.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BARHL2.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEF2D.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYC.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF8.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.THAP5.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TWIST1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.JUNB.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB13.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NEUROG2.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAFK.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BATF.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group32.done'
	shell:
		'touch {output}'
rule graphFP_group33:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.HOXB2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RCOR1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OVOL1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MRPL1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZN333.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ATF3.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IKZF2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PCK2.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HSPA5.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GATA1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA13.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BHLHE23.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMARCC1.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZC3H7A.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MITF.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEF2A.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MNX1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EWSR1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYBL2.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group33.done'
	shell:
		'touch {output}'
rule graphFP_group34:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.CEBPB.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOS.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF711.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.H2AFY.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX4.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC12.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ISL2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.APEX2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NOBOX.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU5F1B.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PO2F1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PITX3.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OTX1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFAC2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNRNPLL.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU2F3.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU5F1P1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SEMA4A.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR2E1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RUNX2.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group34.done'
	shell:
		'touch {output}'
rule graphFP_group35:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.FLI1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMUG1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SRP9.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HES5.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DUSP26.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAX.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HEY2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF7.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EIF5A2.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU3F4.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.C9orf156.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAFB.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA5.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OLIG2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFAP4.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ALX3.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBPL1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMPX.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXC1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AGGF1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group35.done'
	shell:
		'touch {output}'
rule graphFP_group36:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GLYCTK.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SND1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PPARG.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC10.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TIMM8A.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NNT.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RPS4X.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX9.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNRNPA1.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBPJ.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SNRPB2.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BACH1.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARNT.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DDX20.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SP1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB4.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MESP1.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GTF2I.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.VAX2.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF6.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group36.done'
	shell:
		'touch {output}'
rule graphFP_group37:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.PTCD1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXB1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ODC1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.JUN.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PURG.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ELF2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MSRA.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HIC2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GATA5.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRKRIR.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NMI.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX3.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OTX2.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HCLS1.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFE3.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LARP4.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FEZF2.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TEAD3.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LBX1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LUZP1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group37.done'
	shell:
		'touch {output}'
rule graphFP_group38:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.HOXD12.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SNAPC5.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NANOG.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXJ3.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TP73.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.XRCC1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEIS3.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UBE2K.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CNOT6.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SSBP3.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFAC3.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OLIG3.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IL24.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MXI1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TAL1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMAD4.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXC11.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GFI1.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GOT1.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group38.done'
	shell:
		'touch {output}'
rule graphFP_group39:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.DHX36.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RAN.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TRIP10.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LMX1A.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNRNPA0.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXJ2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RUNX3.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FIP1L1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF695.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CXXC1.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.VAX1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZFP3.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZEB1.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FHL2.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PDX1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ESRP1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FEV.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF766.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX2.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.JDP2.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group39.done'
	shell:
		'touch {output}'
rule graphFP_group40:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MIXL1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAF.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXO1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYCN.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FUBP1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ITF2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UQCRB.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LARP1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBX15.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBX4.graphFP.bamcopy10.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ANXA1.graphFP.bamcopy11.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBX1.graphFP.bamcopy12.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF384.graphFP.bamcopy13.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EDN1.graphFP.bamcopy14.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SNAI1.graphFP.bamcopy15.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CREB1.graphFP.bamcopy16.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BSX.graphFP.bamcopy17.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GPD1.graphFP.bamcopy18.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRDX5.graphFP.bamcopy19.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RAB14.graphFP.bamcopy20.done'
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group40.done'
	shell:
		'touch {output}'
rule graphFP_group41:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.TFDP1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX14.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DLX4.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR2E3.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FAM127B.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LINC00471.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEOX2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFAP2A.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCF12.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB18.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group41.done'
	shell:
		'touch {output}'
rule graphFP_group42:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.PDE6H.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CBX3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBED1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXD1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRNP.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX3.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFEB.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARNT2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AHR.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RAB2A.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group42.done'
	shell:
		'touch {output}'
rule graphFP_group43:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.FOXO3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CRX.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZMYND8.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SSX2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KDM5A.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.VSX1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LUZP2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXD4.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RFXANK.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PGAM2.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group43.done'
	shell:
		'touch {output}'
rule graphFP_group44:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.FOSL1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETV4.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXD4.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LMX1B.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TSC22D4.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXD9.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HXD9.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR1I3.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DIABLO.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ESRRA.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group44.done'
	shell:
		'touch {output}'
rule graphFP_group45:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MSX1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GRHL1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAGEF1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB5.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF124.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PDCD11.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TEAD1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYB.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CFL2.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CD59.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group45.done'
	shell:
		'touch {output}'
rule graphFP_group46:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GATA4.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AP1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PDS5A.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NELFB.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXN3.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCF3.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EZR.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF71.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXO6.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXK1.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group46.done'
	shell:
		'touch {output}'
rule graphFP_group47:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.FOXD2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OLIG1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXI1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NANOS1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ACO1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFAP2C.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BOLL.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFIA.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DTL.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POLI.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group47.done'
	shell:
		'touch {output}'
rule graphFP_group48:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.FOXP2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF12.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.INSM1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFIB.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GRHPR.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BRCA1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBMS1.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HIRIP3.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UNCX.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RAXL1.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group48.done'
	shell:
		'touch {output}'
rule graphFP_group49:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.SHOX2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LHX9.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RAX2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DBP.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF250.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MRPL2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ADNP.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ADARB1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TAF4.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TMSB4XP8.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group49.done'
	shell:
		'touch {output}'
rule graphFP_group50:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MEIS2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR1I2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZSWIM1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CREB3L1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MGA.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TAF1A.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.USF2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.OTUD4.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.WDR83.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KDM5D.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group50.done'
	shell:
		'touch {output}'
rule graphFP_group51:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.TULP1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PQBP1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RPP25.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXP3.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DLX3.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXG1.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TIMELESS.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCEAL2.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEX3C.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EEF1D.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group51.done'
	shell:
		'touch {output}'
rule graphFP_group52:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NUCB1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TSN.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MXD4.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ISX.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARFGAP1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TBX5.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DAB2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ISL1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TPI1.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB46.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group52.done'
	shell:
		'touch {output}'
rule graphFP_group53:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.PGR.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GCR.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRGR.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MSX2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HTATIP2.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ESR2.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MBTPS2.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXO4.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZKSCAN8.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNRNPC.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group53.done'
	shell:
		'touch {output}'
rule graphFP_group54:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NKX25.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HHEX.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NUP107.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAFG.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SUCLG1.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOCS4.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXD3.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SHOX.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.USF1.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF76.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group54.done'
	shell:
		'touch {output}'
rule graphFP_group55:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZNF326.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SRBD1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DGCR8.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SNAPC4.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MDM2.graphFP.bamcopy5.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TPPP.graphFP.bamcopy6.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NELFA.graphFP.bamcopy7.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXL1.graphFP.bamcopy8.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GTPBP6.graphFP.bamcopy9.done', 
		'{path}footprints/operations/graphs/{mergedsample}.H1FX.graphFP.bamcopy10.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group55.done'
	shell:
		'touch {output}'
rule graphFP_group56:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MAGEA8.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HP1BP3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU3F2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EBF1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GTF3C5.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group56.done'
	shell:
		'touch {output}'
rule graphFP_group57:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NONO.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CDX2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYOD1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PPP2R3B.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SNRNP70.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group57.done'
	shell:
		'touch {output}'
rule graphFP_group58:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.TOB2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB12.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF6.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ETFB.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DLX2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group58.done'
	shell:
		'touch {output}'
rule graphFP_group59:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.SPIB.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF304.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.VAMP3.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF205.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU5F1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group59.done'
	shell:
		'touch {output}'
rule graphFP_group60:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.EHF.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RAX.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SRY.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MEIS1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TIA1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group60.done'
	shell:
		'touch {output}'
rule graphFP_group61:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RXRA.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMAD3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ESRP2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NRL.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBM42.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group61.done'
	shell:
		'touch {output}'
rule graphFP_group62:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZNF3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CPEB1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CDX1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NF2L1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TGIF2LX.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group62.done'
	shell:
		'touch {output}'
rule graphFP_group63:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.PKNOX2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LEF1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PPARA.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PPP1R10.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SLC18A1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group63.done'
	shell:
		'touch {output}'
rule graphFP_group64:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.DIS3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MAP4K2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.USP39.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GTF2B.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZCCHC14.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group64.done'
	shell:
		'touch {output}'
rule graphFP_group65:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.FOXD3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EXO5.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MSI2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TSNAX.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF696.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group65.done'
	shell:
		'touch {output}'
rule graphFP_group66:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NFATC2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCF4.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TIMM44.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF3.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PIK3C3.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group66.done'
	shell:
		'touch {output}'
rule graphFP_group67:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GTF3C2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX15.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PKM.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RIOK2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF720.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group67.done'
	shell:
		'touch {output}'
rule graphFP_group68:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.DLX5.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PIR.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ECSIT.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IVD.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NCBP2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group68.done'
	shell:
		'touch {output}'
rule graphFP_group69:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.SF3B1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KLF4.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF160.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF354C.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMCR7L.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group69.done'
	shell:
		'touch {output}'
rule graphFP_group70:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.HCFC2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TAF9.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HHAT.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RARB.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RHOXF1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group70.done'
	shell:
		'touch {output}'
rule graphFP_group71:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ARI3A.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ARID3A.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR4A2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF830.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TRMT1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group71.done'
	shell:
		'touch {output}'
rule graphFP_group72:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RFX4.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SRRM3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFIC.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ASPSCR1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ABCF2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group72.done'
	shell:
		'touch {output}'
rule graphFP_group73:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.YEATS4.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF26.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SNAI2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HSF1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA4.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group73.done'
	shell:
		'touch {output}'
rule graphFP_group74:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.AGAP2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB9.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PDLIM5.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ASCC1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SFT2D1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group74.done'
	shell:
		'touch {output}'
rule graphFP_group75:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RARA.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BAX.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRRX1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBFOX2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GIT2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group75.done'
	shell:
		'touch {output}'
rule graphFP_group76:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.CBFA2T2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF706.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RAB7A.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RPS10.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF655.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group76.done'
	shell:
		'touch {output}'
rule graphFP_group77:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NUP133.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CBX7.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAXIP1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NMRAL1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MED30.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group77.done'
	shell:
		'touch {output}'
rule graphFP_group78:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MYEF2.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FEZ1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PAX5.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBM8A.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MSI1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group78.done'
	shell:
		'touch {output}'
rule graphFP_group79:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ANXA11.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RARG.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB25.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.THRA.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RNASEH2C.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group79.done'
	shell:
		'touch {output}'
rule graphFP_group80:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.CAT.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RFC3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RPS6KA5.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX5.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFE2L1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group80.done'
	shell:
		'touch {output}'
rule graphFP_group81:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NXPH3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TCEAL6.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.YY1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PSMA6.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFIL3.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group81.done'
	shell:
		'touch {output}'
rule graphFP_group82:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MZF1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.VDR.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LHX4.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PITX1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ID2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group82.done'
	shell:
		'touch {output}'
rule graphFP_group83:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ETS1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LSM6.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TGIF1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMAD2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UTP18.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group83.done'
	shell:
		'touch {output}'
rule graphFP_group84:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.PLAGL1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POLE3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF503.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PHOX2A.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CLK1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group84.done'
	shell:
		'touch {output}'
rule graphFP_group85:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RAB18.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HNRNPH3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ING3.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ENO1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.U2AF1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group85.done'
	shell:
		'touch {output}'
rule graphFP_group86:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RPL6.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IKZF1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.P4HB.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CCDC25.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NOC2L.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group86.done'
	shell:
		'touch {output}'
rule graphFP_group87:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.HDAC8.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DDX43.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SCAND2P.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CSNK2B.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FOXM1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group87.done'
	shell:
		'touch {output}'
rule graphFP_group88:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.SCMH1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.POU4F3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UGP2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GATA6.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.C19orf40.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group88.done'
	shell:
		'touch {output}'
rule graphFP_group89:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.HIST2H2AB.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SMAP2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TRIM21.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BCL11A.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CSTF2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group89.done'
	shell:
		'touch {output}'
rule graphFP_group90:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.NR2F1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CELF5.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NAP1L1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UBIP1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UBP1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group90.done'
	shell:
		'touch {output}'
rule graphFP_group91:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MYF6.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DLX6.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SPI1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MTHFD1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBBP5.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group91.done'
	shell:
		'touch {output}'
rule graphFP_group92:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MAGED4B.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PHLDA2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TFAM.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LAS1L.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MYLK.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group92.done'
	shell:
		'touch {output}'
rule graphFP_group93:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.AVEN.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZMAT4.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EVX1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFATC4.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PICK1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group93.done'
	shell:
		'touch {output}'
rule graphFP_group94:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.CPSF4.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.VPS4B.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.R3HDM2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.FGF19.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.AFF4.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group94.done'
	shell:
		'touch {output}'
rule graphFP_group95:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.DUSP22.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFATC3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BARX1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HSPA1L.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MSRB3.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group95.done'
	shell:
		'touch {output}'
rule graphFP_group96:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.CEBPG.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZCCHC17.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RNF114.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBM3.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GATA2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group96.done'
	shell:
		'touch {output}'
rule graphFP_group97:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GTPBP1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PLG.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MECP2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CELF6.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GAR1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group97.done'
	shell:
		'touch {output}'
rule graphFP_group98:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GTF2H3.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MORN1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NFIX.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.YWHAZ.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBM17.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group98.done'
	shell:
		'touch {output}'
rule graphFP_group99:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.CKMT1B.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GATA3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX10.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PPP5C.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOD1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group99.done'
	shell:
		'touch {output}'
rule graphFP_group100:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.PTPMT1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HMGA1.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CANX.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SOX13.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NR4A1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group100.done'
	shell:
		'touch {output}'
rule graphFP_group101:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.ZDHHC5.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZMAT2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.WISP2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SF1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF510.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group101.done'
	shell:
		'touch {output}'
rule graphFP_group102:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.GPAM.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF131.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HMG20A.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF671.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RBBP9.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group102.done'
	shell:
		'touch {output}'
rule graphFP_group103:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.UBXN1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.C19orf25.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.STUB1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PHTF1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RUFY3.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group103.done'
	shell:
		'touch {output}'
rule graphFP_group104:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.METTL21B.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CELF4.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.A1CF.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GPANK1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KIF22.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group104.done'
	shell:
		'touch {output}'
rule graphFP_group105:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.HLCS.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RFC2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.IRF1.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.CDK2AP1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GADD45A.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group105.done'
	shell:
		'touch {output}'
rule graphFP_group106:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.DUS3L.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.LHX2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZNF207.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.KIAA0907.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SPATS2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group106.done'
	shell:
		'touch {output}'
rule graphFP_group107:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.MAPK1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.EXOSC3.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DDX53.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZBTB43.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.GLTPD1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group107.done'
	shell:
		'touch {output}'
rule graphFP_group108:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RNF138.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.UBB.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.BAD.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DAZAP1.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SSX3.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group108.done'
	shell:
		'touch {output}'
rule graphFP_group109:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.AKR1A1.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.SPR.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.DDX4.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.RPL35.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ESX1.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group109.done'
	shell:
		'touch {output}'
rule graphFP_group110:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.CBFB.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.PRRX2.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.TROVE2.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZRSR2.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.MCTP2.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group110.done'
	shell:
		'touch {output}'
rule graphFP_group111:
	input:
		'{path}footprints/operations/graphs/{mergedsample}.RBM7.graphFP.bamcopy1.done', 
		'{path}footprints/operations/graphs/{mergedsample}.ZSCAN9.graphFP.bamcopy2.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXA6.graphFP.bamcopy3.done', 
		'{path}footprints/operations/graphs/{mergedsample}.HOXB6.graphFP.bamcopy4.done', 
		'{path}footprints/operations/graphs/{mergedsample}.NA.graphFP.bamcopy5.done', 
	output:
		'{path}footprints/operations/groups/{mergedsample}.graphFP.group111.done'
	shell:
		'touch {output}'
