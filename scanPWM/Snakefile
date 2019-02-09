# You can run the entire set in the terminal with the following command structure:
# for((i=1;i<=65;i+=1)); do snakemake -j 20 group$i; done
#
# removed
# 'T.sites.Rdata', had waaay too many
# alternative proteins:
# 'JUN (var.2).sites.Rdata', 'JUND (var.2).sites.Rdata', 'JUN (var.2).sites.Rdata', 'JUND (var.2).sites.Rdata', 'MZF1(var.2).sites.Rdata', 'RORA(var.2).sites.Rdata',
# 'MZF1(var.2).sites.Rdata', 'JDP2(var.2).sites.Rdata', 'TFAP2C(var.2).sites.Rdata', 'TFAP2C(var.3).sites.Rdata', 'SREBF2(var.2).sites.Rdata', 'TFAP2A(var.3).sites.Rdata',
# 'JDP2(var.2).sites.Rdata', 'RARA(var.2).sites.Rdata', 'RORA(var.2).sites.Rdata',
# 'TFAP2A(var.2).sites.Rdata', 'TFAP2B(var.2).sites.Rdata', 'TFAP2B(var.3).sites.Rdata', 'TFAP2C(var.2).sites.Rdata', 'TFAP2C(var.3).sites.Rdata', 'SREBF2(var.2).sites.Rdata', 'TFAP2A(var.3).sites.Rdata',
# 'TFAP2A(var.2).sites.Rdata', 'TFAP2B(var.2).sites.Rdata', 'TFAP2B(var.3).sites.Rdata', 'RARA(var.2).sites.Rdata',

rule runmotifscan:
		input:
			"logs/all.done.txt"

###################################################################################################################
rule scanPWM:
        output:
            "{gene}.sites.Rdata"
        script:
            "scripts/scanMotifs.R"
##################################################################################################################

rule group1:
	input:
		'TFAP2A.sites.Rdata', 'NFIL3.sites.Rdata', 'HLF.sites.Rdata', 'NHLH1.sites.Rdata', 'MAX.sites.Rdata', 'USF1.sites.Rdata', 'CEBPA.sites.Rdata', 'EBF1.sites.Rdata', 'CEBPB.sites.Rdata', 'FOS.sites.Rdata',
		'TFAP2A.sites.Rdata', 'NFIL3.sites.Rdata', 'HLF.sites.Rdata', 'NHLH1.sites.Rdata', 'MAX.sites.Rdata', 'USF1.sites.Rdata', 'CEBPA.sites.Rdata', 'EBF1.sites.Rdata', 'CEBPB.sites.Rdata', 'FOS.sites.Rdata',
		'SREBF2.sites.Rdata', 'AHR.sites.Rdata', 'TFAP4.sites.Rdata', 'ARNT.sites.Rdata', 'ATF6.sites.Rdata', 'BACH1.sites.Rdata', 'BACH2.sites.Rdata', 'CREB1.sites.Rdata', 'ATF2.sites.Rdata', 'TCF3.sites.Rdata',
		'SREBF2.sites.Rdata', 'AHR.sites.Rdata', 'TFAP4.sites.Rdata', 'ARNT.sites.Rdata', 'ATF6.sites.Rdata', 'BACH1.sites.Rdata', 'BACH2.sites.Rdata', 'CREB1.sites.Rdata', 'ATF2.sites.Rdata', 'TCF3.sites.Rdata',
		'MYC.sites.Rdata', 'MXI1.sites.Rdata', 'BHLHE40.sites.Rdata', 'ARNTL.sites.Rdata', 'ATF4.sites.Rdata', 'ATF7.sites.Rdata', 'BATF3.sites.Rdata', 'BHLHA15.sites.Rdata', 'BHLHE41.sites.Rdata',
		'BHLHE22.sites.Rdata', 'MYC.sites.Rdata', 'MXI1.sites.Rdata', 'BHLHE40.sites.Rdata', 'ARNTL.sites.Rdata', 'ATF4.sites.Rdata', 'ATF7.sites.Rdata', 'BATF3.sites.Rdata', 'BHLHA15.sites.Rdata',
		'BHLHE41.sites.Rdata', 'BHLHE22.sites.Rdata', 'HES7.sites.Rdata', 'HEY1.sites.Rdata', 'HEY2.sites.Rdata', 'ID4.sites.Rdata', 'JDP2.sites.Rdata', 'MAFG.sites.Rdata', 'MESP1.sites.Rdata',
		'MGA.sites.Rdata', 'MLX.sites.Rdata', 'MLXIPL.sites.Rdata', 'HES7.sites.Rdata', 'HEY1.sites.Rdata', 'HEY2.sites.Rdata', 'ID4.sites.Rdata', 'JDP2.sites.Rdata', 'MAFG.sites.Rdata',
		'MESP1.sites.Rdata', 'MGA.sites.Rdata', 'MLX.sites.Rdata', 'MLXIPL.sites.Rdata', 'TFAP2B.sites.Rdata', 'TFE3.sites.Rdata', 'TFEB.sites.Rdata', 'TFEC.sites.Rdata', 'TFAP2D.sites.Rdata',
		'ARID3A.sites.Rdata', 'ARNT2.sites.Rdata', 'ATF1.sites.Rdata', 'ATF5.sites.Rdata', 'CREM.sites.Rdata', 'TFAP2B.sites.Rdata', 'TFE3.sites.Rdata', 'TFEB.sites.Rdata', 'TFEC.sites.Rdata',
		'TFAP2D.sites.Rdata', 'ARID3A.sites.Rdata', 'ARNT2.sites.Rdata', 'ATF1.sites.Rdata', 'ATF5.sites.Rdata', 'CREM.sites.Rdata', 'MAF.sites.Rdata', 'MITF.sites.Rdata', 'MYOG.sites.Rdata',
		'NEUROD1.sites.Rdata', 'NFE2L2.sites.Rdata', 'PTF1A.sites.Rdata', 'TAL1.sites.Rdata', 'TWIST1.sites.Rdata', 'AIRE.sites.Rdata', 'ALX1.sites.Rdata', 'MAF.sites.Rdata', 'MITF.sites.Rdata',
		'MYOG.sites.Rdata', 'NEUROD1.sites.Rdata', 'NFE2L2.sites.Rdata', 'PTF1A.sites.Rdata', 'TAL1.sites.Rdata', 'TWIST1.sites.Rdata', 'AIRE.sites.Rdata', 'ALX1.sites.Rdata', 'ASCL2.sites.Rdata',
		'ATF6A.sites.Rdata', 'ATOH1.sites.Rdata', 'BARH1.sites.Rdata', 'BARH2.sites.Rdata', 'BARX1.sites.Rdata', 'BARX2.sites.Rdata', 'BC11A.sites.Rdata', 'BCL6B.sites.Rdata', 'BCL6.sites.Rdata',
		'ASCL2.sites.Rdata', 'ATF6A.sites.Rdata', 'ATOH1.sites.Rdata', 'BARH1.sites.Rdata', 'BARH2.sites.Rdata', 'BARX1.sites.Rdata', 'BARX2.sites.Rdata', 'BC11A.sites.Rdata', 'BCL6B.sites.Rdata',
		'BCL6.sites.Rdata', 'CDC5L.sites.Rdata', 'CDX1.sites.Rdata', 'CDX2.sites.Rdata', 'CEBPZ.sites.Rdata', 'CENPB.sites.Rdata', 'COE1.sites.Rdata', 'COT1.sites.Rdata', 'COT2.sites.Rdata',
		'CPEB1.sites.Rdata', 'CR3L1.sites.Rdata', 'CDC5L.sites.Rdata', 'CDX1.sites.Rdata', 'CDX2.sites.Rdata', 'CEBPZ.sites.Rdata', 'CENPB.sites.Rdata', 'COE1.sites.Rdata', 'COT1.sites.Rdata',
		'COT2.sites.Rdata', 'CPEB1.sites.Rdata', 'CR3L1.sites.Rdata', 'DLX3.sites.Rdata', 'DLX4.sites.Rdata', 'DLX5.sites.Rdata', 'DLX6.sites.Rdata', 'DMBX1.sites.Rdata', 'DPRX.sites.Rdata',
		'DRGX.sites.Rdata', 'DUXA.sites.Rdata', 'E2F1.sites.Rdata', 'E2F2.sites.Rdata', 'DLX3.sites.Rdata', 'DLX4.sites.Rdata', 'DLX5.sites.Rdata', 'DLX6.sites.Rdata', 'DMBX1.sites.Rdata',
		'DPRX.sites.Rdata', 'DRGX.sites.Rdata', 'DUXA.sites.Rdata', 'E2F1.sites.Rdata', 'E2F2.sites.Rdata', 'EGR4.sites.Rdata', 'EHF.sites.Rdata', 'ELF1.sites.Rdata', 'ELF2.sites.Rdata',
		'ELF3.sites.Rdata', 'ELF5.sites.Rdata', 'ELK1.sites.Rdata', 'ELK3.sites.Rdata', 'ELK4.sites.Rdata', 'EMX1.sites.Rdata', 'EGR4.sites.Rdata', 'EHF.sites.Rdata', 'ELF1.sites.Rdata',
		'ELF2.sites.Rdata', 'ELF3.sites.Rdata', 'ELF5.sites.Rdata', 'ELK1.sites.Rdata', 'ELK3.sites.Rdata', 'ELK4.sites.Rdata', 'EMX1.sites.Rdata'
	output:
		"logs/group1.done.txt"
	shell:
		"touch {output}"

rule group2:
	input:
		'ETS1.sites.Rdata', 'ETS2.sites.Rdata', 'ETV1.sites.Rdata', 'ETV2.sites.Rdata', 'ETV3.sites.Rdata', 'ETV4.sites.Rdata', 'ETV5.sites.Rdata', 'ETV6.sites.Rdata', 'ETV7.sites.Rdata',
		'EVI1.sites.Rdata', 'ETS1.sites.Rdata', 'ETS2.sites.Rdata', 'ETV1.sites.Rdata', 'ETV2.sites.Rdata', 'ETV3.sites.Rdata', 'ETV4.sites.Rdata', 'ETV5.sites.Rdata', 'ETV6.sites.Rdata',
		'ETV7.sites.Rdata', 'EVI1.sites.Rdata', 'FOXD1.sites.Rdata', 'FOXD2.sites.Rdata', 'FOXD3.sites.Rdata', 'FOXF1.sites.Rdata', 'FOXF2.sites.Rdata', 'FOXG1.sites.Rdata', 'FOXH1.sites.Rdata',
		'FOXI1.sites.Rdata', 'FOXJ2.sites.Rdata', 'FOXJ3.sites.Rdata', 'FOXD1.sites.Rdata', 'FOXD2.sites.Rdata', 'FOXD3.sites.Rdata', 'FOXF1.sites.Rdata', 'FOXF2.sites.Rdata', 'FOXG1.sites.Rdata',
		'FOXH1.sites.Rdata', 'FOXI1.sites.Rdata', 'FOXJ2.sites.Rdata', 'FOXJ3.sites.Rdata', 'FUBP1.sites.Rdata', 'GABP1.sites.Rdata', 'GABPA.sites.Rdata', 'GATA1.sites.Rdata', 'GATA2.sites.Rdata',
		'GATA3.sites.Rdata', 'GATA4.sites.Rdata', 'GATA5.sites.Rdata', 'GATA6.sites.Rdata', 'GBX1.sites.Rdata', 'FUBP1.sites.Rdata', 'GABP1.sites.Rdata', 'GABPA.sites.Rdata', 'GATA1.sites.Rdata',
		'GATA2.sites.Rdata', 'GATA3.sites.Rdata', 'GATA4.sites.Rdata', 'GATA5.sites.Rdata', 'GATA6.sites.Rdata', 'GBX1.sites.Rdata', 'GLIS2.sites.Rdata', 'GLIS3.sites.Rdata', 'GMEB2.sites.Rdata',
		'GRHL1.sites.Rdata', 'GSC2.sites.Rdata', 'GSC.sites.Rdata', 'GSX1.sites.Rdata', 'GSX2.sites.Rdata', 'HBP1.sites.Rdata', 'HEN1.sites.Rdata', 'GLIS2.sites.Rdata', 'GLIS3.sites.Rdata',
		'GMEB2.sites.Rdata', 'GRHL1.sites.Rdata', 'GSC2.sites.Rdata', 'GSC.sites.Rdata', 'GSX1.sites.Rdata', 'GSX2.sites.Rdata', 'HBP1.sites.Rdata', 'HEN1.sites.Rdata', 'HMX3.sites.Rdata',
		'HNF1A.sites.Rdata', 'HNF1B.sites.Rdata', 'HNF4A.sites.Rdata', 'HNF4G.sites.Rdata', 'HNF6.sites.Rdata', 'HOMEZ.sites.Rdata', 'HSF1.sites.Rdata', 'HSF2.sites.Rdata', 'HSF4.sites.Rdata',
		'HMX3.sites.Rdata', 'HNF1A.sites.Rdata', 'HNF1B.sites.Rdata', 'HNF4A.sites.Rdata', 'HNF4G.sites.Rdata', 'HNF6.sites.Rdata', 'HOMEZ.sites.Rdata', 'HSF1.sites.Rdata', 'HSF2.sites.Rdata',
		'HSF4.sites.Rdata', 'HXB13.sites.Rdata', 'HXB1.sites.Rdata', 'HXB2.sites.Rdata', 'HXB3.sites.Rdata', 'HXB6.sites.Rdata', 'HXB7.sites.Rdata', 'HXB8.sites.Rdata', 'HXC10.sites.Rdata',
		'HXC11.sites.Rdata', 'HXC12.sites.Rdata', 'HXB13.sites.Rdata', 'HXB1.sites.Rdata', 'HXB2.sites.Rdata', 'HXB3.sites.Rdata', 'HXB6.sites.Rdata', 'HXB7.sites.Rdata', 'HXB8.sites.Rdata',
		'HXC10.sites.Rdata', 'HXC11.sites.Rdata', 'HXC12.sites.Rdata', 'HXD9.sites.Rdata', 'IKZF1.sites.Rdata', 'INSM1.sites.Rdata', 'IRF1.sites.Rdata', 'IRF2.sites.Rdata', 'IRF3.sites.Rdata',
		'IRF4.sites.Rdata', 'IRF5.sites.Rdata', 'IRF7.sites.Rdata', 'IRF8.sites.Rdata', 'HXD9.sites.Rdata', 'IKZF1.sites.Rdata', 'INSM1.sites.Rdata', 'IRF1.sites.Rdata', 'IRF2.sites.Rdata',
		'IRF3.sites.Rdata', 'IRF4.sites.Rdata', 'IRF5.sites.Rdata', 'IRF7.sites.Rdata', 'IRF8.sites.Rdata', 'KLF14.sites.Rdata', 'KLF15.sites.Rdata', 'KLF16.sites.Rdata', 'KLF1.sites.Rdata',
		'KLF3.sites.Rdata', 'KLF4.sites.Rdata', 'KLF6.sites.Rdata', 'KLF8.sites.Rdata', 'LBX2.sites.Rdata', 'LEF1.sites.Rdata', 'KLF14.sites.Rdata', 'KLF15.sites.Rdata', 'KLF16.sites.Rdata',
		'KLF1.sites.Rdata', 'KLF3.sites.Rdata', 'KLF4.sites.Rdata', 'KLF6.sites.Rdata', 'KLF8.sites.Rdata', 'LBX2.sites.Rdata', 'LEF1.sites.Rdata', 'MCR.sites.Rdata', 'MECP2.sites.Rdata',
		'MEF2A.sites.Rdata', 'MEF2B.sites.Rdata', 'MEF2C.sites.Rdata', 'MEF2D.sites.Rdata', 'MEIS1.sites.Rdata', 'MEIS2.sites.Rdata', 'MEIS3.sites.Rdata', 'MEOX1.sites.Rdata', 'MCR.sites.Rdata',
		'MECP2.sites.Rdata', 'MEF2A.sites.Rdata', 'MEF2B.sites.Rdata', 'MEF2C.sites.Rdata', 'MEF2D.sites.Rdata', 'MEIS1.sites.Rdata', 'MEIS2.sites.Rdata', 'MEIS3.sites.Rdata', 'MEOX1.sites.Rdata',
		'MYBB.sites.Rdata', 'MYB.sites.Rdata', 'MZF1.sites.Rdata', 'NANOG.sites.Rdata', 'NDF1.sites.Rdata', 'NDF2.sites.Rdata', 'NF2L1.sites.Rdata', 'NF2L2.sites.Rdata', 'NFAC1.sites.Rdata',
		'NFAC2.sites.Rdata', 'MYBB.sites.Rdata', 'MYB.sites.Rdata', 'MZF1.sites.Rdata', 'NANOG.sites.Rdata', 'NDF1.sites.Rdata', 'NDF2.sites.Rdata', 'NF2L1.sites.Rdata', 'NF2L2.sites.Rdata',
		'NFAC1.sites.Rdata', 'NFAC2.sites.Rdata'
	output:
		"logs/group2.done.txt"
	shell:
		"touch {output}"

rule group3:
	input:
		'NGN2.sites.Rdata', 'NKX21.sites.Rdata', 'NKX22.sites.Rdata', 'NKX23.sites.Rdata', 'NKX25.sites.Rdata', 'NKX28.sites.Rdata', 'NKX31.sites.Rdata', 'NKX32.sites.Rdata', 'NKX61.sites.Rdata',
		'NKX62.sites.Rdata', 'NGN2.sites.Rdata', 'NKX21.sites.Rdata', 'NKX22.sites.Rdata', 'NKX23.sites.Rdata', 'NKX25.sites.Rdata', 'NKX28.sites.Rdata', 'NKX31.sites.Rdata', 'NKX32.sites.Rdata',
		'NKX61.sites.Rdata', 'NKX62.sites.Rdata', 'NR2E1.sites.Rdata', 'NR2E3.sites.Rdata', 'NR2F6.sites.Rdata', 'NR4A1.sites.Rdata', 'NR4A2.sites.Rdata', 'NR4A3.sites.Rdata', 'NR5A2.sites.Rdata',
		'NR6A1.sites.Rdata', 'NRF1.sites.Rdata', 'ONEC2.sites.Rdata', 'NR2E1.sites.Rdata', 'NR2E3.sites.Rdata', 'NR2F6.sites.Rdata', 'NR4A1.sites.Rdata', 'NR4A2.sites.Rdata', 'NR4A3.sites.Rdata',
		'NR5A2.sites.Rdata', 'NR6A1.sites.Rdata', 'NRF1.sites.Rdata', 'ONEC2.sites.Rdata', 'PAX3.sites.Rdata', 'PAX4.sites.Rdata', 'PAX5.sites.Rdata', 'PAX6.sites.Rdata', 'PAX7.sites.Rdata',
		'PAX8.sites.Rdata', 'PBX1.sites.Rdata', 'PBX2.sites.Rdata', 'PBX3.sites.Rdata', 'PDX1.sites.Rdata', 'PAX3.sites.Rdata', 'PAX4.sites.Rdata', 'PAX5.sites.Rdata', 'PAX6.sites.Rdata',
		'PAX7.sites.Rdata', 'PAX8.sites.Rdata', 'PBX1.sites.Rdata', 'PBX2.sites.Rdata', 'PBX3.sites.Rdata', 'PDX1.sites.Rdata', 'PLAL1.sites.Rdata', 'PO2F1.sites.Rdata', 'PO2F2.sites.Rdata',
		'PO2F3.sites.Rdata', 'PO3F1.sites.Rdata', 'PO3F2.sites.Rdata', 'PO3F3.sites.Rdata', 'PO3F4.sites.Rdata', 'PO4F1.sites.Rdata', 'PO4F2.sites.Rdata', 'PLAL1.sites.Rdata', 'PO2F1.sites.Rdata',
		'PO2F2.sites.Rdata', 'PO2F3.sites.Rdata', 'PO3F1.sites.Rdata', 'PO3F2.sites.Rdata', 'PO3F3.sites.Rdata', 'PO3F4.sites.Rdata', 'PO4F1.sites.Rdata', 'PO4F2.sites.Rdata', 'PRGR.sites.Rdata',
		'PROP1.sites.Rdata', 'PROX1.sites.Rdata', 'PRRX1.sites.Rdata', 'PRRX2.sites.Rdata', 'PURA.sites.Rdata', 'RARA.sites.Rdata', 'RARB.sites.Rdata', 'RARG.sites.Rdata', 'RAX2.sites.Rdata',
		'PRGR.sites.Rdata', 'PROP1.sites.Rdata', 'PROX1.sites.Rdata', 'PRRX1.sites.Rdata', 'PRRX2.sites.Rdata', 'PURA.sites.Rdata', 'RARA.sites.Rdata', 'RARB.sites.Rdata', 'RARG.sites.Rdata',
		'RAX2.sites.Rdata', 'RORG.sites.Rdata', 'RREB1.sites.Rdata', 'RUNX1.sites.Rdata', 'RUNX2.sites.Rdata', 'RUNX3.sites.Rdata', 'RXRA.sites.Rdata', 'RXRB.sites.Rdata', 'RXRG.sites.Rdata',
		'RX.sites.Rdata', 'SCRT1.sites.Rdata', 'RORG.sites.Rdata', 'RREB1.sites.Rdata', 'RUNX1.sites.Rdata', 'RUNX2.sites.Rdata', 'RUNX3.sites.Rdata', 'RXRA.sites.Rdata', 'RXRB.sites.Rdata',
		'RXRG.sites.Rdata', 'RX.sites.Rdata', 'SCRT1.sites.Rdata', 'SOX10.sites.Rdata', 'SOX11.sites.Rdata', 'SOX13.sites.Rdata', 'SOX15.sites.Rdata', 'SOX17.sites.Rdata', 'SOX18.sites.Rdata',
		'SOX1.sites.Rdata', 'SOX21.sites.Rdata', 'SOX2.sites.Rdata', 'SOX3.sites.Rdata', 'SOX10.sites.Rdata', 'SOX11.sites.Rdata', 'SOX13.sites.Rdata', 'SOX15.sites.Rdata', 'SOX17.sites.Rdata',
		'SOX18.sites.Rdata', 'SOX1.sites.Rdata', 'SOX21.sites.Rdata', 'SOX2.sites.Rdata', 'SOX3.sites.Rdata', 'SPI1.sites.Rdata', 'SPIB.sites.Rdata', 'SPIC.sites.Rdata', 'SPZ1.sites.Rdata',
		'SRBP1.sites.Rdata', 'SRBP2.sites.Rdata', 'SRF.sites.Rdata', 'SRY.sites.Rdata', 'STA5A.sites.Rdata', 'STA5B.sites.Rdata', 'SPI1.sites.Rdata', 'SPIB.sites.Rdata', 'SPIC.sites.Rdata',
		'SPZ1.sites.Rdata', 'SRBP1.sites.Rdata', 'SRBP2.sites.Rdata', 'SRF.sites.Rdata', 'SRY.sites.Rdata', 'STA5A.sites.Rdata', 'STA5B.sites.Rdata', 'TBX19.sites.Rdata', 'TBX1.sites.Rdata',
		'TBX20.sites.Rdata', 'TBX21.sites.Rdata', 'TBX2.sites.Rdata', 'TBX3.sites.Rdata', 'TBX4.sites.Rdata', 'TBX5.sites.Rdata', 'TCF7.sites.Rdata', 'TEAD1.sites.Rdata', 'TBX19.sites.Rdata',
		'TBX1.sites.Rdata', 'TBX20.sites.Rdata', 'TBX21.sites.Rdata', 'TBX2.sites.Rdata', 'TBX3.sites.Rdata', 'TBX4.sites.Rdata', 'TBX5.sites.Rdata', 'TCF7.sites.Rdata', 'TEAD1.sites.Rdata',
		'TGIF2.sites.Rdata', 'THAP1.sites.Rdata', 'THA.sites.Rdata', 'THB.sites.Rdata', 'TLX1.sites.Rdata', 'TWST1.sites.Rdata', 'TYY1.sites.Rdata', 'TYY2.sites.Rdata', 'UBIP1.sites.Rdata',
		'UNC4.sites.Rdata', 'TGIF2.sites.Rdata', 'THAP1.sites.Rdata', 'THA.sites.Rdata', 'THB.sites.Rdata', 'TLX1.sites.Rdata', 'TWST1.sites.Rdata', 'TYY1.sites.Rdata', 'TYY2.sites.Rdata',
		'UBIP1.sites.Rdata', 'UNC4.sites.Rdata'
	output:
		"logs/group3.done.txt"
	shell:
		"touch {output}"

rule group4:
	input:
		'ZBT49.sites.Rdata', 'ZBT7A.sites.Rdata', 'ZBT7B.sites.Rdata', 'ZBTB4.sites.Rdata', 'ZBTB6.sites.Rdata', 'ZEB1.sites.Rdata', 'ZEP1.sites.Rdata', 'ZEP2.sites.Rdata', 'ZFHX3.sites.Rdata',
		'ZFX.sites.Rdata', 'ZBT49.sites.Rdata', 'ZBT7A.sites.Rdata', 'ZBT7B.sites.Rdata', 'ZBTB4.sites.Rdata', 'ZBTB6.sites.Rdata', 'ZEB1.sites.Rdata', 'ZEP1.sites.Rdata', 'ZEP2.sites.Rdata',
		'ZFHX3.sites.Rdata', 'ZFX.sites.Rdata', 'ZN282.sites.Rdata', 'ZN333.sites.Rdata', 'ZN350.sites.Rdata', 'ZN384.sites.Rdata', 'ZN410.sites.Rdata', 'ZN423.sites.Rdata', 'ZN524.sites.Rdata',
		'ZN589.sites.Rdata', 'ZN639.sites.Rdata', 'ZN652.sites.Rdata', 'ZN282.sites.Rdata', 'ZN333.sites.Rdata', 'ZN350.sites.Rdata', 'ZN384.sites.Rdata', 'ZN410.sites.Rdata', 'ZN423.sites.Rdata',
		'ZN524.sites.Rdata', 'ZN589.sites.Rdata', 'ZN639.sites.Rdata', 'ZN652.sites.Rdata', 'AGGF1.sites.Rdata', 'AKR1A1.sites.Rdata', 'ANXA1.sites.Rdata', 'ANXA11.sites.Rdata', 'APEX2.sites.Rdata',
		'ARFGAP1.sites.Rdata', 'ASCC1.sites.Rdata', 'ASPSCR1.sites.Rdata', 'AVEN.sites.Rdata', 'BAD.sites.Rdata', 'AGGF1.sites.Rdata', 'AKR1A1.sites.Rdata', 'ANXA1.sites.Rdata', 'ANXA11.sites.Rdata',
		'APEX2.sites.Rdata', 'ARFGAP1.sites.Rdata', 'ASCC1.sites.Rdata', 'ASPSCR1.sites.Rdata', 'AVEN.sites.Rdata', 'BAD.sites.Rdata', 'LINC00471.sites.Rdata', 'C9orf156.sites.Rdata', 'CANX.sites.Rdata',
		'CAT.sites.Rdata', 'CBFA2T2.sites.Rdata', 'CBFB.sites.Rdata', 'CBX7.sites.Rdata', 'ZNF830.sites.Rdata', 'CCDC25.sites.Rdata', 'CD59.sites.Rdata', 'LINC00471.sites.Rdata', 'C9orf156.sites.Rdata',
		'CANX.sites.Rdata', 'CAT.sites.Rdata', 'CBFA2T2.sites.Rdata', 'CBFB.sites.Rdata', 'CBX7.sites.Rdata', 'ZNF830.sites.Rdata', 'CCDC25.sites.Rdata', 'CD59.sites.Rdata', 'CSTF2.sites.Rdata',
		'CYB5R1.sites.Rdata', 'CYCS.sites.Rdata', 'DAB2.sites.Rdata', 'DAZAP1.sites.Rdata', 'ASAP3.sites.Rdata', 'DDX20.sites.Rdata', 'DDX4.sites.Rdata', 'DDX43.sites.Rdata', 'DDX53.sites.Rdata',
		'CSTF2.sites.Rdata', 'CYB5R1.sites.Rdata', 'CYCS.sites.Rdata', 'DAB2.sites.Rdata', 'DAZAP1.sites.Rdata', 'ASAP3.sites.Rdata', 'DDX20.sites.Rdata', 'DDX4.sites.Rdata', 'DDX43.sites.Rdata',
		'DDX53.sites.Rdata', 'ECSIT.sites.Rdata', 'EDN1.sites.Rdata', 'EEF1D.sites.Rdata', 'EIF5A2.sites.Rdata', 'ENO1.sites.Rdata', 'ESRRA.sites.Rdata', 'ETFB.sites.Rdata', 'EWSR1.sites.Rdata',
		'EXOSC3.sites.Rdata', 'METTL21B.sites.Rdata', 'ECSIT.sites.Rdata', 'EDN1.sites.Rdata', 'EEF1D.sites.Rdata', 'EIF5A2.sites.Rdata', 'ENO1.sites.Rdata', 'ESRRA.sites.Rdata', 'ETFB.sites.Rdata',
		'EWSR1.sites.Rdata', 'EXOSC3.sites.Rdata', 'METTL21B.sites.Rdata', 'GLYCTK.sites.Rdata', 'GOT1.sites.Rdata', 'GPAM.sites.Rdata', 'GPD1.sites.Rdata', 'GRHPR.sites.Rdata', 'GTF2B.sites.Rdata',
		'GTF2H3.sites.Rdata', 'GTF3C2.sites.Rdata', 'GTF3C5.sites.Rdata', 'GTPBP1.sites.Rdata', 'GLYCTK.sites.Rdata', 'GOT1.sites.Rdata', 'GPAM.sites.Rdata', 'GPD1.sites.Rdata', 'GRHPR.sites.Rdata',
		'GTF2B.sites.Rdata', 'GTF2H3.sites.Rdata', 'GTF3C2.sites.Rdata', 'GTF3C5.sites.Rdata', 'GTPBP1.sites.Rdata', 'HIRIP3.sites.Rdata', 'HIST1H2BN.sites.Rdata', 'HIST2H2AB.sites.Rdata',
		'HIST2H2BE.sites.Rdata', 'HLCS.sites.Rdata', 'HMG20A.sites.Rdata', 'HNRNPA0.sites.Rdata', 'HNRNPA1.sites.Rdata', 'HNRNPC.sites.Rdata', 'HNRNPH3.sites.Rdata', 'HIRIP3.sites.Rdata',
		'HIST1H2BN.sites.Rdata', 'HIST2H2AB.sites.Rdata', 'HIST2H2BE.sites.Rdata', 'HLCS.sites.Rdata', 'HMG20A.sites.Rdata', 'HNRNPA0.sites.Rdata', 'HNRNPA1.sites.Rdata', 'HNRNPC.sites.Rdata',
		'HNRNPH3.sites.Rdata', 'ING3.sites.Rdata', 'IRF6.sites.Rdata', 'IVD.sites.Rdata', 'KDM5A.sites.Rdata', 'KDM5D.sites.Rdata', 'KCNIP1.sites.Rdata', 'KIAA0907.sites.Rdata', 'KIF22.sites.Rdata',
		'LARP1.sites.Rdata', 'LARP4.sites.Rdata', 'ING3.sites.Rdata', 'IRF6.sites.Rdata', 'IVD.sites.Rdata', 'KDM5A.sites.Rdata', 'KDM5D.sites.Rdata', 'KCNIP1.sites.Rdata', 'KIAA0907.sites.Rdata',
		'KIF22.sites.Rdata', 'LARP1.sites.Rdata', 'LARP4.sites.Rdata', 'MAGEF1.sites.Rdata', 'MAGOH.sites.Rdata', 'MAP4K2.sites.Rdata', 'MAPK1.sites.Rdata', 'MBTPS2.sites.Rdata', 'MCTP2.sites.Rdata',
		'MDM2.sites.Rdata', 'MEF2BNB-MEF2B.sites.Rdata', 'GLTPD1.sites.Rdata', 'RBM42.sites.Rdata', 'MAGEF1.sites.Rdata', 'MAGOH.sites.Rdata', 'MAP4K2.sites.Rdata', 'MAPK1.sites.Rdata', 'MBTPS2.sites.Rdata',
		'MCTP2.sites.Rdata', 'MDM2.sites.Rdata', 'MEF2BNB-MEF2B.sites.Rdata', 'GLTPD1.sites.Rdata', 'RBM42.sites.Rdata'
	output:
		"logs/group4.done.txt"
	shell:
		"touch {output}"

rule group5:
	input:
		'MXD4.sites.Rdata', 'MYEF2.sites.Rdata', 'MYLK.sites.Rdata', 'NANOS1.sites.Rdata', 'NAP1L1.sites.Rdata', 'NCALD.sites.Rdata', 'NCBP2.sites.Rdata', 'NFATC3.sites.Rdata', 'NFATC4.sites.Rdata',
		'NFIB.sites.Rdata', 'MXD4.sites.Rdata', 'MYEF2.sites.Rdata', 'MYLK.sites.Rdata', 'NANOS1.sites.Rdata', 'NAP1L1.sites.Rdata', 'NCALD.sites.Rdata', 'NCBP2.sites.Rdata', 'NFATC3.sites.Rdata',
		'NFATC4.sites.Rdata', 'NFIB.sites.Rdata', 'NUCB1.sites.Rdata', 'NUP107.sites.Rdata', 'NUP133.sites.Rdata', 'NXPH3.sites.Rdata', 'ODC1.sites.Rdata', 'OTUD4.sites.Rdata', 'P4HB.sites.Rdata',
		'PAXIP1.sites.Rdata', 'PCK2.sites.Rdata', 'PDCD11.sites.Rdata', 'NUCB1.sites.Rdata', 'NUP107.sites.Rdata', 'NUP133.sites.Rdata', 'NXPH3.sites.Rdata', 'ODC1.sites.Rdata', 'OTUD4.sites.Rdata',
		'P4HB.sites.Rdata', 'PAXIP1.sites.Rdata', 'PCK2.sites.Rdata', 'PDCD11.sites.Rdata', 'PKNOX2.sites.Rdata', 'PLAGL1.sites.Rdata', 'PLG.sites.Rdata', 'POLE3.sites.Rdata', 'POLI.sites.Rdata',
		'POU3F2.sites.Rdata', 'POU4F3.sites.Rdata', 'PPP1R10.sites.Rdata', 'PPP2R3B.sites.Rdata', 'PPP5C.sites.Rdata', 'PKNOX2.sites.Rdata', 'PLAGL1.sites.Rdata', 'PLG.sites.Rdata', 'POLE3.sites.Rdata',
		'POLI.sites.Rdata', 'POU3F2.sites.Rdata', 'POU4F3.sites.Rdata', 'PPP1R10.sites.Rdata', 'PPP2R3B.sites.Rdata', 'PPP5C.sites.Rdata', 'RAB14.sites.Rdata', 'RAB18.sites.Rdata', 'RAB2A.sites.Rdata',
		'RAB7A.sites.Rdata', 'RAN.sites.Rdata', 'RAX.sites.Rdata', 'RBBP5.sites.Rdata', 'RBBP9.sites.Rdata', 'RBM17.sites.Rdata', 'RBM22.sites.Rdata', 'RAB14.sites.Rdata', 'RAB18.sites.Rdata',
		'RAB2A.sites.Rdata', 'RAB7A.sites.Rdata', 'RAN.sites.Rdata', 'RAX.sites.Rdata', 'RBBP5.sites.Rdata', 'RBBP9.sites.Rdata', 'RBM17.sites.Rdata', 'RBM22.sites.Rdata', 'RIOK2.sites.Rdata',
		'MEX3C.sites.Rdata', 'RNASEH2C.sites.Rdata', 'RNF138.sites.Rdata', 'RPL35.sites.Rdata', 'RPL6.sites.Rdata', 'RPP25.sites.Rdata', 'RPS10.sites.Rdata', 'RPS4X.sites.Rdata', 'RPS6KA5.sites.Rdata',
		'RIOK2.sites.Rdata', 'MEX3C.sites.Rdata', 'RNASEH2C.sites.Rdata', 'RNF138.sites.Rdata', 'RPL35.sites.Rdata', 'RPL6.sites.Rdata', 'RPP25.sites.Rdata', 'RPS10.sites.Rdata', 'RPS4X.sites.Rdata',
		'RPS6KA5.sites.Rdata', 'SMAP2.sites.Rdata', 'SMCR7L.sites.Rdata', 'SMPX.sites.Rdata', 'SMUG1.sites.Rdata', 'SNAPC4.sites.Rdata', 'SNAPC5.sites.Rdata', 'SND1.sites.Rdata', 'SNRNP70.sites.Rdata',
		'SNRPB2.sites.Rdata', 'SOCS4.sites.Rdata', 'SMAP2.sites.Rdata', 'SMCR7L.sites.Rdata', 'SMPX.sites.Rdata', 'SMUG1.sites.Rdata', 'SNAPC4.sites.Rdata', 'SNAPC5.sites.Rdata', 'SND1.sites.Rdata',
		'SNRNP70.sites.Rdata', 'SNRPB2.sites.Rdata', 'SOCS4.sites.Rdata', 'STAU2.sites.Rdata', 'STUB1.sites.Rdata', 'SUCLG1.sites.Rdata', 'TAF1A.sites.Rdata', 'TAF9.sites.Rdata', 'TAGLN2.sites.Rdata',
		'TBPL1.sites.Rdata', 'TCEAL2.sites.Rdata', 'TCEAL6.sites.Rdata', 'TFAM.sites.Rdata', 'STAU2.sites.Rdata', 'STUB1.sites.Rdata', 'SUCLG1.sites.Rdata', 'TAF1A.sites.Rdata', 'TAF9.sites.Rdata',
		'TAGLN2.sites.Rdata', 'TBPL1.sites.Rdata', 'TCEAL2.sites.Rdata', 'TCEAL6.sites.Rdata', 'TFAM.sites.Rdata', 'TP73.sites.Rdata', 'TPI1.sites.Rdata', 'TPPP.sites.Rdata', 'TRIM21.sites.Rdata',
		'TRIM69.sites.Rdata', 'TRIP10.sites.Rdata', 'TRMT1.sites.Rdata', 'TROVE2.sites.Rdata', 'TSC22D4.sites.Rdata', 'TSN.sites.Rdata', 'TP73.sites.Rdata', 'TPI1.sites.Rdata', 'TPPP.sites.Rdata',
		'TRIM21.sites.Rdata', 'TRIM69.sites.Rdata', 'TRIP10.sites.Rdata', 'TRMT1.sites.Rdata', 'TROVE2.sites.Rdata', 'TSC22D4.sites.Rdata', 'TSN.sites.Rdata', 'EZR.sites.Rdata', 'VPS4B.sites.Rdata',
		'NELFA.sites.Rdata', 'WISP2.sites.Rdata', 'XG.sites.Rdata', 'XRCC1.sites.Rdata', 'YEATS4.sites.Rdata', 'YWHAE.sites.Rdata', 'YWHAZ.sites.Rdata', 'ZBTB12.sites.Rdata', 'EZR.sites.Rdata',
		'VPS4B.sites.Rdata', 'NELFA.sites.Rdata', 'WISP2.sites.Rdata', 'XG.sites.Rdata', 'XRCC1.sites.Rdata', 'YEATS4.sites.Rdata', 'YWHAE.sites.Rdata', 'YWHAZ.sites.Rdata', 'ZBTB12.sites.Rdata',
		'ZMAT2.sites.Rdata', 'ZMAT4.sites.Rdata', 'ZNF124.sites.Rdata', 'ZNF131.sites.Rdata', 'ZNF160.sites.Rdata', 'ZKSCAN8.sites.Rdata', 'ZSCAN9.sites.Rdata', 'ZNF205.sites.Rdata', 'ZNF207.sites.Rdata',
		'ZBTB18.sites.Rdata', 'ZMAT2.sites.Rdata', 'ZMAT4.sites.Rdata', 'ZNF124.sites.Rdata', 'ZNF131.sites.Rdata', 'ZNF160.sites.Rdata', 'ZKSCAN8.sites.Rdata', 'ZSCAN9.sites.Rdata', 'ZNF205.sites.Rdata',
		'ZNF207.sites.Rdata', 'ZBTB18.sites.Rdata'
	output:
		"logs/group5.done.txt"
	shell:
		"touch {output}"

rule group6:
	input:
		'ZNF655.sites.Rdata', 'ZNF671.sites.Rdata', 'ZNF695.sites.Rdata', 'ZNF706.sites.Rdata', 'ZNF71.sites.Rdata', 'ZNF720.sites.Rdata', 'ZNF76.sites.Rdata', 'ZNF766.sites.Rdata', 'ZRSR2.sites.Rdata',
		'ZSWIM1.sites.Rdata', 'ZNF655.sites.Rdata', 'ZNF671.sites.Rdata', 'ZNF695.sites.Rdata', 'ZNF706.sites.Rdata', 'ZNF71.sites.Rdata', 'ZNF720.sites.Rdata', 'ZNF76.sites.Rdata', 'ZNF766.sites.Rdata',
		'ZRSR2.sites.Rdata', 'ZSWIM1.sites.Rdata', 'TLX1::NFIC.sites.Rdata', 'NKX3-1.sites.Rdata', 'ZNF354C.sites.Rdata', 'EWSR1-FLI1.sites.Rdata', 'RXR::RAR_DR5.sites.Rdata', 'BATF::JUN.sites.Rdata',
		'DUX4.sites.Rdata', 'FOXP1.sites.Rdata', 'TLX1::NFIC.sites.Rdata', 'NKX3-1.sites.Rdata', 'ZNF354C.sites.Rdata', 'EWSR1-FLI1.sites.Rdata', 'RXR::RAR_DR5.sites.Rdata', 'BATF::JUN.sites.Rdata',
		'DUX4.sites.Rdata', 'FOXP1.sites.Rdata', 'TP53.sites.Rdata', 'YY1.sites.Rdata', 'KLF5.sites.Rdata', 'EN1.sites.Rdata', 'MAX::MYC.sites.Rdata', 'PPARG::RXRA.sites.Rdata', 'ZNF143.sites.Rdata',
		'TP53.sites.Rdata', 'YY1.sites.Rdata', 'KLF5.sites.Rdata', 'EN1.sites.Rdata', 'MAX::MYC.sites.Rdata', 'PPARG::RXRA.sites.Rdata', 'ZNF143.sites.Rdata', 'STAT1::STAT2.sites.Rdata', 'DMRT3.sites.Rdata',
		'LBX1.sites.Rdata', 'POU6F1.sites.Rdata', 'BARHL2.sites.Rdata', 'ELF4.sites.Rdata', 'EN2.sites.Rdata', 'HOXA13.sites.Rdata', 'HOXC11.sites.Rdata',  'STAT1::STAT2.sites.Rdata', 'DMRT3.sites.Rdata',
		'LBX1.sites.Rdata', 'POU6F1.sites.Rdata', 'BARHL2.sites.Rdata', 'ELF4.sites.Rdata', 'EN2.sites.Rdata', 'HOXA13.sites.Rdata', 'HOXC11.sites.Rdata', 'SP8.sites.Rdata', 'YY2.sites.Rdata', 'ZBTB7A.sites.Rdata',
		'ZNF410.sites.Rdata', 'ZNF740.sites.Rdata', 'ONECUT2.sites.Rdata', 'ONECUT3.sites.Rdata', 'MYBL1.sites.Rdata', 'MYBL2.sites.Rdata', 'SP8.sites.Rdata', 'YY2.sites.Rdata', 'ZBTB7A.sites.Rdata', 'ZNF410.sites.Rdata',
		'ZNF740.sites.Rdata', 'ONECUT2.sites.Rdata', 'ONECUT3.sites.Rdata', 'MYBL1.sites.Rdata', 'MYBL2.sites.Rdata', 'HOXD12.sites.Rdata', 'BSX.sites.Rdata', 'HMBOX1.sites.Rdata',
		'HOXD12.sites.Rdata', 'BSX.sites.Rdata', 'HMBOX1.sites.Rdata', 'MIZF.sites.Rdata', 'AP1.sites.Rdata', 'HIF1A::ARNT.sites.Rdata', 'HINFP1.sites.Rdata', 'ZNF238.sites.Rdata', 'ZNF306.sites.Rdata',
		'ZNF524.sites.Rdata', 'ZNF75A.sites.Rdata', 'ZNF784.sites.Rdata', 'ZSCAN4.sites.Rdata', 'MIZF.sites.Rdata', 'AP1.sites.Rdata', 'HIF1A::ARNT.sites.Rdata', 'HINFP1.sites.Rdata', 'ZNF238.sites.Rdata',
		'ZNF306.sites.Rdata', 'ZNF524.sites.Rdata', 'ZNF75A.sites.Rdata', 'ZNF784.sites.Rdata', 'ZSCAN4.sites.Rdata', 'HOXD8.sites.Rdata', 'IRX5.sites.Rdata', 'PHOX2B.sites.Rdata', 'RAXL1.sites.Rdata',
		'ESRRG.sites.Rdata', 'THRB.sites.Rdata', 'Trp53.sites.Rdata', 'Trp73.sites.Rdata', 'ZBTB49.sites.Rdata', 'ZNF232.sites.Rdata', 'HOXD8.sites.Rdata', 'IRX5.sites.Rdata', 'PHOX2B.sites.Rdata',
		'RAXL1.sites.Rdata', 'ESRRG.sites.Rdata', 'THRB.sites.Rdata', 'Trp53.sites.Rdata', 'Trp73.sites.Rdata', 'ZBTB49.sites.Rdata', 'ZNF232.sites.Rdata', 'CREB3L2.sites.Rdata', 'DBX2.sites.Rdata',
		'DMC1.sites.Rdata', 'EBF3.sites.Rdata', 'EP300.sites.Rdata', 'EZH2.sites.Rdata', 'FOXJ1.sites.Rdata', 'FOXN1.sites.Rdata', 'GMEB1.sites.Rdata', 'GTF2F1.sites.Rdata', 'CREB3L2.sites.Rdata',
		'DBX2.sites.Rdata', 'DMC1.sites.Rdata', 'EBF3.sites.Rdata', 'EP300.sites.Rdata', 'EZH2.sites.Rdata', 'FOXJ1.sites.Rdata', 'FOXN1.sites.Rdata', 'GMEB1.sites.Rdata', 'GTF2F1.sites.Rdata',
		'HOXA7.sites.Rdata', 'HOXA9.sites.Rdata', 'HOXB1.sites.Rdata', 'HOXB4.sites.Rdata', 'HOXB6.sites.Rdata', 'HOXB7.sites.Rdata', 'HOXB8.sites.Rdata', 'HOXC4.sites.Rdata', 'HOXC5.sites.Rdata',
		'HOXC6.sites.Rdata', 'HOXA7.sites.Rdata', 'HOXA9.sites.Rdata', 'HOXB1.sites.Rdata', 'HOXB4.sites.Rdata', 'HOXB6.sites.Rdata', 'HOXB7.sites.Rdata', 'HOXB8.sites.Rdata', 'HOXC4.sites.Rdata',
		'HOXC5.sites.Rdata', 'HOXC6.sites.Rdata'
	output:
		"logs/group6.done.txt"
	shell:
		"touch {output}"

rule group7:
	input:
		'LHX1.sites.Rdata', 'LHX5.sites.Rdata', 'MECOM.sites.Rdata', 'MTA3.sites.Rdata', 'NFATC2.sites.Rdata', 'NKX1-1.sites.Rdata', 'NKX1-2.sites.Rdata', 'NKX2-1.sites.Rdata', 'NKX2-2.sites.Rdata',
		'NKX2-4.sites.Rdata', 'LHX1.sites.Rdata', 'LHX5.sites.Rdata', 'MECOM.sites.Rdata', 'MTA3.sites.Rdata', 'NFATC2.sites.Rdata', 'NKX1-1.sites.Rdata', 'NKX1-2.sites.Rdata', 'NKX2-1.sites.Rdata',
		'NKX2-2.sites.Rdata', 'NKX2-4.sites.Rdata', 'POU5F1.sites.Rdata', 'PRDM14.sites.Rdata', 'RAD21.sites.Rdata', 'RBPJ.sites.Rdata', 'RCOR1.sites.Rdata', 'RFX7.sites.Rdata', 'RHOXF2.sites.Rdata',
		'RORC.sites.Rdata', 'SIN3A.sites.Rdata', 'SIX1.sites.Rdata', 'POU5F1.sites.Rdata', 'PRDM14.sites.Rdata', 'RAD21.sites.Rdata', 'RBPJ.sites.Rdata', 'RCOR1.sites.Rdata', 'RFX7.sites.Rdata',
		'RHOXF2.sites.Rdata', 'RORC.sites.Rdata', 'SIN3A.sites.Rdata', 'SIX1.sites.Rdata', 'SOX6.sites.Rdata', 'SP100.sites.Rdata', 'STAT5A.sites.Rdata', 'STAT5B.sites.Rdata', 'TAF1.sites.Rdata',
		'TBL1XR1.sites.Rdata', 'TCF21.sites.Rdata', 'TFAP2E.sites.Rdata', 'TFCP2L1.sites.Rdata', 'TLX2.sites.Rdata', 'SOX6.sites.Rdata', 'SP100.sites.Rdata', 'STAT5A.sites.Rdata', 'STAT5B.sites.Rdata',
		'TAF1.sites.Rdata', 'TBL1XR1.sites.Rdata', 'TCF21.sites.Rdata', 'TFAP2E.sites.Rdata', 'TFCP2L1.sites.Rdata', 'TLX2.sites.Rdata', 'ZNF274.sites.Rdata', 'ZNF281.sites.Rdata', 'ZNF333.sites.Rdata',
		'ZNF350.sites.Rdata', 'ZNF35.sites.Rdata', 'ZNF384.sites.Rdata', 'ZNF423.sites.Rdata', 'ZNF652.sites.Rdata', 'ZNF691.sites.Rdata', 'ZNF711.sites.Rdata', 'ZNF274.sites.Rdata', 'ZNF281.sites.Rdata',
		'ZNF333.sites.Rdata', 'ZNF350.sites.Rdata', 'ZNF35.sites.Rdata', 'ZNF384.sites.Rdata', 'ZNF423.sites.Rdata', 'ZNF652.sites.Rdata', 'ZNF691.sites.Rdata', 'ZNF711.sites.Rdata'
	output:
		"logs/group7.done.txt"
	shell:
		"touch {output}"

rule AGGREGATE_groups:
	input:
		"logs/group1.done.txt",
		"logs/group2.done.txt",
		"logs/group3.done.txt",
		"logs/group4.done.txt",
		"logs/group5.done.txt",
		"logs/group6.done.txt",
		"logs/group7.done.txt"
	output:
		"logs/all.done.txt"
	shell:
		"touch {output}"
