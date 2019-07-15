########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################
#
########################################################################################################################################
#### LNCAP #############################################################################################################################
########################################################################################################################################
rule diffpeaks_lncap:
    input:
        expand("lncap/ex01/wt01/operations/footprints/raw/LNCaP-WT-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/wt02/operations/footprints/raw/LNCaP-WT-02-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr01/operations/footprints/raw/LNCaP-CR-01-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr02/operations/footprints/raw/LNCaP-CR-02-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr04/operations/footprints/raw/LNCaP-CR-04-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr05/operations/footprints/raw/LNCaP-CR-05-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr07/operations/footprints/raw/LNCaP-CR-07-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"]),
        expand("lncap/ex01/cr08/operations/footprints/raw/LNCaP-CR-08-REP1.{genename}.rawFPanalysis.done", genename=config["geneNames"])

####
nohup macs2 callpeak -t LNCaP-CR-01-REP1.u.bam -c LNCaP-WT-01-REP1.u.bam -n cr01_wt01_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-02-REP1.u.bam -c LNCaP-WT-01-REP1.u.bam -n cr02_wt01_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-04-REP1.u.bam -c LNCaP-WT-01-REP1.u.bam -n cr04_wt01_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-05-REP1.u.bam -c LNCaP-WT-01-REP1.u.bam -n cr05_wt01_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-07-REP1.u.bam -c LNCaP-WT-01-REP1.u.bam -n cr07_wt01_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-08-REP1.u.bam -c LNCaP-WT-01-REP1.u.bam -n cr08_wt01_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &

####
nohup macs2 callpeak -t LNCaP-CR-01-REP1.u.bam -c LNCaP-WT-02-REP1.u.bam -n cr01_wt02_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-02-REP1.u.bam -c LNCaP-WT-02-REP1.u.bam -n cr02_wt02_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-04-REP1.u.bam -c LNCaP-WT-02-REP1.u.bam -n cr04_wt02_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-05-REP1.u.bam -c LNCaP-WT-02-REP1.u.bam -n cr05_wt02_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-07-REP1.u.bam -c LNCaP-WT-02-REP1.u.bam -n cr07_wt02_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
nohup macs2 callpeak -t LNCaP-CR-08-REP1.u.bam -c LNCaP-WT-02-REP1.u.bam -n cr08_wt02_globalnorm --outdir peaks --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.05 &
