########################################################################################################################################
#### NOTES #############################################################################################################################
########################################################################################################################################
#
########################################################################################################################################
#### LNCAP #############################################################################################################################
########################################################################################################################################
rule diffpeaks_lncap:
    input:
        "pros/lncap/diffpeaks/ctrl-{sample1}.treat-{sample2}_globalnorm_peaks.narrowPeak"
        a="{path}bam/{sample}-REP{repnum}.bam",



##
nohup macs2 callpeak -t wt01/preprocessing/10unique/LNCaP-WT-01-REP1.u.bam -c cr01/preprocessing/10unique/LNCaP-CR-01-REP1.u.bam -n ctrl-wt01.treat-cr01 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt01/preprocessing/10unique/LNCaP-WT-01-REP1.u.bam -c cr02/preprocessing/10unique/LNCaP-CR-02-REP1.u.bam -n ctrl-wt01.treat-cr02 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt01/preprocessing/10unique/LNCaP-WT-01-REP1.u.bam -c cr04/preprocessing/10unique/LNCaP-CR-04-REP1.u.bam -n ctrl-wt01.treat-cr04 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt01/preprocessing/10unique/LNCaP-WT-01-REP1.u.bam -c cr05/preprocessing/10unique/LNCaP-CR-05-REP1.u.bam -n ctrl-wt01.treat-cr05 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt01/preprocessing/10unique/LNCaP-WT-01-REP1.u.bam -c cr07/preprocessing/10unique/LNCaP-CR-07-REP1.u.bam -n ctrl-wt01.treat-cr07 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt01/preprocessing/10unique/LNCaP-WT-01-REP1.u.bam -c cr08/preprocessing/10unique/LNCaP-CR-08-REP1.u.bam -n ctrl-wt01.treat-cr08 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt02/preprocessing/10unique/LNCaP-WT-02-REP1.u.bam -c cr01/preprocessing/10unique/LNCaP-CR-01-REP1.u.bam -n ctrl-wt02.treat-cr01 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt02/preprocessing/10unique/LNCaP-WT-02-REP1.u.bam -c cr02/preprocessing/10unique/LNCaP-CR-02-REP1.u.bam -n ctrl-wt02.treat-cr02 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt02/preprocessing/10unique/LNCaP-WT-02-REP1.u.bam -c cr04/preprocessing/10unique/LNCaP-CR-04-REP1.u.bam -n ctrl-wt02.treat-cr04 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt02/preprocessing/10unique/LNCaP-WT-02-REP1.u.bam -c cr05/preprocessing/10unique/LNCaP-CR-05-REP1.u.bam -n ctrl-wt02.treat-cr05 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt02/preprocessing/10unique/LNCaP-WT-02-REP1.u.bam -c cr07/preprocessing/10unique/LNCaP-CR-07-REP1.u.bam -n ctrl-wt02.treat-cr07 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
nohup macs2 callpeak -t wt02/preprocessing/10unique/LNCaP-WT-02-REP1.u.bam -c cr08/preprocessing/10unique/LNCaP-CR-08-REP1.u.bam -n ctrl-wt02.treat-cr08 --outdir diffpeaks --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 &
