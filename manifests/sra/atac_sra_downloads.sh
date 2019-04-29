#!/bin/bash
#
## For reference information see here: https://edwards.sdsu.edu/research/getting-data-from-the-sra/
## The SRA toolkit needed to run this script can be found here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
## Use this script to download relevant sequence datasets from NCBIs SRA
## For ATACseq data analysis

#### H508 #######################################################################################################

## H508 WGS: https://www.ncbi.nlm.nih.gov/sra/SRX5466690[accn]
bin/fastq-dump --outdir fastq/h508 --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8670690

## H508 RRBS: https://www.ncbi.nlm.nih.gov/sra/SRX5431921[accn]
bin/fastq-dump --outdir fastq/h508 --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8633497	

## H508 RNAseq: https://www.ncbi.nlm.nih.gov/sra/SRX5414591[accn]
bin/fastq-dump --outdir fastq/h508 --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615809

#### LS1034 #######################################################################################################

## LS1034: https://www.ncbi.nlm.nih.gov/sra/SRX5414357[accn]
bin/fastq-dump --outdir fastq/ls1034 --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615396

## LS1034: https://www.ncbi.nlm.nih.gov/sra/ERX183570[accn]
bin/fastq-dump --outdir fastq/ls1034 --gzip --skip-technical --readids --dumpbase --split-files --clip ERR208905

#### MDST8 #######################################################################################################

## MDST8: https://www.ncbi.nlm.nih.gov/sra/SRX5414820[accn]
bin/fastq-dump --outdir fastq/ls1034 --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8615580

#### Hs675T #######################################################################################################

#### SNU-61 #######################################################################################################

## SNU-61: WXS https://www.ncbi.nlm.nih.gov/sra/SRX5418100[accn] Accession number SRR8618987
bin/fastq-dump --outdir fastq --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8618987

## SNU-61: RNAseq https://www.ncbi.nlm.nih.gov/sra/SRX5415026[accn] Accession number SRR8616024
bin/fastq-dump --outdir fastq --gzip --skip-technical --readids --dumpbase --split-files --clip SRR8616024

#### SNU-16 #######################################################################################################



