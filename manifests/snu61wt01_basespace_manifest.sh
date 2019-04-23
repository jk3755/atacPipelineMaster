#!/bin/bash
#
#### DETAILS ####################################################################################################################
#
# Access token: 9e1acf007f014868b78778d705ce1e64
# 
# To find file URLs, log in to basespace, select the run you want to download files from, click "view all" in the samples panel
# Now, click one of the samples. A page open with 8 .fastq.gz files, one for each read in each lane
# Right click each of the files and click "copy link address." Paste the copied links for each file in the table below
#
#### SAMPLES ####################################################################################################################
#
# sample 1 lane 1 read 1 = https://basespace.illumina.com/sample/203218142/files/tree//SNU61-WT-01-S1_S6_L001_R1_001.fastq.gz?id=12822536016
# sample 1 lane 1 read 2 = https://basespace.illumina.com/sample/203218142/files/tree//SNU61-WT-01-S1_S6_L001_R2_001.fastq.gz?id=12822536017
# sample 1 lane 2 read 1 = https://basespace.illumina.com/sample/203218142/files/tree//SNU61-WT-01-S1_S6_L002_R1_001.fastq.gz?id=12822536018
# sample 1 lane 2 read 2 = https://basespace.illumina.com/sample/203218142/files/tree//SNU61-WT-01-S1_S6_L002_R2_001.fastq.gz?id=12822536019
# sample 1 lane 3 read 1 = https://basespace.illumina.com/sample/203218142/files/tree//SNU61-WT-01-S1_S6_L003_R1_001.fastq.gz?id=12822536021
# sample 1 lane 3 read 2 = https://basespace.illumina.com/sample/203218142/files/tree//SNU61-WT-01-S1_S6_L003_R2_001.fastq.gz?id=12822536020
# sample 1 lane 4 read 1 = https://basespace.illumina.com/sample/203218142/files/tree//SNU61-WT-01-S1_S6_L004_R1_001.fastq.gz?id=12822536022
# sample 1 lane 4 read 2 = https://basespace.illumina.com/sample/203218142/files/tree//SNU61-WT-01-S1_S6_L004_R2_001.fastq.gz?id=12822536023
#
# sample 2 lane 1 read 1 = https://basespace.illumina.com/sample/203230029/files/tree//SNU61-WT-01-S2_S4_L001_R1_001.fastq.gz?id=12822520678
# sample 2 lane 1 read 2 = https://basespace.illumina.com/sample/203230029/files/tree//SNU61-WT-01-S2_S4_L001_R2_001.fastq.gz?id=12822520679
# sample 2 lane 2 read 1 = https://basespace.illumina.com/sample/203230029/files/tree//SNU61-WT-01-S2_S4_L002_R1_001.fastq.gz?id=12822520680
# sample 2 lane 2 read 2 = https://basespace.illumina.com/sample/203230029/files/tree//SNU61-WT-01-S2_S4_L002_R2_001.fastq.gz?id=12822520681
# sample 2 lane 3 read 1 = https://basespace.illumina.com/sample/203230029/files/tree//SNU61-WT-01-S2_S4_L003_R1_001.fastq.gz?id=12822520682
# sample 2 lane 3 read 2 = https://basespace.illumina.com/sample/203230029/files/tree//SNU61-WT-01-S2_S4_L003_R2_001.fastq.gz?id=12822520684
# sample 2 lane 4 read 1 = https://basespace.illumina.com/sample/203230029/files/tree//SNU61-WT-01-S2_S4_L004_R1_001.fastq.gz?id=12822520683
# sample 2 lane 4 read 2 = https://basespace.illumina.com/sample/203230029/files/tree//SNU61-WT-01-S2_S4_L004_R2_001.fastq.gz?id=12822520685
#
# sample 3 lane 1 read 1 = https://basespace.illumina.com/sample/203224076/files/tree//SNU61-WT-01-S3_S1_L001_R1_001.fastq.gz?id=12822541273
# sample 3 lane 1 read 2 = https://basespace.illumina.com/sample/203224076/files/tree//SNU61-WT-01-S3_S1_L001_R2_001.fastq.gz?id=12822541274
# sample 3 lane 2 read 1 = https://basespace.illumina.com/sample/203224076/files/tree//SNU61-WT-01-S3_S1_L002_R1_001.fastq.gz?id=12822541275
# sample 3 lane 2 read 2 = https://basespace.illumina.com/sample/203224076/files/tree//SNU61-WT-01-S3_S1_L002_R2_001.fastq.gz?id=12822541277
# sample 3 lane 3 read 1 = https://basespace.illumina.com/sample/203224076/files/tree//SNU61-WT-01-S3_S1_L003_R1_001.fastq.gz?id=12822541276
# sample 3 lane 3 read 2 = https://basespace.illumina.com/sample/203224076/files/tree//SNU61-WT-01-S3_S1_L003_R2_001.fastq.gz?id=12822541278
# sample 3 lane 4 read 1 = https://basespace.illumina.com/sample/203224076/files/tree//SNU61-WT-01-S3_S1_L004_R1_001.fastq.gz?id=12822541280
# sample 3 lane 4 read 2 = https://basespace.illumina.com/sample/203224076/files/tree//SNU61-WT-01-S3_S1_L004_R2_001.fastq.gz?id=12822541279
#
#### INSTRUCTIONS ###############################################################################################################
#
# In the copied link above, look for the file ID at the end of each link, where it says id=XXXXXXXX
# 
# The files can be downloaded with the following shell command:
# wget -O filename 'https://api.basespace.illumina.com/v1pre3/files/{id}/content?access_token={token}'
# replace {id} with each of the 8 file ids, and {token} with the access token above
#
#### IMPORTANT ##################################################################################################################
#
# Launching multiple wget forks at once may cause the downloads to fail, so download the files one at a time
# Use 'wait $!' command, along with sending the wget commands to background (with &) to force the shell to wait for each to finish before starting the next
# $! refers to the PID of the last launched background process
#
#################################################################################################################################
#
# Sample 1
wget -O SNU61-WT-01-S1_S6_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822536016/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S1_S6_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822536017/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S1_S6_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822536018/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S1_S6_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822536019/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S1_S6_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822536021/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S1_S6_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822536020/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S1_S6_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822536022/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S1_S6_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822536023/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
# Sample 2
wget -O SNU61-WT-01-S2_S4_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520678/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S2_S4_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520679/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S2_S4_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520680/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S2_S4_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520681/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S2_S4_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520682/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S2_S4_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520684/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S2_S4_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520683/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S2_S4_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520685/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
# Sample 3
wget -O SNU61-WT-01-S3_S1_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822541273/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S3_S1_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822541274/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S3_S1_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822541275/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S3_S1_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822541277/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S3_S1_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822541276/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S3_S1_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822541278/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S3_S1_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822541280/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O SNU61-WT-01-S3_S1_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822541279/content?access_token=9e1acf007f014868b78778d705ce1e64' &