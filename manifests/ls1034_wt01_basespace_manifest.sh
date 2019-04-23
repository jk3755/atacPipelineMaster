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
# sample 1 lane 1 read 1 = https://basespace.illumina.com/sample/203217141/files/tree//LS1034-WT-01-S1_S5_L001_R1_001.fastq.gz?id=12822520137
# sample 1 lane 1 read 2 = https://basespace.illumina.com/sample/203217141/files/tree//LS1034-WT-01-S1_S5_L001_R2_001.fastq.gz?id=12822520136
# sample 1 lane 2 read 1 = https://basespace.illumina.com/sample/203217141/files/tree//LS1034-WT-01-S1_S5_L002_R1_001.fastq.gz?id=12822520138
# sample 1 lane 2 read 2 = https://basespace.illumina.com/sample/203217141/files/tree//LS1034-WT-01-S1_S5_L002_R2_001.fastq.gz?id=12822520139
# sample 1 lane 3 read 1 = https://basespace.illumina.com/sample/203217141/files/tree//LS1034-WT-01-S1_S5_L003_R1_001.fastq.gz?id=12822520140
# sample 1 lane 3 read 2 = https://basespace.illumina.com/sample/203217141/files/tree//LS1034-WT-01-S1_S5_L003_R2_001.fastq.gz?id=12822520141
# sample 1 lane 4 read 1 = https://basespace.illumina.com/sample/203217141/files/tree//LS1034-WT-01-S1_S5_L004_R1_001.fastq.gz?id=12822520142
# sample 1 lane 4 read 2 = https://basespace.illumina.com/sample/203217141/files/tree//LS1034-WT-01-S1_S5_L004_R2_001.fastq.gz?id=12822520143
#
# sample 2 lane 1 read 1 = https://basespace.illumina.com/sample/203225074/files/tree//LS1034-WT-01-S2_S3_L001_R1_001.fastq.gz?id=12822522475
# sample 2 lane 1 read 2 = https://basespace.illumina.com/sample/203225074/files/tree//LS1034-WT-01-S2_S3_L001_R2_001.fastq.gz?id=12822522476
# sample 2 lane 2 read 1 = https://basespace.illumina.com/sample/203225074/files/tree//LS1034-WT-01-S2_S3_L002_R1_001.fastq.gz?id=12822522477
# sample 2 lane 2 read 2 = https://basespace.illumina.com/sample/203225074/files/tree//LS1034-WT-01-S2_S3_L002_R2_001.fastq.gz?id=12822522478
# sample 2 lane 3 read 1 = https://basespace.illumina.com/sample/203225074/files/tree//LS1034-WT-01-S2_S3_L003_R1_001.fastq.gz?id=12822522479
# sample 2 lane 3 read 2 = https://basespace.illumina.com/sample/203225074/files/tree//LS1034-WT-01-S2_S3_L003_R2_001.fastq.gz?id=12822522480
# sample 2 lane 4 read 1 = https://basespace.illumina.com/sample/203225074/files/tree//LS1034-WT-01-S2_S3_L004_R1_001.fastq.gz?id=12822522481
# sample 2 lane 4 read 2 = https://basespace.illumina.com/sample/203225074/files/tree//LS1034-WT-01-S2_S3_L004_R2_001.fastq.gz?id=12822522482
#
# sample 3 lane 1 read 1 = https://basespace.illumina.com/sample/203227036/files/tree//LS1034-WT-01-S3_S2_L001_R1_001.fastq.gz?id=12822533112
# sample 3 lane 1 read 2 = https://basespace.illumina.com/sample/203227036/files/tree//LS1034-WT-01-S3_S2_L001_R2_001.fastq.gz?id=12822533111
# sample 3 lane 2 read 1 = https://basespace.illumina.com/sample/203227036/files/tree//LS1034-WT-01-S3_S2_L002_R1_001.fastq.gz?id=12822533113
# sample 3 lane 2 read 2 = https://basespace.illumina.com/sample/203227036/files/tree//LS1034-WT-01-S3_S2_L002_R2_001.fastq.gz?id=12822533114
# sample 3 lane 3 read 1 = https://basespace.illumina.com/sample/203227036/files/tree//LS1034-WT-01-S3_S2_L003_R1_001.fastq.gz?id=12822533115
# sample 3 lane 3 read 2 = https://basespace.illumina.com/sample/203227036/files/tree//LS1034-WT-01-S3_S2_L003_R2_001.fastq.gz?id=12822533116
# sample 3 lane 4 read 1 = https://basespace.illumina.com/sample/203227036/files/tree//LS1034-WT-01-S3_S2_L004_R1_001.fastq.gz?id=12822533117
# sample 3 lane 4 read 2 = https://basespace.illumina.com/sample/203227036/files/tree//LS1034-WT-01-S3_S2_L004_R2_001.fastq.gz?id=12822533118
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
wget -O LS1034-WT-01-S1_S5_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520137/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S1_S5_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520136/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S1_S5_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520138/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S1_S5_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520139/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S1_S5_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520140/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S1_S5_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520141/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S1_S5_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520142/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S1_S5_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822520143/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
# Sample 2
wget -O LS1034-WT-01-S2_S3_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822522475/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S2_S3_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822522476/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S2_S3_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822522477/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S2_S3_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822522478/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S2_S3_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822522479/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S2_S3_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822522480/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S2_S3_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822522481/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S2_S3_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822522482/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
# Sample 3
wget -O LS1034-WT-01-S3_S2_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822533112/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S3_S2_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822533111/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S3_S2_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822533113/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S3_S2_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822533114/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S3_S2_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822533115/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S3_S2_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822533116/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S3_S2_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822533117/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O LS1034-WT-01-S3_S2_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/12822533118/content?access_token=9e1acf007f014868b78778d705ce1e64' &