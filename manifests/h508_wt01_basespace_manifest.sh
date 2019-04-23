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
# sample 1 lane 1 read 1 = https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L001_R1_001.fastq.gz?id=10172213598
# sample 1 lane 1 read 2 = https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L001_R2_001.fastq.gz?id=10172213599
# sample 1 lane 2 read 1 = https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L002_R1_001.fastq.gz?id=10172213596
# sample 1 lane 2 read 2 = https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L002_R2_001.fastq.gz?id=10172213597
# sample 1 lane 3 read 1 = https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L003_R1_001.fastq.gz?id=10172213601
# sample 1 lane 3 read 2 = https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L003_R2_001.fastq.gz?id=10172213602
# sample 1 lane 4 read 1 = https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L004_R1_001.fastq.gz?id=10172213600
# sample 1 lane 4 read 2 = https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L004_R2_001.fastq.gz?id=10172213603
#
# sample 2 lane 1 read 1 = https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L001_R1_001.fastq.gz?id=10172226239
# sample 2 lane 1 read 2 = https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L001_R2_001.fastq.gz?id=10172226240
# sample 2 lane 2 read 1 = https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L002_R1_001.fastq.gz?id=10172226242
# sample 2 lane 2 read 2 = https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L002_R2_001.fastq.gz?id=10172226241
# sample 2 lane 3 read 1 = https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L003_R1_001.fastq.gz?id=10172226243
# sample 2 lane 3 read 2 = https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L003_R2_001.fastq.gz?id=10172226244
# sample 2 lane 4 read 1 = https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L004_R1_001.fastq.gz?id=10172226245
# sample 2 lane 4 read 2 = https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L004_R2_001.fastq.gz?id=10172226246
#
# sample 3 lane 1 read 1 = https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L001_R1_001.fastq.gz?id=10172224329
# sample 3 lane 1 read 2 = https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L001_R2_001.fastq.gz?id=10172224328
# sample 3 lane 2 read 1 = https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L002_R1_001.fastq.gz?id=10172224331
# sample 3 lane 2 read 2 = https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L002_R2_001.fastq.gz?id=10172224330
# sample 3 lane 3 read 1 = https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L003_R1_001.fastq.gz?id=10172224332
# sample 3 lane 3 read 2 = https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L003_R2_001.fastq.gz?id=10172224334
# sample 3 lane 4 read 1 = https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L004_R1_001.fastq.gz?id=10172224333
# sample 3 lane 4 read 2 = https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L004_R2_001.fastq.gz?id=10172224335
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
wget -O H508-1_S3_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213598/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213599/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213596/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213597/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213601/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213602/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213600/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213603/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
# Sample 2
wget -O H508-2_S2_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226239/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226240/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226242/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226241/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226243/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226244/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226245/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226246/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
# Sample 3
wget -O H508-3_S1_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224329/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224328/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224331/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224330/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224332/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224334/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224333/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224335/content?access_token=9e1acf007f014868b78778d705ce1e64' &