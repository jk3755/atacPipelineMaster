####################################################################################################
##### METHYLATION ANALYSIS WITH BICYCLE PIPELINE ###################################################
####################################################################################################
# Bicycle pipeline manual can be found here:
# http://www.sing-group.org/bicycle/manual.html#faq
#
# hg38 reference sequence in .fa (fasta) format can be found here:
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
#
###################################################################################################
#
# bicycle will want to create the PROJECT directory itself, so don't do it beforehand
# -p specifies the project directory
# -r specifies the reference sequence directory (.fa format)
# -f specifies the reads data directory (.fastq format)
#
# CREATE THE PROJECT FIRST
# bash /home/ubuntu2/atac/programs/bicycle/cmd/bicycle create-project -p /home/ubuntu2/atac/h508/wt02a/rrbs/project -r /home/ubuntu2/atac/h508/wt02a/rrbs/reference -f /home/ubuntu2/atac/h508/wt02a/rrbs/reads
#
# CREATE THE BISULFATION REFERENCE
# bash /home/ubuntu2/atac/programs/bicycle/cmd/bicycle reference-bisulfitation -p /home/ubuntu2/atac/h508/wt02a/rrbs/project
#
# CREATE THE REFERENCE INDEX
# This may take several hours, even with 20 threads
# -v specifies the bowtie version to use
# -t specifies the number of bowtie2 threads to use
# bash /home/ubuntu2/atac/programs/bicycle/cmd/bicycle reference-index -p /home/ubuntu2/atac/h508/wt02a/rrbs/project -v 2 -t 20
#
# ALIGN THE READS
# This will take some time, depending on how big the input .fastq is
# -t specifies alignment threads
# -v specifies bowtie version
## ERRORS MAY BE THROWN IF YOU USE THE WRONG Q SCORES, TRY SPECIFYING WITH -Q@ OPTION ##
# bash /home/ubuntu2/atac/programs/bicycle/cmd/bicycle align -p /home/ubuntu2/atac/h508/wt02a/rrbs/project -t 10 -v 2
# bash /home/ubuntu2/atac/programs/bicycle/cmd/bicycle align -p /home/ubuntu2/atac/h508/wt02a/rrbs/project -t 10 -v 2 -q2 phred33
#
# ANALYZE METHYLATION
# Will sort, convert to bam, and index first, methylation analysis takes a few minutes
# -n specifies the number of threads
# -r specifies ignore non-correctly bisulfite-converted reads
# -a specifies remove ambiguous (aligned to two strands) reads
# -o only use uniquely mapping reads
# -f sets FDR threshold (default = 0.01)
# -c remove clonal reads
# bash /home/ubuntu2/atac/programs/bicycle/cmd/bicycle analyze-methylation -p /home/ubuntu2/atac/h508/wt02a/rrbs/project -n 20 -f 0.05


#rule generate_motifData:
#    output:
#        "sites/motifData.Rdata"
#    script:
#        "scripts/scanPWM/generateMotifData.R"
