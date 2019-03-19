# rule threesample_plotcorrspearman:
#     # parameters:
#     # -b input bam files
#     # -o output file name
#     # -bs set the bin size used for comparison, default is 10000 bp
#     # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
#     # -p set the number of computing processors to use
#     # -v verbose mode
#     input:
#         a="{sample1}/{wt1}{num1}/preprocessing/12all/{s1}.bam",
#         b="{sample1}/{wt1}{num1}/preprocessing/12all/{s2}.bam",
#         c="{sample1}/{wt1}{num1}/preprocessing/12all/{s3}.bam"
#     output:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.corrTest"
#     shell:
#         "multiBamSummary bins -b {input.a} {input.b} {input.c} -o {output} -bs 10000 -p 20 -v"

# rule threesample_makecorrheatmap:
#     input:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.corrTest"
#     output:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.heatmap.svg"
#     shell:
#         "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"

# rule ninesample_plotcorrspearman:
#     # parameters:
#     # -b input bam files
#     # -o output file name
#     # -bs set the bin size used for comparison, default is 10000 bp
#     # -r to reduce computation time, a specific region of genome can be set, format: chr1:10000:20000
#     # -p set the number of computing processors to use
#     # -v verbose mode
#     # For processing nine sample, restrict analysis to chr1, or computation will take forever
#     input:
#         a="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
#         b="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
#         c="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam",
#         d="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
#         e="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
#         f="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam",
#         g="{sample1}/{wt1}{num1}/preprocessing/10unique/{s1}-REP1.u.bam",
#         h="{sample1}/{wt1}{num1}/preprocessing/10unique/{s2}-REP2.u.bam",
#         i="{sample1}/{wt1}{num1}/preprocessing/10unique/{s3}-REP3.u.bam"
#     output:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.REPS.spearman.corrTest"
#     shell:
#         "multiBamSummary bins -b {input.a} {input.b} {input.c} {input.d} {input.e} {input.f} {input.g} {input.h} {input.i} -o {output} -bs 10000 -p 20 -v -r chr1"

# rule ninesample_makecorrheatmap:
#     input:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.REPS.spearman.corrTest"
#     output:
#         "xsample_analysis/correlation/{sample1}-{wt1}-{num1}.{sample2}-{wt2}-{num2}.{sample3}-{wt3}-{num3}.spearman.heatmap.svg"
#     shell:
#         "plotCorrelation -in {input} -c spearman -p heatmap -o {output} --plotNumbers"