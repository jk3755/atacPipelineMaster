



## Build a null model simulation similar to the one you use to identify actual footprints
## Start with a total signal =100, end at toal signal = 1000
## Do 1000 iterations for each value, distribute signal randomly, calculate mean at the motif locus
## Will need to do this for different values of motif length, start from 6-14 bp
## This might be computationally expensive, do save these models for future use



## Once the model is built, convert the bam files into a vector for each chromosome with signal at each bp position
## This will significantly speed up the processing



## Now, at each interval of 200 bp + each simulated motif size, calculate the total signal, signal in the motif, and mean
## Compare these values to the corresponding null model. Given a predetermined p-value cutoff, decide if a footprint is present or not
## Procees base pair by base pair across the entire set of accessible (peak) regions in the dataset


## Once putative footprints have been identified, motif enrichment analysis can be performed. You can use known motifs as a positive control


## Compare the results to ARACNe networks to see if you can recover or validate specific interactions