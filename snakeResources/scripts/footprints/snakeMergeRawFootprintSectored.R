
##
cat("Loading libraries", "\n")
suppressMessages(library(GenomicRanges))

##
cat("Setting snakemake vars", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
sitesPath <- snakemake@input[[3]]
#
outPath <- snakemake@output[[1]]
#
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

## Scan for all available input files
sectorPath <- paste0(dirPath, "/footprints/data/temp")
cat("Scanning for sector footprint files at path:", sectorPath, "\n")
sectorFiles <- Sys.glob(file.path(sectorPath, paste0(sampleName, "-REP", sampleRep, ".", geneName, "*rawFootprintData.Rdata")))
cat("Found sector footprint files:", sectorFiles, "\n")

## Load the files and store in a list
numFiles <- length(sectorFiles)
cat("Loading", numFiles, "sectored files", "\n")

## 
footprintDataList <- list()
##
for (a in 1:numFiles){
  ##
  load(sectorFiles[a])
  ##
  footprintDataList[[a]] <- footprintData
} # end for (a in 1:numFiles)

## Determine the number of motifs
numMotifs <- length(footprintDataList[[1]])
cat("Found", numMotifs, "motifs", "\n")

####
cat("Determining number of motifs...", "\n")
load(sitespath)
num_motifs <- length(bindingSites)
cat("Found ", num_motifs, " motifs", "\n")

for (x in 1:num_motifs){
  
  signalpath <- paste0(dirpath, "footprints/data/merged/", sample, ".", gene, ".", "motif", x, ".merged.Rdata")
  cat("Output path for signal object: ", signalpath, "\n")
  
  if (file.exists(signalpath) == TRUE){
    
    cat("Merged file already exists, skipping...", "\n")
    next
    
  } else {
    
    cat("Merged file not found, processing...", "\n")
    
    ## Load each chromosome
    ## Note that, because some chromosomes may have been skipped due to finding to errors
    ## Will need to check that the file exists before loading it and run an error catching loop
    cat("Loading data by chromosome...", "\n")
    chr_names <- paste0("chr", c(1:22, "X", "Y"))
    found_chr <- c()
    
    #
    for (b in chr_names){
      
      cat("Checking for file for", b, "\n")
      
      com <- paste0(b, "_in <- gsub('", b, ".done.bychr.txt', paste0('motif', x, '.", b ,".Rdata'), ", b, "_input)")
      eval(parse(text = com))
      
      com <- paste0(b, "_in <- gsub('operations', 'data/bychr', ", b, "_in)")
      eval(parse(text = com))
      
      com <- paste0("curfile <- '", b, "_in'")
      eval(parse(text = com))
      
      ## check if the file for the current chr was output or not
      if (file.exists(get(curfile)) == FALSE){
        
        cat("No file found for", b, "skipping", "\n")
        next
        
      } else {
        
        #
        cat("Found file for", b, "loading...", "\n")
        found_chr <- c(found_chr, b)
        com <- paste0("load(", curfile, ")")
        eval(parse(text = com))
        com <- paste0("sigs_", b, " <- sigs")
        eval(parse(text = com))
        
      } # if (file.exists(get(curfile)) == FALSE)
    } # end for (b in chr_names)
    
    cat("Chromosome files found: ", found_chr, "\n")
    ## Perform the merge
    cat("Merging data by chromosome...", "\n")
    merged_signal <- list()
    merge_names_plus <- paste0("sigs_", found_chr, "[['signal']][['+']]")
    merge_names_minus <- paste0("sigs_", found_chr, "[['signal']][['-']]")
    
    mplus <- paste(merge_names_plus[1:length(merge_names_plus)], collapse = ",")
    mminus <- paste(merge_names_minus[1:length(merge_names_minus)], collapse = ",")
    
    nplus <- gsub((paste0(merge_names_plus[length(merge_names_plus)], ",")), paste0("sigs_", found_chr[(length(merge_names_plus))], "[['signal']][['+']])"), mplus)
    nminus <- gsub((paste0(merge_names_minus[length(merge_names_minus)], ",")), paste0("sigs_", found_chr[(length(merge_names_minus))], "[['signal']][['+']])"), mminus)
    
    com <- paste0("merged_signal$'+' <- rbind(", nplus, ")")
    eval(parse(text = com))
    com <- paste0("merged_signal$'-' <- rbind(", nminus, ")")
    eval(parse(text = com))
    
    merge_names_sites <- paste0("sigs_", found_chr, "[['bindingSites']]")
    msites <- paste(merge_names_sites[1:length(merge_names_sites)], collapse = ",")
    nsites <- gsub((paste0(merge_names_sites[length(merge_names_sites)], ",")), paste0("sigs_", found_chr[(length(merge_names_sites))], "[['signal']][['+']])"), msites)
    
    com <- paste0("merged_sites <- c(", nsites, ")")
    eval(parse(text = com))
    
    #
    merged_signals <- list()
    merged_signals$"signal" <- merged_signal
    merged_signals$"bindingSites" <- merged_sites
    
    ##
    save(merged_signals, file = signalpath)
    
    
  }
} # end for (x in 1:num_motifs)

cat("Finished merging!", "\n")
file.create(output)








