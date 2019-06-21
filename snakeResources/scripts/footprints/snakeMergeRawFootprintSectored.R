
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

#### Merge the sectored data
cat("Merging sectored data", "\n")

##
footPrintData <- list()

for (b in 1:numMotifs){
  
  ##
  tempList <- list()
  
  ## PWM
  com <- paste0("tempList$PWM <- footprintDataList[[1]][['motif", b, "']][['PWM']]")
  eval(parse(text = com))
  
  ## Motif Width
  com <- paste0("tempList$motifWidth <- footprintDataList[[1]][['motif", b, "']][['motifWidth']]")
  eval(parse(text = com))
  
  ## allSites
  com <- paste0("mergedSites <- c(footprintDataList[[1]][['motif", b, "']][['allSites']], footprintDataList[[2]][['motif", b, "']][['allSites']],
                footprintDataList[[3]][['motif", b, "']][['allSites']], footprintDataList[[4]][['motif", b, "']][['allSites']],
                footprintDataList[[5]][['motif", b, "']][['allSites']], footprintDataList[[6]][['motif", b, "']][['allSites']],
                footprintDataList[[7]][['motif", b, "']][['allSites']], footprintDataList[[8]][['motif", b, "']][['allSites']],
                footprintDataList[[9]][['motif", b, "']][['allSites']], footprintDataList[[10]][['motif", b, "']][['allSites']],
                footprintDataList[[11]][['motif", b, "']][['allSites']], footprintDataList[[12]][['motif", b, "']][['allSites']],
                footprintDataList[[13]][['motif", b, "']][['allSites']], footprintDataList[[14]][['motif", b, "']][['allSites']],
                footprintDataList[[15]][['motif", b, "']][['allSites']], footprintDataList[[16]][['motif", b, "']][['allSites']],
                footprintDataList[[17]][['motif", b, "']][['allSites']], footprintDataList[[18]][['motif", b, "']][['allSites']],
                footprintDataList[[19]][['motif", b, "']][['allSites']], footprintDataList[[20]][['motif", b, "']][['allSites']])")
  eval(parse(text = com))
  ##
  com <- paste0("tempList$allSites <- mergedSites")
  eval(parse(text = com))
  
  ####
  com <- paste0("footprintData$motif", b, " <- tempList")
  eval(parse(text = com))
  
} # end for (b in 1:numMotifs)


###############
x <- c(footprintDataList[[1]][["motif1"]][["allSites"]], footprintDataList[[2]][["motif1"]][["allSites"]], footprintDataList[[3]][["motif1"]][["allSites"]],
       footprintDataList[[4]][["motif1"]][["allSites"]], footprintDataList[[5]][["motif1"]][["allSites"]], footprintDataList[[6]][["motif1"]][["allSites"]],
       footprintDataList[[7]][["motif1"]][["allSites"]], footprintDataList[[8]][["motif1"]][["allSites"]], footprintDataList[[9]][["motif1"]][["allSites"]],
       footprintDataList[[10]][["motif1"]][["allSites"]], footprintDataList[[11]][["motif1"]][["allSites"]], footprintDataList[[12]][["motif1"]][["allSites"]],
       footprintDataList[[13]][["motif1"]][["allSites"]], footprintDataList[[14]][["motif1"]][["allSites"]], footprintDataList[[15]][["motif1"]][["allSites"]],
       footprintDataList[[16]][["motif1"]][["allSites"]], footprintDataList[[17]][["motif1"]][["allSites"]], footprintDataList[[18]][["motif1"]][["allSites"]],
       footprintDataList[[19]][["motif1"]][["allSites"]], footprintDataList[[20]][["motif1"]][["allSites"]])


g1 <- footprintDataList[[1]][["motif1"]][["allSites"]]
g2 <- footprintDataList[[2]][["motif1"]][["allSites"]]

g3 <- c(g1, g2)




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








