########################################################################################################################
#### Raw FP Analysis Groups 1-40, Normal Processing ####################################################################
########################################################################################################################

## Path for the size-ordered gene name list
namePath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\snakeResources\\scripts\\footprints\\names\\bindingSitesSizeOrdered.txt"
## Load the size-ordered gene names
orderedNames <- readLines(namePath)
numGenes <- length(orderedNames)

## Output path for the text file
outPath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\snakeResources\\scripts\\footprints\\names\\rawFPnames.txt"

## A vector of strings to hold the file line text for output
writeText <- c()

##
group1Names <- c()

####
while (a <= 40){
  
  strings <- c()
  a <- 1 # count for genes
  b <- 1 # string index (group)
  
  #### Groups 1-40
  while (b <= 40){
    
    ##
    c <- a
    d <- a+1
    e <- a+2
    f <- a+3
    g <- a+4
    h <- a+5
    i <- a+6
    j <- a+7
    k <- a+8
    l <- a+9
    m <- a+10
    n <- a+11
    o <- a+12
    p <- a+13
    q <- a+14
    r <- a+15
    s <- a+16
    t <- a+17
    u <- a+18
    v <- a+19
    
    ##
    tmp1 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[c], ".rawFPanalysis.bamcopy1.done", "', ")
    tmp2 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[d], ".rawFPanalysis.bamcopy2.done", "', ")
    tmp3 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[e], ".rawFPanalysis.bamcopy3.done", "', ")
    tmp4 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[f], ".rawFPanalysis.bamcopy4.done", "', ")
    tmp5 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[g], ".rawFPanalysis.bamcopy5.done", "', ")
    tmp6 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[h], ".rawFPanalysis.bamcopy6.done", "', ")
    tmp7 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[i], ".rawFPanalysis.bamcopy7.done", "', ")
    tmp8 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[j], ".rawFPanalysis.bamcopy8.done", "', ")
    tmp9 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[k], ".rawFPanalysis.bamcopy9.done", "', ")
    tmp10 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[l], ".rawFPanalysis.bamcopy10.done", "', ")
    tmp11 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[m], ".rawFPanalysis.bamcopy11.done", "', ")
    tmp12 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[n], ".rawFPanalysis.bamcopy12.done", "', ")
    tmp13 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[o], ".rawFPanalysis.bamcopy13.done", "', ")
    tmp14 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[p], ".rawFPanalysis.bamcopy14.done", "', ")
    tmp15 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[q], ".rawFPanalysis.bamcopy15.done", "', ")
    tmp16 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[r], ".rawFPanalysis.bamcopy16.done", "', ")
    tmp17 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[s], ".rawFPanalysis.bamcopy17.done", "', ")
    tmp18 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[t], ".rawFPanalysis.bamcopy18.done", "', ")
    tmp19 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[u], ".rawFPanalysis.bamcopy19.done", "', ")
    tmp20 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[v], ".rawFPanalysis.bamcopy20.done", "'")
    
    ##
    strings[b] <- paste0(
      "rule rawFPanalysis_group",
      b,
      ":\n",
      "\tinput:\n\t\t",
      tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
      tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t",
      tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t",
      tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20, "\n\t",
      "output:\n\t\t",
      "'{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group", b, ".done'\n",
      "\tshell:\n\t\t",
      "'touch {output}'")
    
    ##
    a <- a+20
    b <- b+1
  
  
  
} # end while (a <= 40)

  
  
##
b <- 1
while (b <= 200){
  tmp <- c()
  for (c in 1:20){
    group1Names[b] <- paste0("\t\t'{path}operations/footprints/raw/{sample}.", orderedNames[b], ".rawFPanalysis.bamcopy", c, ".done", "',")
    b <- (b + 1)
  } # end for (c in 1:20)
} # end while (b <= 200)

##
writeText[1] <- "rule rawFPanalysis_group1:"
writeText[2] <- "\tinput:"
writeText[3:202] <- group1Names
writeText[203] <- "\toutput:"
writeText[204] <- "\t\t'{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group1.done'"
writeText[205] <- "\tshell:"
writeText[206] <- "\t\t'touch {output}'"




#### Write the file
write.table(
  writeText,
  file = outPath,
  quote = FALSE,
  sep = ",",
  eol = "\n",
  row.names = FALSE,
  col.names = FALSE)


########################################################################################################################
#### Raw FP Analysis Groups 1-4, Normal Processing #####################################################################
########################################################################################################################
#### The first 800 genes are small enough to run with the basi script
#### They can be organized into 4 groups of 200

## Path for the size-ordered gene name list
namePath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\snakeResources\\scripts\\footprints\\names\\bindingSitesSizeOrdered.txt"
## Load the size-ordered gene names
orderedNames <- readLines(namePath)
numGenes <- length(orderedNames)

## Output path for the text file
outPath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\snakeResources\\scripts\\footprints\\names\\rawFPnames.txt"

## A vector of strings to hold the file line text for output
writeText <- c()

#### Group 1, Genes 1-200
group1Names <- c()
b <- 1
while (b <= 200){
  tmp <- c()
  for (c in 1:20){
    group1Names[b] <- paste0("\t\t'{path}operations/footprints/raw/{sample}.", orderedNames[b], ".rawFPanalysis.bamcopy", c, ".done", "',")
    b <- (b + 1)
  } # end for (c in 1:20)
} # end while (b <= 200)

##
writeText[1] <- "rule rawFPanalysis_group1:"
writeText[2] <- "\tinput:"
writeText[3:202] <- group1Names
writeText[203] <- "\toutput:"
writeText[204] <- "\t\t'{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group1.done'"
writeText[205] <- "\tshell:"
writeText[206] <- "\t\t'touch {output}'"

#### Group 2, Genes 201-400
group2Names <- c()
b <- 201
while (b <= 400){
  tmp <- c()
  for (c in 1:20){
    group2Names[b-200] <- paste0("\t\t'{path}operations/footprints/raw/{sample}.", orderedNames[b], ".rawFPanalysis.bamcopy", c, ".done", "',")
    b <- (b + 1)
  } # end for (c in 1:20)
} # end while (b <= 200)

##
writeText[207] <- "rule rawFPanalysis_group2:"
writeText[208] <- "\tinput:"
writeText[209:408] <- group2Names
writeText[409] <- "\toutput:"
writeText[410] <- "\t\t'{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group2.done'"
writeText[411] <- "\tshell:"
writeText[412] <- "\t\t'touch {output}'"

#### Group 3, Genes 401-600
group3Names <- c()
b <- 401
while (b <= 600){
  tmp <- c()
  for (c in 1:20){
    group3Names[b-400] <- paste0("\t\t'{path}operations/footprints/raw/{sample}.", orderedNames[b], ".rawFPanalysis.bamcopy", c, ".done", "',")
    b <- (b + 1)
  } # end for (c in 1:20)
} # end while (b <= 200)

##
writeText[413] <- "rule rawFPanalysis_group3:"
writeText[414] <- "\tinput:"
writeText[415:614] <- group3Names
writeText[615] <- "\toutput:"
writeText[616] <- "\t\t'{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group3.done'"
writeText[617] <- "\tshell:"
writeText[618] <- "\t\t'touch {output}'"

#### Group 4, Genes 601-800
group4Names <- c()
b <- 601
while (b <= 800){
  tmp <- c()
  for (c in 1:20){
    group4Names[b-600] <- paste0("\t\t'{path}operations/footprints/raw/{sample}.", orderedNames[b], ".rawFPanalysis.bamcopy", c, ".done", "',")
    b <- (b + 1)
  } # end for (c in 1:20)
} # end while (b <= 200)

##
writeText[619] <- "rule rawFPanalysis_group4:"
writeText[620] <- "\tinput:"
writeText[621:820] <- group4Names
writeText[821] <- "\toutput:"
writeText[822] <- "\t\t'{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group4.done'"
writeText[823] <- "\tshell:"
writeText[824] <- "\t\t'touch {output}'"
  

#### Write the file
write.table(
            writeText,
            file = outPath,
            quote = FALSE,
            sep = ",",
            eol = "\n",
            row.names = FALSE,
            col.names = FALSE)


########################################################################################################################
#### Raw FP Analysis Groups 1-4, Normal Processing #####################################################################
########################################################################################################################