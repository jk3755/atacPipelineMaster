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