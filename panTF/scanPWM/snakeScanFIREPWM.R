
#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("Biostrings", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates = TRUE)

#### Library loading
cat("Loading libraries...", "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))
suppressMessages(library(MotifDb))


##
load("C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\panTF\\scanPWM\\coad_mr_fire_motifs.RData")


#### PFM ####
## ADNP
ADNP <- list()
ADNP$"motif1" <- COAD_FIRE_PFM[["ADNP_motif1_fire_PFM"]]
## TAF4
TAF4 <- list()
TAF4$"motif1" <- COAD_FIRE_PFM[["TAF4_motif1_fire_PFM"]]
TAF4$"motif2" <- COAD_FIRE_PFM[["TAF4_motif2_fire_PFM"]]
## ZMYND8
ZMYND8 <- list()
ZMYND8$"motif1" <- COAD_FIRE_PFM[["ZMYND8_motif1_fire_PFM"]]
## ZNF696
ZNF696 <- list()
ZNF696$"motif1" <- COAD_FIRE_PFM[["ZNF696_motif1_fire_PFM"]]
ZNF696$"motif2" <- COAD_FIRE_PFM[["ZNF696_motif2_fire_PFM"]]
ZNF696$"motif3" <- COAD_FIRE_PFM[["ZNF696_motif3_fire_PFM"]]

#### PWM ####
ADNP <- list()
ADNP$"motif1" <- COAD_FIRE_PWM[["ADNP_motif1_fire_PWM"]]
## TAF4
TAF4 <- list()
TAF4$"motif1" <- COAD_FIRE_PWM[["TAF4_motif1_fire_PWM"]]
TAF4$"motif2" <- COAD_FIRE_PWM[["TAF4_motif2_fire_PWM"]]
## ZMYND8
ZMYND8 <- list()
ZMYND8$"motif1" <- COAD_FIRE_PWM[["ZMYND8_motif1_fire_PWM"]]
## ZNF696
ZNF696 <- list()
ZNF696$"motif1" <- COAD_FIRE_PWM[["ZNF696_motif1_fire_PWM"]]
ZNF696$"motif2" <- COAD_FIRE_PWM[["ZNF696_motif2_fire_PWM"]]
ZNF696$"motif3" <- COAD_FIRE_PWM[["ZNF696_motif3_fire_PWM"]]


#### ADNP
outPath <- "C:\\Users\\jsk33\\Desktop\\ADNP.bindingSites.Rdata"

## Set parameters
genome <- Hsapiens
score <- "99%"
bindingSites <- list()

## Scan the genome for matches to each unique motif
PWM <- ADNP[["motif1"]]
rownames(PWM) <- DNA_BASES
tempSites1 <- list()
sites <- matchPWM(PWM, genome, min.score=score)
tempSites$PWM <- PWM
tempSites$sites <- sites

bindingSites$"motif1" <- tempSites
save(bindingSites, file = outPath)


#### TAF4
outPath <- "C:\\Users\\jsk33\\Desktop\\TAF4.bindingSites.Rdata"

## Set parameters
genome <- Hsapiens
score <- "99%"
bindingSites <- list()

## Scan the genome for matches to each unique motif
PWM <- TAF4[["motif1"]]
rownames(PWM) <- DNA_BASES
tempSites <- list()
sites <- matchPWM(PWM, genome, min.score=score)
tempSites$PWM <- PWM
tempSites$sites <- sites
bindingSites$"motif1" <- tempSites

## Scan the genome for matches to each unique motif
PWM <- TAF4[["motif2"]]
rownames(PWM) <- DNA_BASES
tempSites <- list()
sites <- matchPWM(PWM, genome, min.score=score)
tempSites$PWM <- PWM
tempSites$sites <- sites
bindingSites$"motif2" <- tempSites

save(bindingSites, file = outPath)


#### ZMYND8
outPath <- "C:\\Users\\jsk33\\Desktop\\ZMYND8.bindingSites.Rdata"

## Set parameters
genome <- Hsapiens
score <- "99%"
bindingSites <- list()

## Scan the genome for matches to each unique motif
PWM <- ZMYND8[["motif1"]]
rownames(PWM) <- DNA_BASES
tempSites <- list()
sites <- matchPWM(PWM, genome, min.score=score)
tempSites$PWM <- PWM
tempSites$sites <- sites
bindingSites$"motif1" <- tempSites

save(bindingSites, file = outPath)

#### ZNF696
outPath <- "C:\\Users\\jsk33\\Desktop\\ZNF696.bindingSites.Rdata"

## Set parameters
genome <- Hsapiens
score <- "99%"
bindingSites <- list()

## Scan the genome for matches to each unique motif
PWM <- ZNF696[["motif1"]]
rownames(PWM) <- DNA_BASES
tempSites <- list()
sites <- matchPWM(PWM, genome, min.score=score)
tempSites$PWM <- PWM
tempSites$sites <- sites
bindingSites$"motif1" <- tempSites

## Scan the genome for matches to each unique motif
PWM <- ZNF696[["motif2"]]
rownames(PWM) <- DNA_BASES
tempSites <- list()
sites <- matchPWM(PWM, genome, min.score=score)
tempSites$PWM <- PWM
tempSites$sites <- sites
bindingSites$"motif2" <- tempSites

## Scan the genome for matches to each unique motif
PWM <- ZNF696[["motif3"]]
rownames(PWM) <- DNA_BASES
tempSites <- list()
sites <- matchPWM(PWM, genome, min.score=score)
tempSites$PWM <- PWM
tempSites$sites <- sites
bindingSites$"motif3" <- tempSites

save(bindingSites, file = outPath)