################################################################################
## Run a pre-step of BSLMM in GEMMA
################################################################################

## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(coda)
library(Rcpp)
library(bigmemory)

library(moiR)
load_all('.')

####************************************************************************####
#### Genomematrix and fam #####

# not necessary # x<-attach.big.matrix("databig/example.desc")
fam<-read.table("databig/example.fam",header=T)
y=fam[,6]
colnames(fam)[6]

####************************************************************************####

run_gemma(plinkbase = "example",
          plinkfolder = paste0("databig/",colnames(fam)[6]),
          out = colnames(fam)[6],type = "bslmm",
          dryrun = F)


####************************************************************************####
#### RUN ALL

for(i in 7:ncol(fam)){

run_gemma(plinkbase = "example",
          plinkfolder = paste0("databig/",colnames(fam)[i]),
          out = colnames(fam)[i],type = "bslmm",
          dryrun = F)


}

