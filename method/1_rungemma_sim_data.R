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
#### fam random marix#####
fam<-read.table("../databig/simexample.fam",header=T)

mysub<-sample(10,7:ncol(fam))

for(i in mysub){
  plinkfolder<-paste0("../databig/",colnames(fam)[i])
  if(!file.exists(paste0(plinkfolder,".param.txt"))){
    run_gemma(plinkbase = "example",
              plinkfolder = plinkfolder,
              out = colnames(fam)[i],
              type = "bslmm",
              dryrun = T)
  }
}

####************************************************************************####
#### fam LD random marix#####
ldfam<-read.table("../databig/ldsimexample.fam",header=T)

mysub<-sample(10,7:ncol(ldfam))

for(i in mysub){
  plinkfolder<-paste0("../databig/",colnames(ldfam)[i])
  if(!file.exists(paste0(plinkfolder,".param.txt"))){
    run_gemma(plinkbase = "example",
              plinkfolder = plinkfolder,
              out = colnames(ldfam)[i],
              type = "bslmm",
              dryrun = T)
  }
}

####************************************************************************####
# Run predictions?

pred_gemma()
