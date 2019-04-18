################################################################################
## Run nap MCMC inference over simulated traits
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

x<-attach.big.matrix("databig/example.desc")
fam<-read.table("databig/simexample.fam",header=T)
pheno<-colnames(fam)[6]
pheno<-"b0.01_a0.01_p0_svar0.01_epi0.9_mod1_h20.96"
y<-fam[,pheno]

hyp<-.read_gemma(folder="output",name=pheno, what = "heritability")
gammas<-.read_gemma(folder="output",name=pheno, what = "bslmm")
bv<-.read_gemma(folder="output",name=pheno, what = "bv") %>% fn
cor(bv,y)
plot(bv,y)


betas<-.read_gemma(folder="output",name=pheno, what = "lm")
mycols<- which(gammas$gamma >= sort(gammas$gamma,decreasing = T)[100] )
s<-gammas[mycols,"effect"]
hist(s)
mycols2<- which(betas$beta >= sort(betas$beta,decreasing = T)[100] )
s<-betas[mycols2,"effect"]
hist(s)

length(setdiff(mycols,mycols2)) /length(unique(c(mycols,mycols2)))
