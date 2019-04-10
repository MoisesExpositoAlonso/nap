## Load packages
library(devtools)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)

library(coda)
library(Rcpp)
library(MCMCglmm)

library(moiR)
load_all('.')
# library(nap)
library(bigmemory)


################################################################################
#### SIMULATE GENOME MATRIX ####

n=500
m=10

XBMsimulate<-function(backingpath="databig",filename="example",n,m,force=F){
  # Create a filebacked big matrix in R, and fill it
  # with 2s (homozygous for minor) at a MAF simulated from runif 0-0.5
  # To read matrix, attach.big.matrix("descriptionfile.desc")
  if(force){
    system(paste("rm",paste0(backingpath,"/",filename,".bk")))
  }
  x <- filebacked.big.matrix(n, m, type='int', init= 0,
                             backingpath = backingpath,
                             backingfile=paste0(filename,".bk"),
                             descriptorfile=paste0(filename,".desc")
                            )
  x
  BMsimulate(x@address)
  return(x)
}


start_time <- Sys.time()
x<-XBMsimulate(n=1e3,m=1e4,force=T)
end_time <- Sys.time()
end_time - start_time

################################################################################
#### SIMULATE ARRAY OF GENOTYPES ####

as<-c(0.01,0.5)
bs<-c(0.01,0.5)
p<-c(0,0.5)
mu=1
svars=c(0.1,0.5)
ss=0
epi<-c()
FITs<-
FITmode=3

replicates=1

  maf=mafsim(m)

X=A <- Xsim(n,m,maf)

s= ssim(m,svar)
Ey=wsim(X,s,mode=FITmode,epi = epi,mu = mu)
y=sampleW(Ey,a,b,p,rep = replicates)
h=sort(rep.int(1:n,replicates))
title <- ggdraw() + draw_label(paste("h2= ",format(var(Ey) / var(y),digits=2)),fontface = 'bold')
plot_grid(title,
          plot_grid(qplot(Ey),qplot(y)),
          ncol = 1, rel_heights = c(0.1, 1)
         )

# x<-napMCMCwrap()

d<-napMCMCcompare(100)
d
