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


hyp<-.read_gemma(folder="output",name=pheno, what = "heritability")
betas<-.read_gemma(folder="output",name=pheno, what = "lm")
gammas<-.read_gemma(folder="output",name=pheno, what = "bslmm")
bv<-.read_gemma(folder="output",name=pheno, what = "bv") %>% fn
length(bv)

mycols<- which(gammas$gamma >= sort(gammas$gamma,decreasing = T)[100] )
s<-gammas[mycols,"effect"]
hist(s)


mycols2<- which(betas$beta >= sort(betas$beta,decreasing = T)[100] )
s<-betas[mycols2,"effect"]
hist(s)

length(setdiff(mycols,mycols2)) /length(unique(c(mycols,mycols2)))


qplot(ssimC(100,svar=0.01), gammas$effect)


####************************************************************************####
####************************************************************************####

#### RUN EXAMPLE
m=500
n=1500
fitmod=2
napMCMC(y = y[1:n],
    h = 1:n,
    A = x@address,
    mycols = 1:m, # Only effect SNPs
    myrows = 1:n,
    s=rnorm(m,0,0.01),
    iterations=5e4,
    verbose = F ,
    b=0.5, bmin=0.01, bmax=0.01,
    a=0.5, amin=0.01, amax=0.01,
    p=0.2,pmin=0.2,pmax=0.2,
    svar=0.25,svarmin=0.25,svarmax=0.25,
    epi=1.1,
    # epimin=0.8,epimax=1.2,
    FITmode = 1, Smode = 1, PRImode = 1
    # mumin=0.5,mumax=5
    ) ->
  res

res$accuracy
posteriorplot(res$likelihood)
cbind(parnames(),res$par)

plot( main="true",
     y,
     wCBM(A = x@address,
          s=ssimC(1,0.05)[1:500],
                mycols=1:500,myrows=1:nrow(x), mode=fitmod))
plot(main="mGWA",
     y,
     wCBM(A = x@address,
          s=BMs(x@address,y,mycols=1:500,myrows=1:nrow(x)),
                mycols=1:500,myrows=1:nrow(x), mode=1))
plot(main="nap",
     y,
     res$w)

scorplot(
    res$shat,
    ssimC(m,svar=0.05)[1:m]
    )
scorplot(
    BMs(x@address,y,mycols=1:m,myrows=1:nrow(x)),
     ssimC(1,0.05)[1:m]
    )


# gwam<-lm(y ~ x[])
# gwa %>% summary()

LIKELIHOOD(y,(1:nrow(x)) -1,
           wCBM(A = x@address,s=BMs(x@address,y,mycols=1:500,myrows=1:nrow(x)),
                mycols=1:500,myrows=1:nrow(x), mode=1),
            b=0.5,a=0.5,p=0.5,mu=1,epi=1,
            verbose=T,
           printall=T
            )
LIKELIHOOD(y,(1:nrow(x)) -1,
           wCBM(A = x@address,s=ssim(500,svar=0.05)[1:500],
                mycols=1:500,myrows=1:nrow(x), mode=1),
            b=0.5,a=0.5,p=0.5,mu=1,epi=1,
            verbose=T
            )
# LIKELIHOOD(y,(1:nrow(x)) -1,
#            wCBM(A = x@address,s=res$shat,
#                 mycols=1:500,myrows=1:nrow(x), mode=1),
#             b=0.1,a=0.1,p=0,mu=1,epi=1,
#             verbose=T
#             )
