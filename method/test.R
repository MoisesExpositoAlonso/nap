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

library(moiR)
load_all('.')
library(bigmemory)

####************************************************************************####
#### Genomematrix and fam #####

x<-attach.big.matrix("databig/example.desc")
fam<-read.table("databig/example.fam",header=T)
head(fam[,1:6])
fam[is.na(fam)]<-0

sinit<- BMs(x@address,y,mycols=1:500,myrows=1:nrow(x))
sinit[abs(sinit)> sd(sinit,na.rm = T) *3]<-0
sinit[is.na(sinit)]<-0

y=fam[,6]
length(y)

# Analyze only effect snps
nap(y = y,
    h = 1:nrow(x),
    A = x@address,
    m = 1:500, # Only effect SNPs
    n = 1:nrow(x),
    s=sinit,
    iterations=10000,
    verbose = F
    ) ->
  res


accuracies(y,wCBM(A = x@address,
                s=sinit,
                # s=rnorm(500,0,0.01),
                mycols=1:500,myrows=1:nrow(x), mode=1))


plot(sinit,
     ssimC(500,svar = 0.05)[1:500])
plot(rnorm(500,0,0.01),
     ssimC(500,svar = 0.05)[1:500])


LIKELIHOOD(y,(1:nrow(x)) -1,
           wCBM(A = x@address,
                s=sinit,
                # s=rnorm(500,0,0.01),
                mycols=1:500,myrows=1:nrow(x), mode=1),
            b=0.2,a=0.2,p=0,mu=1,epi=1,
            verbose=T
           )

LIKELIHOOD(y,(1:nrow(x)) -1,
           # wCBM(A = x@address,s=ssimC(500,0.05),mycols=1:500,myrows=1:nrow(x), mode=1),
           # wCBM(A = x@address,s=BMridge(x[],y),mycols=1:500,myrows=1:nrow(x), mode=1),
           # wCBM(A = x@address,s=rnorm(500,mean = 0,0.01),mycols=1:500,myrows=1:nrow(x), mode=1),
           # wCBM(A = x@address,s=exp(rnorm(500,mean = 0,0.01))-1,mycols=1:500,myrows=1:nrow(x), mode=1),
           # wCBM(A = x@address,s=BMs(x@address,y,mycols=1:500,myrows=1:nrow(x))* 5.793e-01 ,
           #      mycols=1:500,myrows=1:nrow(x), mode=1),
            b=0.1,a=0.1,p=0,mu=1,epi=1,
            verbose=T
           )

LIKELIHOOD(y,(1:nrow(x)) -1,
           # wCBM(A = x@address,s=ssimC(500,0.05),mycols=1:500,myrows=1:nrow(x), mode=1),
           # wCBM(A = x@address,s=BMridge(x[],y),mycols=1:500,myrows=1:nrow(x), mode=1),
           # wCBM(A = x@address,s=rnorm(500,mean = 0,0.01),mycols=1:500,myrows=1:nrow(x), mode=1),
           # wCBM(A = x@address,s=exp(rnorm(500,mean = 0,0.01))-1,mycols=1:500,myrows=1:nrow(x), mode=1),
           # wCBM(A = x@address,s=BMs(x@address,y,mycols=1:500,myrows=1:nrow(x))* 5.793e-01 ,
           #      mycols=1:500,myrows=1:nrow(x), mode=1),
            b=0.1,a=0.1,p=0,mu=1,epi=1,
            verbose=T
           )

plot(BMs(x@address,y,mycols=1:500,myrows=1:nrow(x)),
     ssimC(500,svar = 0.05)[1:500])

lm(ssimC(500,svar = 0.05)[1:500] ~
     BMs(x@address,y,mycols=1:500,myrows=1:nrow(x))
     ) %>% summary


hist(BMs(x@address,y,mycols=1:500,myrows=1:nrow(x)))
(BMs(x@address,y,mycols=1:500,myrows=1:nrow(x))/4) %>% hist

hshist(ssimC(500,svar = 0.05)[1:500])

hist(fam[,6])

# library(microbenchmark)
# microbenchmark(
# wC(x[],
#       s=rnorm(10000,0,0.1),
#       mode=1,epi=1
#     ),
# wCBM2(A = x@address,
#       s=rnorm(10000,0,0.1),
#       mycols = 5:10, # Only effect SNPs
#       myrows = 5:10,
#       mode=1,epi=1
#     ),
# wCBM(A = x@address,
#       s=rnorm(500,0,0.1),
#       mycols = (1:500), # Only effect SNPs
#       myrows = (1:nrow(x)),
#       mode=1,epi=1
#     )
# )

