################################################################################
## Run napML inference over simulated traits
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

G<-attach.big.matrix("databig/example.desc")
fam<-read.table("databig/simexample.fam",header=T)
pheno<-colnames(fam)[6]
pheno<-"b0.01_a0.01_p0_svar0.01_epi0.9_mod1_h20.96"
y<-fam[,pheno]
h=(1:nrow(G))

hyp<-.read_gemma(folder="output",name=pheno, what = "heritability")
h2=hyp[1]
gammas<-.read_gemma(folder="output",name=pheno, what = "bslmm")
bv<-.read_gemma(folder="output",name=pheno, what = "bv") %>% fn

pheno

####************************************************************************####
#### Try once more optimization ####
totsnps=500
m=(1:totsnps)
s<-rep(0,ncol(G))
s[gammas$pos]<-gammas$effect
s<-s[1:totsnps]

start_time <- Sys.time()
r<-napML(y = y,
        h = h,
        m = m,
        A = G,
        s=s)
end_time <- Sys.time()
end_time - start_time


cor(r$par[1:500],s)
cor(r$par[1:500],ssimC(100,svar=0.01)[1:500])


####
plotresults<-plot_grid(labels=c("ML","BSLMM"),ncol=1,
        indplot(y,wCBM(A = G@address,s=r$par[1:totsnps],
                        mycols=1:totsnps,myrows=1:nrow(G), mode=1))$pind  ,

        indplot(y,wCBM(A = G@address,s=gammas$effect,
                        mycols=1:totsnps,myrows=1:nrow(G), mode=1))$pind
        )
save_plot(paste0("figs/",pheno,"_indplot.pdf"),
          plotresults,base_height = 2*5,base_width = 5)

plotresults2<-plot_grid(labels=c("ML","BSLMM"), ncol=1,
                        scorplot(ssimC(100,svar=0.01)[1:totsnps],
                                 r$par[1:totsnps])$psel         ,
                        scorplot(ssimC(100,svar=0.01)[1:totsnps],
                                 gammas$effect[gammas$pos][1:totsnps])$psel
              )
save_plot(paste0("figs/",pheno,"_splot.pdf"),
          plotresults2,base_height = 2*5,base_width = 5)


