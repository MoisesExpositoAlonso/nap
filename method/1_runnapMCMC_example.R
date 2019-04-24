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

G<-attach.big.matrix("databig/example.desc")
fam<-read.table("databig/simexample.fam",header=T)
pheno<-colnames(fam)[6]
pheno<-"b0.01_a0.01_p0_svar0.01_epi0.9_mod1_h20.96"
y<-fam[,pheno]
h=(1:nrow(G))-1

hyp<-.read_gemma(folder="output",name=pheno, what = "heritability")
h2=hyp[1]
gammas<-.read_gemma(folder="output",name=pheno, what = "bslmm")
bv<-.read_gemma(folder="output",name=pheno, what = "bv") %>% fn


pheno

####************************************************************************####
#### Try once more optimization ####
# m=1000
s<-rep(0,ncol(G))
s[gammas$pos]<-gammas$effect

start_time <- Sys.time()
r<-napML(y = y,
        h = h,
        A = G,
        s=gammas$effect)
end_time <- Sys.time()
end_time - start_time



plotresults<-plot_grid(labels=c("ML","BSLMM"),ncol=1,
        indplot(y,wCBM(A = G@address,s=r$par[1:ncol(G)],
                        mycols=1:ncol(G),myrows=1:nrow(G), mode=1))$pind  ,

        indplot(y,wCBM(A = G@address,s=gammas$effect,
                        mycols=1:ncol(G),myrows=1:nrow(G), mode=1))$pind
        )
save_plot(paste0("figs/",pheno,"_indplot.pdf"),
          plotresults,base_height = 2*5,base_width = 5)

plotresults2<-plot_grid(labels=c("ML","BSLMM"), ncol=1,
                        scorplot(ssimC(100,svar=0.01),r$par[1:ncol(G)])$psel         ,
                        scorplot(ssimC(100,svar=0.01)[gammas$pos],
                                 gammas$effect)$psel
              )
save_plot("exampleworkedApril18_selection.pdf",plotresults2,base_height = 2*5,base_width = 5)


#### TRY MCMC ####

fitmod=2
napMCMC(y = y,
        h = 1:n,
        A = x@address,
        mycols = 1:m,
        myrows = 1:n,
        s=s,
        iterations=1e4,
        verbose = F ,
        bw=0.1#,
        # updateratio=0.9, # set to 1 for only S optim
        # bw=0.0001#,
        # b=0.1,
        # # bmin=0.1,bmax=0.1,
        # a=0.1,
        # p=p,
        # svar=0.05,
        # epi=1,epimin=0.9,epimax=1.1,
        # FITmode = fitmod,
        # Smode = 2, PRImode = 2
        ) ->
      res

posteriorplot(res$posterior)
#
# cor(data.frame(res2$shat,res$shat,s))
# # accuracies(y,wCBM(A = x@address,s=gammas$effect[1:500]
         # mycols=1:500,myrows=1:nrow(x), mode=1))


# ####************************************************************************####
# #### RUN EXAMPLE  -- GOOD
# m=500 # to limit only to effect SNPs.
# n=1500
# sstart=rep(0,ncol(x))
# sstart[gammas$pos]<-gammas$effect
# s=sstart[1:m]
#
# a= (1-h2) * var(y)
# p=length(y[y==0])/length(y)
#
# fitmod=2
# napMCMC(y = y,
#         h = 1:n,
#         A = x@address,
#         mycols = 1:m,
#         myrows = 1:n,
#         s=s,
#         iterations=1e5,
#         verbose = F ,
#         epi=1.1
#         # epimin=0.9,epimax=1.1,
#         # updateratio=0.9, # set to 1 for only S optim
#         # bw=0.0001#,
#         # b=0.1,
#         # # bmin=0.1,bmax=0.1,
#         # a=0.1,
#         # p=p,
#         # svar=0.05,
#         # FITmode = fitmod,
#         # Smode = 2, PRImode = 2
#         ) ->
#       res
# res$accuracy
# accuracies(y,wCBM(A = x@address,s=gammas$effect[1:500],
#                 mycols=1:500,myrows=1:nrow(x), mode=1))
# ####************************************************************************####
# #
# # cbind(parnames(),res$par)
# #
# # indplot(y,wCBM(A = x@address,s=gammas$effect[1:500],
# #                 mycols=1:500,myrows=1:nrow(x), mode=1))
# # indplot(y,res$w)
# #
# # scorplot(ssimC(100,svar=0.01)[1:500],res2$shat)
# # scorplot(ssimC(100,svar=0.01)[1:500],res$shat)
# # scorplot(ssimC(100,svar=0.01)[gammas$pos][1:500],
# #          gammas$effect[1:500])
#
#
# ####************************************************************************####
# #### erase at some point ####
# # LIKELIHOOD(y,(1:nrow(x)) -1,
# #            wCBM(A = x@address,s=s,
# #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# #             b=0.5,a=0.5,p=0.5,mu=1,epi=1,
# #             verbose=T,
# #            printall=T
# #             )
# # LIKELIHOOD(y,(1:nrow(x)) -1,
# #            wCBM(A = x@address,s=res$shat,
# #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# #             b=0.1,a=0.1,p=0,mu=1,epi=1,
# #             verbose=T
# #             )
# # LIKELIHOOD(y,(1:nrow(x)) -1,
# #            wCBM(A = x@address,s=ssim(500,svar=0.01)[1:500],
# #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# #             b=0.1,a=0.1,p=0,mu=1,epi=1,
# #             verbose=T
# #             )
#
#
#
# # res$accuracy
# # posteriorplot(res$likelihood)
# # cbind(parnames(),res$par)
# #
# # plot( main="true",
# #      y,
# #      wCBM(A = x@address,
# #           s=ssimC(1,0.05)[1:500],
# #                 mycols=1:500,myrows=1:nrow(x), mode=fitmod))
# # plot(main="mGWA",
# #      y,
# #      wCBM(A = x@address,
# #           s=BMs(x@address,y,mycols=1:500,myrows=1:nrow(x)),
# #                 mycols=1:500,myrows=1:nrow(x), mode=1))
# # plot(main="nap",
# #      y,
# #      res$w)
# #
# # scorplot(
# #     res$shat,
# #     ssimC(m,svar=0.05)[1:m]
# #     )
# # scorplot(
# #     BMs(x@address,y,mycols=1:m,myrows=1:nrow(x)),
# #      ssimC(1,0.05)[1:m]
# #     )
#
#
# # gwam<-lm(y ~ x[])
# # gwa %>% summary()
#
# # LIKELIHOOD(y,(1:nrow(x)) -1,
# #            wCBM(A = x@address,s=BMs(x@address,y,mycols=1:500,myrows=1:nrow(x)),
# #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# #             b=0.5,a=0.5,p=0.5,mu=1,epi=1,
# #             verbose=T,
# #            printall=T
# #             )
# # LIKELIHOOD(y,(1:nrow(x)) -1,
# #            wCBM(A = x@address,s=ssim(500,svar=0.05)[1:500],
# #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# #             b=0.5,a=0.5,p=0.5,mu=1,epi=1,
# #             verbose=T
# #             )
# # LIKELIHOOD(y,(1:nrow(x)) -1,
# #            wCBM(A = x@address,s=res$shat,
# #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# #             b=0.1,a=0.1,p=0,mu=1,epi=1,
# #             verbose=T
# #             )
