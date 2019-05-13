################################################################################
## Run nap MCMC inference over simulated traits
################################################################################

## Load packages
library(BB)
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(coda)
library(Rcpp)
library(bigmemory)

library(moiR)

setwd("~/nap")
load_all('napML')

####************************************************************************####
## Arguments

library('argparser')
cat("\n")
p <- arg_parser("Run Non-Additive Polygenic ML model:")
p <- add_argument(p, "--p", help="path to .param.txt")
p <- add_argument(p, "--f", help="path to .fam")
p <- add_argument(p, "--g", help="path to .desc")
p <- add_argument(p, "--l", help="number of SNPs to analyze", default=500)
p <- add_argument(p, "--q", help="fraction of known effect SNPs (if 0, it guesses)", default=1)
p <- add_argument(p, "--n", help="fraction of individuals to analyze", default=0.95)
p <- add_argument(p, "--e", help="global epistasis power", default=1)
p <- add_argument(p, "--m", help="fitness mode: 1 for additive, 2 for multiplicative", default=2)
p <- add_argument(p, "--i", help="maximum optimization iterations", default=100)
p <- add_argument(p, "--d", help="dry run?", default=FALSE)
argv<-parse_args(p)

if(argv$d){
  print(argv)
  stop("User specified this is a dry run")
}

if(is.na(argv$p)){
  # argv$p<-"output/b0.5_a0.5_p0.2_svar0.1_epi1.2_mod1_h20.74.param.txt"
  # argv$f<-"databig/b0.5_a0.5_p0.2_svar0.1_epi1.2_mod1_h20.74/example.fam"
# pheno<-"b0.01_a0.01_p0_svar0.01_epi0.9_mod1_h20.96" # old
# pheno<-"b0.01_a0.01_p0_svar0.01_epi1_mod1_h20.99"
# pheno<-"b0.01_a0.01_p0_svar0.01_epi1_mod2_h21" # from reworked example
  argv$p<-"output/b0.01_a0.01_p0_svar0.01_epi1_mod2_h21.param.txt"
  argv$f<-"databig/b0.01_a0.01_p0_svar0.01_epi1_mod2_h21/example.fam"
  argv$g<-"databig/example.desc"
}


####************************************************************************####
#### Filenames
bname<-pheno<-basename(argv$p)
dname<-gsub(x=argv$p,pattern = bname,replacement = "")
bname<-gsub(x=bname,pattern=".param.txt",replacement="")

finalbase=paste0(dname,"/",bname,".results.",
            ifelse(argv$m==1,"a.","m."),
            paste0("e",argv$e),
            ifelse(argv$q==1,paste0("all",argv$l),paste0("frac",argv$q,"-",argv$l))
          )
finalfile<-paste0(finalbase,".tsv")
finallog<-paste0(finalbase,".log")
finalrda<-paste0(finalbase,".rda")
finaltsv<-paste0(finalbase,".tsv")
finaliplot<-paste0(finalbase,"_indplot.pdf")
finalsplot<-paste0(finalbase,"_splot.pdf")
successfile<-paste0(finalbase,".success")

#if(all(file.exists(c(finalfile,finallog,finalrda)) )){
#  stop("Result files already exist, want to override?")
#}
#
####************************************************************************####
cat("Loading data  \n")
# Read genome
G<-A<-attach.big.matrix(argv$g)

# Read phenotyes
pheno=bname
y<-read.table(argv$f)[,6]

hall<-h<-1:nrow(G)
htrain=sample(1:nrow(G),size = round(argv$n *nrow(G)),replace = F) %>% sort
htest=(1:nrow(G))[-htrain]
ytrain<-y[htrain]
ytest<-y[-htrain]

# y<-fam[,pheno]
# hist(y)
# h<-h_<-1:nrow(G)


####************************************************************************####
#### BSLMM starting point
cat("Loading gwa for starting conditions  \n")
gammas<-.read_gemma(folder=dname,name=bname, what = "bslmm")
bslmm<-rep(0,ncol(G))
bslmm[gammas$pos]<-gammas$effect

# Define SNPs to analyse
if(argv$q < 1){ # trick to start -500
  frac<-round((1-argv$q)*500)
  m=(frac:(frac+argv$l))
}else if(argv$q ==0){ # if known fraction is 0, let the software guess
  m<-which(rank(max(abs(bslmm))-abs(bslmm)) <= argv$l )
}else{
  m=(1:argv$l)
}

# Propose starting point of selection coefficients based on BSLMM
s<-rep(0,length(m))
s<-bslmm[m]


# gammas<-.read_gemma(folder="output",name=pheno, what = "bslmm")
#
# m=500 # to limit only to effect SNPs.
# n=1500
# sstart=rep(0,ncol(x))
# sstart[gammas$pos]<-gammas$effect
# s=sstart[1:m]
# spos<-m_<-1:m
# hist(s)

####************************************************************************####
#### optimization ####
cat("Optimization  \n")
w = wCBM(A = A@address,
          s=s,
          mycols =m,
          myrows = h,
          mode=1,
          epi=1
          )
cor(w,y)
plot(w,y)
r<-napSPG(y = ytrain,
        m_ =m,
        h_ = htrain,
        A = G,
        s=s,
        mod=1,
        epi=1,
        iter=500
        )
res<-parseoptim(r)
sinf<-res$s
####************************************************************************####
#### analyze results ####
cat("Results  \n")
winf=wCBM(A = G@address,
                      s=sinf,
                      mycols=m,
                      myrows= hall,
                      mode=argv$m,
                      epi=argv$e
          )
wtest=wCBM(A = G@address,
                      s=sinf,
                      mycols=m,
                      myrows= htest,
                      mode=argv$e,
                      epi=argv$m
          )
wgwa=wCBM(A = G@address,
                      s=bslmm,
                      mycols=1:ncol(G),
                      myrows= hall,
                      mode=1,
                      epi=1
          )

####************************************************************************####
#### save ####
saveRDS(file = finalrda,object = list(res,winf,wgwa,bslmm))
####************************************************************************####


#### plots ####
plotresults<-plot_grid(labels=c("ML","BSLMM"),ncol=1,
        indplot(y,winf)$pind  ,
        indplot(y,wgwa)$pind
        )
plotresults

realsvar<- ifelse(grepl("svar0.01", argv$p),0.01,0.1)
sreal=ssimC(m,svar=realsvar)
plotresults2<-plot_grid(labels=c("ML","BSLMM"), ncol=1,
                        scorplot(sreal,sinf)$psel         ,
                        scorplot(sreal,s)$psel
              )
plotresults2
save_plot(finaliplot,
          plotresults,base_height = 2*5,base_width = 5)
save_plot(finalsplot,
          plotresults2,base_height = 2*5,base_width = 5)

####************************************************************************####
accuracyresults(
  pheno=bname,
  mod=argv$m,
  epi=argv$e,
  wgwa=wgwa,winf=winf,wtest=wtest,y=y,ytest=ytest,
  finaltsv=finaltsv
)

####************************************************************************####
# sink() # end sinc

sink(successfile)
cat(res$AIC)
sink()



# winf = wCBM(A = A@address,
#             s=res$s,
#             mycols =m,
#             myrows =h,
#             mode=1,
#             epi=1
#             )
# cor(winf,y)
# plot(winf,y)
# plot(wCBM(A = A@address,
#       s=sstart[1:m],
#       mycols =spos,
#       myrows = h,
#       mode=1,
#       epi=1
#       ),y)
#
# cor(d$s,s)
# cor(d$s,ssimC(m_,svar=0.01))
# plot(sstart[1:m],ssimC(m_,svar=0.01))
# plot(d$s,ssimC(m_,svar=0.01))
# plot(d2$s,ssimC(m_,svar=0.01))
#
# cor(r$par[1:500],s)
# cor(r$par[1:500],ssimC(100,svar=0.01)[1:500])
#
# sinf<-fn(unlist(r)[1:500])
# cor(sinf,s)
# cor(sinf,ssimC(100,svar=0.01)[1:500])

# 21.22
#
# r<-readRDS("exampleoptim.rda")
#
#
# plotresults<-plot_grid(labels=c("ML","BSLMM","MCMC"),ncol=1,
#         indplot(y,wCBM(A = G@address,s=r$par[1:500],
#                         mycols=1:500,myrows=1:nrow(G), mode=1))$pind  ,
#
#         indplot(y,wCBM(A = G@address,s=gammas$effect[1:500],
#                         mycols=1:500,myrows=1:nrow(G), mode=1))$pind   ,
#         indplot(y,res$w)$pind
#         )
# save_plot("exampleworkedApril18_individuals.pdf",plotresults,base_height = 3*5,base_width = 5)
#
# plotresults2<-plot_grid(labels=c("ML","BSLMM","MCMC"), ncol=1,
#                         scorplot(ssimC(100,svar=0.01)[1:500],r$par[1:500])$psel         ,
#                         scorplot(ssimC(100,svar=0.01)[gammas$pos][1:500],
#                                  gammas$effect[1:500])$psel                      ,
#                         scorplot(ssimC(100,svar=0.01)[1:500],res$shat)$psel
#               )
# save_plot("exampleworkedApril18_selection.pdf",plotresults2,base_height = 3*5,base_width = 5)
#
#
# #### TRY MCMC ####
#
# fitmod=2
# napMCMC(y = y,
#         h = 1:n,
#         A = x@address,
#         mycols = 1:m,
#         myrows = 1:n,
#         s=s,
#         iterations=1e4,
#         verbose = F ,
#         bw=0.1#,
#         # updateratio=0.9, # set to 1 for only S optim
#         # bw=0.0001#,
#         # b=0.1,
#         # # bmin=0.1,bmax=0.1,
#         # a=0.1,
#         # p=p,
#         # svar=0.05,
#         # epi=1,epimin=0.9,epimax=1.1,
#         # FITmode = fitmod,
#         # Smode = 2, PRImode = 2
#         ) ->
#       res
#
# posteriorplot(res$posterior)
# #
# # cor(data.frame(res2$shat,res$shat,s))
# # # accuracies(y,wCBM(A = x@address,s=gammas$effect[1:500]
#          # mycols=1:500,myrows=1:nrow(x), mode=1))
#
#
# # ####************************************************************************####
# # #### RUN EXAMPLE  -- GOOD
# # m=500 # to limit only to effect SNPs.
# # n=1500
# # sstart=rep(0,ncol(x))
# # sstart[gammas$pos]<-gammas$effect
# # s=sstart[1:m]
# #
# # a= (1-h2) * var(y)
# # p=length(y[y==0])/length(y)
# #
# # fitmod=2
# # napMCMC(y = y,
# #         h = 1:n,
# #         A = x@address,
# #         mycols = 1:m,
# #         myrows = 1:n,
# #         s=s,
# #         iterations=1e5,
# #         verbose = F ,
# #         epi=1.1
# #         # epimin=0.9,epimax=1.1,
# #         # updateratio=0.9, # set to 1 for only S optim
# #         # bw=0.0001#,
# #         # b=0.1,
# #         # # bmin=0.1,bmax=0.1,
# #         # a=0.1,
# #         # p=p,
# #         # svar=0.05,
# #         # FITmode = fitmod,
# #         # Smode = 2, PRImode = 2
# #         ) ->
# #       res
# # res$accuracy
# # accuracies(y,wCBM(A = x@address,s=gammas$effect[1:500],
# #                 mycols=1:500,myrows=1:nrow(x), mode=1))
# # ####************************************************************************####
# # #
# # # cbind(parnames(),res$par)
# # #
# # # indplot(y,wCBM(A = x@address,s=gammas$effect[1:500],
# # #                 mycols=1:500,myrows=1:nrow(x), mode=1))
# # # indplot(y,res$w)
# # #
# # # scorplot(ssimC(100,svar=0.01)[1:500],res2$shat)
# # # scorplot(ssimC(100,svar=0.01)[1:500],res$shat)
# # # scorplot(ssimC(100,svar=0.01)[gammas$pos][1:500],
# # #          gammas$effect[1:500])
# #
# #
# # ####************************************************************************####
# # #### erase at some point ####
# # # LIKELIHOOD(y,(1:nrow(x)) -1,
# # #            wCBM(A = x@address,s=s,
# # #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# # #             b=0.5,a=0.5,p=0.5,mu=1,epi=1,
# # #             verbose=T,
# # #            printall=T
# # #             )
# # # LIKELIHOOD(y,(1:nrow(x)) -1,
# # #            wCBM(A = x@address,s=res$shat,
# # #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# # #             b=0.1,a=0.1,p=0,mu=1,epi=1,
# # #             verbose=T
# # #             )
# # # LIKELIHOOD(y,(1:nrow(x)) -1,
# # #            wCBM(A = x@address,s=ssim(500,svar=0.01)[1:500],
# # #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# # #             b=0.1,a=0.1,p=0,mu=1,epi=1,
# # #             verbose=T
# # #             )
# #
# #
# #
# # # res$accuracy
# # # posteriorplot(res$likelihood)
# # # cbind(parnames(),res$par)
# # #
# # # plot( main="true",
# # #      y,
# # #      wCBM(A = x@address,
# # #           s=ssimC(1,0.05)[1:500],
# # #                 mycols=1:500,myrows=1:nrow(x), mode=fitmod))
# # # plot(main="mGWA",
# # #      y,
# # #      wCBM(A = x@address,
# # #           s=BMs(x@address,y,mycols=1:500,myrows=1:nrow(x)),
# # #                 mycols=1:500,myrows=1:nrow(x), mode=1))
# # # plot(main="nap",
# # #      y,
# # #      res$w)
# # #
# # # scorplot(
# # #     res$shat,
# # #     ssimC(m,svar=0.05)[1:m]
# # #     )
# # # scorplot(
# # #     BMs(x@address,y,mycols=1:m,myrows=1:nrow(x)),
# # #      ssimC(1,0.05)[1:m]
# # #     )
# #
# #
# # # gwam<-lm(y ~ x[])
# # # gwa %>% summary()
# #
# # # LIKELIHOOD(y,(1:nrow(x)) -1,
# # #            wCBM(A = x@address,s=BMs(x@address,y,mycols=1:500,myrows=1:nrow(x)),
# # #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# # #             b=0.5,a=0.5,p=0.5,mu=1,epi=1,
# # #             verbose=T,
# # #            printall=T
# # #             )
# # # LIKELIHOOD(y,(1:nrow(x)) -1,
# # #            wCBM(A = x@address,s=ssim(500,svar=0.05)[1:500],
# # #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# # #             b=0.5,a=0.5,p=0.5,mu=1,epi=1,
# # #             verbose=T
# # #             )
# # # LIKELIHOOD(y,(1:nrow(x)) -1,
# # #            wCBM(A = x@address,s=res$shat,
# # #                 mycols=1:500,myrows=1:nrow(x), mode=1),
# # #             b=0.1,a=0.1,p=0,mu=1,epi=1,
# # #             verbose=T
# # #             )
