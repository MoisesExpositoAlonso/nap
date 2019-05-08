################################################################################
## Run napML inference over simulated traits
################################################################################
message("Loading packages required for NAP")
## Load packages
library(latex2exp);
library(optimx)
library(BB)
library(cowplot);
library(devtools);
library(dplyr);
library(bigmemory);
library(Rcpp);
devtools::load_all('napML')
sourceCpp("napML/src/ML.cpp")
source("napML/R/napML.R")

setwd("~/nap")

####************************************************************************####
## Arguments

library('argparser')
cat("\n")
p <- arg_parser("Run Non-Additive Polygenic ML model:")
p <- add_argument(p, "--p", help="path to .param.txt")
p <- add_argument(p, "--f", help="path to .fam")
p <- add_argument(p, "--g", help="path to .desc")
p <- add_argument(p, "--l", help="number of SNPs to analyze", default=500)
p <- add_argument(p, "--q", help="fraction of known effect SNPs", default=1)
p <- add_argument(p, "--n", help="fraction of individuals to analyze", default=0.9)
p <- add_argument(p, "--e", help="global epistasis power", default=1)
p <- add_argument(p, "--m", help="fitness mode: 1 for additive, 2 for multiplicative", default=1)
p <- add_argument(p, "--d", help="dry run?", default=FALSE)
argv<-parse_args(p)

if(argv$d){
  print(argv)
  stop("User specified this is a dry run")
}

if(is.na(argv$p)){
  argv$p<-"output/b0.5_a0.5_p0.2_svar0.1_epi1.2_mod1_h20.74.param.txt"
  argv$f<-"databig/b0.5_a0.5_p0.2_svar0.1_epi1.2_mod1_h20.74/example.fam"
  argv$g<-"databig/example.desc"
}
#
# argv$l<- 10

bname<-basename(argv$p)
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

#if(all(file.exists(c(finalfile,finallog,finalrda)) )){
#  stop("Result files already exist, want to override?")
#}
#
####************************************************************************####
#### Data #####
message("Loading data and defining starting conditions")
# Read genome
G<-attach.big.matrix(argv$g)

# Read phenotyes
y<-read.table(argv$f)[,6]

# Read run BSLMM
gammas<-.read_gemma(folder=dname,name=bname, what = "bslmm")
bslmm<-rep(0,ncol(G))
bslmm[gammas$pos]<-gammas$effect


# Define training set of individuals
hall<-1:nrow(G)
h=sample(1:nrow(G),size = round(argv$n *nrow(G)),replace = F) %>% sort
htest=(1:nrow(G))[-h]
ytrain<-y[h]
ytest<-y[-h]

# Define SNPs to analyse
if(argv$q < 1){ # trick to start -500
  frac<-round((1-argv$q)*500)
  m=(frac:(frac+argv$l))
}else{
  m=(1:argv$l)
}

# Propose starting point of selection coefficients based on BSLMM
s<-rep(0,length(m))
s<-bslmm[m]
# s[abs(s)>1]<- 0
# s[abs(s)< -1]<- 0
s_range<-sd(s)*3
s_range<-diff(range(s))
slow<-s-s_range
shigh<-s+s_range



####************************************************************************####
#### Inference ####
redo<-lm(y ~G[,m])
redo$coefficients
s<-redo$coefficients[-1]
s[is.na(s)]<-0
mycor<-cor(G[,m])
table(mycor==1)

w=wCBM(A = G@address,
          s=s,
          # s=rep(0,length(m)),
          mycols =m,
          myrows =hall,
          mode=1,
          epi=1)
plot(w,y)
likelihoodR(y = y,
            w = w,
            b=0.5,
            a=0.5,
            p=0) * (-1)
lik.nap(y=y,
        h_=hall,
        m_=m,
        A=G,
        par=c(s,0.5,0.5,0),
        mod=1,
        e=1
        )

parstart<-list(
                "s"=s,
                "b"=1,
                "a"=1,
                "p"=0
                ) %>% unlist
parlow<- list(
              # "s"=slow,
              "s"=rep(-0.5, length(s)),
              "b"=0.01,
              "a"=0.01,
              "p"=0
              ) %>% unlist
parhigh<- list(
              # "s"=shigh,
              "s"=rep(0.5, length(s)),
              "b"=1,
              "a"=1,
              "p"=1
              )%>% unlist
r<-spg(
        fn=lik.nap,
        y=y,
        h_=hall,
        m_=m,
        A=G,
        mod=argv$m,
        e=argv$e,
        # control=list(maxit=100),
        control=list(maxit=20),
        par = parstart ,
        lower = parlow,
        upper=parhigh
)
r<-optimx(fn=lik.nap,
          y=y,
          h_=hall,
          m_=m,
          A=G,
          mod=argv$m,
          e=argv$e,
          par = parstart ,
          lower = parlow,
          upper=parhigh,
          control=list(
                      trace=1,
                      maxit=20
                          ),
          method = "L-BFGS-B"
          )
r



####************************************************************************####
#### Get inferences
#### Get inferences
realsvar<- ifelse(grepl("svar0.01", argv$p),0.01,0.1)
sreal=ssimC(m,svar=realsvar)
sinf<-moiR::fn(r[1,1:length(s)])
sinf<-moiR::fn(r$par[1:length(s)])
winf=wCBM(A = G@address,
                      s=sinf,
                      mycols=m,
                      myrows= hall,
                      mode=argv$m,
                      epi=argv$e
          )
wgwa=wCBM(A = G@address,
                      # s=gammas$effect,
                      # mycols=mapping,
                      s=bslmm,
                      mycols=1:ncol(G),
                      # s=s,
                      # mycols=m,
                      # s=bslmm[m],
                      # s=alphas[m],
                      # s=betas[m],
                      myrows= hall,
                      mode=1,
                      epi=1
          )

saveRDS(file = finalrda,object = list(r,winf,sinf,wgwa,bslmm))

#### plot

plotresults<-plot_grid(labels=c("ML","BSLMM"),ncol=1,
        indplot(y,winf)$pind  ,
        indplot(y,wgwa)$pind
        )
plotresults2<-plot_grid(labels=c("ML","BSLMM"), ncol=1,
                        scorplot(sreal,sinf)$psel         ,
                        scorplot(sreal,s)$psel
              )

save_plot(paste0("figs/",bname,"_indplot","_l",argv$l,"_m",argv$m,"_e",argv$e ,".pdf"),
          plotresults,base_height = 2*5,base_width = 5)

save_plot(paste0("figs/",bname,"_splot","_l",argv$l,"_m",argv$m,"_e",argv$e ,".pdf"),
          plotresults2,base_height = 2*5,base_width = 5)



# #### Get inferences
# sinf<-r$par[1:length(s)]
# winf=wCBM(A = G@address,
#                       s=sinf,
#                       mycols=m, myrows= hall,
#                       mode=argv$m,
#                       epi=argv$e
#           )
# winftrain<-wCBM(A = G@address,
#                       s=sinf,
#                       mycols=m, myrows= h,
#                       mode=argv$m,
#                       epi=argv$e
#           )
# wgwatrain=wCBM(A = G@address,
#                       s=bslmm,
#                       mycols=1:ncol(G), myrows= h,
#                       mode=1)
#
# wgwa=wCBM(A = G@address,
#                       s=bslmm,
#                       mycols=1:ncol(G), myrows= hall,
#                       mode=1)
# realsvar<- ifelse(grepl("svar0.01", argv$p),0.01,0.1)
# sreal=ssimC(m,svar=realsvar)
#
# # plot(sinf,sreal)
# # plot(winftrain,ytest)
# # plot(winf,ytest)
# # plot(wgwa,ytest)
# # plot(s,sreal)
# plot(wgwatrain,ytrain)


####************************************************************************####
# Accuracies and write to file
# message("Writing results")
# resu<-rbind(
#   cbind("method"="nap_s",rbind(accuracies(sreal,sinf)) ),
#   # cbind("method"="nap_i",rbind(accuracies(ytest,winf)) ),
#   cbind("method"="nap_i",rbind(accuracies(y,winf)) ),
#   cbind("method"="gwa_s",rbind(accuracies(sreal,bslmm[m])) ),
#   cbind("method"="gwa_i",rbind(accuracies(y, wgwa)) )
# )
# resu
# write.tsv(file=finalfile,resu  )
# sink(file=finallog)
#   cat(r$value);cat("\n")
#   cat(end_time-start_time);cat("\n")
#   r$par[-c(1:length(s))]
# sink()



# ####************************************************************************####
# ####************************************************************************####


####************************************************************************####
# #### Optimization ####
# start_time <- Sys.time()
# print(paste("Starting optimiztion on", start_time))
# w=wCBM(A = G@address,
#           s=s,
#           mycols =m,
#           myrows =hall,
#           mode=1,
#           epi=1)
# # plot(w,y)
# likelihoodR(y = y,
#             w = w,
#             b=0.5,
#             a=0.5,
#             p=0)
#
# parstart<-list(
#                 "s"=s,
#                 "b"=0.5,
#                 "a"=0.5,
#                 "p"=0
#                 ) %>% unlist
# parlow<- list(
#               "s"=slow,
#               "b"=0.01,
#               "a"=0.01,
#               "p"=0
#               ) %>% unlist
# parhigh<- list(
#               "s"=shigh,
#               "b"=1,
#               "a"=1,
#               "p"=1
#               )%>% unlist
# myscale<-c( rep(median(abs(s)), length(s) ), rep(0.1,3))
# r<-optimx(fn=lik.nap,
#       y=y,
#       h_=hall,
#       m_=m,
#       A=G,
#       mod=argv$m,
#       e=argv$e,
#       par = parstart ,
#       lower = parlow,
#       upper=parhigh,
#       control=list(
#                   trace=1,
#                   fnscale=myscale,
#                   maxit=20
#                       ),
#       method = "spg"
#       # method ="BFGS"
#       # method = "L-BFGS-B"
#       # method = "CG"
#       )
# r
# end_time<-Sys.time()
# message("Ending optimiztion ")
# print(end_time - start_time)
