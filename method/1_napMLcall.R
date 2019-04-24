################################################################################
## Run napML inference over simulated traits
################################################################################
message("Loading packages required for NAP")
## Load packages
library(devtools)
library(dplyr)
library(bigmemory)
library(Rcpp)
library(moiR)
load_all('.')

####************************************************************************####
## Arguments

library('argparser')
cat("\n")
p <- arg_parser("Run Non-Additive Polygenic ML model:")
p <- add_argument(p, "--p", help="phenotype file name prefix of _.param.txt")
p <- add_argument(p, "--l", help="number of SNPs to analyze", default=500)
p <- add_argument(p, "--f", help="fraction of known effect SNPs", default=1)
p <- add_argument(p, "--n", help="fraction of individuals to analyze", default=0.9)
p <- add_argument(p, "--e", help="global epistasis power", default=1)
p <- add_argument(p, "--m", help="fitness mode: 1 for additive, 2 for multiplicative", default=1)
p <- add_argument(p, "--o", help="output folder", default="output")
p <- add_argument(p, "--d", help="dry run?", default=FALSE)
p <- add_argument(p, "--x", help="force run and override files?", default=FALSE)
argv<-parse_args(p)

if(argv$d){
  print(argv)
  stop("User specified this is a dry run")
}

finalbase=paste0(argv$o,"/",argv$p,".results.",
            ifelse(argv$m==1,"a.","m."),
            paste0("e",argv$e),
            ifelse(argv$f==1,paste0("all",argv$l),paste0("frac",argv$f,"-",argv$l))
          )
finalfile<paste0(finalbase,".tsv")
finallog<paste0(finalbase,".log")
finalrda<paste0(finalbase,".rda")

if(all(file.exists(paste0(argv$o,"/",c(finalfile,finallog,finalrda)) ))){
  stop("Result files already exist, want to override? use --x TRUE")
}

####************************************************************************####
#### Data #####
message("Loading data and defining starting conditions")
# Read genome
G<-attach.big.matrix("databig/example.desc")

# Read phenotyes
fam<-read.table("databig/simexample.fam",header=T)
pheno<-colnames(fam)[6]
pheno<-argv$p
y<-fam[,pheno]

# Read run BSLMM
gammas<-.read_gemma(folder=argv$o,name=pheno, what = "bslmm")
bslmm<-rep(0,ncol(G))
bslmm[gammas$pos]<-gammas$effect

# Define training set of individuals
h=sample(1:nrow(G),size = round(argv$n *nrow(G)),replace = F) %>% sort
htest=(1:nrow(G))[-h]
ytrain<-y[h]
ytest<-y[-h]

# Define SNPs to analyse
## the first 500 are effect SNPs, the rest 9500 are neutral
totsnps=argv$l
if(argv$f < 1){ # trick to start -500
  frac<-round((1-argv$f)*500)
  m=(frac:(frac+totsnps))
}else{
  m=(1:totsnps)
}

# Propose starting point of selection coefficients based on BSLMM
s<-rep(0,length(m))
s<-bslmm[m]
s[abs(s)>1]<- 0
s[abs(s)< -1]<- 0
s_range<-range(s) %>% diff


####************************************************************************####
#### Optimization ####
message("Starting optimiztion")
start_time <- Sys.time()
r<-napML(y = ytrain,
        h = h,
        m = m,
        A = G,
        mod=argv$m,
        e=argv$e,
        s=s,
        slow=s-s_range,
        shigh=s+s_range)
end_time <- Sys.time()
end_time - start_time

#### Get inferences
sinf<-r$par[1:length(s)]
winf=wCBM(A = G@address,
                      s=sinf,
                      mycols=m, myrows= htest,
                      mode=argv$m)
wgwa=wCBM(A = G@address,
                      s=bslmm,
                      mycols=1:ncol(G), myrows= htest,
                      mode=1)
realsvar<- ifelse(grepl("svar0.01", argv$p),0.01,0.1)
sreal=ssimC(m,svar=realsvar)

save(file = finalrda,list = list(r,winf,sfin))

####************************************************************************####
# Accuracies and write to file
message("Writing results")
resu<-rbind(
  c("method"="nap_s",rbind(accuracies(sreal,sinf)) ),
  c("method"="gwa_s",rbind(accuracies(sreal,bslmm[m])) ),
  c("method"="nap_i",rbind(accuracies(ytest,winf)) ),
  c("method"="gwa_s",rbind(accuracies(ytest, wgwa)) )
)
write.tsv(file=finalfile,resu  )
sink(file=finallog)
  cat(r$value)
  cat("\n")
  r$par[-c(1:length(s))]
sink()



# ####
# plotresults<-plot_grid(labels=c("ML","BSLMM"),ncol=1,
#         indplot(y,wCBM(A = G@address,s=r$par[1:totsnps],
#                         mycols=1:totsnps,myrows=1:nrow(G), mode=1))$pind  ,
#
#         indplot(y,wCBM(A = G@address,s=gammas$effect,
#                         mycols=1:totsnps,myrows=1:nrow(G), mode=1))$pind
#         )
# save_plot(paste0("figs/",pheno,"_indplot.pdf"),
#           plotresults,base_height = 2*5,base_width = 5)
#
# plotresults2<-plot_grid(labels=c("ML","BSLMM"), ncol=1,
#                         scorplot(ssimC(totsnps,svar=0.01),
#                                  r$par[1:totsnps])$psel         ,
#                         scorplot(ssimC(totsnps,svar=0.01),
#                                  bslmm[1:totsnps])$psel
#               )
# save_plot(paste0("figs/",pheno,"_splot.pdf"),
#           plotresults2,base_height = 2*5,base_width = 5)
#
