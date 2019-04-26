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
  argv$p<-"output/b0.01_a0.01_p0_svar0.1_epi1.1_mod2_h21.param.txt"
  argv$f<-"databig/b0.01_a0.01_p0_svar0.1_epi1.1_mod2_h21/example.fam"
  argv$g<-"databig/example.desc"
}


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

if(all(file.exists(c(finalfile,finallog,finalrda)) )){
  stop("Result files already exist, want to override?")
}

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
timeittook<-end_time - start_time


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
  cat(r$value);cat("\n")
  cat(timeittook);cat("\n")
  r$par[-c(1:length(s))]
sink()
