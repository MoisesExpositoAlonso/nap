################################################################################
## Run napML inference over simulated traits
################################################################################
message("Loading packages required for NAP")
## Load packages
loadpackages<-function(){
library(devtools);
library(dplyr);
library(bigmemory);
library(Rcpp);
library(moiR);
}
suppressWarnings( loadpackages() )
devtools::load_all('napML')
# sourceCpp("napML/src/ML.cpp")
#library('nap')

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


####************************************************************************####
#### Optimization ####
message("Starting optimiztion")
start_time <- Sys.time()
Sys.time()
r<-napML(
        y = y,
        h = hall,
        # y = ytrain,
        # h = htrain,
        m = m,
        A = G,
        mod=argv$m,
        e=argv$e,
        s=s,
        slow=s-s_range,
        shigh=s+s_range)
Sys.time()
end_time <- Sys.time()
end_time - start_time
r


#### Get inferences
sinf<-r$par[1:length(s)]
winf=wCBM(A = G@address,
                      s=sinf,
                      mycols=m, myrows= htest,
                      mode=argv$m,
                      epi=argv$e
          )
winftrain<-wCBM(A = G@address,
                      s=sinf,
                      mycols=m, myrows= h,
                      mode=argv$m,
                      epi=argv$e
          )
wgwatrain=wCBM(A = G@address,
                      s=bslmm,
                      mycols=1:ncol(G), myrows= h,
                      mode=1)

wgwa=wCBM(A = G@address,
                      s=bslmm,
                      mycols=1:ncol(G), myrows= htest,
                      mode=1)
realsvar<- ifelse(grepl("svar0.01", argv$p),0.01,0.1)
sreal=ssimC(m,svar=realsvar)

plot(sinf,sreal)
plot(winftrain,ytrain)
plot(winf,ytest)
plot(wgwa,ytest)
plot(wgwatrain,ytrain)

saveRDS(file = finalrda,object = list(r,winf,sinf))

####************************************************************************####
# Accuracies and write to file
message("Writing results")
resu<-rbind(
  cbind("method"="nap_s",rbind(accuracies(sreal,sinf)) ),
  cbind("method"="nap_i",rbind(accuracies(ytest,winf)) ),
  cbind("method"="gwa_s",rbind(accuracies(sreal,bslmm[m])) ),
  cbind("method"="gwa_i",rbind(accuracies(ytest, wgwa)) )
)
resu
write.tsv(file=finalfile,resu  )
sink(file=finallog)
  cat(r$value);cat("\n")
  cat(timeittook);cat("\n")
  r$par[-c(1:length(s))]
sink()



# ####************************************************************************####
# ####************************************************************************####
#### DEBUG ####
par = unlist(list(
                  # "s"=rep(0,length(s)),
                  "s"=s,
                  # "s"=sreal,
                  "b"=1,
                  "a"=1,
                  "p"=0.5#,
                  ))
w=wCBM(A = G@address,
          s=par[1:length(m)],
          mycols =m,
          myrows =h,
          mode=argv$m,
          epi=argv$e)
hist(w)
LIKELIHOOD(y = ytrain,
            w = w,
            b=par[length(m)+1],
            a=par[length(m)+2],
            p=par[length(m)+3],
            mu=1,
            epi=argv$e
            )



####************************************************************************####
####************************************************************************####
