library(bigmemory)
library(data.table)
library(devtools)
library(dplyr)
library(moiR)
load_all('.')


set.seed(1)

################################################################################
#### SIMULATE GENOME MATRIX ####
totsnps=1e4
totinds=1.5e3
# x<-XBMsimulate(n=1.5e3,m=1e4,force=T)
# MAPwrite(x,path="../databig/example")
# RAW012write(x,path="../databig/example")
# PEDwrite(x,path="../databig/example")
# BEDmake(path="../databig/example")

#### Reatach
x<-attach.big.matrix("../databig/example.desc")

################################################################################
#### SIMULATE ARRAY OF GENOTYPES ####

as<-c(0.01,0.5)
bs<-c(0.01,0.5)
ps<-c(0,0.2)
mu=1
svars=c(0.01,0.1)
ss=0
epi<-c(0.8,1,1.2)
FITs<-c(1,2)
replicates=1 # assuming one replicate

#### grid of simulation
mysim<-expand.grid(bs,as,ps,svars,epi,FITs)
colnames(mysim)<-c("b","a","p","svar","epi","mod")
head(mysim)
dim(mysim)

#### effect SNPs
## generate some distributions, that I can compare later
# s_svar01<- c(exp(rnorm(500,0,0.01))-1 , rep(0,1e4 -500))
# hist(s_svar01[s_svar01!=0])
# ssaveC(s_svar01,"databig/s_svar01.txt") # all simulations will have same S but scaled to a Svar


#### Run simulations of phenotypes and create FAM
d<-matrix(ncol=nrow(mysim),nrow=length(inds) ) %>% data.frame
i=1
h2s<-c()
for( i in 1:nrow(mysim)){
  l=mysim[i,]
  s = ssimC(1:(totsnps),fn(l["svar"]));
  # all(wC(x[],s,1,1,1) == wCBM(x@address,s,1:totsnps-1,1:totinds-1,1,1,1))
  w=wCBM(x@address,
         s,
         mycols = 1:totsnps,myrows= 1:totinds,
         mode=fn(l["mod"]),
         epi=fn(l["epi"])
         )
  y=sampleWC(w,
             b=fn(l["b"]),
             a=fn(l["a"]),
             p=fn(l["p"]),
             rep=1)
  d[,i]<-y
  h2<-format(var(w,na.rm = T) / var(y,na.rm = T),digits=2)
  h2s[i]<-h2
  colnames(d)[i]<- paste0( collapse="_",
                           c(paste0(  colnames(mysim), (mysim[i,])),
                           paste0("h2",h2)
                          )
                        )
}

head(d)
d[is.na(d)]<-0

#### write simulations array (and add h2)
mysim$h2<-h2s
write.csv(file = "../databig/simulationgrid.tsv",mysim)

#### write general fam
structure<-data.frame(FID=1:totinds,IID=1:totinds,PAT=0,MAT=0,SEX=0)
dfam<-cbind(structure,d)
dim(dfam)
head(dfam[,1:6])
write.table(row.names = F,col.names = T, quote = F,
            file="../databig/simexample.fam",
            dfam
            )


#### write fams per folder
for( i in 1:ncol(d)){
  message(colnames(d)[i])
  system(paste0("mkdir ../databig/",colnames(d)[i]))
  write.table(row.names = F,col.names = F, quote = F,
            file=paste0("../databig/",colnames(d)[i],"/example.fam"),
            cbind(structure,data.frame(PHENO=d[,i]))
            )
  system(paste0("ln ../databig/example.bed ../databig/",colnames(d)[i],"/example.bed"))
  system(paste0("ln ../databig/example.bim ../databig/",colnames(d)[i],"/example.bim"))
}

