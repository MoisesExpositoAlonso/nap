library(Rcpp)
library(bigmemory)
library(data.table)
library(devtools)
library(dplyr)
library(moiR)
load_all('.')

set.seed(1)

################################################################################
#### SIMULATE GENOME MATRIX ####

x<-attach.big.matrix("../databig/genome.desc") # this was read from 515g.012
totinds=nrow(x)
# Choose SNPs only chr1
totsnps=1e4
map<-fread("../databig/515g.map")
# do not run # snpsort<-sample(which(map[,1]==1),totsnps)
# do not run #  write.table(quote = F,row.names = F,col.names = F,file="databig/snp_sort.txt",snpsort)
snpsort<-moiR::fn(read.table("../databig/snp_sort.txt",header = F))


################################################################################
#### SIMULATE ARRAY OF GENOTYPES ####
# same parameter ranges as random LD matrix

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


#### Run simulations of phenotypes and create FAM
inds<-1:nrow(x) # all individuals

d<-matrix(ncol=nrow(mysim),nrow=length(inds) ) %>% data.frame
h2s<-c()
for( i in 1:nrow(mysim)){
  l=mysim[i,]
  s = c(ssimC(1:totsnps,fn(l["svar"])))
  w=wCBM(x@address,s,
         mycols=snpsort,
         myrows=1:totinds,
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
                           c( "ld",
                              paste0(colnames(mysim), (mysim[i,])),
                              paste0("h2",h2)
                            )
                        )
}
head(d)
d[is.na(d)]<-0

#### write simulations array (and add h2)
mysim$h2<-h2s
write.csv(file = "../databig/ldsimulationgrid.tsv",mysim)


#### write general fam
structure<-data.frame(FID=1:totinds,IID=1:totinds,PAT=0,MAT=0,SEX=0)
dfam<-cbind(structure,d)
dim(dfam)
head(dfam[,1:6])
write.table(row.names = F,col.names = T, quote = F,
            file="../databig/ldsimexample.fam",
            dfam
            )
# d<-read.table("../databig/ldsimexample.fam",header=T)[,-c(1:5)]

#### write fams per folder
for( i in 1:ncol(d)){
  message(colnames(d)[i])
  system(paste0("mkdir ../databig/",colnames(d)[i]))
  write.table(row.names = F,col.names = F, quote = F,
            file=paste0("../databig/",colnames(d)[i],"/ldexample.fam"),
            cbind(structure,data.frame(PHENO=d[,i]))
            )
  system(paste0("ln ../databig/515g.bed ../databig/",colnames(d)[i],"/ldexample.bed"))
  system(paste0("ln ../databig/515g.bim ../databig/",colnames(d)[i],"/ldexample.bim"))
}
####### IMPORTANT 515genome matrix!
