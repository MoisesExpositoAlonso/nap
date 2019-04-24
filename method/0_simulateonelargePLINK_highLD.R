library(bigmemory)
library(data.table)
library(devtools)
library(dplyr)
library(moiR)
load_all('.')

set.seed(1)

################################################################################
#### SIMULATE GENOME MATRIX ####

x<-attach.big.matrix("databig/genomes.desc")
totinds=nrow(x)
# Choose SNPs only chr1
totsnps=1e4
map<-fread("databig/515g.map")
# do not run # snpsort<-sample(which(map[,1]==1),totsnps)
# do not run #  write.table(quote = F,row.names = F,col.names = F,file="databig/snp_sort.txt",snpsort)
snpsort<-moiR::fn(read.table("databig/snp_sort.txt",header = F))


################################################################################
#### SIMULATE ARRAY OF GENOTYPES ####
# same parameter ranges as random LD matrix
mysim<-read.csv(file = "databig/simulationgrid.tsv")[,-1] # remove row names

#### Run simulations of phenotypes and create FAM
inds<-1:nrow(x) # all individuals

# d<-matrix(ncol=nrow(mysim),nrow=length(inds) ) %>% data.frame
# for( i in 1:nrow(mysim)){
#   l=mysim[i,]
#   s = c(ssimC(1:10000,fn(l["svar"])))
#   w=wCBM(x@address,s,
#          mycols=snpsort,
#          myrows=1:totinds,
#          mode=fn(l["mod"]),
#          epi=fn(l["epi"])
#          )
#   y=sampleWC(w,
#              b=fn(l["b"]),
#              a=fn(l["a"]),
#              p=fn(l["p"]),
#              rep=1)
#   d[,i]<-y
#   h2<-format(var(w,na.rm = T) / var(y,na.rm = T),digits=2)
#   colnames(d)[i]<- paste0( collapse="_",
#                            c( "ld",
#                               paste0(colnames(mysim), (mysim[i,])),
#                               paste0("h2",h2)
#                             )
#                         )
# }
# head(d)
# d[is.na(d)]<-0

#### write general fam
structure<-data.frame(FID=1:totinds,IID=1:totinds,PAT=0,MAT=0,SEX=0)
# dfam<-cbind(structure,d)
# dim(dfam)
# head(dfam[,1:6])
# write.table(row.names = F,col.names = T, quote = F,
#             file="databig/ldsimexample.fam",
#             dfam
#             )
dfam<-read.table("databig/ldsimexample.fam",header = T)
d<-dfam[,-c(1:5)]

#### write fams per folder
for( i in 1:ncol(d)){
  message(colnames(d)[i])
  system(paste0("mkdir databig/",colnames(d)[i]))
  write.table(row.names = F,col.names = F, quote = F,
            file=paste0("databig/",colnames(d)[i],"/ldexample.fam"),
            cbind(structure,data.frame(PHENO=d[,i]))
            )
  system(paste0("ln databig/example.bed databig/",colnames(d)[i],"/ldexample.bed"))
  system(paste0("ln databig/example.bim databig/",colnames(d)[i],"/ldexample.bim"))
}

#### check things are not crazy
#
# pdf("databig/histograms.pdf")
# for( i in 1:ncol(d)){
#   hist(d[,i],breaks=50,main=
#          paste(paste(  colnames(mysim), (mysim[i,])), collapse="  "))
# }
# dev.off()


# maxes<-apply(d[,1:ncol(d)],2,function(i) {max(i)}) %>% fn
# table(maxes>50)
# toy<-mysim
# mysim$max <- maxes
# View(mysim)
