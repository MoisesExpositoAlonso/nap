### Assemble figure 1
## depends on files
# figure_fitness_pca.R
# figure_fitness_distribution.R
# figure_kondrasov.R
# figure_variance_mean_fitness.

library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(bigmemory)
library(moiR)
library(nap)

####************************************************************************####
#### PCA ####
fam<-read.table("../databig/simexample.fam",header=T)
fam<-fam[,6:ncol(fam)]
for(i in 1:ncol(fam)){
  fam[,i]<-sort(fam[,i])
  fam[sample(100,1:nrow(fam)),i]<- 0##### Add some zeroes
}
fam[,ncol(fam)] %>% plot
siminfo<-read.csv("../databig/simulationgrid.tsv",header=T)
colnames(fam)<-paste0("fam",1:ncol(fam))
head(fam)
tail(fam)

summary(fam)

### Load simulated fitness with LD
lfam<-read.table("../databig/ldsimexample.fam",header=T)
lfam<-lfam[,6:ncol(lfam)]
for(i in 1:ncol(lfam)){
  lfam[,i]<-sort(lfam[,i])
  lfam[sample(50,1:nrow(lfam)),i]<- 0##### Add some zeroes
}
lfam[,ncol(lfam)] %>% plot
ldsiminfo<-read.csv("../databig/simulationgrid.tsv",header=T)
colnames(lfam)<-paste0("lfam",1:ncol(lfam))
head(lfam)
tail(lfam)

summary(lfam)



### Load real fitness
# fit<-sapply(paste0("../databig/rFitness_",fieldcodes(),"/515g.fam"), FUN=function(x)read.table(x)[,6])
# for(i in 1:ncol(fit)){
#   fit[,i]<-sort(fit[,i])
#   fit[,i][ fit[,i]< 0 ] <-0
#   fit[,i][ fit[,i] == -9 ] <-0
# }
# summary(fit)
fit<-read.csv("../databig/ARAPHENO Exposito-Alonso 2019 fitness field experiments and climate of origin of 515 accessions.csv",
              header=T)[c(3:10)] # these are relative fitness
fit[is.na(fit)]<-0
for(i in 1:ncol(fit)){
  fit[,i]<-sort(fit[,i])
  fit[,i][ fit[,i]< 0 ] <-0
  fit[,i][ fit[,i] == -9 ] <-0
}
head(fit)
dim(fit)


#### Additional simultion

G<-attach.big.matrix("../databig/example.desc")
b=0.5
a=0.5
p=0.4
svar=0.1
epi=1.2
mod=2

s = c(exp(rnorm(500,0,svar))-1 )
hist(s)
# s = ssimC(1:500,svar=svar)
w=wCBM(G@address,
       s,
       mycols=1:500,
       myrows=1:515,
       mode=mod,
       epi=epi
       )
y=sampleWC(w,
           b=b,
           a=a,
           p=p,
           rep=1)

#### Run pca
dat<-cbind(fit,fam[sort(sample(1:nrow(fam),nrow(fit))),], lfam , sort(y))
pca<-prcomp(t(dat),scale. = T)
# pca<-prcomp(t(fit))
pca$sdev[1:5] / sum(pca$sdev)


to<-data.frame(pc1=pca$x[,1],pc2=pca$x[,2],
                   co=c(rep(1,ncol(fit)),rep(2,ncol(fam)),rep(3,ncol(lfam)) ,4 ))

p<-ggplot(to) +
  geom_point(aes(x=pc1,y=pc2,col=factor(co)))+
  scale_color_manual("reality,simulation",values=c("orange","black","grey","red")) +
  ylab(paste("PC2", format(pca$sdev[2]/sum(pca$sdev),digits= 2)))+
  xlab(paste("PC1", format(pca$sdev[1]/sum(pca$sdev),digits= 2)))
p
