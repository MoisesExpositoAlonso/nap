################################################################################
## Run napML inference over simulated traits
################################################################################

## Load packages
library(devtools)
library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(coda)
library(cowplot)
library(RColorBrewer)

library(moiR)
library(nap)

####************************************************************************####
### Load simulated
fam<-read.table("databig/simexample.fam",header=T)
fam<-fam[,6:ncol(fam)]
for(i in 1:ncol(fam)){
  fam[,i]<-sort(fam[,i])
}
fam[,ncol(fam)] %>% plot
siminfo<-read.csv("databig/simulationgrid.tsv",header=T)
colnames(fam)<-paste0("fam",1:ncol(fam))
head(fam)
tail(fam)

summary(fam)

### Load simulated fitness with LD
lfam<-read.table("databig/ldsimexample.fam",header=T)
lfam<-lfam[,6:ncol(lfam)]
for(i in 1:ncol(lfam)){
  lfam[,i]<-sort(lfam[,i])
}
lfam[,ncol(lfam)] %>% plot
ldsiminfo<-read.csv("databig/simulationgrid.tsv",header=T)
colnames(lfam)<-paste0("lfam",1:ncol(lfam))
head(lfam)
tail(lfam)

summary(lfam)

### Load real fitness
fit<-sapply(paste0("databig/rFitness_",fieldcodes(),"/515g.fam"), FUN=function(x)read.table(x)[,6])
for(i in 1:ncol(fit)){
  fit[,i]<-sort(fit[,i])
  fit[,i][ fit[,i]< 0 ] <-0
  fit[,i][ fit[,i] == -9 ] <-0
}

summary(fit)

####************************************************************************####
#Run pca
dat<-cbind(fit,fam[sort(sample(1:nrow(fam),nrow(fit))),] )
pca<-prcomp(t(dat),scale. = T)
# pca<-prcomp(t(fit))
pca$sdev[1:5] / sum(pca$sdev)

plot(pca$x[-c(1:8),1]~pca$x[-c(1:8),2])
points(pca$x[c(1:8),1]~pca$x[c(1:8),2],col="orange",cex=5)

p<-ggplot(data.frame(pc1=pca$x[,1],pc2=pca$x[,2],
                  co=c(rep(1,ncol(x)),rep(2,ncol(fam))) ) ) +
  geom_point(aes(x=pc1,y=pc2,col=factor(co)))+
  scale_color_manual("reality,simulation",values=c("orange","black")) +
  ylab(paste("PC2", format(pca$sdev[2]/sum(pca$sdev),digits= 2)))+
  xlab(paste("PC1", format(pca$sdev[1]/sum(pca$sdev),digits= 2)))
p
save_plot("figs/fitness-real-simulated-pca.pdf",p)


qplot(dat[,7] , dat[,86])

plot(dat[,86])

####************************************************************************####
### Color them based on epistasis

