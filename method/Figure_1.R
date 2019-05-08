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
fam<-read.table("databig/simexample.fam",header=T)
fam<-fam[,6:ncol(fam)]
for(i in 1:ncol(fam)){
  fam[,i]<-sort(fam[,i])
  fam[sample(100,1:nrow(fam)),i]<- 0##### Add some zeroes
}
fam[,ncol(fam)] %>% plot
siminfo<-read.csv("databig/simulationgrid.tsv",header=T)
colnames(fam)<-paste0("fam",1:ncol(fam))

### Load simulated fitness with LD
lfam<-read.table("databig/ldsimexample.fam",header=T)
lfam<-lfam[,6:ncol(lfam)]
for(i in 1:ncol(lfam)){
  lfam[,i]<-sort(lfam[,i])
  lfam[sample(50,1:nrow(lfam)),i]<- 0##### Add some zeroes
}
lfam[,ncol(lfam)] %>% plot
ldsiminfo<-read.csv("databig/simulationgrid.tsv",header=T)
colnames(lfam)<-paste0("lfam",1:ncol(lfam))


### Load real fitness ## 513 accessions because in Nature paper some excluded
# fit<-sapply(paste0("databig/rFitness_",fieldcodes(),"/515g.fam"), FUN=function(x)read.table(x)[,6])
# for(i in 1:ncol(fit)){
#   fit[,i]<-sort(fit[,i])
#   fit[,i][ fit[,i]< 0 ] <-0
#   fit[,i][ fit[,i] == -9 ] <-0
# }
# summary(fit)
fit<-read.csv("databig/ARAPHENO Exposito-Alonso 2019 fitness field experiments and climate of origin of 515 accessions.csv",
              header=T)[c(3:10)] # these are relative fitness
fit[is.na(fit)]<-0
for(i in 1:ncol(fit)){
  fit[,i]<-sort(fit[,i])
  fit[,i][ fit[,i]< 0 ] <-0
  fit[,i][ fit[,i] == -9 ] <-0
}
head(fit)
dim(fit)


#### Run pca
dat<-cbind(fit,fam[sort(sample(1:nrow(fam),nrow(fit))),], lfam)
pca<-prcomp(t(dat),scale. = F)
pca$sdev[1:5] / sum(pca$sdev)


to<-data.frame(pc1=pca$x[,1],pc2=pca$x[,2],
                   co=c(1,1,2,2,1,1,2,2,
                        rep(3,ncol(fam)),rep(4,ncol(lfam))))
plot(to$pc1 ~to$pc2)
p<-ggplot(to) +
  geom_point(aes(x=pc1,y=pc2,col=factor(co)))+
  scale_color_manual("reality,simulation",
                     values=c(
                              watercolors()[[2]],watercolors()[[1]],
                              transparent("black"),transparent("grey"))) +
  ylab(paste("PC2", format(pca$sdev[2]/sum(pca$sdev),digits= 2)))+
  xlab(paste("PC1", format(pca$sdev[1]/sum(pca$sdev),digits= 2)))
p

# save_plot("../figs/figure_pcafitness.pdf",base_height = 3,base_width = 5,p)

####************************************************************************####
#### DISTRIBUTION ####

### Load real fitness


fitplot<-fit %>%
        ggplot(.) +
        geom_histogram(aes(x = rFitness_mhi),bins = 60, fill=transparent(watercolors()[2] ,0.8))+
        geom_histogram(aes(x = rFitness_mlp),bins = 60, fill=transparent(watercolors()[1] ,0.8))+
        ylab('Number of genotypes') +
        xlab('Mean-scaled fitness (survival x # seeds)')+
        geom_vline(xintercept = 1,lty='dashed')
fitplot

fitplot<-fam %>%
        ggplot(.) +
        geom_histogram(aes(x = fam1),bins = 60, fill=transparent("grey",0.8))+
        geom_histogram(aes(x = fam10),bins = 60, fill=transparent("black" ,0.8))+
        ylab('Number of genotypes') +
        xlab('Mean-scaled fitness')+
        geom_vline(xintercept = 1,lty='dashed')
fitplot


save_plot("figs/figure_distribution.pdf",base_height = 3,base_width = 5,fitplot)

fitplot<-fit %>%
        ggplot(.) +
        geom_histogram(data=sample_n(fam,515),aes(x = fam1),bins = 60, fill=transparent("black",0.8))+
        geom_histogram(data=sample_n(fam,515),aes(x = fam10),bins = 60, fill=transparent("grey" ,0.8))+
        geom_histogram(aes(x = rFitness_mhi),bins = 60, fill=transparent(watercolors()[2] ,0.8))+
        geom_histogram(aes(x = rFitness_mlp),bins = 60, fill=transparent(watercolors()[1] ,0.8))+
        ylab('Number of genotypes') +
        xlab('Mean-scaled fitness')+
        geom_vline(xintercept = 1,lty='dashed')
fitplot

# save_plot("figs/figure_distribution_realsimulated.pdf",base_height = 3,base_width = 5,fitplot)

####************************************************************************####
#### DISTRIBUTION ####

# ## In file data-cleaning/Yfitness.R
# Y<-readRDS('data/Y.rds')
#
# # pdf("figs/mean_to_sd_ratio.pdf",useDingbats = FALSE,width = 5,height = 5)
# yforvariance=Y
# toplot<-data.frame(y=sqrt(Vy(yforvariance$oFitness,yforvariance$row)),
#        x=My(yforvariance$oFitness,yforvariance$row)
#        )
# p1<-ggdotscolor(toplot$y,
# 						x=toplot$x) %>%
# 	addggregression(se=F) + ylim(c(0,2))+
# 	ylab("Standard Deviation of fitness")+
# 	xlab("Mean fitness")
#
# yforvariance=filter(Y, oFitness!=0)
# toplot<-data.frame(y=sqrt(Vy(yforvariance$oFitness,yforvariance$row)),
#        x=My(yforvariance$oFitness,yforvariance$row)
#        )
# p2<-ggdotscolor(toplot$y,
# 						x=toplot$x) %>%
#   addggregression(se=F) + ylim(c(0,2))+
# 	ylab("Standard Deviation of fitness")+
# 	xlab("Mean fitness")
# # dev.off()

####************************************************************************####
#### KONDRASHOV ####

g<-attach.big.matrix("databig/example.desc")
s<-ssimC(1:100,0.1)
nmut<-apply(g[,1:100],1,function(i){sum(i* s/(abs(s))) })
nmut<-nmut+ abs(min(nmut))

# additive
# y1<-fam[ ,5+dplyr::filter(siminfo, b==0.01, a==0.01, p==0, epi==1, mod==1,svar==0.1)[1,1] ]
# y1%>% hist
# qplot(y1 , x=nmut)
#
# # multiplicative
# y2<-fam[ ,5+dplyr::filter(siminfo, b==0.01, a==0.01, p==0, epi==1, mod==2,svar==0.1)[1,1] ]
# qplot(y2 , x=nmut)
#
# d<-data.frame(a=y1,m=y2,nmut=nmut)
# ggplot(d) +
#   geom_point(aes(y=y1, x=nmut),col=transparent("black" ,0.8))+
#   stat_smooth(aes(y=y1, x=nmut), formula = y~x,method="lm", col=transparent("black" ,0.8))+
#   geom_point(aes(y=y2, x=nmut),col=transparent("grey" ,0.8))+
#   stat_smooth(aes(y=y2, x=nmut), formula = y~exp(x),method="lm",col=transparent("grey" ,0.8))
#
# ggplot(d) +
#   geom_point(aes(y=log(y1), x=nmut),col=transparent("black" ,0.8))+
#   stat_smooth(aes(y=log(y1), x=nmut), formula = y~x,method="lm", col=transparent("black" ,0.8))+
#   geom_point(aes(y=log(y2), x=nmut),col=transparent("grey" ,0.8))+
#   stat_smooth(aes(y=log(y2), x=nmut), formula = y~(x),method="lm",col=transparent("grey" ,0.8))


hist(nmut)
y1<-1 - nmut* 0.01
y2<-1 - exp(nmut* 0.01)
y3<-1 - exp(nmut* 0.01)^1.2
y4<-1 - exp(nmut* 0.01)^1.8

y1<- wCBM(g@address,mycols =1:100,myrows=1:1500 , s=rep(0.01,100), mode=1)
y2<- wCBM(g@address,mycols =1:100,myrows=1:1500,s=rep(0.01,100), mode=2)
y3<- wCBM(g@address,mycols =1:100,myrows=1:1500,s=rep(0.01,100), mode=1,epi=0.6)
y4<- wCBM(g@address,mycols =1:100,myrows=1:1500,s=rep(0.01,100), mode=1,epi=1.4)
y5<- wCBM(g@address,mycols =1:100,myrows=1:1500,s=rep(0.01,100), mode=2,epi=0.6)
y6<- wCBM(g@address,mycols =1:100,myrows=1:1500,s=rep(0.01,100), mode=2,epi=1.4)
nmut<-apply(g[,1:100],1,sum)
# nmut<-wCBM(g[,1:100],mycols =1:100,myrows=1:1500 , s=rep(0.5,100), mode=1)
plot(y1~nmut)
plot(y2~nmut)
plot(y3~nmut)
plot(y4~nmut)

p<-ggplot() +
  geom_point(aes(y=(y1)+rnorm(100,0,0.001), x=nmut),col=transparent("black" ,0.8))+
  # stat_smooth(aes(y=(y1), x=nmut), formula = y~x,method="lm", col=transparent("black" ,0.8),se=F)+

  geom_point(aes(y=(y2)+rnorm(100,0,0.001), x=nmut),col=transparent("darkgrey" ,0.8))+
  # stat_smooth(aes(y=(y2), x=nmut), formula = y~(x),method="lm",col=transparent("darkgrey" ,0.8),se=F)+

  geom_point(aes(y=(y3)+rnorm(100,0,0.001), x=nmut),col=transparent("grey" ,0.8),shape=2)+
  # stat_smooth(aes(y=(y3), x=nmut), formula = y~(x),method="lm",col=transparent("darkgrey" ,0.8),se=F,lty=2)+

  geom_point(aes(y=(y4)+rnorm(100,0,0.001), x=nmut),col=transparent("grey" ,0.8),shape=3)+
  # stat_smooth(aes(y=(y4), x=nmut), formula = y~(x),method="lm",col=transparent("darkgrey" ,0.8),se=F,lty=3)+

  geom_point(aes(y=(y5)+rnorm(100,0,0.001), x=nmut),col=transparent("grey" ,0.8),shape=4)+

  geom_point(aes(y=(y6)+rnorm(100,0,0.001), x=nmut),col=transparent("grey" ,0.8),shape=5)+
  xlab("# of effect mutations") +ylab("individual relative fitness")
p

p+ coord_trans(y="log")

save_plot("figs/figure_kondrashov.pdf",base_height = 3,base_width = 3.5,p)

####************************************************************************####
#### simple gaussian distribution ####
normplot<-ggplot(data.frame(y=rnorm(1e6))) + geom_density(aes(y))
lognormplot<-ggplot(data.frame(y=rlnorm(1e6,sdlog=0.8))) + geom_density(aes(y))
lognormplot
save_plot("figs/gaussian.pdf",base_height = 3,base_width = 3.5,normplot)
save_plot("figs/lgaussian.pdf",base_height = 3,base_width = 3.5,lognormplot)

