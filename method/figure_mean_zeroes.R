Y<-readRDS('data/Y.rds')

hist(Y$rFitness, breaks=1000)
hist(Y$rFitness, breaks=100)


mymeans<-tapply(Y$Fitness, Y$sample.ID,mean)
myzer<-tapply(Y$Fitness, Y$sample.ID,countzer)
myzersim<-sapply(mymeans,FUN = function(x){sum(rnorm(5,x,sd=x)<=0 )})
myzersim<-sapply(mymeans,FUN = function(x){sum(rnorm(5,x,sd=0.5*x)<=0 )})
myzersim<-sapply(mymeans,FUN = function(x){sum(rnorm(5,x,sd=mean(mymeans))<=0 )})

qplot(y=mymeans, x=myzer, col=transparent("purple")) + geom_point(aes(y=mymeans,x=myzersim),col=transparent("orange"),shape=2)
cor(mymeans,myzer,method="spearman")
cor(mymeans,myzersim,method="spearman")
