load_all('.')

####************************************************************************####
#### estimate marginal s

x<-attach.big.matrix("databig/example.desc")
y<-read.table("databig/example.fam",header=T)[,6]
pheno<-colnames(read.table("databig/example.fam",header=T))[6]
y<-read.table("databig/example.fam",header=T)[,80]
pheno<-colnames(read.table("databig/example.fam",header=T))[80]
qplot(y,main=pheno)


s<- BMs(x@address,y,1:ncol(x),1:nrow(x))
qplot(s)


qplot( y,
       wCBM(x@address,s,1:ncol(x),1:nrow(x), mode = 1)
       )

yinf=(wCBM(x@address,s,1:ncol(x),1:nrow(x), mode = 1) )
mlin<-lm(  y ~ yinf )%>% summary
mpol<-lm(  y ~ yinf + I(yinf^2) )  %>% summary
mpol2<-lm(  y[y!=0] ~ yinf[y!=0] + I(yinf[y!=0]^2) )  %>% summary
plot(y[y!=0] ,x= yinf[y!=0])

qplot( y,
      wCBM(x@address,s,1:ncol(x),1:nrow(x), mode = 2)
       )
