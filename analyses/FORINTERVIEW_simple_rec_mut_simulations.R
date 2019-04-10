library(dplyr)
library(cowplot)
library(ggplot2)
library(Rcpp)

sourceCpp("multitools.cpp")

runsim<-function(
                  ind= 100,
                  p=500,
                  # gm=matrix(ncol=p,nrow=n,0),
                  gm=matrix(ncol=p,nrow=ind,sample(c(-1,+1),size = p*ind,replace = T, prob=c(0.9,0.1))),
                  No=5,
                  s=PropoS(p,0.05),
                  # s<-rep(0,p),
                  absfit=0.8,
                  mu= 7e-9,
                  outrate=1-0.98,
                  rectot= 9,
                  gens=5,
                  probmut=c(1,1)
                  ){

    N=rep(No,ind)
    Ntot=c()
    Ntot=rep(0,gens)

    for (i in 1:gens){
      Ntot[i]<-sum(N)
      message("number of genotypes:", ind)
      message("population size:", Ntot[i])

      if(Ntot[i]>1500){
        Ntot[i:gens]<-1500
        message("Population rescued!")
        break
      }
      if(Ntot[i]==0 | ind<=1){ #
        message("population dead!")
        # break
      }else{
          # mutation loop
          for(j in 1:(Ntot[i]*p*mu)){
            indi<-sample(1:ind,1)
            snp<-sample(1:p,1)
            vari<-sample(c(-1,+1),1,prob = probmut,replace=F)
            # message(indi, " ",snp)
            gm[indi,snp] <- (vari)
          }
          # recombination loop
          if(
            (outrate * Ntot[i]) >= 1
            ){
              for( j in 1:ceiling(outrate * Ntot[i])){
                pare=sample(1:ind,size = 2, replace = F)
                breakingpoints<-c(1,sort(sample(p-1,size = rectot,replace = F)),p+1)
                chunks<-sapply(1:(length(breakingpoints)-1),function(l){
                  seq(breakingpoints[l],breakingpoints[l+1]-1)
                })
                pares<-c(matrix(nrow=rectot+1,ncol=1,pare))
                newgenome<-c()
                for( l in 1:(rectot+1) ){
                  newgenome<-append(newgenome,gm[pares[l],unlist(chunks[l])] )
                }
                if(length(newgenome)==p+1) {newgenome<-newgenome[1:p]}
                # print(dim(gm))
                # print(length(newgenome))
                gm<-rbind(gm, matrix(newgenome,nrow=1))
                N<-c(N,1)
                ind=length(N)
              }
          }
          # reproduction loop
          w=wC(gm, s, 1)*absfit
          N_new<-sapply(1:length(N), function(u){
            sum(rpois(w[u],N[u]))
          })
          if(all(N_new<1)){N_new[1]<-1} # check so no breaks
          gm<-gm[which(N_new>0),]
          N<-N[which(N_new>0)]
          ind=length(N)
        }# pop not dead


    }# end loop

    return(list(N=Ntot,div=nrow(gm)))
}

res<-runsim(gens=10,absfit = 0.8)

# wC(res[[2]],s,3)


reps=100
gens=10
manyres<-lapply(1:reps,function(h) moiR::fn(runsim(gens=gens,absfit = 0.8)[[1]]) ) %>% do.call(cbind,.)
manyres_<-data.frame(N=c(manyres),g=1:gens,rep=sort(rep(1:reps,gens)))

# recboos<-lapply(1:reps,function(h) moiR::fn(runsim(rectot = 9*5,outrate = 0.9,mu= 7e-8,gens=gens,probmut=c(1,2),absfit = 0.8)[[1]]) ) %>% do.call(cbind,.)
recboos2<-lapply(1:reps,function(h) moiR::fn(runsim(rectot = 9*5,outrate = 0.95,mu= 7e-8,gens=gens,probmut=c(1,2),absfit = 0.8)[[1]]) ) %>% do.call(cbind,.)
# recboos_<-data.frame(N=c(recboos),g=1:gens,rep=sort(rep(1:reps,gens)))
recboos_<-data.frame(N=c(recboos2),g=1:gens,rep=sort(rep(1:reps,gens)))

ggplot() +
 stat_summary(data=manyres_,aes(x=g,y=N))+
 stat_summary(data=recboos_,aes(x=g,y=N),col="#e0903a")+
  ylab("Population size")+
  xlab("Generations")+
  scale_x_continuous(breaks = c(1:10))

p<-ggplot() +
 # geom_point(data=manyres_,aes(x=g,y=N),color="#e0903a", alpha=0.5)+
 # geom_point(data=recboos_,aes(x=g,y=N),col="#90c643",alpha=0.5)+
 geom_line(data=manyres_,aes(x=g,y=N,group=rep),color="#4286f4", alpha=0.1)+
 geom_line(data=recboos_,aes(x=g,y=N,group=rep),col="#90c643",alpha=0.1)+
 stat_smooth(data=manyres_,aes(x=g,y=N),col="#4286f4",se=F,lwd=2,method="loess",span=0.8)+
 stat_smooth(data=recboos_,aes(x=g,y=N),col="#90c643",se=F,lwd=2,method="loess",span=0.8)+
 # stat_smooth(data=manyres_,aes(x=g,y=N),col="#e0903a",se=F,lwd=2, method="lm",formula=y~500+poly(x,4))+
 # stat_smooth(data=recboos_,aes(x=g,y=N),col="#90c643",se=F,lwd=2, method="lm",formula=y~500+poly(x,4))+
 # stat_summary(data=manyres_,aes(x=g,y=N),col="#e0903a",geom="line")+
 # stat_summary(data=recboos_,aes(x=g,y=N),col="#90c643",geom="line")+
 # stat_summary(data=manyres_,aes(x=g,y=N),col="#e0903a")+
 # stat_summary(data=recboos_,aes(x=g,y=N),col="#90c643")+
  ylab("Population size")+
  xlab("Generations")+
  scale_x_continuous(breaks = c(1:10))
p

p<-ggplot() +
 # geom_point(data=manyres_,aes(x=g,y=N),color="#e0903a", alpha=0.5)+
 # geom_point(data=recboos_,aes(x=g,y=N),col="#90c643",alpha=0.5)+
 geom_line(data=manyres_,aes(x=g,y=N,group=rep),color="grey", alpha=0.1)+
 geom_line(data=recboos_,aes(x=g,y=N,group=rep),col="black",alpha=0.1)+
 stat_smooth(data=manyres_,aes(x=g,y=N),col="grey40",se=F,lwd=2,method="loess",span=0.8, lty="dotted")+
 stat_smooth(data=recboos_,aes(x=g,y=N),col="black",se=F,lwd=2,method="loess",span=0.8, lty="dashed")+
 # stat_smooth(data=manyres_,aes(x=g,y=N),col="#e0903a",se=F,lwd=2, method="lm",formula=y~500+poly(x,4))+
 # stat_smooth(data=recboos_,aes(x=g,y=N),col="#90c643",se=F,lwd=2, method="lm",formula=y~500+poly(x,4))+
 # stat_summary(data=manyres_,aes(x=g,y=N),col="#e0903a",geom="line")+
 # stat_summary(data=recboos_,aes(x=g,y=N),col="#90c643",geom="line")+
 # stat_summary(data=manyres_,aes(x=g,y=N),col="#e0903a")+
 # stat_summary(data=recboos_,aes(x=g,y=N),col="#90c643")+
  ylab("Population size")+
  xlab("Generations")+
  scale_x_continuous(breaks = c(1:10))
p
# save_plot("figs/for_future_research.pdf",p,base_height = 3,base_width = 4)
# save_plot("figs/for_future_research2.pdf",p,base_height = 3,base_width = 4)
save_plot("figs/for_future_research3.pdf",p,base_height = 2.8,base_width = 3.1)


(tail(head(manyres,5),1) >= 500) %>% table
(tail(head(recboos2,5),1) >= 500) %>% table


# ggplot() +
#  geom_line(data=manyres_,aes(x=g,y=N,group=rep))+
#  geom_line(data=recboos_,aes(x=g,y=N,group=rep),col="#e0903a")+
#   ylab("Population size")+
#   xlab("Generations")
#
#
