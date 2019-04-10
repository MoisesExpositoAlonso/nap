parnames<-function() return(c("b","a","p","mu","epi","svar","ss"))

parametersummary<-function(parchain,pnames=c("b","a","p","mu","epi","svar","ss")){
  require(coda)
  r<-as.mcmc(parchain)
  colnames(r) <- pnames
  summary(r)
}

paramplot<-function(parchain,pnames=c("b","a","p","mu","epi","svar","ss"),
                    truevalues=NULL){

  parchain<-data.frame(parchain)
  colnames(parchain)<-pnames

  plotlist<-list()

  for(i in pnames){

    p1<-qplot(parchain[,i],
            x=1:nrow(parchain),
            # ylim=c(0,1),
            geom='line',xlab='iterations',ylab=i) #+geom_hline(yintercept = sdat[[i]], col='grey',lty='dashed')   ,
    p2<-qplot(parchain[,i],geom="density",
            xlab=i,
            # xlim=c(0,1),
            fill=I(transparent("black"))) #+geom_vline(xintercept = sdat[[i]], lty="dashed",col='grey')

    if(!is.null(truevalues)){
      p2+geom_vline(xintercept = truevalues[i],color="red")
    }

    plotlist[[i]]<-plot_grid(p1,p2)
  }

  bigpanel<-plot_grid(plotlist=plotlist,ncol=1)

  return(bigpanel)
}

