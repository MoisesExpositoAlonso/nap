

scorplot <-function(s,sinf,sinf_range=NULL){

    ## Bias and accuracy at the selsection level
    lmobj<-summary(lm(s~sinf))
    # lmobj<-summary(lm(sinf ~s))
    accuracy<-lmobj$r.squared %>% round(.,digits=3)
    if(dim(lmobj$coefficients)[1] ==1){
      bias<-"na"
    }else{
      bias<-coefficients(lmobj)[2,1]  %>% round(.,digits=3)
      bias2<-coefficients(lmobj)[2,2]  %>% round(.,digits=3)
    }

    if( is.null(sinf_range)) sinf_range =cbind(sinf,sinf)

    ## PLot selection coefficients
  psel<-qplot(y=s,x=sinf,
              xlim=c(range(c(s,sinf_range))),
              ylim=c(range(c(s,sinf_range)))
              ) +
  # geom_segment(aes(x=s,xend=s,y=sinf_range[,1],yend=sinf_range[,2]) )+
      geom_abline(intercept = 0,slope = 1,lty="dotted")+
      xlab("Inferred selection coefficients")+
      ylab("True selection coefficients")+
      ggtitle(TeX(paste("$R^2$ = ",accuracy, ", $\\beta = $",bias,", $\\a = $",bias2 )))
    # print(psel)
    if(!is.null(sinf_range)) {
      psel<-psel+geom_segment(aes(x=s,xend=s,y=sinf_range[,1],yend=sinf_range[,2]) )
    }
  return(list(psel=psel,accuracy=accuracy,bias=bias))
}

ssummary<-function(parchain){
  require(coda)
  r<-as.mcmc(parchain)
  colnames(r) <- paste0("SNP",1:dim(parchain)[2])
  summary(r)
}


shist <- function(parchain){

  if(ncol(parchain)>5) parchain=parchain[,1:5 ]

  parchain<-data.frame(parchain)
  colnames(parchain)<- paste0("SNP",1:dim(parchain)[2])

  pnames<-paste0("SNP",1:dim(parchain)[2])

  plotlist<-list()

  for(i in pnames){
    panel<-plot_grid(
      qplot(parchain[,i],
            x=1:nrow(parchain),
            # ylim=c(0,1),
            geom='line',xlab='iterations',ylab=i), #+geom_hline(yintercept = sdat[[i]], col='grey',lty='dashed')   ,
      qplot(parchain[,i],geom="density",
            xlab=i,
            # xlim=c(0,1),
            fill=I(transparent("black"))) #+geom_vline(xintercept = sdat[[i]], lty="dashed",col='grey')

    )
    plotlist[[i]]<-panel
  }

  bigpanel<-plot_grid(plotlist=plotlist,ncol=1)
  return(bigpanel)
}
