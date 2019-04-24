plotNAP<-function(x){
  p1<-paramplot(x$parchain)
  p0<-posteriorplot
  if(ncol(x$chain)>5) toplot=x$chain[,1:5]
  else toplot=x$chain
  p2<-shist(x$chain)

  plot_grid(
            plot_grid(p0,p2,ncol=1),
            p1,
            ncol=2
  )
}

posteriorplot<-function(post,mylab="Posterior"){
  qplot(post,
            x=1:nrow(post),
            geom='line',xlab='iterations',ylab=mylab) #+
      # geom_vline(xintercept=(1:nrow(x$accept))[x$accept==1], col=transparent("darkgreen",0.2))
}
