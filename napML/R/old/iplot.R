

indplot<-function(y,h,inf){
  ## Prepare data
  toplot<-data.frame(y=y, x=inf[h])

  ## Bias and accuracy at the individual level
  lmobj2<-summary(lm(toplot$y~toplot$x))
  accuracy2<-lmobj2$r.squared %>% round(.,digits=3)

  if(dim(lmobj2$coefficients)[1] ==1){
    bias2<-"na"
  }else{
    bias2<-coefficients(lmobj2)[2,1]  %>% round(.,digits=3)
    bias3<-coefficients(lmobj2)[2,2]  %>% round(.,digits=3)
  }

  pind<-ggplot(data=toplot,aes(y=y,x=x)) +
    geom_point(col="grey") +
    stat_smooth(aes(y=y , x=x), se=F,
                method="glm",lty="dashed",col="black")+
          ylim(range(c(toplot$x,toplot$y))) +
          xlim(range(c(toplot$x,toplot$y)))+
          xlab("Inferred individual fitness")+
          ylab("True individual fitness")+
      geom_abline(intercept = 0,slope = 1,lty="dotted")+
      geom_hline(yintercept = 0,lty="dotted")+
      geom_vline(xintercept = 0,lty="dotted")+
      ggtitle(TeX(paste("$R^2$ = ",accuracy2, ", $\\beta = $",bias2,", $\\a = $",bias3 )))
}

# iplotreal<-function(true, inf){
#
#   toplot<-data.frame(x=true, y=inf)
#   toplot$xnozero <- toplot$x
#   toplot$xnozero[toplot$xnozero==0 ] <- NA
#
#   return(iplot_(toplot))
# }
#
# iplot_<-function(toplot){
#     ## Bias and accuracy at the individual level
#     lmobj2<-summary(lm(toplot$y~toplot$xnozero))
#     # lmobj2<-summary(lm(toplot$xnozero ~toplot$y))
#     accuracy2<-lmobj2$r.squared %>% round(.,digits=3)
#     if(dim(lmobj2$coefficients)[1] ==1){
#     bias2<-"na"
#     }else{
#     bias2<-coefficients(lmobj2)[2,1]  %>% round(.,digits=3)
#     }
#
#
#     pind<-ggplot(data=toplot,aes(y=y,x=x)) +
#       geom_point(col="grey") +
#       stat_smooth(aes(y=y , x=xnozero), se=F,
#                   method="glm",lty="dashed",col="black")+
#             ylim(range(c(toplot$x,toplot$y))) +
#             xlim(range(c(toplot$x,toplot$y)))+
#             ylab("Inferred individual fitness")+
#             xlab("True individual fitness")+
#         geom_abline(intercept = 0,slope = 1,lty="dotted")+
#         geom_hline(yintercept = 0,lty="dotted")+
#         geom_vline(xintercept = 0,lty="dotted")+
#         ggtitle(TeX(paste("$R^2$ = ",accuracy2, ", $\\beta = $",bias2)))
#   return(list(pind=pind,accuracy2=accuracy2,bias2=bias2))
# }
