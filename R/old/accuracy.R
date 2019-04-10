saccuracy <-function(s,sinf){
  ## Bias and accuracy at the individual level
  lmobj<-summary(lm(s~sinf))
  accuracy<-lmobj$r.squared %>% round(.,digits=3)

  if(dim(lmobj$coefficients)[1] ==1){
    bias2<-"na"
  }else{
    bias1<-coefficients(lmobj)[2,1]  %>% round(.,digits=3)
    bias2<-coefficients(lmobj)[2,2]  %>% round(.,digits=3)
  }
  return(c("R2"=accuracy,
              "b"=bias1,
              "a"=bias2)
         )
}

iaccuracy<-function(y,h,inf){
    ## Prepare data
    yinf=inf[h]

  ## Bias and accuracy at the selsection level
    lmobj<-summary(lm(y~yinf))
    accuracy<-lmobj$r.squared %>% round(.,digits=3)
    if(dim(lmobj$coefficients)[1] ==1){
      bias<-"na"
    }else{
      bias1<-coefficients(lmobj)[2,1]  %>% round(.,digits=3)
      bias2<-coefficients(lmobj)[2,2]  %>% round(.,digits=3)
    }

  return(c("R2"=accuracy,
            "b"=bias1,
            "a"=bias2)
       )
}


