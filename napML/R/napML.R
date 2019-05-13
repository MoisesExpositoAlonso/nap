####************************************************************************####
# May 12 functions reworked from successful example April 28 and
# wihout box constraints to improve BFGS algorithm


likelihoodR<-function(y, w, b,a, p){
  mymin<-  -(.Machine$double.xmax / (length(y)))

  tmp<-ifelse(y==0,
              p  + (1-p) *  pnorm(0,w,a+w*b,TRUE,FALSE),
              (1-p) * dnorm(y,w,a+w*b,FALSE)
             )
  tmp<-log(tmp)
  tmp[is.infinite(tmp) & tmp<0 ] <- mymin
  tmp[is.infinite(tmp) & tmp>0 ] <- -(mymin)
  tmp[is.na(tmp)]<-mymin
  LL<-sum(tmp)

  if(is.infinite(LL) & LL<0) LL<-mymin
  if(is.na(LL)) LL<-mymin
  return(LL)
}

####************************************************************************####
lik.nap.free<-function(y,h_,m_,A,par,mod,epi){
    w = wCBM(A = A@address,
      s= istr(par[1:length(m_)]),
      mycols =m_,
      myrows = h_,
      mode=mod,
      epi=epi
      )
   cor(w,y)
  -likelihoodR(
              y = y,
              w,
              b=iptr(par[length(m_)+1]),
              a=iptr(par[length(m_)+2]),
              p=iptr(par[length(m_)+3])
              )
}
####************************************************************************####
napSPG<-function(y,h_,m_,A,s,mod=1,epi=1,iter=20){
  parstart<-par<-list(
                  "s"=str(s),
                  "b"=ptr(0.2),
                  "a"=ptr(0.2),
                  "p"=ptr(0.2)
                  ) %>% unlist
start_time <- Sys.time()
cat("Spectral Projected Gradient SPG \n")
cat("Starting",iter,"iterations ... \n")
r<-spg(
        fn= lik.nap.free,
        y=y,
        A=A,
        h_=h_,
        m_=m_,
        mod=mod,
        epi=epi,
        par = parstart,
        control = list(trace=1,maxit=iter)
        )
diftime<-Sys.time() - start_time
cat("Finished after",as.numeric(diftime, units = "mins"), "mins \n")
return(r)
}

napML<-function(y,h_,m_,A,s,mod=1,epi=1,iter=20){
  parstart<-par<-list(
                  "s"=str(s),
                  "b"=ptr(0.2),
                  "a"=ptr(0.2),
                  "p"=ptr(0.2)
                  ) %>% unlist
start_time <- Sys.time()
cat("Broyden-Fletcher-Goldfarb-Shanno BFGS \n")
cat("Starting",iter,"iterations ... \n")
r<-optim(
        fn= lik.nap.free,
        y=y,
        A=A,
        h_=h_,
        m_=m_,
        mod=mod,
        epi=epi,
        par = parstart,
        method = "BFGS",
        control = list(trace=1,maxit=iter)
        )
diftime<-Sys.time() - start_time
cat("Finished after",as.numeric(diftime, units = "mins"), "mins \n")
return(r)
}
####************************************************************************####
ptr<-function(p) sapply(p,function(i) 1+(-log((1-i)/(i))))
iptr<-function(h) sapply(h,function(i) (1/(1+exp(1-i))))
str<-function(s) sapply(s,function(i) 1+ log(1+i))
istr<-function(x) sapply(x,function(i) exp(i-1) -1 )


####************************************************************************####
parseoptim<-function(r,hpar=3){
  s<-r$par[1:(length(r$par)-hpar)]
  b=r$par[length(s)+1]
  a=r$par[length(s)+2]
  p=r$par[length(s)+3]
  # e=r$par[length(s)+4]
  AIC= 2* (length(r$par)) -2 * (-r$value)
  return(list(
    s=istr(s),
    b=iptr(b),
    a=iptr(a),
    p=iptr(p),
    # e=e,
    AIC=AIC
  ))
}


accuracyresults<-function(pheno,mod,epi,wgwa,winf,wtest,y,ytest,finaltsv){
  d<-data.frame(phenotype=pheno,
              method=c("gwa","nap","nap_t"),
              mod=c(NA,mod,mod),
              epi=c(NA,epi,epi)
              )

  d_<-rbind(
    accuracies(wgwa,y),
    accuracies(winf,y),
    accuracies(wtest,ytest)
  )
  d<-data.frame(cbind(d,d_))
  write.table(row.names = F, quote = F, sep = "\t",
              x=sapply(d,as.character),
              file=finaltsv
              )
  return(d)
}


####************************************************************************####


#
# ####************************************************************************####
# #### example ####
# pois example
# x<-rpois(1000,10)
# # x <- rep(obs, freq)
# plot(table(x), main="Count data")
# lklh.poisson <- function(x, lambda) lambda^x/factorial(x) * exp(-lambda)
# log.lklh.poisson <- function(x, lambda){
#                      -sum(x * log(lambda) - log(factorial(x)) - lambda)
# }
# optim(par = 2, log.lklh.poisson, x = x, method = "BFGS")
#


# ###### Optimization
# napoptim<-function(s,lik.nap,y,hall,m,G,mod,epi,maxit = 500){
#   parstart<-list(
#                   "s"=s,
#                   "b"=0.5,
#                   "a"=0.5,
#                   "p"=0.1
#                   ) %>% unlist
#   parlow<- list(
#               # "s"=s-sd(s)*2,
#               "s"=rep(-0.99, length(s)),
#                 "b"=0.0,
#                 "a"=0.0,
#                 "p"=0
#                 ) %>% unlist
#   parhigh<- list(
#               # "s"=s+sd(s)*2,
#               "s"=rep(0.99, length(s)),
#                 "b"=1,
#                 "a"=1,
#                 "p"=1
#                 )%>% unlist
#
#   cat("Likelihood with starting parameters ")
#   cat(lik.nap(y=y,
#           h_=hall,
#           m_=m,
#           A=G,
#           par=c(s,parstart["b"],parstart["a"],parstart["p"]),
#           mod=mod,
#           e=parstart["e"]
#           )
#   )
#   cat("\n")
#   cat("Low-storage Broyden-Fletcher-Goldfarb-Shanno L-BFGS ... \n")
#   r <- tryCatch(
#                   # lbfgs(
#                   #     x0 = parstart ,
#                   optim(
#                       par = parstart ,
#                       fn=lik.nap,
#                       y=y,
#                       h_=hall,
#                       m_=m,
#                       A=G,
#                       mod=mod,
#                       e=epi,
#                       lower = parlow,
#                       upper=parhigh,
#                       control=list(maxit=maxit,trace=1),
#                       method="L-BFGS-B"
#                       )
#                    ,
#                    error = function(e){list(value = NA)}
#                    )
#   if(is.na(r$value)){
#     cat("Spectral Projected Gradient method SPG ... \n")
#     r <- tryCatch(
#                     # spg( # wrappter for spg
#                     BBoptim( # wrappter for spg
#                         par = parstart ,
#                         fn=lik.nap,
#                         y=y,
#                         h_=hall,
#                         m_=m,
#                         A=G,
#                         mod=mod,
#                         e=epi,
#                         lower = parlow,
#                         upper=parhigh,
#                         control=list(maxit=maxit)
#                         )
#                      ,
#                      error = function(e){list(value = NA)}
#                      )
#   }
#
#   # results<-(c(r1$value,r2$value))
#   # if(!all(is.na(results))){
#   #   decidemode<-which(results == min(results,na.rm = T))
#   #   cat(paste("Best model is ", c("spg","lbfgsb")[decidemode],"\n"))
#   #   r<-list(r1,r2)[head(decidemode,1)]
#   # }else{
#   #   r<-list(value = NA)
#   # }
#   return(r)
# }
#
# ###### R implementation
# # likelihoodR<-function(y, w, b,a, p){
# #   mymin<-  -(.Machine$double.xmax / (length(y)))
# #
# #   tmp<-ifelse(y==0,
# #               p  + (1-p) *  pnorm(0,w,a+w*b,TRUE,FALSE),
# #               (1-p) * dnorm(y,w,a+w*b,FALSE)
# #              )
# #   tmp<-log(tmp)
# #   tmp[is.infinite(tmp) & tmp<0 ]<- mymin
# #   tmp[is.infinite(tmp) & tmp>0 ]<- -(mymin)
# #   tmp[is.na(tmp)]<-mymin
# #   LL<-sum(tmp)
# #
# #   if(is.infinite(LL)) LL<-log(.Machine$double.xmin)
# #   if(is.na(LL)) LL<-log(.Machine$double.xmin)
# #   return(LL)
# # }
# likelihoodR<-function(y, w, b,a, p){
#   mymin<-  -(.Machine$double.xmax / (length(y)))
#
#   tmp<-ifelse(y==0,
#               p  + (1-p) *  pnorm(0,w,a+w*b,TRUE,FALSE),
#               (1-p) * dnorm(y,w,a+w*b,FALSE)
#              )
#   tmp<-log(tmp)
#   tmp[is.infinite(tmp) & tmp<0 ] <- mymin
#   tmp[is.infinite(tmp) & tmp>0 ] <- -(mymin)
#   tmp[is.na(tmp)]<-mymin
#   LL<-sum(tmp)
#
#   if(is.infinite(LL) & LL<0) LL<-mymin
#   if(is.na(LL)) LL<-mymin
#   return(LL)
# }
#
# lik.nap<-function(y,h_,m_,A,par,mod,e,debug=F){
#   h=h_
#   m=m_
#   w=wCBM(A = A@address,
#             s=par[1:length(m)],
#             mycols =m,
#             myrows =h,
#             mode=mod,
#             epi=e)
#   LL<-likelihoodR(y = y,
#             w = w,
#             b=par[length(m)+1],
#             a=par[length(m)+2],
#             p=par[length(m)+3]
#             )
#   if(debug) print(LL)
#   if(debug) print(par[seq(length(m)+1,length(par))])
#   if(debug) print(par[1:5])
#
#   return( - LL) # DO NOT FORGET MINUS
# }
#
