napMCMCcompare<-function(iterations=10000){
add<-napMCMCwrap(iterations)
mul<-napMCMCwrap(iterations)
imul<-napMCMCwrap(iterations)


res<-list(
      cbind(data.frame(method=c("add","mul","imul","gwa"),par="s"),
            rbind(
                  saccuracy(s,add$shat),
                  saccuracy(s,mul$shat),
                  saccuracy(s,imul$shat),
                  saccuracy(s,gwa(y,h,A))
                  )
      ),
      cbind(data.frame(method=c("add","mul","imul","gwa","true"),par="i"),
      rbind(
            iaccuracy(y,h,add$w),
            iaccuracy(y,h,mul$w),
            iaccuracy(y,h,imul$w),
            iaccuracy(y,h,A[] %*% BMridge(X,My(y,h))+1),
            iaccuracy(y,h,wC(A[],s,FITmode,epi = epi ,mu= mu))
            )
      )
)

return(res)
}


napMCMCwrap<-function(iterations=10000){
x<-napMCMC(y,h,X,m=1:m,n=1:n,
           rnorm(m,0,0.1),
           test=F,
           verbose=F,
           iterations=iterations,
           FITmode = FITmode,PRImode = 1,LIKmode = 2,
           bw=0.001)

x2<-napMCMC(y,h,X,m=1:m,n=1:n,
           x$shat,
           test=F,
           verbose=F,
           iterations=iterations,
           FITmode = FITmode,PRImode = 1,LIKmode = 2,
           bw=0.05)

x3<-napMCMC(y,h,X,m=1:m,n=1:n,
           x2$shat,
           test=F,
           verbose=F,
           iterations=iterations,
           FITmode = FITmode,PRImode = 1,LIKmode = 2,
           bw=0.001)

chosen<-which(c(x$final_likelihood, x2$final_likelihood, x3$final_likelihood) == max(c(x$final_likelihood, x2$final_likelihood, x3$final_likelihood)))

if(chosen==1){
  return(x)
}else if(chosen==2){
  return(x2)
}else if(chosen==3){
  return(x2)
}else{
  return(x3)
}
}


gwa<-function(y, h, A){
  sinf<-BMridge(A,My(y,h),lambda = 1)
  return(sinf)
}

napMCMC <- function(y, h, A, m, n, s=rnorm(m,0,0.1),
                    b = 0.5, bmin = 0, bmax = 1,
                    a = 0.1, amin = 0, amax = 1,
                    p = 0.5, pmin = 0, pmax = 1,
                    mu = 1, mumin = 0, mumax = 50,
                    epi = 1, epimin = 1, epimax = 1,
                    svar = 0.1, svarmin = 0.001, svarmax = 1,
                    ss = 0.0, ssmin = 0, ssmax = 1,
                    smin = -0.98039215686274505668, smax = 10,
                    iterations = 1e4L, burnin = 0.1,
                    Smode = 1L, LIKmode = 2L, PRImode = 1L, FITmode = 2L,
                    verbose = FALSE, test = FALSE,
                    updateratio = -9.0, bw = 0.01,
                    file2sink = "output.log", sink2file = FALSE,
                    RTEST=F) {

 # Calculate likelihood based on original proposal
 prop=LIKELIHOOD(y,h-1,wC(A,s,FITmode),b,a,p,mu,epi,verbose = F)


 ################################################################################
 # Propose start of S
 lambdas=c(0.0001,0.001,0.01,0.1,1,10,100)
 prob=c()
 sgwa<-c()
 for(l in 1:length(lambdas)){
   # snew
   prob[l]=LIKELIHOOD(y,h-1,wC(A,BMridge(A,My(y,h),lambda = lambdas[l]),FITmode),
              b,a,p,mu,epi,verbose = T)
 }
prob
 ################################################################################
 # Decide start
prob[is.na(prob)] <- Inf
if(is.na(prop)) prop= Inf
 if(all(is.infinite(c(prob))) & is.infinite(prop)){
   stop("Could not find a good starting point for vector s, all proposals are -inf")
 }

if(all(is.infinite(c(prob))) & !is.infinite(prop)){
   message("Proposed start is better than rGWA inferred")
  sstart=s
 }

if(!all(is.infinite(c(prob))) & is.infinite(prop)){
   message("rGWA start for vector s is better than proposed start")
  sstart=BMridge(A,My(y,h),lambda=lambdas[which(prob==max(prob))])
}

if(!all(is.infinite(c(prob))) & !is.infinite(prop)){
  if(max(prob[!is.infinite(prob)],na.rm = T)> prop){
   message("rGWA start for vector s is better than proposed start ;)")
   sstart=BMridge(A,My(y,h),lambda=lambdas[which(prob==max(prob[!is.infinite(prob)]))[1] ])
  }else{
   message("Proposed start is better than rGWA inferred ;)")
   sstart=s
  }
}


################################################################################
 # Try to infer hyperparameters
  b_hat=as.numeric(lm(tapply(y[y!=0],h[y!=0],sd) ~ tapply(y[y!=0],h[y!=0],mean) )$coef[2])
  a_hat=as.numeric(lm(tapply(y[y!=0],h[y!=0],sd) ~ tapply(y[y!=0],h[y!=0],mean) )$coef[1])
  p_hat=length(y[y==0])/length(y)
  mu_hat=mean(y[y!=0])
  pre1=LIKELIHOOD(y,h-1,wC(A,s,FITmode),b_hat,a_hat,p_hat,mu_hat,epi,verbose = F)

  # Common defaults
  b_def=0.1
  a_def=0.1
  p_def=0
  mu_def=1
  pre2=LIKELIHOOD(y,h-1,wC(A,s,FITmode),b_def,a_def,p_def,mu_def,epi,verbose = F)

  if(pre1>prop | pre2>prop){
    message("Inferred starting hyperparameters are better than proposed start")
    if(pre2>pre1){
      b_hat=b_def;
      a_hat=a_def
      p_hat=p_def
      mu_hat=mu_def
    }
  }else{
    message("Proposed hyperparameters are better than the inferred")
    b_hat=b
    a_hat=a
    p_hat=p
    mu_hat=mu
  }
  prop=LIKELIHOOD(y,h-1,wC(A,s,FITmode),b_hat,a_hat,p_hat,mu_hat,epi,verbose = F)

 # Ranges
  err=0.5 # only allow to vary within 5% error during MCMC
    a_hatmin=a_hat*(1-err); if(a_hatmin<0) a_hatmin=0
    a_hatmax=a_hat*(1+err)
    b_hatmin=b_hat*(1-err); if(b_hatmin<0) b_hatmin=0
    b_hatmax=b_hat*(1+err)
    p_hatmin=0;
    p_hatmax=p_hat*(1+err); if(p_hatmax>1) p_hatmax=1
    # mu_hatmin=mu_hat
    # mu_hatmax=mu_hat
    mu_hatmin=mu_hat*(1-err); if(mu_hatmin<0) mu_hatmin=0
    mu_hatmax=mu_hat*(1+err);

################################################################################
  message("The starting likelihood is: ",LIKELIHOOD(y,h-1,
                            wC(A,sstart,FITmode),
                            b_hat,a_hat,p_hat,mu_hat,epi,verbose = T))
 ################################################################################

 # message("First Likelihood ",prob)
 # message("First s ",BMridge(A,My(y,h)))
 # message("First w ",wC(A,BMridge(A,My(y,h)),FITmode)[1])

  message("Calling napMCMC C++ function")
  # Run MCMC
  if(RTEST==F){
    res<-.Call('_nap_napMCMC', PACKAGE = 'nap', y, h, A, m, n, s=sstart,
          b=b_hat, bmin=b_hatmin, bmax=b_hatmax,
          a=a_hat, amin=a_hatmin, amax=a_hatmax,
          p=p_hat, pmin=p_hatmin, pmax=p_hatmax,
          mu=mu_hat, mumin=mu_hatmin, mumax=p_hatmax,
          epi, epimin, epimax,
          svar, svarmin, svarmax,
          ss, ssmin, ssmax,
          smin, smax,
          iterations, burnin, Smode,
          LIKmode, PRImode, FITmode,
          verbose, test, updateratio, bw, file2sink, sink2file)
    # Add names
    rownames(res$par)<-parnames()
  }else{
    res<-.Call('_nap_test_napMCMC', PACKAGE = 'nap', y, h, A, m, n, s=sstart,
          b=b_hat, bmin=b_hatmin, bmax=b_hatmax,
          a=a_hat, amin=a_hatmin, amax=a_hatmax,
          p=p_hat, pmin=p_hatmin, pmax=p_hatmax,
          mu=mu_hat, mumin=mu_hatmin, mumax=p_hatmax,
          epi, epimin, epimax,
          svar, svarmin, svarmax,
          ss, ssmin, ssmax,
          smin, smax,
          iterations, burnin, Smode,
          LIKmode, PRImode, FITmode,
          verbose, test, updateratio, bw, file2sink, sink2file)
  }
  return(res)
}
