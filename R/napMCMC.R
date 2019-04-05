napMCMC <- function(y, h, A, m, n, s,
                    b = 0.5, bmin = 0, bmax = 1,
                    a = 0.1, amin = 0, amax = 1,
                    p = 0.5, pmin = 0, pmax = 1,
                    mu = 1, mumin = 0, mumax = 50,
                    epi = 1, epimin = 1, epimax = 1,
                    svar = 0.1, svarmin = 0.001, svarmax = 1,
                    ss = 0.1, ssmin = 0, ssmax = 1,
                    smin = -0.98039215686274505668, smax = 10,
                    iterations = 1e4L, burnin = 0.1,
                    Smode = 2L, LIKmode = 3L, PRImode = 2L, FITmode = 3L,
                    verbose = FALSE, test = FALSE,
                    updateratio = -9.0, bw = 0.01,
                    file2sink = "output.log", sink2file = FALSE) {

 # Infer the mean and variance of population fitness
  err=0.25 # only allow to vary within 5% error during MCMC
  a_hat=as.numeric(lm(tapply(y[y!=0],h[y!=0],sd) ~ tapply(y[y!=0],h[y!=0],mean) )$coef[1])
  a_hat=0.1
    a_hatmin=a_hat*(1-err); if(a_hatmin<0) a_hatmin=0
    a_hatmax=a_hat*(1+err)
  b_hat=as.numeric(lm(tapply(y[y!=0],h[y!=0],sd) ~ tapply(y[y!=0],h[y!=0],mean) )$coef[2])
  b_hat=0.1
    b_hatmin=b_hat*(1-err); if(b_hatmin<0) b_hatmin=0
    b_hatmax=b_hat*(1+err)

  # Infer the proportion of zeros and mean of the population
  p_hat=length(y[y==0])/length(y)
    p_hatmin=0;
    p_hatmax=p_hat*(1+err); if(p_hatmax>1) p_hatmax=1
  mu_hat=mean(y[y!=0])
  # mu_hat=1
    mu_hatmin=mu_hat
    mu_hatmax=mu_hat
    # mu_hatmin=mu_hat*(1-err); if(mu_hatmin<0) mu_hatmin=0
    # mu_hatmax=mu_hat*(1+err);


 prob=NA
 while(is.na(prob) | is.infinite(prob)){
   prob=LIKELIHOOD(y,h-1,wC(A,BMridge(A,My(y,h)),FITmode),
            b_hat,a_hat,p_hat,mu_hat,epi,verbose = T)
 }
 message("First Likelihood ",prob)
 message("First s ",BMridge(A,My(y,h)))
 message("First w ",wC(A,BMridge(A,My(y,h)),FITmode)[1])


  # Run MCMC
  res<-.Call('_nap_napMCMC', PACKAGE = 'nap', y, h, A, m, n, s=BMridge(A,My(y,h)),
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
  return(res)
}


test_napMCMC <- function(y, h, A, m, n, s,
                    b = 0.5, bmin = 0, bmax = 1,
                    a = 0.1, amin = 0, amax = 1,
                    p = 0.5, pmin = 0, pmax = 1,
                    mu = 1, mumin = 0, mumax = 50,
                    epi = 1, epimin = 1, epimax = 1,
                    svar = 0.1, svarmin = 0.001, svarmax = 1,
                    ss = 0.1, ssmin = 0, ssmax = 1,
                    smin = -0.98039215686274505668, smax = 10,
                    iterations = 1e4L, burnin = 0.1,
                    Smode = 2L, LIKmode = 3L, PRImode = 2L, FITmode = 3L,
                    verbose = FALSE, test = FALSE,
                    updateratio = -9.0, bw = 0.01,
                    file2sink = "output.log", sink2file = FALSE) {

  # Infer the mean and variance of population fitness
  err=0.25 # only allow to vary within 5% error during MCMC
  a_hat=as.numeric(lm(tapply(y[y!=0],h[y!=0],sd) ~ tapply(y[y!=0],h[y!=0],mean) )$coef[1])
  a_hat=0.1
    a_hatmin=a_hat*(1-err); if(a_hatmin<0) a_hatmin=0
    a_hatmax=a_hat*(1+err)
  b_hat=as.numeric(lm(tapply(y[y!=0],h[y!=0],sd) ~ tapply(y[y!=0],h[y!=0],mean) )$coef[2])
  b_hat=0.1
    b_hatmin=b_hat*(1-err); if(b_hatmin<0) b_hatmin=0
    b_hatmax=b_hat*(1+err)

  # Infer the proportion of zeros and mean of the population
  p_hat=length(y[y==0])/length(y)
    p_hatmin=0;
    p_hatmax=p_hat*(1+err); if(p_hatmax>1) p_hatmax=1
  mu_hat=mean(y[y!=0])
  mu_hat=1
    mu_hatmin=mu_hat
    mu_hatmax=mu_hat
    # mu_hatmin=mu_hat*(1-err); if(mu_hatmin<0) mu_hatmin=0
    # mu_hatmax=mu_hat*(1+err);

 prob=NA
 while(is.na(prob) | is.infinite(prob)){
   prob=LIKELIHOOD(y,h-1,wC(A,BMridge(A,My(y,h)),FITmode),
            b_hat,a_hat,p_hat,mu_hat,epi,verbose = T)
 }
 message("First Likelihood ")
 message("First s ")
 BMridge(A,My(y,h))
 message("First w ",wC(A,BMridge(A,My(y,h)),FITmode)[1])


  # Run MCMC
  res<-.Call('_nap_test_napMCMC', PACKAGE = 'nap', y, h, A, m, n, s=BMridge(A,My(y,h)),
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
  return(res)
}
