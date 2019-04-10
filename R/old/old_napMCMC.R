

#####**********************************************************************#####
napMCMCR<-function (y=rnorm(n,1,1),
                   h=1:100,
                   A=Xsim(100,500,mafsim(500)),
                   s=ssim(500,0.1),
                   m=1:500,
                   n=1:100,
                  b = 0.5, bmin = 0, bmax = 1,
                  a = 0.4, amin = 0, amax = 1,
                  p = 0.5,
                  mu = 1, mumin = 0,mumax = 10,
                  epi = 1, epimin = 1, epimax = 1,
                  svar = .1, svarmin = 0, svarmax = 2,
                  ss = .1, ssmin = 0, ssmax = 1,

                  bw =  ifelse(Fitnessmode ==3L, 1, 0.1),
                  nupdates = 1L,

                  min = ifelse(Fitnessmode ==3L, 0, -1),
                  max = ifelse(Fitnessmode ==3L, +2.5, 1),

                  iterations = 10000,
                  TEST = FALSE,
                  verbose = FALSE,
                  debug = FALSE,
                  Fitnessmode = 1L,
                  Priormode = 1L,
                  Proposalmode = ifelse(Fitnessmode ==3L, 2L, 1L),
                  file2sink = "output.log", sink2file = FALSE){

      # parnames.push_back("b");
      # parnames.push_back("a");
      # parnames.push_back("p");
      # parnames.push_back("mu");
      # parnames.push_back("epi");
      # parnames.push_back("svar");
      # parnames.push_back("ss");

  # If starting value epi is not 1, do not try to estimate
  if(epi!=1 & (epimin !=1 | epimax!=1) ){
    epimin=0.5
    epimax=2
  }

  # if(mu!=1 & (mumin !=1 | mumin!=1)){
  #   mu=mean(y[y!=0])
  #   mumin=0.1
  #   mumax=quantile(y[y!=0],probs = 0.75)
  # # }

  # Guess better starting values
  # message('Generate starting values from marginal GWA')
  # sprior<-BMgwa(X0=A[n,m],yraw=My(y,h),type=1)$coefficients
  # rank(sprior,ties.method = 'first')
  # sproposal= PropoS(length(m), svar=svar)
  # sstart=  sort(sproposal)[rank(sprior,ties.method = 'first')]
  # sstart=  sprior
  sstart=PropoS(length(m), svar=svar)

  message('Running napMCMC')
  napMCMCC(y, h,
          # A@address,
          A[n,m],
          sstart, m, n,
    b, bmin, bmax, a, amin, amax, p, mu, mumin, mumax, epi, epimin,
    epimax, svar,svarmin,svarmax,ss,ssmin,ssmax, bw, nupdates, min, max, iterations, TEST, verbose,
    debug, fn(Fitnessmode), fn(Priormode), fn(Proposalmode), file2sink, sink2file)

  # .Primitive(".Call")(<pointer: 0x12df0adc0>, y, h, A, s, m, n,
  #   b, bmin, bmax, a, amin, amax, p, mu, mumin, mumax, epi, epimin,
  #   epimax, bw, nupdates, min, max, iterations, TEST, verbose,
  #   debug, Fitnessmode, Priormode, Proposalmode, file2sink, sink2file)
}

Wgo<- function(X, s, mode, w=NULL,epi = 1, ref = 1) {
  if(is.null(w)){
  #   .Call('_nap_W_go', PACKAGE = 'nap', X, s, mode, epi, ref)
    W_go (X, s, mode, epi, ref)

  }else{
  #   .Call('_nap_W_update', PACKAGE = 'nap', X, s, mode, w, epi, ref)
    W_update (X, s, mode, w, epi, ref)

  }
}


#####**********************************************************************#####
#### Test fitness ####

testFITNESS<-function(){
  message("Simulating 100 ind and 100 SNPs")
  X<-matrix(ncol=100,nrow=100,sample(c(-1,1), 100*100,replace=T))
  message("selection coefficients uniform" )
  s = exp(rnorm(100,0,0.1)) - 1
  # s = rnorm(100,0,0.5)
  # s[s<0]<-0

  p<-plot_grid(
    nrow=2,
    labels=c('100 selection coefficients','100 individuals'),
    qplot(s,geom='histogram',xlab='s',bins=30),
    plot_grid(
                        qplot(Wgo(X,s,1),geom='histogram',bins=30,xlab='w Multiplicative'),
                        qplot(Wgo(X,s,2),geom='histogram',bins=30,xlab='w Additive'),
                        qplot(Wgo(X,s,3),geom='histogram',bins=30,xlab='w Inverse multiplicative'),
                        ncol=3
                      )
  )
  suppressWarnings(p)
}
# testFITNESS()


#####**********************************************************************#####
test_GPROPOSALR<-function(){
  x<-t( test_GPROPOSAL() )
  colnames(x)<-c('b','a','p','mu','epi','svar','ss')
  x<-data.frame(x)
  p<-ggplot(tidyr::gather(x), aes(value)) +
    xlab('')+
    geom_histogram(bins=30) +
    facet_wrap(~key, scales = 'free_x') +
    ggtitle('Proposals of hyperparameters')
  return(p)
}
#test_GPROPOSALR()


#####**********************************************************************#####

test_SPROPOSALR<-function(mode=2){
  x<-(PropoS(mode=mode))
  # x<-as.mcmc(x)
  # plot(x)

  # colnames(x)<-c('b','a','p','mu','epi','svar','ss')
  x<-data.frame(x)
  p<-ggplot(tidyr::gather(x), aes(value)) +
    xlab('')+
    geom_histogram(bins=30) +
    facet_wrap(~key, scales = 'free_x') +
    ggtitle('Proposals of S')

  if(mode==2){
    p2<-ggplot(tidyr::gather(log(1+x)), aes(value)) +
    xlab('')+
    geom_histogram(bins=30) +
    facet_wrap(~key, scales = 'free_x')+
          ggtitle('in the space log(1+s)')
  p<-plot_grid(p,p2,nrow=2)
  }
  return(p)
}

#####**********************************************************************#####

test_PRIORR<-function(min=-1,
                      max=10,
                      cent=0,
                      svar=0.5,
                      sparsity=0.1,
                      mod=2){

  x<-PropoS(mode=mod)

  qplot(x[,1],xlab='s')

  y<-test_PRIOR(s=as.numeric(x[,1]),
               min= min,
               max = max,
               mean= cent,
               svar=svar,
               sparsity=sparsity) # this returns 4 values

  return(y)
}

#####**********************************************************************#####


test_LIKELIHOODR<-function(

                            test=FALSE){
  s<-PropoS(mode=2)
  # s = exp(rnorm(100,0,0.1)) - 1
  message("Simulating 100 ind and 100 SNPs")
  X<-matrix(ncol=100,nrow=100,sample(c(-1,1), 100*100,replace=T))
  message("selection coefficients uniform" )

  test_LIKELIHOOD( s=s

  )
  # p<-plot_grid(
  #   nrow=2,
  #   labels=c('100 selection coefficients','100 individuals'),
  #   qplot(s,geom='histogram',xlab='s',bins=30),
  #   plot_grid(
  #                       qplot(Wgo(X,s,1),geom='histogram',bins=30,xlab='w Multiplicative'),
  #                       qplot(Wgo(X,s,2),geom='histogram',bins=30,xlab='w Additive'),
  #                       qplot(Wgo(X,s,3),geom='histogram',bins=30,xlab='w Inverse multiplicative'),
  #                       ncol=3
  #                     )
  # )

}

