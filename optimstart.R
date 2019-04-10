# #
# # napMCMC <- function(y, h, A, m, n, s=rnorm(m,0,0.1),
# #                     b = 0.5, bmin = 0, bmax = 1,
# #                     a = 0.1, amin = 0, amax = 1,
# #                     p = 0.5, pmin = 0, pmax = 1,
# #                     mu = 1, mumin = 0, mumax = 50,
# #                     epi = 1, epimin = 1, epimax = 1,
# #                     svar = 0.1, svarmin = 0.001, svarmax = 1,
# #                     ss = 0.0, ssmin = 0, ssmax = 1,
# #                     smin = -0.98039215686274505668, smax = 10,
# #                     iterations = 1e4L, burnin = 0.1,
# #                     Smode = 1L, LIKmode = 2L, PRImode = 1L, FITmode = 2L,
# #                     verbose = FALSE, test = FALSE,
# #                     updateratio = -9.0, bw = 0.01,
# #                     file2sink = "output.log", sink2file = FALSE,
# #                     RTEST=F) {
#
#  # Calculate likelihood based on original proposal
#  prop=LIKELIHOOD(y,h-1,wC(A,s,FITmode),b,a,p,mu,epi,verbose = F)
#
#
#  ################################################################################
#  # Propose start of S
#  lambdas=c(0.0001,0.001,0.01,0.1,1,10,100)
#  prob=c()
#  for(l in 1:length(lambdas)){
#    prob[l]=LIKELIHOOD(y,h-1,wC(A,BMridge(A,My(y,h),lambda = lambdas[l]),2),
#               b,a,p,mu,1)
#  }
#
#
# ################################################################################
#
#
# # pois example
# x<-rpois(1000,10)
# # x <- rep(obs, freq)
# plot(table(x), main="Count data")
# lklh.poisson <- function(x, lambda) lambda^x/factorial(x) * exp(-lambda)
# log.lklh.poisson <- function(x, lambda){
#                      -sum(x * log(lambda) - log(factorial(x)) - lambda)
# }
# optim(par = 2, log.lklh.poisson, x = x, method = "BFGS")
#
# #my model
# x<-list("y"=y,"h"=h-1,"A"=X)
#
# parstart<-list(
#                 # "lambda"=10,
#                 "b"=0.1,
#                 "a"=0.1,
#                 "p"=0.1,
#                 "mu"=1)
#
# lik.nap<-function(x,par){
#   -LIKELIHOOD(x[[1]],x[[2]],
#               wC(x[[3]],BMridge(x[[3]],My(x[[1]],x[[2]]),lambda = 10 ),2),
#               par[[1]],par[[2]],par[[3]],par[[4]],1)
# }
#
# optim(par = parstart, lik.nap, x = x, method = "L-BFGS-B")
#
