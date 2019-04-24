lik.nap.add<-function(y,h,A,par){
  -LIKELIHOOD(y = y,
              hs = h,
              w = wCBM(A = A,
                  s=par[1:500],
                  mycols =1:500,
                  myrows = 1:1500,
                  mode=1),
              b=par[500+1],
              a=par[500+2],
              p=par[500+3],
              epi=par[500+4],
              mu=1
              )
}

napML<-function(yy,hh,AA,s){
  parstart<-list(
                  "s"=s,
                  "b"=0.1,
                  "a"=0.1,
                  "p"=0.1,
                  "epi"=1
                  ) %>% unlist
  parlow<- list(
                "s"=rep(-0.98039215686274505668,500),
                "b"=0.001,
                "a"=0.001,
                "p"=0.0,
                "epi"=0.5
                ) %>% unlist
  parhigh<- list(
                "s"=rep(10,500),
                "b"=1,
                "a"=1,
                "p"=1,
                "epi"=1.5
                ) %>% unlist
  optim(fn= lik.nap.add, 
        y=yy,
        h=hh,
        A=AA,
        par = parstart, 
        lower = parlow,
        upper=parhigh,
        method = "L-BFGS-B"
        )
}


# lik.nap.add(y,h,A,parstart)


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
