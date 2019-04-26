#' Likelihood function of NAP ML model
#'
#' @param y
#' @param h
#' @param m
#' @param A
#' @param par
#' @param mod
#' @param e
#'
#' @return
#' @export
#'
#' @examples
lik.nap<-function(y,h,m,A,par,mod,e){
  # print(par[seq(length(m)+1,length(par))])
  # print(summary(par[1:length(m)]))
  # print(paste(par[1:length(m)],collapse=","))
  w=wCBM(A = A@address,
            s=par[1:length(m)],
            mycols =m,
            myrows =h,
            mode=mod)
  -LIKELIHOOD(y = y,
              w = w,
              b=par[length(m)+1],
              a=par[length(m)+2],
              p=par[length(m)+3],
              mu=1,
              epi=e
              # epi=par[length(m)+4]
              )
}

#' napML call
#'
#' @param y
#' @param h
#' @param m
#' @param A
#' @param mod
#' @param e
#' @param s
#' @param slow
#' @param shigh
#'
#' @return
#' @export
#'
#' @examples
napML<-function(y,h,m,A,mod,e,s,slow=rep(-0.5,length(s)),shigh= rep(0.5,length(s))){
  parstart<-list(
                  "s"=s,
                  "b"=0.1,
                  "a"=0.1,
                  "p"=0.1#,
                  # "mu"=1,
                  # "epi"=1
                  ) %>% unlist
  parlow<- list(
                # "s"=rep(-0.98039215686274505668,length(m)),
                # "s"=rep(-0.5,length(m)),
                "s"=slow,
                "b"=0.01,
                "a"=0.01,
                "p"=0.0#,
                # "mu"=0.99,
                # "epi"=0.8
                ) %>% unlist
  parhigh<- list(
                # "s"=rep(10,length(m)),
                # "s"=rep(0.5,length(m)),
                "s"=shigh,
                "b"=1,
                "a"=1,
                "p"=0.2#,
                # "mu"=1.02,
                # "epi"=1.2
                ) %>% unlist
  optim(fn= lik.nap,
        y=y,
        h=h,
        m=m,
        A=A,
        mod=mod,
        e=e,
        par = parstart,
        lower = parlow,
        upper=parhigh,
        control=list(pgtol=0,maxit=500,factr=1e8),
        method = "L-BFGS-B"
        )
}
