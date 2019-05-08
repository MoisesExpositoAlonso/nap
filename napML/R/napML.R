###### R implementation
likelihoodR<-function(y, w, b,a, p){
  # mymin<-.Machine$double.xmin / (length(y)*10)
  mymin<-.Machine$double.xmin
  logmymin<-log(mymin)
  tmp<-ifelse(y==0,
              p  + (1-p) *  pnorm(0,w,a+w*b,TRUE,FALSE),
              (1-p) * dnorm(y,w,a+w*b,FALSE)
             )
  tmp[is.infinite(tmp)]<-mymin
  LL<-sum(log(tmp))
  if(is.infinite(LL)) LL<-logmymin
  return(LL)
}

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
lik.nap<-function(y,h_,m_,A,par,mod,e,debug=F){
  h=h_
  m=m_
  w=wCBM(A = A@address,
            s=par[1:length(m)],
            mycols =m,
            myrows =h,
            mode=mod)
  LL<-likelihoodR(y = y,
            w = w,
            b=par[length(m)+1],
            a=par[length(m)+2],
            p=par[length(m)+3]
            )
  if(debug) print(LL)
  if(debug) print(par[seq(length(m)+1,length(par))])
  if(debug) print(par[1:5])

  return( - LL) # DO NOT FORGET MINUS
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
napML<-function(y,h,m,A,mod,e,s,slow=rep(-0.5,length(s)),shigh= rep(0.5,length(s)),debug=F){
  parstart<-list(
                  "s"=s,
                  "b"=0.5,
                  "a"=1,
                  "p"=0.5#,
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
                "p"=0.3#,
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
        control=list(
                     trace=1,
                     pgtol=0,
                     maxit=500,
                     factr=1e-12),
        method = "L-BFGS-B"
        )
}


###### R implementation
llmix<-function(y, w, b,a, p){
  tmp<-ifelse(y==0,
              p  + (1-p) *  pnorm(0,w,a+w*b,TRUE,FALSE),
              (1-p) * dnorm(y,w,a+w*b,FALSE)
             )
  return(sum(log(tmp)))
}
#
# llmix(y = ytrain,
#             w = w,
#             b=par[length(m)+1],
#             a=par[length(m)+2],
#             p=par[length(m)+3])
