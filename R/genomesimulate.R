XBMsimulate<-function(backingpath="databig",filename="example",n,m,force=F){
  # Create a filebacked big matrix in R, and fill it
  # with 2s (homozygous for minor) at a MAF simulated from runif 0-0.5
  # To read matrix, attach.big.matrix("descriptionfile.desc")
  if(force){
    system(paste("rm",paste0(backingpath,"/",filename,".bk")))
  }
  x <- filebacked.big.matrix(n, m, type='double', init= 0,
                             backingpath = backingpath,
                             backingfile=paste0(filename,".bk"),
                             descriptorfile=paste0(filename,".desc")
                            )
  x
  BMsimulate(x@address)
  return(x)
}

## Need to also produce a Map and BIM files

## Need to produce a 012, raw, and PED files
