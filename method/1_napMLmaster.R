################################################################################

################################################################################
## Find BSLMM GEMMA runs and send napML runs

runs<-list.files(path = "output/",pattern="param")
runs<-gsub(pattern = ".param.txt",replacement = "",runs)
runs

## Parameters
es<-c(0.9,1,1.1)
mods<-c(1,2)
snps<-c(500,2000,1000)
frac<-c(1,1,0.75)
snpfrac<-list(c(snps[1],frac[1]),
                c(snps[2],frac[2]),
                c(snps[3],frac[3])
                )

#### Send runs
# background<-" "3*
background<-" &"

message("Calling NAP runs")
for(file in runs[1:2]){
  for(e in es){
    for(mod in mods){
      for(loci in snpfrac){
      message("file", file)
      message("mode", mod)
      message("epistasis", e)
      message(loci)
    system(paste("nice -n 19 Rscript method/1_napMLcall.R",
                 "--p", file,
                 "--e", e,
                 "--m", mod,
                 "--l", loci[1],
                 "--f", loci[2],
                  background
                 )
           )
      }
    }
  }
}
