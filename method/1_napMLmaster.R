################################################################################
## Find BSLMM GEMMA runs and send napML runs
pathdata<-"databig/"
exname<-"example"
pname<-list.files(path = "output/",pattern="param",full=T)
pname<-pname[which(grepl("epi1.2",pname))]
pname

bname<-gsub(pattern = ".param.txt",replacement = "",basename(pname))
fname<-paste0(pathdata,bname,"/",exname,".fam")
gname<-paste0(pathdata,exname,".desc")

## Parameters
es<-c(0.8,1,1.2)
mods<-c(1,2)
snps<-c(500,2000,1000)
frac<-c(1,1,0.75)
snpfrac<-list(c(snps[1],frac[1]),
                c(snps[2],frac[2]),
                c(snps[3],frac[3])
                )

#### Send runs
# background<-" "
background<-" &"

message("Calling NAP runs")

for(i in 1){
  for(e in es[1]){
    for(mod in mods[1]){
      for(loci in snpfrac[1]){

      message("par name: ", pname[i])
      message("fam name: ", fname[i])
      message("g name: ", gname[i])
      message("mode: ", mod)
      message("epistasis: ", e)
      message("loci: ",loci)
      command<-
           paste("nice -n 19 Rscript method/1_napMLcall.R",
                 "--p", pname[i],
                 "--f", fname[i],
                 "--g", gname[i],
                 "--e", e,
                 "--m", mod,
                 "--l", loci[1],
                 "--q", loci[2],
                  background
                 )
      message("command: ",command)
      system(command)
      }
    }
  }
}
