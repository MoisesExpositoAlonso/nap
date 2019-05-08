################################################################################
## Find BSLMM GEMMA runs and send napML runs
pathdata<-"databig/"


exname<-"example"
pname<-list.files(path = "output/",pattern="param",full=T)

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
  for(e in es){
    for(mod in mods){
      for(loci in snpfrac){

      message("par name: ", pname[i])
      message("fam name: ", fname[i])
      message("g name: ", gname)
      message("mode: ", mod)
      message("epistasis: ", e)
      message("loci: ",loci)
      command<-
           paste("nice -n 19 Rscript method/1_napMLcall.R",
                 "--p", pname[i],
                 "--f", fname[i],
                 "--g", gname,
                 "--e", e,
                 "--m", mod,
                 "--l", loci[1],
                 "--q", loci[2],
                  background
                 )
      message("command: ",command)
##      system(command)
      }
    }
  }
}


###

rdaname<-list.files(path = "outputremote/",pattern=".rda",full=T)

r<-readRDS(rdaname[4])

hist(r[[2]])
r[[1]]$value

r<-readRDS("outputremote/b0.5_a0.5_p0.2_svar0.1_epi1.2_mod1_h20.74.results.a.e1all2000.rda")
r<-readRDS("outputremote/b0.01_a0.01_p0.2_svar0.01_epi0.8_mod2_h20.12.results.a.e1all500.rda")
