load_all('.')

set.seed(1)

################################################################################
#### SIMULATE GENOME MATRIX ####

totsnps=1e4
totinds=1.5e3
#### Simulate matrix filebacked
# start_time <- Sys.time()
# x<-XBMsimulate(n=1.5e3,m=1e4,force=T)
# end_time <- Sys.time()
# end_time - start_time

#### Write bim map
dmap<-data.frame(CHR=1,SNP=paste(sep="_",1,1:totsnps),
                 X=0,POS=1:totsnps,R="C",A="A")
write.table(row.names = F,col.names = F, quote = F,
            file="databig/example.map",
            dmap[,1:4] # the last columns are for a .bim file, but gets reconstructed based on map and ped using plink later
            )
# #### write 012 PED
# BMwrite012(x@address,"databig/example.012") # only needed for my own purposes
BMwritePED(x@address,"databig/example.ped.tmp") # to convert to bed
# dum<-data.frame(FID=1:totinds,IID=1:totinds,PAT=0,MAT=0,SEX=0,PHENO=-9)
# write.table(row.names = F,col.names = F, quote = F,
#             file="databig/dummy.fam",
#             dum
#             )
# system("cut -d ' ' -f-6 dummy.fam > example.ped.tmp2")
# system("paste databig/example.ped.tmp2 databig/example.ped.tmp > databig/example.ped")
# system("plink -file databig/example --make-bed")
# system("plink --noweb --file databig/example --make-bed --out databig/example")

#### Reatach
x<-attach.big.matrix("databig/example.desc")

################################################################################
#### SIMULATE ARRAY OF GENOTYPES ####

as<-c(0.01,0.5)
bs<-c(0.01,0.5)
ps<-c(0,0.2)
mu=1
svars=c(0.01,0.1)
ss=0
epi<-c(0.9,1,1.1)
FITs<-c(1,2)
replicates=1 # assuming one replicate

#### grid of simulation
mysim<-expand.grid(bs,as,ps,svars,epi,FITs)
colnames(mysim)<-c("b","a","p","svar","epi","mod")
head(mysim)
dim(mysim)
write.csv(file = "databig/simulationgrid.tsv",mysim)

#### effect SNPs
inds<-1:1.5e3 # all individuals
snps<-1:500 # only a fraction of SNPs hve effect

## generate some distributions, that I can compare later
# s_svar01<- c(exp(rnorm(500,0,0.01))-1 , rep(0,1e4 -500))
# hist(s_svar01[s_svar01!=0])
# ssaveC(s_svar01,"databig/s_svar01.txt") # all simulations will have same S but scaled to a Svar


#### Run simulations of phenotypes and create FAM
d<-matrix(ncol=nrow(mysim),nrow=length(inds) ) %>% data.frame
i=1
for( i in 1:nrow(mysim)){
  l=mysim[i,]
  s = ssimC(1:length(totsnps),fn(l["svar"]));
  # all(wC(x[],s,1,1,1) == wCBM(x@address,s,1:totsnps-1,1:totinds-1,1,1,1))
  w=wCBM(x@address,s,1:totsnps,1:totinds,
         mode=fn(l["mod"]),
         epi=fn(l["epi"])
         )
  y=sampleWC(w,
             b=fn(l["b"]),
             a=fn(l["a"]),
             p=fn(l["p"]),
             rep=1)
  d[,i]<-y
  h2<-format(var(w,na.rm = T) / var(y,na.rm = T),digits=2)
  colnames(d)[i]<- paste0( collapse="_",
                           c(paste0(  colnames(mysim), (mysim[i,])),
                           paste0("h2",h2)
                          )
                        )
}

head(d)
d[is.na(d)]<-0

#### write general fam
structure<-data.frame(FID=1:totinds,IID=1:totinds,PAT=0,MAT=0,SEX=0)
dfam<-cbind(structure,d)
dim(dfam)
head(dfam[,1:6])
write.table(row.names = F,col.names = T, quote = F,
            file="databig/simexample.fam",
            dfam
            )
#### write fams per folder
for( i in 1:ncol(d)){
  message(colnames(d)[i])
  system(paste0("mkdir databig/",colnames(d)[i]))
  write.table(row.names = F,col.names = F, quote = F,
            file=paste0("databig/",colnames(d)[i],"/example.fam"),
            cbind(structure,data.frame(PHENO=d[,i]))
            )
  system(paste0("ln databig/example.bed databig/",colnames(d)[i],"/example.bed"))
  system(paste0("ln databig/example.bim databig/",colnames(d)[i],"/example.bim"))
}

#### check things are not crazy

pdf("databig/histograms.pdf")
for( i in 1:ncol(d)){
  hist(d[,i],breaks=50,main=
         paste(paste(  colnames(mysim), (mysim[i,])), collapse="  "))
}
dev.off()


# maxes<-apply(d[,1:ncol(d)],2,function(i) {max(i)}) %>% fn
# table(maxes>50)
# toy<-mysim
# mysim$max <- maxes
# View(mysim)
