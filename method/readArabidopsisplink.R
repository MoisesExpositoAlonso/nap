library(Rcpp)
library(bigmemory)
library(devtools)
load_all('.')

####************************************************************************####
#### reading genome ####
# g<-readplink(g012file = "databig/515g.012",
#           backingfile = "genomes",
#           backingpath = "databig/")

x<-attach.big.matrix("databig/genomes.desc")
dim(x)
x[1:10,1:10]
## Recode if necessary
BMrecode(x@address,2,1)
x[1:10,1:10]
