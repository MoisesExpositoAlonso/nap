################################################################################
################################################################################
## Extract fitness averages per genotype for GWA
################################################################################
################################################################################

devtools::load_all('.')

## load data sets

data(dry)
head(data)

################################################################################
## (0) ## Calculate relative fitness
################################################################################

thevars=c('Survival_fruit','Seeds', 'Fitness')
thevars=c('Fitness')


dryr<-dry

 Survivalrelative<-
   apply(dryr[, grep('Survival_fruit', colnames(dryr))],
         2, relative)
 colnames(Survivalrelative)<-paste0("r",colnames(Survivalrelative))

 Seedrelative<-
   apply(dryr[, grep('Seeds', colnames(dryr))],
         2, relative)
 colnames(Seedrelative)<-paste0("r",colnames(Seedrelative))

Fitnessrelative<-
  apply(dryr[, grep('Fitness', colnames(dryr))],
        2,
        relative
        ) ## log

colnames(Fitnessrelative)<-paste0("r",colnames(Fitnessrelative))


dryr<-cbind(dryr, Survivalrelative, Seedrelative, Fitnessrelative)
#dryr<-cbind(dryr, Fitnessrelative)

colnames(dryr)

################################################################################
## Prepare the plink files to run the gwa
################################################################################

colnames(dry)


thevars=c("rFitness","rSeeds","rSurvival_fruit")

codes=fieldcodes()

tooutput=expand.grid(thevars, codes)
print(tooutput)

lapply(1:nrow(tooutput),
        FUN = function(r){
        row=fc(tooutput[r,])
        genplink(data = dryr,
                 phenotype = paste(row[1],row[2],sep="_"),
                 out=paste(row[1],row[2],sep="_"),
                 cleardir = TRUE )
        }
)


# lapply(fieldcodes(rep='i'),function(code) subset_2plink(data = dry,'Survival_flowering', code,dorelative = TRUE))
# lapply(fieldcodes(rep='i'),function(code) subset_2plink(data = dry,'Survival_fruit', code,dorelative = TRUE))
# lapply(fieldcodes(rep='i'),function(code) subset_2plink(data = dry,'log10Seeds', code,dorelative = TRUE))
# lapply(fieldcodes(rep='i'),function(code) subset_2plink(data = dry,'rlog10Fitness', code,dorelative=FALSE,naval=0) )
# lapply(fieldcodes(rep='i'),function(code) subset_2plink(data = dry,'Fitness', code,dorelative=FALSE,naval=0) )
  # the idea of having naval in the fitness is that it accounts for survival as well


################################################################################
## Run GWA BSLMM for a phenotype
################################################################################
## Define which want to be run

thevars<- c("rFitness")
thevars<- c("rSeeds","rSurvival_fruit")
codes=fieldcodes()

expand.grid(thevars,codes) %>%
  dplyr::mutate(combined=paste(Var1, Var2, sep="_")) %>%
  dplyr::select(combined)%>%
  as.matrix %>% c() -> torun
torun

torun<-torun[1]

## Send the runs

for(out in torun){
  run_gemma(out =out ,background=T,type = 'lm')
}

for(out in torun){
  run_gemma(out =out ,background=T,type = 'bslmm')
}


# for(out in torun){
#   run_gemma(out =out ,background=T,type = 'lmm')
# }

trialXX<-as.matrix(read.table('outputold/extinctionfreq_lmm.cXX.txt',header=F))
trialXX[1:5,1:5]
diag(trialXX)  %>% sum
trialXX[1,] %>% sum


####************************************************************************####
# Above was 2018 paper # 
# Below NAP project to compare predictions #
####************************************************************************####
### Run predictions

pred_gemma
