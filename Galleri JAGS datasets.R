## NHS-Galleri EE: Stat models for exploratory analyses of sharing
## CCGA data only (v1)
## S Dias, May 2024

##
## make data for JAGS models
##

wd <- getwd()
print(wd)

##
library("R2jags")
library("R2WinBUGS")
library("epitools")
#packageVersion("epitools")  #  ‘0.5.10.1’

##
## CCGA data
##
## removed missing stage and adjusted totals (including correcting N for "Other"), ignore unstageable tumours
# txt file contains simplified dataset 
CCGAv1 <- read.table("CCGAdatav1.txt", header=TRUE)  # Read matrix data
CCGAv1    # Data as matrix

##
## Create data 
##
# for "ones trick" - add z[j,k] 
# 'vec1' defines sharing cancer types (1:22), 'vecn' non-sharing cancer types (23:25)
# exchangeability across stages defined for mean in 'stcm' and het in 'stcv'
# Match with model variable names

##
## Basic data: separate means and het across all stages --> "dataJ.txt"
# Data for Base Model and Model 1A
data <- list(nCa=nrow(CCGAv1),
             s=cbind(as.matrix(CCGAv1[,5:8])),
             S=cbind(as.matrix(CCGAv1[,1:4])),
             z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
             vec1=1:22, vecn=23:25,
             stcm=1:4, nm=4,       # define separate means for all stages
             stcv=1:4, nv=4        # define separate het for all stages
             ) 
bugs.data(data, data.file = "data.txt") # make BUGS data file "data.txt"
bugs2jags("data.txt", "dataJ.txt")      # "dataJ.txt" has data for analysis in JAGS

##
## Basic data: separate means but common het across all stages --> "dataJb.txt"
# Data for Model 1B
datab <- list(nCa=nrow(CCGAv1),
             s=cbind(as.matrix(CCGAv1[,5:8])),
             S=cbind(as.matrix(CCGAv1[,1:4])),
             z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
             vec1=1:22, vecn=23:25,
             stcm=1:4, nm=4,       # define separate means for all stages
             stcv=rep(1,4), nv=1   # define common het for all stages
             ) 
bugs.data(datab, data.file = "datab.txt") # make BUGS data file "datab.txt"
bugs2jags("datab.txt", "dataJb.txt")      # "dataJb.txt" has data for analysis in JAGS


##
## Data for model with varying exchangeability across stages in ves1 
## for cancer types in vec1
##
## Stage 4 exchangeable across all cancer types; other stages independent
# Data for Model 2
data4b <- list(nCa=nrow(CCGAv1),
               s=cbind(as.matrix(CCGAv1[,5:8])),
               S=cbind(as.matrix(CCGAv1[,1:4])),
               z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
               vec1=1:22, vecn=23:25, 
               ves1=4, ves2=1:3,     # exch stages for cancer types in vec1
               stcm=1:4, nm=4,       # define separate means for all stages
               stcv=1:4, nv=4        # define separate het for all stages
) 
bugs.data(data4b, data.file = "data4b.txt") # make BUGS data file "data4b.txt"
bugs2jags("data4b.txt", "dataJ4b.txt")      # "dataJ4b.txt" has data for analysis in JAGS


## Stage 4 exchangeable across all cancer types except "high sens"; other stages independent
# Data for Model 4
# define vector of high sensitivity cancer types
hig_sens <-c(2,3,7)
vec <- 1:22  # define vector of defined cancer types
vecNHS <- vec[! vec %in% hig_sens]  # define vector of cancer types to be exchangeable
print(vecNHS)

data4e <- list(nCa=nrow(CCGAv1),
               s=cbind(as.matrix(CCGAv1[,5:8])),
               S=cbind(as.matrix(CCGAv1[,1:4])),
               z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
               vec1=vecNHS, vecn=c(hig_sens, 23:25), 
               ves1=4, ves2=1:3,     # exch stages for cancer types in vec1
               stcm=1:4, nm=4,       # define separate means for all stages
               stcv=1:4, nv=4        # define separate het for all stages
) 
bugs.data(data4e, data.file = "data4e.txt") # make BUGS data file "data4e.txt"
bugs2jags("data4e.txt", "dataJ4e.txt")      # "dataJ4e.txt" has data for analysis in JAGS

##
## Low sensitivity cancers independent: bladder, kidney, thyroid (vec1/vecn adjusted)
## separate means and het across all stages
# Data for Model 3
# define vector of low sensitivity cancer types
low_sens <- c(1,5,13)
vec <- 1:22  # define vector of defined cancer types
vecNLS <- vec[! vec %in% low_sens]  # define vector of cancer types to be exchangeable
print(vecNLS)
# make data for JAGS
datae <- list(nCa=nrow(CCGAv1),
             s=cbind(as.matrix(CCGAv1[,5:8])),
             S=cbind(as.matrix(CCGAv1[,1:4])),
             z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
             vec1=vecNLS, vecn=c(low_sens,23:25), # low sens cancer types independent
             stcm=1:4, nm=4,       # define separate means for all stages
             stcv=1:4, nv=4        # define separate het for all stages
) 
bugs.data(datae, data.file = "datae.txt") # make BUGS data file "datae.txt"
bugs2jags("datae.txt", "dataJe.txt")      # "dataJe.txt" has data for analysis in JAGS

##
## Data for EX-NEX models
## (1a) mixture by cancer type and stage; separate het, beta(1,1) prior for mixture
# Data for Model 5
dataEX1a <- list(nCa=nrow(CCGAv1),
                 s=cbind(as.matrix(CCGAv1[,5:8])),
                 S=cbind(as.matrix(CCGAv1[,1:4])),
                 z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
                 aMix=matrix(rep(1,22*4), nrow = 22, ncol = 4), # parameters for beta prior for EX
                 bMix=matrix(rep(1,22*4), nrow = 22, ncol = 4), 
                 vec1=1:22, vecn=23:25,
                 stcm=1:4, nm=4,       # define separate means for all stages
                 stcv=1:4, nv=4        # define separate het for all stages
) 
bugs.data(dataEX1a, data.file = "dataEX1a.txt") # make BUGS data file "dataEX1a.txt"
bugs2jags("dataEX1a.txt", "dataJEX1a.txt")      # "dataJEX1a.txt" has data for analysis in JAGS

##
## (2a) mixture by cancer type (all stages); separate het, beta(1,1) prior for mixture
# Data for Model 6
dataEX2a <- list(nCa=nrow(CCGAv1),
                 s=cbind(as.matrix(CCGAv1[,5:8])),
                 S=cbind(as.matrix(CCGAv1[,1:4])),
                 z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
                 aMix=rep(1,22), bMix=rep(1,22), # parameters for beta prior for EX
                 vec1=1:22, vecn=23:25,
                 stcm=1:4, nm=4,       # define separate means for all stages
                 stcv=1:4, nv=4        # define separate het for all stages
) 
bugs.data(dataEX2a, data.file = "dataEX2a.txt") # make BUGS data file "dataEX2a.txt"
bugs2jags("dataEX2a.txt", "dataJEX2a.txt")      # "dataJEXb.txt" has data for analysis in JAGS

##
##
## Read Cancer Groups for class models
##
# Empirically informed Groups
EmpGroups <- read.table("EmpiricalGroups.txt", header=FALSE)

# Hubbell Dwell Groups
HubbGroups <- read.table("Hubbell Groups.txt", header=FALSE)

# Dai 5yr survival Groups
DaiGroups <- read.table("Dai Groups.txt", header=FALSE)
# make j=11 (Prostate) 'no class' as single element
DaiGroups[which(DaiGroups==5),1] <- NA

# Sasieni stage shift Groups
SGGroups <- read.table("SasieniGroups.txt", header=FALSE)

# cluster analysis based Groups
CluGroups <- read.table("ClusterGroups.txt", header=FALSE)

##
## Data for Class models, Empirical groups
# Data for Model 7
dataEmp <- list(nCa=nrow(CCGAv1),
                s=cbind(as.matrix(CCGAv1[,5:8])),
                S=cbind(as.matrix(CCGAv1[,1:4])), 
                z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
                vec1=which(!is.na(EmpGroups)), vecn=which(is.na(EmpGroups)),  # class model cancers
                D=EmpGroups[,1], ncl=max(EmpGroups[,1], na.rm = TRUE)
) 
bugs.data(dataEmp, data.file = "dataEmp.txt") # make BUGS data file "dataEmp.txt"
bugs2jags("dataEmp.txt", "dataJEmp.txt")      # "dataJEmp.txt" has data for analysis in JAGS

##
## Data for Class models, Hubbell groups
# Sensitivity analysis
dataHG <- list(nCa=nrow(CCGAv1),
               s=cbind(as.matrix(CCGAv1[,5:8])),
               S=cbind(as.matrix(CCGAv1[,1:4])), 
               z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
               vec1=which(!is.na(HubbGroups)), vecn=which(is.na(HubbGroups)),  # class model cancers
               D=HubbGroups[,1], ncl=max(HubbGroups[,1], na.rm = TRUE)
) 
bugs.data(dataHG, data.file = "dataHG.txt") # make BUGS data file "dataHG.txt"
bugs2jags("dataHG.txt", "dataJHG.txt")      # "dataJHG.txt" has data for analysis in JAGS

##
## Data for Class models, Dai groups
# Sensitivity analysis
dataDG <- list(nCa=nrow(CCGAv1),
               s=cbind(as.matrix(CCGAv1[,5:8])),
               S=cbind(as.matrix(CCGAv1[,1:4])), 
               z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
               vec1=which(!is.na(DaiGroups)), vecn=which(is.na(DaiGroups)),  # class model cancers
               D=DaiGroups[,1], ncl=max(DaiGroups[,1], na.rm = TRUE)
) 
bugs.data(dataDG, data.file = "dataDG.txt") # make BUGS data file "dataDG.txt"
bugs2jags("dataDG.txt", "dataJDG.txt")      # "dataJDG.txt" has data for analysis in JAGS

##
## Data for Class models, Sasieni groups
# Sensitivity analysis
dataSG <- list(nCa=nrow(CCGAv1),
               s=cbind(as.matrix(CCGAv1[,5:8])),
               S=cbind(as.matrix(CCGAv1[,1:4])), 
               z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
               vec1=which(!is.na(SGGroups)), vecn=which(is.na(SGGroups)),  # class model cancers
               D=SGGroups[,1], ncl=max(SGGroups[,1], na.rm = TRUE)
) 
bugs.data(dataSG, data.file = "dataSG.txt") # make BUGS data file "dataSG.txt"
bugs2jags("dataSG.txt", "dataJSG.txt")      # "dataJSG.txt" has data for analysis in JAGS

##
## Data for Class models, Cluster based groups
# Sensitivity analysis
dataClu <- list(nCa=nrow(CCGAv1),
                s=cbind(as.matrix(CCGAv1[,5:8])),
                S=cbind(as.matrix(CCGAv1[,1:4])), 
                z=matrix(rep(1,nrow(CCGAv1)*4), nrow = nrow(CCGAv1), ncol = 4), # for "ones trick"
                vec1=which(!is.na(CluGroups)), vecn=which(is.na(CluGroups)),  # class model cancers
                D=CluGroups[,1], ncl=max(CluGroups[,1], na.rm = TRUE)
) 
bugs.data(dataClu, data.file = "dataClu.txt") # make BUGS data file "dataClu.txt"
bugs2jags("dataClu.txt", "dataJClu.txt")      # "dataJClu.txt" has data for analysis in JAGS

##
## calculate point estimate and 95%CI for probabilities for sensitivities in data
dataP <- list(cbind(p = rep(NA,25), low = rep(NA,25), upp = rep(NA,25)), 
              cbind(p = rep(NA,25), low = rep(NA,25), upp = rep(NA,25)),
              cbind(p = rep(NA,25), low = rep(NA,25), upp = rep(NA,25)),
              cbind(p = rep(NA,25), low = rep(NA,25), upp = rep(NA,25))
)
for(k in 1:4){
  for(j in 1:nrow(data$S)){
    if(data$S[j,k] !=0){
      dataP[[k]][j,"p"] <- binom.exact(data$s[j,k], data$S[j,k])$p
      dataP[[k]][j,"low"] <- binom.exact(data$s[j,k], data$S[j,k])$low
      dataP[[k]][j,"upp"] <- binom.exact(data$s[j,k], data$S[j,k])$upper
    }
  }
}

print(dataP)

##
## END
save.image("Modelsv1.RData")