## NHS-Galleri EE: Stat models for exploratory analyses of sharing
## CCGA data only (v1)
## S Dias, May 2024

##
## JAGS models and outputs
##

wd <- getwd()
print(wd)

##
library("R2jags")
# packageVersion("R2jags")    # check package version: 0.8.5
# jags.version()      # version: 4.3.1
library("R2WinBUGS")
# packageVersion("R2WinBUGS")    # check package version: 2.1.22.1

########################
##
## BASE MODEL
##
## Run JAGS model 3-0: p constrained to increase with stage, no sharing
## "ones trick" to constrain mu's to be increasing
## dataJ.txt --> data not used: vec1, vecn, stcm, nm, stcv, nv
set.seed(123)
SModel30v1 <- jags(data = "dataJ.txt", 
                   parameters.to.save = c("p",                          # inference
                                          "dev", "resdev", "totresdev", # model fit
                                          "mu"                          # check convergence 
                   ),
                   model.file = "SensModel3-0_v1.txt", n.chains = 2,
                   n.iter = 150000, n.burnin = 50000, n.thin = 1
                   )
##
## check convergence of mu
cbind(pos=which(SModel30v1$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel30v1$BUGSoutput$summary[which(SModel30v1$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M30v1mcmc <- as.mcmc(SModel30v1)  # make MCMC object
sum30v1 <- summary(M30v1mcmc, quantiles = c(0.025, 0.5, 0.975)) # summary stats

## save/remove large objects
rm(SModel30v1)
save(M30v1mcmc, sum30v1, file = "Model30v1.RData")
rm(M30v1mcmc, sum30v1)
#load("Model30v1.RData")


########################
##
## MODEL 1A
##
## Run JAGS model 3-3G: p constrained to increase with stage
# "ones trick" to constrain mu's to be increasing
# exchangeable across cancer types in 'vec1' and across stages defined in 'stcm' 
# variance common across stages defined in 'stcv' 
# cancer types in 'vecn' independent
## dataJ.txt --> (a) stcm=1:4, nm=4, stcv=1:4, nv=4
set.seed(123)
SModel33v1a <- jags(data = "dataJ.txt", 
                   parameters.to.save = c("p","M", "SD",                # inference
                                          "dev", "resdev", "totresdev", # model fit
                                          "mu"                          # check convergence 
                   ),
                   model.file = "SensModel3-3Gvec_v1.txt", n.chains = 2,
                   n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel33v1a$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel33v1a$BUGSoutput$summary[which(SModel33v1a$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M33v1amcmc <- as.mcmc(SModel33v1a)
sum33v1a <- summary(M33v1amcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel33v1a)
save(M33v1amcmc, sum33v1a, file = "Model33v1a.RData")
rm(M33v1amcmc, sum33v1a)
#load("Model33v1a.RData")


##########################
##
## MODEL 1B
##
## dataJb.txt --> (b) stcm=1:4, nm=4, stcv=rep(1,4), nv=1
set.seed(123)
SModel33v1b <- jags(data = "dataJb.txt", 
                    parameters.to.save = c("p","M", "SD",                # inference
                                           "dev", "resdev", "totresdev", # model fit
                                           "mu"                          # check convergence 
                    ),
                    model.file = "SensModel3-3Gvec_v1.txt", n.chains = 2,
                    n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel33v1b$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel33v1b$BUGSoutput$summary[which(SModel33v1b$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M33v1bmcmc <- as.mcmc(SModel33v1b)
sum33v1b <- summary(M33v1bmcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel33v1b)
save(M33v1bmcmc, sum33v1b, file = "Model33v1b.RData")
rm(M33v1bmcmc, sum33v1b)
#load("Model33v1b.RData")


########################
##
## MODEL 2
##
## Run JAGS model 4-1: p constrained to increase with stage
# "ones trick" to constrain mu's to be increasing
# exchangeable across Ca types in 'vec1' for stages defined in 'ves1' 
# assumptions of  common mean and variance across exch stages defined in 'stcm' and 'stcv' 
# cancer types in 'vecn' independent

##
## dataJ4b.txt --> Stage 4 exchangeable across all cancer types; other stages independent
##stcm=1:4, nm=4, stcv=1:4, nv=4
set.seed(123)
SModel41v1b <- jags(data = "dataJ4b.txt", 
                    parameters.to.save = c("p","M", "SD",                # inference
                                           "dev", "resdev", "totresdev", # model fit
                                           "mu"                          # check convergence 
                    ),
                    model.file = "SensModel4-1vec_v1.txt", n.chains = 2,
                    n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel41v1b$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel41v1b$BUGSoutput$summary[which(SModel41v1b$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M41v1bmcmc <- as.mcmc(SModel41v1b)
sum41v1b <- summary(M41v1bmcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel41v1b)
save(M41v1bmcmc, sum41v1b, file = "Model41v1b.RData")
rm(M41v1bmcmc, sum41v1b)
#load("Model41v1b.RData")


########################
##
## MODEL 4
##
## dataJ4e.txt --> (e) Stage 4 exchangeable across all cancer types except "high sens"; 
## other stages independent
##stcm=1:4, nm=4, stcv=1:4, nv=4
set.seed(123)
SModel41v1e <- jags(data = "dataJ4e.txt", 
                    parameters.to.save = c("p","M", "SD",                # inference
                                           "dev", "resdev", "totresdev", # model fit
                                           "mu"                          # check convergence 
                    ),
                    model.file = "SensModel4-1vec_v1.txt", n.chains = 2,
                    n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel41v1e$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel41v1e$BUGSoutput$summary[which(SModel41v1e$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M41v1emcmc <- as.mcmc(SModel41v1e)
sum41v1e <- summary(M41v1emcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel41v1e)
save(M41v1emcmc, sum41v1e, file = "Model41v1e.RData")
rm(M41v1emcmc, sum41v1e)
#load("Model41v1e.RData")


########################
##
## MODEL 3
##
## dataJe.txt --> (e) low sensitivity cancers independent
## stcm=1:4, nm=4, stcv=1:4, nv=4
set.seed(123)
SModel33v1e <- jags(data = "dataJe.txt", 
                    parameters.to.save = c("p","M", "SD",                # inference
                                           "dev", "resdev", "totresdev", # model fit
                                           "mu"                          # check convergence 
                    ),
                    model.file = "SensModel3-3Gvec_v1.txt", n.chains = 2,
                    n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel33v1e$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel33v1e$BUGSoutput$summary[which(SModel33v1e$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M33v1emcmc <- as.mcmc(SModel33v1e)
sum33v1e <- summary(M33v1emcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel33v1e)
save(M33v1emcmc, sum33v1e, file = "Model33v1e.RData")
rm(M33v1emcmc, sum33v1e)
#load("Model33v1e.RData")


########################
##
## MODEL 5
##
## Run JAGS model 5-1: p constrained to increase with stage
# "ones trick" to constrain mu's to be increasing
# EX-NEX model across types and stages for cancer types vec1
# Means common across stages defined in 'stcm' 
# variance common across stages defined in 'stcv' 
# cancer types in 'vecn' independent
##
## dataJEX1a.txt --> (a) mixture by cancer type and stage
# stcm=1:4, nm=4, stcv=1:4, nv=4, aMix=1, bMix=1
set.seed(123)
SModel51v1a <- jags(data = "dataJEX1a.txt", 
                    parameters.to.save = c("p","M", "SD",                # inference
                                           "pick", "pMix", "theta",      # EXNEX
                                           "dev", "resdev", "totresdev", # model fit
                                           "mu"                          # check convergence 
                                           ),
                    model.file = "SensModel5-1vec_v1.txt", n.chains = 2,
                    n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel51v1a$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel51v1a$BUGSoutput$summary[which(SModel51v1a$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M51v1amcmc <- as.mcmc(SModel51v1a)
sum51v1a <- summary(M51v1amcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel51v1a)
save(M51v1amcmc, sum51v1a, file = "Model51v1a.RData")
rm(M51v1amcmc, sum51v1a)
#load("Model51v1a.RData")


########################
##
## MODEL 6
##
## Run JAGS model 5-2: p constrained to increase with stage
# "ones trick" to constrain mu's to be increasing
# EX-NEX model across types for cancer types vec1
# Means common across stages defined in 'stcm' 
# variance common across stages defined in 'stcv' 
# cancer types in 'vecn' independent

##
## dataJEX2a.txt --> (a) mixture by cancer type only
# stcm=1:4, nm=4, stcv=1:4, nv=4, aMix=1, bMix=1
set.seed(123)
SModel52v1a <- jags(data = "dataJEX2a.txt", 
                    parameters.to.save = c("p","M", "SD",                # inference
                                           "pick", "pMix", "theta",      # EXNEX
                                           "dev", "resdev", "totresdev", # model fit
                                           "mu"                          # check convergence 
                    ),
                    model.file = "SensModel5-2vec_v1.txt", n.chains = 2,
                    n.iter = 200000, n.burnin = 100000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel52v1a$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel52v1a$BUGSoutput$summary[which(SModel52v1a$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

# ***convergence issues on mus and thetas***

M52v1amcmc <- as.mcmc(SModel52v1a)
sum52v1a <- summary(M52v1amcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel52v1a)
save(M52v1amcmc, sum52v1a, file = "Model52v1a.RData")
rm(M52v1amcmc, sum52v1a)
#load("Model52v1a.RData")


########################
##
## MODEL 7
##
## Run JAGS model 3-4G: p constrained to increase with stage  CLASS MODEL
# "ones trick" to constrain mu's to be increasing
# cancer types in 'vec1' exchangeable across Ca types, defined by "class" 
# (common heterogeneity across stages within class)
# cancer types in 'vec2' independent

##
## dataJEmp.txt --> (e) Empirically defined Groups
set.seed(123)
SModel34v1Emp <- jags(data = "dataJEmp.txt", 
                      parameters.to.save = c("p","m", "sd.m",              # inference
                                             "dev", "resdev", "totresdev", # model fit
                                             "mu"                          # check convergence 
                      ),
                      model.file = "SensModel3-4Gvec_v1.txt", n.chains = 2,
                      n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel34v1Emp$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel34v1Emp$BUGSoutput$summary[which(SModel34v1Emp$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M34v1Empmcmc <- as.mcmc(SModel34v1Emp)
sum34v1Emp <- summary(M34v1Empmcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel34v1Emp)
save(M34v1Empmcmc, sum34v1Emp, file = "Model34v1Emp.RData")
rm(M34v1Empmcmc, sum34v1Emp)
#load("Model34v1Emp.RData")

########################
## SENSITIVITY ANALYSES
########################
##
## dataJHG.txt -->  (a) Hubbell Dwell Groups
set.seed(123)
SModel34v1HG <- jags(data = "dataJHG.txt", 
                     parameters.to.save = c("p","m", "sd.m",              # inference
                                            "dev", "resdev", "totresdev", # model fit
                                            "mu"                          # check convergence 
                     ),
                     model.file = "SensModel3-4Gvec_v1.txt", n.chains = 2,
                     n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel34v1HG$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel34v1HG$BUGSoutput$summary[which(SModel34v1HG$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M34v1HGmcmc <- as.mcmc(SModel34v1HG)
sum34v1HG <- summary(M34v1HGmcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel34v1HG)
save(M34v1HGmcmc, sum34v1HG, file = "Model34v1HG.RData")
rm(M34v1HGmcmc, sum34v1HG)
#load("Model34v1HG.RData")


##
## dataJDG.txt --> (b) Dai 5 yr survival Groups
set.seed(123)
SModel34v1DG <- jags(data = "dataJDG.txt", 
                     parameters.to.save = c("p","m", "sd.m",              # inference
                                            "dev", "resdev", "totresdev", # model fit
                                            "mu"                          # check convergence 
                     ),
                     model.file = "SensModel3-4Gvec_v1.txt", n.chains = 2,
                     n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel34v1DG$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel34v1DG$BUGSoutput$summary[which(SModel34v1DG$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M34v1DGmcmc <- as.mcmc(SModel34v1DG)
sum34v1DG <- summary(M34v1DGmcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel34v1DG)
save(M34v1DGmcmc, sum34v1DG, file = "Model34v1DG.RData")
rm(M34v1DGmcmc, sum34v1DG)
#load("Model34v1DG.RData")


##
## dataJSG.txt --> (c) Sasieni Groups
set.seed(123)
SModel34v1SG <- jags(data = "dataJSG.txt", 
                     parameters.to.save = c("p","m", "sd.m",              # inference
                                            "dev", "resdev", "totresdev", # model fit
                                            "mu"                          # check convergence 
                     ),
                     model.file = "SensModel3-4Gvec_v1.txt", n.chains = 2,
                     n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel34v1SG$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel34v1SG$BUGSoutput$summary[which(SModel34v1SG$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M34v1SGmcmc <- as.mcmc(SModel34v1SG)
sum34v1SG <- summary(M34v1SGmcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel34v1SG)
save(M34v1SGmcmc, sum34v1SG, file = "Model34v1SG.RData")
rm(M34v1SGmcmc, sum34v1SG)
#load("Model34v1SG.RData")


##
## dataJSG.txt --> (d) Cluster Groups
set.seed(123)
SModel34v1Clu <- jags(data = "dataJClu.txt", 
                      parameters.to.save = c("p","m", "sd.m",              # inference
                                             "dev", "resdev", "totresdev", # model fit
                                             "mu"                          # check convergence 
                      ),
                      model.file = "SensModel3-4Gvec_v1.txt", n.chains = 2,
                      n.iter = 150000, n.burnin = 50000, n.thin = 1)
##
## check convergence of mu
cbind(pos=which(SModel34v1Clu$BUGSoutput$summary[,"Rhat"] > 1.01), 
      Rhat=SModel34v1Clu$BUGSoutput$summary[which(SModel34v1Clu$BUGSoutput$summary[,"Rhat"] > 1.01),"Rhat"])

M34v1Clumcmc <- as.mcmc(SModel34v1Clu)
sum34v1Clu <- summary(M34v1Clumcmc, quantiles = c(0.025, 0.5, 0.975))

## save/remove large objects
rm(SModel34v1Clu)
save(M34v1Clumcmc, sum34v1Clu, file = "Model34v1Clu.RData")
rm(M34v1Clumcmc, sum34v1Clu)
#load("Model34v1Clu.RData")


##
## END
save.image("Modelsv1.RData")