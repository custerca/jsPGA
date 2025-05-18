library(rstan)
library(tidyverse)
library(kableExtra)
library(scoringRules)


##########################################################################################################################
# This needs to be run for each basis vector count group of model realizations

# folder name containing all saved RDS files containing stan output
M=32
# vector of all file names within folder
x <- list.files(paste0("output/sim/M",M))
# grabbing model realization number of all length(x) files
simnum <- sapply(strsplit(str_remove(x,".rds"),"_"),"[[",2)

# reading in each of the RDS files
out_list <- list()
for(i in 1:length(x)){
  name <- paste0("iter",simnum[i])
  out_list[[name]] <- readRDS(paste0("output/sim/M",M,"/",x[i]))
  print(i)
}
# merging all stan objects (each model realization posterior) into single object (single posterior for analysis)
out <- sflist2stanfit(out_list)
rm(out_list)

# saving single stan object containing all model psoteriors
saveRDS(out,file=paste0("output/sim/simpost_M",M,".rds"))
##########################################################################################################################
# Calculating precision and accuracy measures for parameter estiamtes
# naming for each of the basis vector groups of models
modnm <- c("M16","M32","M64","M128")

# quick function to summarize parameters from posterior chains
chainfun <- function(x){ 
  
  my975 <- function(a) quantile(a,0.975)
  my025 <- function(a) quantile(a,0.025)
  tmpf <- function(y,func){
    switch(as.character(length(dim(y))),
           "1" = eval(call(func,y)),
           "2" = eval(call("apply",y,2,func)),
           "3" = eval(call("apply",y,c(2,3),func)))
  }
  return(list(
    avg = tmpf(x,func="mean"),
    median = tmpf(x,func="median"),
    sd = tmpf(x,func="sd"),
    u95 = tmpf(x,func="my975"),
    l95 = tmpf(x,func="my025")
  )
  )
}

# read in each merged stan object (from code above)
out16 <- readRDS("output/sim/simpost_M16.rds")
out32 <- readRDS("output/sim/simpost_M32.rds")
out64 <- readRDS("output/sim/simpost_M64.rds")
out128 <- readRDS("output/sim/simpost_M128.rds")

# read in simulated "TRUE" data
dat <- readRDS("data/simdat.rds")
n_species <- length(dat$levs)
n_covs <- dat$P

################## extracting chains and selecting only parameters of interest ##################

####### M=16
chains16 <- rstan::extract(out16,permuted=TRUE,
                           pars=c("beta_0","beta","theta","phi","tau","A","Sigma_A"))

# beta_0
dimnames(chains16$beta_0)[[2]] <- dat$levs
# beta
dimnames(chains16$beta)[[2]] <- dat$covs
dimnames(chains16$beta)[[3]] <- dat$levs
# theta
dimnames(chains16$theta)[[2]] <- dat$levs
dimnames(chains16$theta)[[3]] <- c("TN","GN") #Preliminary guess...need to confirm
# A
dimnames(chains16$A)[[2]] <- paste0("A",1:16)
dimnames(chains16$A)[[3]] <- dat$levs
#Sigma
dimnames(chains16$Sigma_A)[[2]] <- dat$levs
dimnames(chains16$Sigma_A)[[3]] <- dat$levs
# tau
dimnames(chains16$tau)[[2]] <- dat$levs

# Combining betas
chains16$allbeta = array(dim = dim(chains16$beta) + c(0,1,0))
for(i in 1:n_species){
  chains16$allbeta[,,i] <- cbind(beta0=chains16$beta_0[,i],chains16$beta[,,i])
}
dimnames(chains16$allbeta)[[2]] <- c("beta0",dat$covs)
dimnames(chains16$allbeta)[[3]] <- dat$levs


# Full covariance matrix of A: TSigmaT
chains16$TST = array(dim=dim(chains16$Sigma_A))
for(j in 1:dim(chains16$A)[[1]]){
  chains16$TST[j,,] <- diag(chains16$tau[j,]) %*% chains16$Sigma_A[j,,] %*% diag(chains16$tau[j,])
}
dimnames(chains16$TST)[[2]] <- dat$levs
dimnames(chains16$TST)[[3]] <- dat$levs


####### M=32
chains32 <- rstan::extract(out32,permuted=TRUE,
                           pars=c("beta_0","beta","theta","phi","tau","A","Sigma_A"))

# beta_0
dimnames(chains32$beta_0)[[2]] <- dat$levs
# beta
dimnames(chains32$beta)[[2]] <- dat$covs
dimnames(chains32$beta)[[3]] <- dat$levs
# theta
dimnames(chains32$theta)[[2]] <- dat$levs
dimnames(chains32$theta)[[3]] <- c("TN","GN") #Preliminary guess...need to confirm
# A
dimnames(chains32$A)[[2]] <- paste0("A",1:32)
dimnames(chains32$A)[[3]] <- dat$levs
#Sigma
dimnames(chains32$Sigma_A)[[2]] <- dat$levs
dimnames(chains32$Sigma_A)[[3]] <- dat$levs
# tau
dimnames(chains32$tau)[[2]] <- dat$levs

# Combining betas
chains32$allbeta = array(dim = dim(chains32$beta) + c(0,1,0))
for(i in 1:n_species){
  chains32$allbeta[,,i] <- cbind(beta0=chains32$beta_0[,i],chains32$beta[,,i])
}
dimnames(chains32$allbeta)[[2]] <- c("beta0",dat$covs)
dimnames(chains32$allbeta)[[3]] <- dat$levs


# Full covariance matrix of A: TSigmaT
chains32$TST = array(dim=dim(chains32$Sigma_A))
for(j in 1:dim(chains32$A)[[1]]){
  chains32$TST[j,,] <- diag(chains32$tau[j,]) %*% chains32$Sigma_A[j,,] %*% diag(chains32$tau[j,])
}
dimnames(chains32$TST)[[2]] <- dat$levs
dimnames(chains32$TST)[[3]] <- dat$levs





####### M=64
chains64 <- rstan::extract(out64,permuted=TRUE,
                           pars=c("beta_0","beta","theta","phi","tau","A","Sigma_A"))
# beta_0
dimnames(chains64$beta_0)[[2]] <- dat$levs
# beta
dimnames(chains64$beta)[[2]] <- dat$covs
dimnames(chains64$beta)[[3]] <- dat$levs
# theta
dimnames(chains64$theta)[[2]] <- dat$levs
dimnames(chains64$theta)[[3]] <- c("TN","GN") #Preliminary guess...need to confirm
# A
dimnames(chains64$A)[[2]] <- paste0("A",1:64)
dimnames(chains64$A)[[3]] <- dat$levs
#Sigma
dimnames(chains64$Sigma_A)[[2]] <- dat$levs
dimnames(chains64$Sigma_A)[[3]] <- dat$levs
# tau
dimnames(chains64$tau)[[2]] <- dat$levs

# Combining betas
chains64$allbeta = array(dim = dim(chains64$beta) + c(0,1,0))
for(i in 1:n_species){
  chains64$allbeta[,,i] <- cbind(beta0=chains64$beta_0[,i],chains64$beta[,,i])
}
dimnames(chains64$allbeta)[[2]] <- c("beta0",dat$covs)
dimnames(chains64$allbeta)[[3]] <- dat$levs


# Full covariance matrix of A: TSigmaT
chains64$TST = array(dim=dim(chains64$Sigma_A))
for(j in 1:dim(chains64$A)[[1]]){
  chains64$TST[j,,] <- diag(chains64$tau[j,]) %*% chains64$Sigma_A[j,,] %*% diag(chains64$tau[j,])
}
dimnames(chains64$TST)[[2]] <- dat$levs
dimnames(chains64$TST)[[3]] <- dat$levs



####### M=128
chains128 <- rstan::extract(out128,permuted=TRUE,
                            pars=c("beta_0","beta","theta","phi","tau","A","Sigma_A"))
# beta_0
dimnames(chains128$beta_0)[[2]] <- dat$levs
# beta
dimnames(chains128$beta)[[2]] <- dat$covs
dimnames(chains128$beta)[[3]] <- dat$levs
# theta
dimnames(chains128$theta)[[2]] <- dat$levs
dimnames(chains128$theta)[[3]] <- c("TN","GN") #Preliminary guess...need to confirm
# A
dimnames(chains128$A)[[2]] <- paste0("A",1:128)
dimnames(chains128$A)[[3]] <- dat$levs
#Sigma
dimnames(chains128$Sigma_A)[[2]] <- dat$levs
dimnames(chains128$Sigma_A)[[3]] <- dat$levs
# tau
dimnames(chains128$tau)[[2]] <- dat$levs

# Combining betas
chains128$allbeta = array(dim = dim(chains128$beta) + c(0,1,0))
for(i in 1:n_species){
  chains128$allbeta[,,i] <- cbind(beta0=chains128$beta_0[,i],chains128$beta[,,i])
}
dimnames(chains128$allbeta)[[2]] <- c("beta0",dat$covs)
dimnames(chains128$allbeta)[[3]] <- dat$levs


# Full covariance matrix of A: TSigmaT
chains128$TST = array(dim=dim(chains128$Sigma_A))
for(j in 1:dim(chains128$A)[[1]]){
  chains128$TST[j,,] <- diag(chains128$tau[j,]) %*% chains128$Sigma_A[j,,] %*% diag(chains128$tau[j,])
}
dimnames(chains128$TST)[[2]] <- dat$levs
dimnames(chains128$TST)[[3]] <- dat$levs



# applying custom summary function to each chain
chaintab16 <- lapply(chains16,chainfun)
chaintab32 <- lapply(chains32,chainfun)
chaintab64 <- lapply(chains64,chainfun)
chaintab128 <- lapply(chains128,chainfun)

# Avg time to complete model (from seconds to hours)
mean(rowSums(get_elapsed_time(out16)))/(60^2)
mean(rowSums(get_elapsed_time(out32)))/(60^2)
mean(rowSums(get_elapsed_time(out64)))/(60^2)
mean(rowSums(get_elapsed_time(out128)))/(60^2)

# large objects no longer needed
rm(out32,out64,out16,out128)

# Model summaries for each parameter
parms <- c("beta","theta","w","TST")
# MSE
bigtab <- matrix(nrow=length(parms),ncol=length(modnm))
colnames(bigtab) <- modnm
rownames(bigtab) <- parms

# CRPS
bigtabCRPS <- matrix(nrow=length(parms),ncol=length(modnm))
colnames(bigtabCRPS) <- modnm
rownames(bigtabCRPS) <- parms


chain_list <- list(M16=chains16,M32=chains32,M64=chains64,M128=chains128)
####################################################### beta #######################################################
realBeta <- t(cbind(dat$REALbeta0,dat$REALbeta1))


betaRMSE16=array(dim=c(6,5))
betaCRPS16=array(dim=c(6,5))

betaRMSE32=array(dim=c(6,5))
betaCRPS32=array(dim=c(6,5))

betaRMSE64=array(dim=c(6,5))
betaCRPS64=array(dim=c(6,5))

betaRMSE128=array(dim=c(6,5))
betaCRPS128=array(dim=c(6,5))

for(i in 1:6){
  for(j in 1:5){
    betaRMSE16[i,j] = sqrt(mean(((chains16$allbeta[,i,j] - realBeta[i,j])^2)))
    betaCRPS16[i,j] = crps_sample(y=realBeta[i,j],dat=chains16$allbeta[,i,j])
    
    betaRMSE32[i,j] = sqrt(mean(((chains32$allbeta[,i,j] - realBeta[i,j])^2)))
    betaCRPS32[i,j] = crps_sample(y=realBeta[i,j],dat=chains32$allbeta[,i,j])
    
    betaRMSE64[i,j] = sqrt(mean(((chains64$allbeta[,i,j] - realBeta[i,j])^2)))
    betaCRPS64[i,j] = crps_sample(y=realBeta[i,j],dat=chains64$allbeta[,i,j])
    
    betaRMSE128[i,j] = sqrt(mean(((chains128$allbeta[,i,j] - realBeta[i,j])^2)))
    betaCRPS128[i,j] = crps_sample(y=realBeta[i,j],dat=chains128$allbeta[,i,j])
  }
}

bigtab[1,] <- c(mean(betaRMSE16),mean(betaRMSE32),mean(betaRMSE64),mean(betaRMSE128))
bigtabCRPS[1,] <- c(mean(betaCRPS16),mean(betaCRPS32),mean(betaCRPS64),mean(betaCRPS128))


############################################# theta #######################################################
realtheta = dat$REALtheta[,1]


thetaRMSE16=c()
thetaCRPS16=c()

thetaRMSE32=c()
thetaCRPS32=c()

thetaRMSE64=c()
thetaCRPS64=c()

thetaRMSE128=c()
thetaCRPS128=c()

for(j in 1:5){
  thetaRMSE16[j] = sqrt(mean(((chains16$theta[,j,1] - realtheta[j])^2)))
  thetaCRPS16[j] = crps_sample(y=realtheta[j],dat=chains16$theta[,j,1])
  
  thetaRMSE32[j] = sqrt(mean(((chains32$theta[,j,1] - realtheta[j])^2)))
  thetaCRPS32[j] = crps_sample(y=realtheta[j],dat=chains32$theta[,j,1])
  
  thetaRMSE64[j] = sqrt(mean(((chains64$theta[,j,1] - realtheta[j])^2)))
  thetaCRPS64[j] = crps_sample(y=realtheta[j],dat=chains64$theta[,j,1])
  
  thetaRMSE128[j] = sqrt(mean(((chains128$theta[,j,1] - realtheta[j])^2)))
  thetaCRPS128[j] = crps_sample(y=realtheta[j],dat=chains128$theta[,j,1])
}



bigtab[2,] <- c(mean(thetaRMSE16),mean(thetaRMSE32),mean(thetaRMSE64),mean(thetaRMSE128))
bigtabCRPS[2,] <- c(mean(thetaCRPS16),mean(thetaCRPS32),mean(thetaCRPS64),mean(thetaCRPS128))

############################################### w (Psi * A) #######################################################
colnames(dat$Upsi) <- paste0("psi",1:500)
Psimat <- as_tibble(dat$Upsi) %>%
  mutate(lakeid=dat$lake_index) %>%
  select(lakeid,everything()) %>%
  distinct() %>%
  select(-lakeid) %>%
  as.matrix()

realOmega = Psimat %*% dat$REALA


w16 = array(dim=c(20000,500,5))
for(i in 1:20000){
  w16[i,,] = Psimat[,1:16] %*% chains16$A[i,,]
}

w32 = array(dim=c(20000,500,5))
for(i in 1:20000){
  w32[i,,] = Psimat[,1:32] %*% chains32$A[i,,]
}

w64 = array(dim=c(20000,500,5))
for(i in 1:20000){
  w64[i,,] = Psimat[,1:64] %*% chains64$A[i,,]
}

w128 = array(dim=c(20000,500,5))
for(i in 1:20000){
  w128[i,,] = Psimat[,1:128] %*% chains128$A[i,,]
}


wRMSE16=array(dim=c(500,5))
wCRPS16=array(dim=c(500,5))

wRMSE32=array(dim=c(500,5))
wCRPS32=array(dim=c(500,5))

wRMSE64=array(dim=c(500,5))
wCRPS64=array(dim=c(500,5))

wRMSE128=array(dim=c(500,5))
wCRPS128=array(dim=c(500,5))

for(i in 1:500){
  for(j in 1:5){
    wRMSE16[i,j] = sqrt(mean(((w16[,i,j] - realOmega[i,j])^2)))
    wCRPS16[i,j] = crps_sample(y=realOmega[i,j],dat=w16[,i,j])
    
    wRMSE32[i,j] = sqrt(mean(((w32[,i,j] - realOmega[i,j])^2)))
    wCRPS32[i,j] = crps_sample(y=realOmega[i,j],dat=w32[,i,j])
    
    wRMSE64[i,j] = sqrt(mean(((w64[,i,j] - realOmega[i,j])^2)))
    wCRPS64[i,j] = crps_sample(y=realOmega[i,j],dat=w64[,i,j])
    
    wRMSE128[i,j] = sqrt(mean(((w128[,i,j] - realOmega[i,j])^2)))
    wCRPS128[i,j] = crps_sample(y=realOmega[i,j],dat=w128[,i,j])
  }
}

bigtab[3,] <- c(mean(wRMSE16),mean(wRMSE32),mean(wRMSE64),mean(wRMSE128))
bigtabCRPS[3,] <- c(mean(wCRPS16),mean(wCRPS32),mean(wCRPS64),mean(wCRPS128))


############################################# TST (covariance) #######################################################
realTST <- diag(sqrt(dat$REALtau)) %*% dat$REALspeccorr %*% diag(sqrt(dat$REALtau)) 


TSTRMSE16=array(dim=c(5,5))
TSTCRPS16=array(dim=c(5,5))

TSTRMSE32=array(dim=c(5,5))
TSTCRPS32=array(dim=c(5,5))

TSTRMSE64=array(dim=c(5,5))
TSTCRPS64=array(dim=c(5,5))

TSTRMSE128=array(dim=c(5,5))
TSTCRPS128=array(dim=c(5,5))

for(i in 1:5){
  for(j in 1:5){
    TSTRMSE16[i,j] = sqrt(mean(((chains16$TST[,i,j] - realTST[i,j])^2)))
    TSTCRPS16[i,j] = crps_sample(y=realOmega[i,j],dat=w16[,i,j])
    
    TSTRMSE32[i,j] = sqrt(mean(((chains32$TST[,i,j] - realTST[i,j])^2)))
    TSTCRPS32[i,j] = crps_sample(y=realOmega[i,j],dat=w32[,i,j])
    
    TSTRMSE64[i,j] = sqrt(mean(((chains64$TST[,i,j] - realTST[i,j])^2)))
    TSTCRPS64[i,j] = crps_sample(y=realOmega[i,j],dat=w64[,i,j])
    
    TSTRMSE128[i,j] = sqrt(mean(((chains128$TST[,i,j] - realTST[i,j])^2)))
    TSTCRPS128[i,j] = crps_sample(y=realOmega[i,j],dat=w128[,i,j])

  }
}

bigtab[4,] <- c(mean(TSTRMSE16),mean(TSTRMSE32),mean(TSTRMSE64),mean(TSTRMSE128))
bigtabCRPS[4,] <- c(mean(TSTCRPS16),mean(TSTCRPS32),mean(TSTCRPS64),mean(TSTCRPS128))

tab_ltx <- bind_rows(mutate(as_tibble(bigtab,rownames="par"),st="RMSE"),
                     mutate(as_tibble(bigtabCRPS,rownames="par"),st="CRPS"))  %>%
  mutate(par = factor(par,levels=parms)) %>%
  arrange(par) %>%
  select(st,M16,M32,M64,M128,par) 

# Manually coding "scientific notation" latex table for manuscript
# beta RMSE
tab_ltx[1,2:5] <- tab_ltx[1,2:5]*10
tab_ltx[1,1] <- paste0(tab_ltx[1,1],"*")
# beta CRPS
tab_ltx[2,2:5] <- tab_ltx[2,2:5]*10
tab_ltx[2,1] <- paste0(tab_ltx[2,1],"*")
# theta RMSE
tab_ltx[3,2:5] <- tab_ltx[3,2:5]*1000
tab_ltx[3,1] <- paste0(tab_ltx[3,1],"***")
# theta CRPS
tab_ltx[4,2:5] <- tab_ltx[4,2:5]*1000
tab_ltx[4,1] <- paste0(tab_ltx[4,1],"***")
# w RMSE
tab_ltx[5,2:5] <- tab_ltx[5,2:5]*10
tab_ltx[5,1] <- paste0(tab_ltx[5,1],"*")
# w CRPS
tab_ltx[6,2:5] <- tab_ltx[6,2:5]*10
tab_ltx[6,1] <- paste0(tab_ltx[6,1],"*")
# TST RMSE
tab_ltx[7,2:5] <- tab_ltx[7,2:5]*10
tab_ltx[7,1] <- paste0(tab_ltx[7,1],"*")
# TST CRPS
tab_ltx[8,2:5] <- tab_ltx[8,2:5]*10
tab_ltx[8,1] <- paste0(tab_ltx[8,1],"*")

tab_ltx %>%
  select(-par) %>%
  rowwise() %>%
  mutate(across(.cols=c(M16,M32,M64,M128),
                .fns=function(x) sprintf("%.3f",x))) %>%
  knitr::kable(format = "latex", linesep="",booktabs=TRUE) %>%
  pack_rows(index=c("$\\boldsymbol{\\\\beta}$"=2,
                    "$\\boldsymbol{\\\\theta}$"=2,
                    "$\\boldsymbol{\\\\omega}$"=2,
                    "$\\boldsymbol{\\\\Omega}$"=2),
            escape = FALSE,bold = "FALSE")







