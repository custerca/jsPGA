library(tidyverse)
library(fields)
library(Matrix)
library(lubridate)
library(stringr)
library(data.table)

# Set appropriate working directory
# setwd()

# Thermal performance curve function 
TRC = function(temp, CTmax, Topt, sigma){
  
  trc = case_when(
    temp <= Topt ~ exp(-((temp - Topt)/(2*sigma))^2), 
    CTmax >= temp & temp > Topt ~ 1 - ((temp - Topt)/(Topt - CTmax))^2,
    temp > CTmax ~ 0)
  
  return(trc)
  
}

# Reading in full fish data 
fishfull <- readRDS("data/MNfishz.rds")

# selecting names of environmental covariates
covs <- names(select(fishfull,ends_with(".z"))) %>% `[`(.!= "total.for.z")

# Reading in additional MN lake data information for coordinates of unsampled lakes, used for creating spatial basis vectors
MNlakeraw <- read.csv("data/mn_lakeinformation_allopenwat/mn_lakeinformation_1acplus.csv")
MNlakes <- MNlakeraw %>%
  filter(!is.na(dowlknum)) %>%
  mutate(DOW=paste0("mndow_",dowlknum)) %>%
  select(DOW,NHD_ID,lon=Lon,lat=Lat) %>%
  distinct()

# Correcting an incorrect coordinate
# RAW - lat: 34.55589, lon: -43.99598
# Replacing with coordinates from fishfull (from LAGOS)
# unique(fishfull[fishfull$DOW=="mndow_11030500"]$lon)
# unique(fishfull[fishfull$DOW=="mndow_11030500"]$lat)

MNlakes[MNlakes$DOW=="mndow_11030500",]$lon <- -94.35107
MNlakes[MNlakes$DOW=="mndow_11030500",]$lat <- 46.44632

# For consistency, filtering to lakes included in lake coordinate data, we use that data frame for spatial basis vectors
fishMN <- fishfull[,`:=`(GN=ifelse(GEAR=="GN",1,0),TN=ifelse(GEAR=="TN",1,0))][DOW %in% MNlakes$DOW]

# Number of species
n_species <- 5
set.seed(1)
# just arbitrarily selecting species
specsamp <- sample(c(unique(fishMN$COMMON_NAME)),size=n_species,replace=FALSE)

# randomly sampling n lakes
n_lakes <- 500
dowsamp <- sample(unique(fishMN$DOW),size = n_lakes,replace=FALSE)

# filtering down to desired number of species, changing data types
fishraw <- fishMN[COMMON_NAME %in% specsamp
                    ][DOW %in% dowsamp
                      ][,`:=`(DOW=factor(DOW),
                              COMMON_NAME=factor(COMMON_NAME))]

# removing catch data for simulation
fishraw$TOTAL_CATCH <- NULL

# the most inefficient way of changing back to not using real species names for simulation
fishraw$COMMON_NAME <- as.factor(paste0("Species",as.numeric(fishraw$COMMON_NAME)))

setorder(fishraw,DOW,SURVEYDATE,COMMON_NAME,GN)

MNlakes <- filter(MNlakes,DOW %in% dowsamp)
# Begin simulating data ---------------------------------------------------------------

dat = list()

# Number of gear types
n_gears <- 2

# Number of predictor variables
n_pred <- length(covs)

# Randomly generating parameter values for each species
set.seed(1)
beta_0 <- runif(n_species,-1,1)

beta1 <- matrix(ncol=n_pred,nrow=n_species)
set.seed(1)
for(a in 1:n_pred){
  beta1[,a] <- runif(n_species,-2,2)
}

dat$K <- n_species
dat$levs <- levels(fishraw$COMMON_NAME)
dat$lake_index <- (fishraw %>% filter(COMMON_NAME == dat$levs[1]))$DOW
dat$lake_id <- levels(fishraw$DOW)
dat$n_lakes <- length(levels(fishraw$DOW))
# making A ----------------------------------------------------------------
# simulating species variances
set.seed(1)
tau <- runif(5,0.5,2)
sr_tau <- diag(sqrt(tau))

# simulating species correlations
set.seed(1)
spec_corr <- cov2cor(LaplacesDemon::rinvwishart(n_species+10, diag(rep(1, n_species))))

# simulating A
A <- matrix(NA,
            nrow=dat$n_lakes,
            ncol=n_species)

set.seed(1)
A <- rmvnorm(dat$n_lakes,rep(0,n_species),sr_tau %*% spec_corr %*% sr_tau)


# create data from data ---------------------------------------------------

dat$E <- fishraw %>% 
  filter(EFFORT != 0) %>%
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  select(DOW, SURVEYDATE, COMMON_NAME, GN, TN, EFFORT) %>% 
  pivot_wider(names_from = COMMON_NAME, values_from = EFFORT) %>% 
  mutate_at(dat$levs[1], list(G1 = ~ TN *., G2 = ~ GN *.)) %>% 
  select(G1, G2) %>% 
  as.matrix()

dat$X <- fishraw %>% 
  filter(EFFORT != 0) %>%
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  filter(COMMON_NAME == dat$levs[1]) %>% 
  select(all_of(covs)) %>%
  as.matrix()

dat$covs <- fishraw %>% 
  filter(EFFORT != 0) %>%
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  filter(COMMON_NAME == dat$levs[1]) %>% 
  select(all_of(covs)) %>%
  names()

dat$temp <- fishraw %>% 
  filter(EFFORT != 0) %>% 
  filter(COMMON_NAME == dat$levs[1]) %>% 
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  select(july5yr) %>% 
  pull()

dat$J <- n_gears # number of gears, hard coded for now
dat$P <- dim(dat$X)[2] # number of variables

dat$alpha <- rep(1, dat$J) # Dirichlet prior parameter

##### create spatial basis functions using svd

locs <- MNlakes %>%
  arrange(DOW) %>%
  distinct(DOW, .keep_all = T) %>% 
  select(lon, lat) %>% 
  as.matrix()

dmat <- rdist.earth(locs, miles = F)
diag(dmat) <- 0

C = fields::Matern(dmat,range = 100,smoothness = 2.5)

dcomp <- svd(C)
saveRDS(dcomp,"data/dcomp_sim500.rds")
# dcomp <- readRDS("data/dcomp_sim500.rds")

dow_order <- MNlakes %>%
  arrange(DOW) %>%
  distinct(DOW, .keep_all = T) %>% 
  select(DOW)

decomp_df <- as_tibble(dcomp$u %*% diag(sqrt(dcomp$d)), 
                       .name_repair = ~paste0("comp", 1:dat$n_lakes)) %>% 
  mutate(lon = locs[,1], 
         lat = locs[,2],
         DOW = dow_order$DOW)


Upsi <- fishraw %>% 
  filter(EFFORT != 0) %>%
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  filter(COMMON_NAME == dat$levs[1]) %>% 
  left_join(decomp_df, by = "DOW") %>% 
  select(starts_with("comp")) %>% 
  as.matrix() %>% 
  unname()

dat$Upsi <- Upsi

##### simulating total catch
nobs = dim(dat$E)[1]
Etilde = matrix(NA, nrow = nobs, ncol = n_species)
lambda = matrix(NA, nrow = nobs, ncol = n_species)
trc = matrix(NA, nrow = nobs, ncol = n_species)
Y = matrix(NA, nrow = nobs, ncol = n_species)

set.seed(1)
theta = runif(n_species, 0, 1)
theta = cbind(theta,1-theta)

# effort
for(k in 1:n_species){
  Etilde[,k] = dat$E %*% theta[k,] * n_gears
}

### TRC
ctmax_prior <- c()
topt_prior <-  c()
j=1
set.seed(1)
while(j <= n_species){
  ctmp <- runif(1,26,40) #mean(ctmax_lit)
  toptp <- runif(1,15,30) #mean(topt_lit)
  
  if((ctmp-toptp) < 1) next
  
  ctmax_prior[j] <- ctmp
  topt_prior[j] <- toptp
  
  j=j+1L
  
}

# simulating uncertainty in literature derived values
set.seed(1)
ctmax_sd <- runif(n_species,1,4) #sd(ctmax_lit)
set.seed(1)
topt_sd <- runif(n_species,1,4) #sd(topt_lit)

# CTmin to derive sigma
set.seed(1)
mean_ctmin <- runif(n_species,0,3) #mean(ctmin_lit)

# sigma: Derived sd parameter given topt and ctmin (sd = (topt-tmin)/4)
sigma_lit <- (topt_prior - mean_ctmin) / 4


# TRC
for(k in 1:n_species){
  trc[,k] = TRC(dat$temp, ctmax_prior[k], topt_prior[k], sigma_lit[k])
}

# intercept
B0 = matrix(rep(beta_0, nobs), ncol = n_species, byrow = T)

# lambda
lambda = exp(B0 + (dat$X %*% t(beta1)) + (Upsi %*% A) + log(trc))

# simulate data
# arbitrary value for phi (overdispersion parameter)
phi = 3
mu = Etilde * lambda
set.seed(1)
for(k in 1:n_species){
  Y[,k] = rnbinom(nobs, size = phi, mu = mu[,k])
}

dat$Y <- Y
dat$y_vec = c(dat$Y)
dat$N = dim(dat$Y)[1] # number of observations


# Adding max temp where each species was observed
colnames(Y) <- dat$levs
mxtmp <- as_tibble(Y) %>%
  mutate(temp=dat$temp) %>%
  pivot_longer(-temp,names_to = "COMMON_NAME",values_to = "TOTAL_CATCH") %>%
  filter(TOTAL_CATCH > 0) %>%
  group_by(COMMON_NAME) %>%
  summarise(mxtmp=max(temp)) %>%
  select(mxtmp) %>%
  pull()


# Number of iterations
nsim <- 500
ctmax_mat <- matrix(NA,nrow=nsim,ncol=n_species)
topt_mat <- matrix(NA,nrow=nsim,ncol=n_species)

ctmax_mat[1,] <- ctmax_prior
topt_mat[1,] <- topt_prior


# loop
j=2L

set.seed(1)
while(j <= nsim){
  
  ctmax_vec <- rnorm(dat$K, mean = ctmax_prior, sd = ctmax_sd)
  topt_vec <- rnorm(dat$K, mean = topt_prior, sd = topt_sd)
  
  if(any(ctmax_vec - topt_vec < 1)) next
  # Testing how it looks if I make sure the fixed CTmax is never less than the observed max
  if(any(ctmax_vec < mxtmp)) next
  
  #if(sum(round(ctmax_vec) == round(topt_vec)) > 0) next
  
  ctmax_mat[j,] <- ctmax_vec
  topt_mat[j,] <- topt_vec
  j=j+1L
  
}

dat$ctmax <- ctmax_mat
dat$topt <- topt_mat
dat$sigma <- sigma_lit

# saving true values
dat$REALbeta0 <- beta_0
dat$REALbeta1 <- beta1
dat$REALspeccorr <- spec_corr
dat$REALA <- A
dat$REALtau <- tau
dat$REALtheta <- theta
dat$REALctmax <- ctmax_prior
dat$REALtopt <- topt_prior
dat$lambda <- lambda
dat$therm_dist <- tibble(ctmax_mean=ctmax_prior,
                         ctmax_sd=ctmax_sd,
                         topt_mean=topt_prior,
                         topt_sd=topt_sd)


saveRDS(dat,file=paste0("data/simdat.rds"))
