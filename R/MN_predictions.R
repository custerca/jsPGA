library(rstan)
library(tidyverse)
library(data.table)

#setwd()



############################################  Predicting thermal performance and relative abundance  #############################################

TRC = function(temp, CTmax, Topt, sigma){
  trc = case_when(
    temp <= Topt ~ exp(-((temp - Topt)/(2*sigma))^2), 
    CTmax >= temp & temp > Topt ~ 1 - ((temp - Topt)/(Topt - CTmax))^2,
    temp > CTmax ~ 0)
  return(trc)
}

MN_fullcc <- readRDS("data/MNfullcc.rds")

out <- readRDS("output/MN/MNpost_M16.rds")
dat <- readRDS("data/MNdat_M16.rds")

chains <- rstan::extract(out,pars=c("beta_0","beta","theta","phi","tau","A","Sigma_A"))

# beta_0
dimnames(chains$beta_0)[[2]] <- dat$levs
# beta
dimnames(chains$beta)[[2]] <- dat$covs
dimnames(chains$beta)[[3]] <- dat$levs
# theta
dimnames(chains$theta)[[2]] <- dat$levs
dimnames(chains$theta)[[3]] <- c("TN","GN") #Preliminary guess...need to confirm
# A
dimnames(chains$A)[[2]] <- paste0("A",1:dat$M)
dimnames(chains$A)[[3]] <- dat$levs
#Sigma
dimnames(chains$Sigma_A)[[2]] <- dat$levs
dimnames(chains$Sigma_A)[[3]] <- dat$levs

dow_order_fin <- MN_fullcc %>%
  select(DOW) %>%
  distinct() %>%
  arrange(DOW)

# Psi
Xmat <- MN_fullcc %>%
  mutate(b0=1) %>%
  select(DOW,b0,secchi.z,lakearea.z,total.dev.z,total.ag.z,elevation.z) %>%
  distinct() %>%
  right_join(dow_order_fin,by = join_by(DOW)) %>%
  select(-DOW) %>%
  as.matrix()

Psimat <- select(MN_fullcc,DOW,starts_with("comp")) %>%
  distinct() %>%
  right_join(dow_order_fin,by = join_by(DOW)) %>%
  select(-DOW) %>%
  as.matrix()


# number of unique lakes
n_lakes <- nrow(dow_order_fin)

# Creating beta matrix for predictions
beta_0 = chains$beta_0
dim(beta_0)
#20000 x 8

beta1 = chains$beta
dim(beta1)
#20000 x 5 x8

beta_mat = array(dim=c(20000,6,8))
beta_mat[,1,] <- beta_0
beta_mat[,2:6,] <- beta1

#number of iterations across all chains
n_samps = dim(beta1)[1]

# A - basis coefficient matrix
Achain <- chains$A
dim(Achain)
# 20000 x 16 x 8

# Create chain indicator for each mcmc iteration because chains are 'stacked' on top of one another
n_chains <- length(out@sim$samples) # number of chains ran
c_ind <- sort(rep(1:n_chains, n_samps/n_chains))

# removing larger object as no longer needed
rm(out)

t_mods <- unique(MN_fullcc$GCM)
temp_list <- list()

temp_list[["LC"]] <- lapply(t_mods,function(x) 
  MN_fullcc %>% filter(GCM==x) %>% mutate(temp=NLDAS_2000+LC) %>% arrange(DOW) %>% select(temp) %>% pull()
)


names(temp_list[["LC"]]) <- t_mods

# We only discuss late-century predictions in paper
# temp_list[["MC"]] <- lapply(t_mods,function(x) 
#   MN_fullcc %>% filter(GCM==x) %>% mutate(temp=NLDAS_2000+MC) %>% arrange(DOW) %>% select(temp) %>% pull()
# )
# 
# names(temp_list[["MC"]]) <- t_mods



# Predictions for late century
lamdf <- array(dim=c(n_samps,n_lakes,dat$K,6))
scadf <- lamdf

dimnames(lamdf)[[2]] <- dow_order_fin$DOW 
dimnames(scadf)[[2]] <- dow_order_fin$DOW
dimnames(lamdf)[[3]] <- dat$levs 
dimnames(scadf)[[3]] <- dat$levs
dimnames(lamdf)[[4]] <- t_mods
dimnames(scadf)[[4]] <- t_mods


x = Sys.time()
for(z in 1:6){
  for(i in 1:n_samps){
    # predicting values of P(T)it
    for(j in 1:dat$K){
      scadf[i,,j,z] = TRC(temp_list$LC[[z]],dat$ctmax[c_ind[i],j],dat$topt[c_ind[i],j],dat$sigma[j])
    }

    lamdf[i,,,z] = exp((Xmat %*% beta_mat[i,,]) + (Psimat %*% Achain[i,,])) * scadf[i,,,z]
    
    if(((i %% 100)==0)|i==1) print(paste0(t_mods[z],"-",i))
  }
}
Sys.time() - x

x = Sys.time()
saveRDS(list(lamdf,scadf),file="output/MNpreds_LC.rds")
Sys.time() - x


# Predictions for "Current"
tempNLDAS <- MN_fullcc %>% filter(GCM==t_mods[1]) %>% mutate(temp=NLDAS_2021) %>% arrange(DOW) %>% select(temp) %>% pull()

lamdf <- array(dim=c(n_samps,n_lakes,dat$K))
scadf <- lamdf

dimnames(lamdf)[[2]] <- dow_order_fin$DOW 
dimnames(scadf)[[2]] <- dow_order_fin$DOW
dimnames(lamdf)[[3]] <- dat$levs 
dimnames(scadf)[[3]] <- dat$levs


x = Sys.time()

  for(i in 1:n_samps){
    # predicting values of P(T)it
    for(j in 1:dat$K){
      scadf[i,,j] = TRC(tempNLDAS,dat$ctmax[c_ind[i],j],dat$topt[c_ind[i],j],dat$sigma[j])
    }
    
    lamdf[i,,] = exp((Xmat %*% beta_mat[i,,]) + (Psimat %*% Achain[i,,])) * scadf[i,,]
    
    if(((i %% 100)==0)|i==1) print(paste0("NLDAS-",i))
  }

Sys.time() - x

x = Sys.time()
saveRDS(list(lamdf,scadf),file="output/MNpreds_Current.rds")
Sys.time() - x

# Predictions for "Historical"
temp2k <- MN_fullcc %>% filter(GCM==t_mods[1]) %>% mutate(temp=NLDAS_2000) %>% arrange(DOW) %>% select(temp) %>% pull()

lamdf <- array(dim=c(n_samps,n_lakes,dat$K))
scadf <- lamdf

dimnames(lamdf)[[2]] <- dow_order_fin$DOW 
dimnames(scadf)[[2]] <- dow_order_fin$DOW
dimnames(lamdf)[[3]] <- dat$levs 
dimnames(scadf)[[3]] <- dat$levs


x = Sys.time()

for(i in 1:n_samps){
  # predicting values of P(T)it
  for(j in 1:dat$K){
    scadf[i,,j] = TRC(temp2k,dat$ctmax[c_ind[i],j],dat$topt[c_ind[i],j],dat$sigma[j])
  }
  
  lamdf[i,,] = exp((Xmat %*% beta_mat[i,,]) + (Psimat %*% Achain[i,,])) * scadf[i,,]
  
  if(((i %% 100)==0)|i==1) print(paste0("NLDAS-",i))
}

Sys.time() - x

x = Sys.time()
saveRDS(list(lamdf,scadf),file="output/MNpreds_2k.rds")
Sys.time() - x


############################################  Summarizing predictions  #############################################
MN_fullcc <- readRDS("data/MNfullcc.rds")
MN_df_LC <- MN_fullcc %>%
  select(-MC) %>%
  pivot_wider(names_from = "GCM",values_from = "LC")

# MN_df_MC <- MN_fullcc %>%
#   select(-LC) %>%
#   pivot_wider(names_from = "GCM",values_from = "MC")


# Late century
lamscaLC <- readRDS(file="output/MNpreds_LC.rds")
lam_LC <- lamscaLC[[1]]
sca_LC <- lamscaLC[[2]]
rm(lamscaLC)
dim(lam_LC)
# Averaging predicted relative abundance across all 6 GCMs for each lake
lammean_LC <- apply(lam_LC,c(2,3),mean)
# Mean predicted relative abundance for each GCM at each lake
lammean6_LC <- apply(lam_LC,c(2,3,4),mean)
saveRDS(lammean_LC,file="output/lammean_LC.rds")
saveRDS(lammean6_LC,file="output/lammean6_LC.rds")

# Average percent extinction across all 6 GCMs for each lake
pct_extinct_LC <- apply(sca_LC,c(2,3),function(x) mean(x==0))
saveRDS(pct_extinct_LC,file="output/pct_extinct_LC.rds")
pct_extinct_LC6 <- apply(sca_LC,c(2,3,4),function(x) mean(x==0))
saveRDS(pct_extinct_LC6,file="output/pct_extinct_LC6.rds")

# # Mid century
# lamscaMC <- readRDS(file="output/MNpreds_MC_2k.rds")
# lam_MC <- lamscaMC[[1]]
# sca_MC <- lamscaMC[[2]]
# rm(lamscaMC)
# dim(lam_MC)
# lammean_MC <- apply(lam_MC,c(2,3),mean)
# #MN_lammean_MC <- bind_cols(select(MN_df_MC,DOW,ends_with(".z")),lammean_MC)
# lammean6_MC <- apply(lam_MC,c(2,3,4),mean)
# saveRDS(lammean_MC,file="output/lammean_MC_2k.rds")
# saveRDS(lammean6_MC,file="output/lammean6_MC_2k.rds")
# pct_extinct_MC <- apply(sca_MC,c(2,3),function(x) mean(x==0))
# saveRDS(pct_extinct_MC,file="output/pct_extinct_MC_2k.rds")
# pct_extinct_MC6 <- apply(sca_MC,c(2,3,4),function(x) mean(x==0))
# saveRDS(pct_extinct_MC6,file="output/pct_extinct_MC6_2k.rds")
#

# Remove everything for memory
rm(list=ls())

# Current
lamscaNLDAS <- readRDS(file="output/MNpreds_Current.rds")
lam_NLDAS <- lamscaNLDAS[[1]]
sca_NLDAS <- lamscaNLDAS[[2]]
rm(lamscaNLDAS)
dim(lam_NLDAS)
lammean_NLDAS <- apply(lam_NLDAS,c(2,3),mean)
#MN_lammean_NLDAS <- bind_cols(select(MN_df_LC,DOW,ends_with(".z")),lammean_NLDAS)
saveRDS(lammean_NLDAS,file="output/lammean_Current.rds")

# Remove everything for memory
rm(list=ls())

# "Historical" NLDAS
lamsca2k <- readRDS(file="output/MNpreds_2k.rds")
lam_2k <- lamsca2k[[1]]
sca_2k <- lamsca2k[[2]]
rm(lamsca2k)
dim(lam_2k)
lammean_2k <- apply(lam_2k,c(2,3),mean)
#MN_lammean_NLDAS <- bind_cols(select(MN_df_LC,DOW,ends_with(".z")),lammean_NLDAS)
saveRDS(lammean_2k,file="output/lammean_2k.rds")

####################################################

# 2021 NLDAS predictions
lm_curr <- readRDS(file="output/lammean_Current.rds")

# 2000 NLDAS predictions
lm_2k <- readRDS(file="output/lammean_2k.rds")

# Late century GCM
lammean_LC <- readRDS(file="output/lammean_LC.rds")

## Mid century GCM
#lammean_MC <- readRDS(file="output/lammean_MC_2k.rds")


# to get correct order of DOW
MN_fullcc <- readRDS("data/MNfullcc.rds") %>%
  filter(DOW %in% rownames(lm_2k))

MN_df_LC <- MN_fullcc %>%
  select(-MC) %>%
  pivot_wider(names_from = "GCM",values_from = "LC")

# MN_df_MC <- MN_fullcc %>%
#   select(-LC) %>%
#   pivot_wider(names_from = "GCM",values_from = "MC")

# coordinates
MN_locs <- select(MN_df_LC,DOW,lon,lat)

# predicted extinctions
pct_extinct_LC <- readRDS(file="output/pct_extinct_LC.rds")
#pct_extinct_LC6 <- readRDS(file="output/pct_extinct_LC6.rds")

#pct_extinct_MC <- readRDS(file="output/pct_extinct_MC_2k.rds")

# percent change in relative abundance Late Century vs Current (2021)
lampc_LC <- (lammean_LC-lm_curr)/lm_curr * 100


# Joining data frames
MN_lm_curr <- select(MN_fullcc,DOW,ends_with(".z")) %>%
  left_join(as_tibble(lm_curr,rownames="DOW"),by="DOW")

MN_lammean_LC <- select(MN_fullcc,DOW,ends_with(".z")) %>%
  left_join(as_tibble(lammean_LC,rownames="DOW"),by="DOW")  

# Creating spatial objects
library(sf)
lampc_LC_sf <- left_join(MN_locs,as_tibble(lampc_LC,rownames="DOW"),by="DOW") %>%
  pivot_longer(-c(DOW,lon,lat),names_to = "species",values_to = "pc") %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(pc) %>%
  filter(!(species=="cisco" & pc > 50))

usa <-st_as_sf(maps::map("state",fill=TRUE,plot=FALSE))

# Plot of percent change for each species across Minnesota
# Late century vs Current (2021)
ggplot() +
  geom_sf(data=filter(usa,ID=="minnesota")) + 
  geom_sf(data=lampc_LC_sf,aes(color=pc)) +
  facet_wrap(~species,nrow=2) +
  scale_color_viridis_c(na.value = "#5E5D61") +
  scale_shape_manual(values=c("Y"=4,"N"=16),guide=NULL) +
  theme_bw() +
  labs(color="% change") +
  theme(text=element_text(size=18))


# Summary table of percent change
lampc_LC_sf %>%
  st_drop_geometry() %>%
  group_by(species) %>%
  summarise(avg=mean(pc),
            mdn=median(pc),
            l95=quantile(pc,0.25),
            u95=quantile(pc,0.75),
            lt0=mean(pc<0)) %>%
  mutate(tab=paste0(round(mdn,2)," (",round(l95,2),", ",round(u95,2),")")) %>%
  select(species,avg,tab,lt0) %>%
  knitr::kable(format = "latex", digits=2,linesep="",booktabs=TRUE)



# Reading in data used to fit jsPGA
MNfish <- readRDS("data/MNfishz.rds") %>%
  filter(DOW %in% MN_fullcc$DOW)


# Filtering to only lakes with catch > 0 for each species
MNfish0 <- MNfish %>%
  filter(TOTAL_CATCH>0,COMMON_NAME %in% colnames(lm_curr)) %>%
  select(DOW,species=COMMON_NAME) %>%
  distinct()


# Spatial object of predicted percent change in relative abundance within sampled distributions of each species
lampc_0_LC <- inner_join(MNfish0,MN_locs,by="DOW") %>%
  left_join(as_tibble(lampc_LC,rownames="DOW") %>%
              pivot_longer(-c(DOW),names_to = "species",values_to = "pc"),
            by=c("DOW","species")) %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(pc)

summary(pivot_wider(st_drop_geometry(lampc_0_LC),names_from = "species",values_from = "pc"))

st_drop_geometry(lampc_0_LC) %>%
  group_by(species) %>%
  summarise(n())

# Map of predicted percent change in relative abundance within sampled distributions of each species
ggplot() +
  geom_sf(data=filter(usa,ID=="minnesota")) + 
  geom_sf(data=lampc_0_LC,aes(color=pc)) +
  facet_wrap(~species,nrow=2) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(color="% change") +
  theme(text=element_text(size=18))

# Distribution of predicted percent change in relative abundance within sampled distributions of each species
st_drop_geometry(lampc_0_LC) %>%
  ggplot(aes(x=pc)) +
  geom_density() +
  facet_wrap(~species,nrow=2) +
  theme_bw() +
  labs(color="% change",y="Density",x="Percent change (%)") +
  theme(text=element_text(size=18))

# Number of lakes with an predicted increase in relative abundance within sampled distributions of each species
st_drop_geometry(lampc_0_LC) %>%
  filter(!is.na(pc)) %>%
  group_by(species) %>%
  summarise(N=sum(pc>0),
            avg=mean(pc>0))

# Adding coordinates to probability of extinction data frame
ext_probs_LC <- as_tibble(pct_extinct_LC,rownames="DOW") %>%
  left_join(MN_locs,by="DOW") %>%
  mutate(tp="Late-Century")

# summary of predicted probability of extinction within sampled distributions of each species
inner_join(MNfish0,
           pivot_longer(ext_probs_LC,-c(DOW,lon,lat,tp),names_to = "species",values_to = "prob"),
           by=c("DOW","species")) %>%
  select(-c(lon,lat,tp)) %>%
  pivot_wider(names_from = 'species',values_from = 'prob',values_fill = NA) %>%
  summary()

# Creating spatial object of predicted probability of extinction within sampled distributions of each species
ext_prob0 <- inner_join(MNfish0,
                        pivot_longer(ext_probs_LC,-c(DOW,lon,lat,tp),names_to = "species",values_to = "prob"),
                        by=c("DOW","species")) %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(prob) 

# Map of predicted probability of extinction within sampled distributions of each species
ggplot() +
  geom_sf(data=filter(usa,ID=="minnesota")) + 
  geom_sf(data=ext_prob0,
          aes(color=prob),size=3) +
  facet_wrap(~species,nrow=2) +
  scale_color_viridis_c() +
  theme(text=element_text(size=14)) +
  theme_bw() +
  labs(color="Probability \nof extinction")


# Cisco (coldwater species) have significantly higher probabilities of extinction
# Focal species of case study for predicting extinctions

# Sampled distribution of cisco
cisco_DOW <- MNfish %>%
  filter(COMMON_NAME=="cisco",TOTAL_CATCH>0) %>%
  select(DOW) %>%
  distinct()

# Creating spatial object of cisco predicted probability of extinction within sampled distribution
cisco_probs <- ext_probs_LC %>%
  select(DOW,lon,lat,cisco,tp) %>%
  filter(DOW %in% cisco_DOW$DOW)  %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(cisco)

# Map of cisco predicted probability of extinction within sampled distribution
ggplot() +
  geom_sf(data=filter(usa,ID=="minnesota")) + 
  geom_sf(data=cisco_probs,
          aes(fill=cisco),shape=21,size=4) +
  scale_fill_viridis_c() +
  theme(text=element_text(size=14)) +
  theme_bw() +
  labs(fill="Probability \nof extinction")

summary(st_drop_geometry(cisco_probs))

cisco_probs %>%
  st_drop_geometry() %>%
  summarise(avg=mean(cisco),
            tot=sum(cisco > 0.9),
            N=n())


# Predicted increases in lake surface temperature by GCM model
LCchange <- MN_locs %>%
  left_join(select(MN_fullcc,DOW,GCM,LC),by="DOW") %>%
  st_as_sf(coords = c("lon","lat"),crs=st_crs(4269)) %>%
  arrange(LC)

ggplot() +
  geom_sf(data=filter(usa,ID=="minnesota")) + 
  geom_sf(data=LCchange,aes(color=LC),size=3) +
  facet_wrap(~GCM,nrow=2) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(color="Increase °C") +
  theme(text=element_text(size=18))


# Late century GCM6
lammean6_LC <- readRDS(file="output/lammean6_LC.rds")

lc_gcm6 <- array(dim=c(8,6))

lc_gcm6 <- matrix(paste0(
  round(apply(lammean6_LC,c(2,3),mean),2), " (",
  round(apply(lammean6_LC,c(2,3),quantile,0.025),2),", ",
  round(apply(lammean6_LC,c(2,3),quantile,0.975),2),")"
  ), nrow=8)

dimnames(lc_gcm6)[[1]] <- dimnames(lammean6_LC)[[2]]
dimnames(lc_gcm6)[[2]] <- dimnames(lammean6_LC)[[3]]












