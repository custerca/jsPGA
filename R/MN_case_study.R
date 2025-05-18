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

fishfull <- readRDS("data/MNfishz.rds")[COMMON_NAME %in% c("black crappie",
                                                           "bluegill",
                                                           "cisco",
                                                           "largemouth bass",
                                                           "northern pike",
                                                           "smallmouth bass",
                                                           "walleye",
                                                           "yellow perch")]

# spatial basis vectors
MN_sbv <- readRDS("data/MNsbv_M16.rds") %>%
  distinct()

fishraw <- fishfull[,`:=`(GN=ifelse(GEAR=="GN",1,0),TN=ifelse(GEAR=="TN",1,0))][DOW %in% MN_sbv$DOW]

fishraw$DOW <- factor(fishraw$DOW)
fishraw$COMMON_NAME <- factor(fishraw$COMMON_NAME)

setorder(fishraw,DOW,SURVEYDATE,COMMON_NAME,GN)

# Begin data manipulation for stan (and general return to dplyr) ---------------------------------------------------------------

dat = list()

# number of basis functions
M=16

# create data from data ---------------------------------------------------

dat$K <- n_distinct(fishraw$COMMON_NAME)
dat$levs <- levels(fishraw$COMMON_NAME)
dat$lake_index <- (fishraw %>% filter(COMMON_NAME == dat$levs[1]))$DOW
dat$lake_id <- levels(fishraw$DOW)
dat$n_lakes <- length(levels(fishraw$DOW))

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
  select(secchi.z,lakearea.z,total.dev.z,total.ag.z,elevation.z) %>% 
  as.matrix()

dat$covs <- fishraw %>% 
  filter(EFFORT != 0) %>%
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  filter(COMMON_NAME == dat$levs[1]) %>% 
  select(secchi.z,lakearea.z,total.dev.z,total.ag.z,elevation.z) %>% 
  names()

dat$temp <- fishraw %>% 
  filter(EFFORT != 0) %>% 
  filter(COMMON_NAME == dat$levs[1]) %>% 
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  select(july5yr) %>% 
  pull()

dat$J <- 2 # number of gears, hard coded for now
dat$P <- dim(dat$X)[2] # number of variables

dat$alpha <- rep(1, dat$J) # Dirichlet prior parameter


dat$Y <- fishraw %>% 
  filter(EFFORT != 0) %>% 
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  select(DOW, SURVEYDATE, COMMON_NAME, GN, TOTAL_CATCH) %>% 
  pivot_wider(names_from = COMMON_NAME, values_from = TOTAL_CATCH) %>% 
  select(-c(DOW, SURVEYDATE, GN)) %>% 
  as.matrix()

dat$y_vec = c(dat$Y)
dat$N = dim(dat$Y)[1] # number of observations


dat$Psi <- fishraw %>% 
  filter(EFFORT != 0) %>%
  arrange(DOW, SURVEYDATE, COMMON_NAME, GN) %>% 
  filter(COMMON_NAME == dat$levs[1]) %>% 
  left_join(MN_sbv, by = "DOW") %>% 
  select(starts_with("comp")) %>% 
  as.matrix() %>% 
  unname()

dat$M <- M

# Literature values

# CTmax
ctmax <- list()
# # black bullhead
# ctmax[["black bullhead"]] <- c(38.1,37.5,35:39)

# black crappie
ctmax[["black crappie"]] <- c(34.9,36:40,38:40)

# bluegill
ctmax[["bluegill"]] <- c(35.8,39.1,36.6,37.5,37.9,37.5,41.4,
                         35.6,38.5,37.9,36.8,38.8,40.4,36.7,
                         38,40.9,36.3,37.4,40.9,37,39.6,41.4,38.3)

# # brown bullhead
# ctmax[["brown bullhead"]] <- c(37.1,37.8,37.5,38)

# cisco species
ctmax[["cisco"]] <- c(26.2)

# largemouth bass
ctmax[["largemouth bass"]] <- c(37.80,36.70,40.10,36.30,35.40,36.70,38.50)

# northern pike
ctmax[["northern pike"]] <- c(33.6,33.25,30.8)

# smallmouth bass
ctmax[["smallmouth bass"]] <- c(36.90,36.30)

# walleye
ctmax[["walleye"]] <- c(34.40,34.80,35.00,34.30)

# # white sucker
# ctmax[["white sucker"]] <- c(34.9,31.6,35.1,36.1,32.7)

# # yellow bullhead
# ctmax[["yellow bullhead"]] <- c(36.4,38,37.9,35)

# yellow perch
ctmax[["yellow perch"]] <- c(33.50,35.00,33.4,34)

ctmax_mean <- sapply(ctmax,mean)
ctmax_sd <- sapply(ctmax,sd)
ctmax_sd["cisco"] <- mean(ctmax_sd,na.rm=TRUE)


# Topt
topt <- list()
# # black bullhead
# topt[["black bullhead"]] <- c(18:29,23:24)

# black crappie
topt[["black crappie"]] <- c(22:25)

# bluegill
topt[["bluegill"]] <- c(30,30,30,31,30.1,29,30,31,24,27,27,31.2)

# # brown bullhead
# topt[["brown bullhead"]] <- c(32,28.2,29.9)

# cisco species
topt[["cisco"]] <- c(18.1, 13:18) # Ty's code, C. artedii

# largemouth bass
topt[["largemouth bass"]] <- c(30.69,25.00,26.00,27.00,28.00,27.00,
                               30.00,23.90,27.00,25:30)

# northern pike
topt[["northern pike"]] <- c(19,21,20.9,19.8,26,18:25,19:21,21,26)

# smallmouth bass
topt[["smallmouth bass"]] <- c(28.00,25.00,26.00,25.00,29.00,27.00)

# walleye
topt[["walleye"]] <- c(22.1,25.2,22,22.6,22)

# # white sucker
# topt[["white sucker"]] <- c(27,24,26.9,26,24,26.9,24,24,24)

# # yellow bullhead
# topt[["yellow bullhead"]] <- c(28.3,28.8,27.6)

# yellow perch
topt[["yellow perch"]] <- c(22.5,23,24.2,23,28,23:24,29,23,26:30,24.7)

topt_mean <- sapply(topt,mean)
topt_sd <- sapply(topt,sd)


# CTmin
ctmin <- list()
# # black bullhead
# ctmin[["black bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT from brown bullhead

# black crappie
ctmin[["black crappie"]] <- c(3, 0.5, 0.1, rep(1,11),0.017, #CTmin from bluegill
                              3,3,5,7,10,11,15,6,11 # LILT from bluegill
)

# bluegill
ctmin[["bluegill"]] <-  c(3, 0.5, 0.1, rep(1,11),0.017, #CTmin
                          3,3,5,7,10,11,15,6,11 # LILT
)

# # brown bullhead
# ctmin[["brown bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT

# cisco species
ctmin[["cisco"]] <- c(0.3, # CTmin
                      0,0.5,3,4.7 # LILT
)

# largemouth bass
ctmin[["largemouth bass"]] <- c(3.2,7.3, 10.7, # CTmin
                                5,7,11,5.5,11.8,10 # LILT
)

# northern pike
ctmin[["northern pike"]] <- c(0.1,5,3) #LILT

# smallmouth bass
ctmin[["smallmouth bass"]] <- c(2,4,4,7,10,10,2,4,7,10,10.1,1.6) # LILT

# walleye
ctmin[["walleye"]] <- c(0.1,5,3) #LILT from NP

# 

# # white sucker
# ctmin[["white sucker"]] <- c(2,3,6,6,2.5,6.6,4.8,6.1,4.8) #LILT
# 

# # yellow bullhead
# ctmin[["yellow bullhead"]] <- c(-1, 1, 1.3, 3.7, 7, 0.5, 4, 6.8) # LILT from brown bullhead

# yellow perch
ctmin[["yellow perch"]] <- c(1.1, #CTmin
                             4,1.1,3.7,6.8 #LILT
)


ctmin_mean <- sapply(ctmin,mean)

# sigma
sigma_trc <- (topt_mean - ctmin_mean)/4

# Number of iterations
nsim <- 500
ctmax_mat <- matrix(NA,nrow=nsim,ncol=dat$K)
topt_mat <- matrix(NA,nrow=nsim,ncol=dat$K)

mxtmp <- fishraw %>%
  filter(TOTAL_CATCH > 0) %>%
  group_by(COMMON_NAME) %>%
  summarise(mxtmp=max(july5yr)) %>%
  select(mxtmp) %>%
  pull()


# container to to hold stan object from each model
j=1L
set.seed(1)
while(j <= nsim){
  
  ctmax_vec <- rnorm(dat$K, mean = ctmax_mean, sd = ctmax_sd)
  topt_vec <- rnorm(dat$K, mean = topt_mean, sd = topt_sd)
  
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
dat$sigma <- sigma_trc

saveRDS(dat,file=paste0("data/MNdat_M",M,".rds"))

################################# TRC Curves ################################# 
# Thermal performance function 
TRC = function(temp, CTmax, Topt, sigma){
  
  trc = case_when(
    temp <= Topt ~ exp(-((temp - Topt)/(2*sigma))^2), 
    CTmax >= temp & temp > Topt ~ 1 - ((temp - Topt)/(Topt - CTmax))^2,
    temp > CTmax ~ 0)
  
  return(trc)
  
}

colnames(dat$ctmax) <- dat$levs
ctmax_df <- as_tibble(rbind(ctmax_mean,dat$ctmax)) %>%
  mutate(iter=row_number())

ctmax_l <- ctmax_df %>%
  pivot_longer(-iter,
               names_to = "Spp",
               values_to = "ctmax")


colnames(dat$topt) <- dat$levs
topt_df <- as_tibble(rbind(topt_mean,dat$topt)) %>%
  mutate(iter=row_number())

topt_l <- topt_df %>%
  pivot_longer(-iter,
               names_to = "Spp",
               values_to = "topt")

sigma_df <- tibble(Spp=dat$levs,
                    sigma=sigma_trc)


trc_plot <- tibble(iter=1:100) %>%
  mutate(t0=list(seq(1,43,by=0.5)),
         Spp=list(dat$levs)) %>%
  unnest(cols=c(t0)) %>%
  unnest(cols=Spp) %>%
  left_join(ctmax_l,by=c("iter","Spp")) %>%
  left_join(topt_l,by=c("iter","Spp")) %>%
  left_join(sigma_df,by="Spp") %>%
  rowwise() %>%
  mutate(trc=TRC(temp = t0,CTmax = ctmax,Topt = topt,sigma = sigma)) %>%
  ungroup() %>%
  mutate(iter=paste0("iter",iter))


ggplot() +
  geom_line(data=trc_plot[trc_plot$iter != "iter1",],
            aes(x=t0,y=trc,
                color=Spp,
                group=interaction(iter,Spp)
            ),
            alpha=0.2,size=0.4) +
  geom_line(data=trc_plot[trc_plot$iter == "iter1",],
            aes(x=t0,y=trc,
                color=Spp,
                group=interaction(iter,Spp)
            ),
            alpha=2,size=2) +
  facet_wrap(~Spp,nrow=2) +
  scale_color_viridis_d(guide=guide_legend(override.aes = list(linewidth = 10))) +
  labs(x="Temperature (\u00B0C)",y="Thermal Performance Scalar") +
  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.grid = element_blank(),
        strip.text = element_text(size=20),
        legend.position = "none")

