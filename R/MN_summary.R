library(rstan)
library(tidyverse)

M="16"
x <- list.files(paste0("output/MN/M",M))
simnum <- sapply(strsplit(str_remove(x,".rds"),"_"),"[[",2)

out_list <- list()
for(i in 1:length(x)){
  name <- paste0("iter",simnum[i])
  out_list[[name]] <- readRDS(paste0("output/MN/M",M,"/",x[i]))
  print(i)
}

out <- sflist2stanfit(out_list)
rm(out_list)

saveRDS(out,file=paste0("output/MN/MNpost_M",M,".rds"))


###################################################################################
library(rstan)
library(tidyverse)
out <- readRDS("output/MN/MNpost_M16.rds")
dat <- readRDS("data/MNdat_M16.rds")

mean(rowSums(get_elapsed_time(out)))/(60^2)

chains <- rstan::extract(out,pars=c("beta_0","beta","theta","phi","tau","A","Sigma_A"))
names(chains)
lapply(chains,dim)
rm(out)
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
# tau
dimnames(chains$tau)[[2]] <- dat$levs


# Combining betas
chains$allbeta = array(dim = dim(chains$beta) + c(0,1,0))
for(i in 1:dat$K){
  chains$allbeta[,,i] <- cbind(beta0=chains$beta_0[,i],chains$beta[,,i])
}
dimnames(chains$allbeta)[[2]] <- c("beta0",dat$covs)
dimnames(chains$allbeta)[[3]] <- dat$levs

# Full covariance matrix of A: TSigmaT
chains$TST = array(dim=dim(chains$Sigma_A))
for(j in 1:dim(chains$A)[[1]]){
  chains$TST[j,,] <- diag(chains$tau[j,]) %*% chains$Sigma_A[j,,] %*% diag(chains$tau[j,])
}
dimnames(chains$TST)[[2]] <- dat$levs
dimnames(chains$TST)[[3]] <- dat$levs

################################# Tables #################################

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
    avg = tmpf(x,func='mean'),
    median = tmpf(x,func="median"),
    sd = tmpf(x,func="sd"),
    u95 = tmpf(x,func="my975"),
    l95 = tmpf(x,func="my025")
  )
  )
}

chaintab <- lapply(chains,chainfun)

###### chaintab figures

# beta

betagg <- bind_rows(as_tibble(chaintab$beta$median,rownames="Prmtr") %>%
                      pivot_longer(-Prmtr,names_to="Species") %>%
                      mutate(St="median"),
                    as_tibble(chaintab$beta$l95,rownames="Prmtr") %>%
                      pivot_longer(-Prmtr,names_to="Species") %>%
                      mutate(St="l95"),
                    as_tibble(chaintab$beta$u95,rownames="Prmtr") %>%
                      pivot_longer(-Prmtr,names_to="Species") %>%
                      mutate(St="u95")
) %>%
  pivot_wider(names_from=St,values_from=value) %>%
  mutate(Prmtr=factor(Prmtr,levels=c("elevation.z","lakearea.z","secchi.z","total.ag.z","total.dev.z"),
                      labels=c("Elevation","Lake area","Secchi","Agriculture","Developed")),
         zero=ifelse(l95<0 & u95>0,"0","1"))

ggplot(betagg) +
  geom_hline(yintercept=0,linetype=2) +
  geom_point(aes(x=Species,y=median,color=zero,shape=zero),size=3) +
  geom_errorbar(aes(x=Species,ymin=l95,ymax=u95,color=zero)) +
  facet_grid(cols=vars(Prmtr),scales = "free_y") +
  scale_color_manual(values = c("0"='red',"1"='blue')) +
  scale_shape_manual(values = c("0"=17,"1"=19)) +
  labs(y=expression(beta)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=30),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=14,angle=55,hjust=1),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=18),
        legend.position = "none")



# theta NOT USED IN MANUSCRIPT
# thetagg <- bind_rows(as_tibble(chaintab$theta$median,rownames="Species") %>%
#                        pivot_longer(-Species,names_to="Gear") %>%
#                        mutate(St="median"),
#                      as_tibble(chaintab$theta$l95,rownames="Species") %>%
#                        pivot_longer(-Species,names_to="Gear") %>%
#                        mutate(St="l95"),
#                      as_tibble(chaintab$theta$u95,rownames="Species") %>%
#                        pivot_longer(-Species,names_to="Gear") %>%
#                        mutate(St="u95")
# ) %>%
#   pivot_wider(names_from=St,values_from=value)
# 
# ggplot(thetagg) +
#   geom_point(aes(x=Gear,y=median)) +
#   geom_errorbar(aes(x=Gear,ymin=l95,ymax=u95)) +
#   geom_hline(yintercept=0.5,linetype=2) +
#   facet_grid(cols=vars(Species)) +
#   labs(y=expression(theta)) +
#   theme_bw() +
#   theme(axis.title.y = element_text(size=24),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size=20),
#         axis.text.x = element_text(size=14,angle=90,hjust=1,vjust=0.5),
#         panel.grid.major.x = element_blank(),
#         strip.text = element_text(size=14),
#         legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black"))


# Covariance parameters


Sigmagg <- bind_rows(as_tibble(chaintab$Sigma_A$median,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="median"),
                     as_tibble(chaintab$Sigma_A$l95,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="l95"),
                     as_tibble(chaintab$Sigma_A$u95,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="u95")) %>%
  pivot_wider(names_from=St,values_from=value) %>%
  mutate(median=ifelse(l95<0 & u95>0,round(median,2),paste0(round(median,2),"*"))) %>%
  mutate(txt=ifelse(Sp1==Sp2,"",
                    paste0(median,
                           "\n (",
                           round(l95,2),
                           ",",
                           round(u95,2),
                           ")")))

taugg <- bind_rows(as_tibble(chaintab$tau$median,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="median"),
                   as_tibble(chaintab$tau$l95,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="l95"),
                   as_tibble(chaintab$tau$u95,rownames="Sp1") %>%
                     pivot_longer(-Sp1,names_to="Sp2") %>%
                     mutate(St="u95")) %>%
  select(-Sp2) %>%
  pivot_wider(names_from=St,values_from=value) %>%
  mutate(txt=paste0(round(median,2),"\n (",round(l95,2),",",round(u95,2),")"))

tstgg <- bind_rows(as_tibble(chaintab$TST$median,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="median"),
                     as_tibble(chaintab$TST$l95,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="l95"),
                     as_tibble(chaintab$TST$u95,rownames="Sp1") %>%
                       pivot_longer(-Sp1,names_to="Sp2") %>%
                       mutate(St="u95")) %>%
  pivot_wider(names_from=St,values_from=value) %>%
  mutate(median=ifelse(l95<0 & u95>0,round(median,2),paste0(round(median,2),"*"))) %>%
  mutate(txt=ifelse(Sp1==Sp2,"",
                    paste0(median,
                           "\n (",
                           round(l95,2),
                           ",",
                           round(u95,2),
                           ")")))

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat,diag=TRUE)]<- NA
  return(cormat)
}

sigcormed <- get_upper_tri(chaintab$Sigma_A$median)

taumed <- diag(chaintab$tau$median)
taumed[lower.tri(taumed,diag=FALSE)] = NA
taumed[upper.tri(taumed,diag=FALSE)] <- NA
dimnames(taumed) = dimnames(sigcormed)

tstmed <- t(get_upper_tri(chaintab$TST$median))

melted_sigcormed <- reshape2::melt(sigcormed, na.rm = TRUE) %>%
  rename(Sp1=Var1,Sp2=Var2) %>%
  mutate(Sp1=factor(Sp1),
         Sp2=factor(Sp2)) %>%
  left_join(select(Sigmagg,Sp1,Sp2,txt),
            by=c("Sp1","Sp2"))

melted_taumed <- reshape2::melt(taumed, na.rm = TRUE) %>%
  rename(Sp1=Var1,Sp2=Var2) %>%
  mutate(Sp1=factor(Sp1),
         Sp2=factor(Sp2)) %>%
  left_join(select(taugg,Sp1,txt),
            by=c("Sp1"))

melted_tstmed <- reshape2::melt(tstmed, na.rm = TRUE) %>%
  rename(Sp1=Var1,Sp2=Var2) %>%
  mutate(Sp1=factor(Sp1),
         Sp2=factor(Sp2)) %>%
  left_join(select(tstgg,Sp1,Sp2,txt),
            by=c("Sp1","Sp2"))


# Heatmap
ggplot() +
  geom_tile(data = melted_sigcormed, 
            aes(Sp2, Sp1, fill = value),
            color = "white") +
  geom_tile(data=melted_tstmed,aes(Sp2,Sp1),fill="white",color='black') +
  geom_tile(data=melted_taumed,aes(Sp2,Sp1),fill="gray",color="black",linewidth=1) +
  geom_text(data = melted_sigcormed, aes(Sp2, Sp1,label=txt),size=5) +
  geom_text(data = melted_taumed, aes(Sp2, Sp1,label=txt),size=5) +
  geom_text(data = melted_tstmed, aes(Sp2, Sp1,label=txt),size=5) +
  scale_fill_gradient2(low="red",
                       mid="white",
                       high="blue",
                       limits = c(-1,1),) +
  theme_void() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 20, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title = element_blank(),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18)) +
  coord_fixed() +
  labs(fill="Species \ndependency")





#########################################################################



theta_med <-t(chaintab$theta$median)
theta_l95 <- t(chaintab$theta$l95)
theta_u95 <- t(chaintab$theta$u95)
thetamat <- matrix(paste0(theta_med," (",theta_l95,", ",theta_u95,")"),
                   nrow=nrow(theta_med),
                   dimnames=dimnames(theta_med))

beta_med <- chaintab$allbeta$median
beta_l95 <- chaintab$allbeta$l95
beta_u95 <- chaintab$allbeta$u95
betamat <- matrix(paste0(beta_med," (",beta_l95,", ",beta_u95,")"),
                  nrow=nrow(beta_med),
                  dimnames=dimnames(beta_med))

sum(!(beta_l95<0 & beta_u95>0))
mean(!(beta_l95<0 & beta_u95>0))

A_med <- chaintab$A$median
A_l95 <- chaintab$A$l95
A_u95 <- chaintab$A$u95
Amat <- matrix(paste0(A_med," (",A_l95,", ",A_u95,")"),
                 nrow=nrow(A_med),
                 dimnames=dimnames(A_med))

#Species summaries
summary((A_med))

#Overlapping zero
colSums(((A_l95) < 0 & (A_u95) > 0))
colMeans(((A_l95) < 0 & (A_u95) > 0))



sigma_mat <- matrix(paste0(round(chaintab$Sigma_A$median,2)," (",
                           round(chaintab$Sigma_A$l95,2),", ",
                           round(chaintab$Sigma_A$u95,2),")"),
                    nrow=nrow(chaintab$Sigma_A$median),
                    dimnames=dimnames(chaintab$Sigma_A$median))

sigma_mat[upper.tri(sigma_mat,diag = FALSE)] <- "-"
diag(sigma_mat) <- 1

tau_med <- chaintab$tau$median
tau_l95 <- chaintab$tau$l95
tau_u95 <- chaintab$tau$u95
taumat <- matrix(paste0(tau_med," (",tau_l95,", ",tau_u95,")"),
                 nrow=1)
dimnames(taumat)[[1]] <- "tau"
dimnames(taumat)[[2]] <- names(tau_med)


as_tibble(rbind(beta_0mat,
                betamat,
                thetamat,
                Amat,
                taumat,
                sigma_mat),
          rownames="Parameter") %>%
  knitr::kable(format = "latex", digits=2,linesep="",
               booktabs=TRUE,longtable=TRUE,caption = "My caption") %>%
  kableExtra::kable_styling(latex_options=c("repeat_header")) %>%
  kableExtra::landscape()













































