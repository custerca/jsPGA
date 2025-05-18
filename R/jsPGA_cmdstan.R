library(cmdstanr)

# This code works by bringing in the jobid from the Roar supercomputer
# Each model realization is indexed by jobid to use a different randomly sampled value for CTmax and Topt
# Output from each model is saved and then merged into single posterior within analysis scripts

# jobid = as.integer(Sys.getenv("INPUT"))
# print(jobid)

# Read in list object (simulated or case study)
dat <- readRDS("data/simdat_M16_TC.rds")

# Grab raondomly sampled values for thermal performance parameters
dat$CTmax_prior = dat$ctmax[jobid,]
dat$Topt_prior = dat$topt[jobid,]

# cmdstanr doesn't accept non-numerical data
# easier (i.e., lazier) to just NULL here than switch all previous code (originally used rstan package)
dat$ctmax <- NULL
dat$topt <- NULL
dat$levs <- NULL
dat$lake_id <- NULL
dat$covs <- NULL

# Initial values function for parameters that were problematic
init_fun <- list(theta=matrix(rep(0.5,dat$K * dat$J),nrow=dat$K),
                 recip_phi=1)

#identifying stan file path
stanfile <- "stan/jspga.stan"

#creating stan model
mod <- cmdstan_model(stanfile)

# sampling from posterior
fit <- mod$sample(data = dat, 
                  seed=1,
                  init = list(init_fun),
                  chains = 1,
                  thin=5,
                  #max_treedepth = 12,
                  iter_warmup = 100,
                  iter_sampling = 100,
                  refresh = 10)

# converting cmdstanr output into stanfit object to utilize familiar functions for posterior merging and analysis
stancmd <- rstan::read_stan_csv(fit$output_files())

# saving stanfit object
saveRDS(stancmd,file=paste0("output/Sim/M16/simoutM16_",jobid,".rds"))
