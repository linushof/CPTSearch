# Preparation -------------------------------------------------------------

# load pkgs
pacman::p_load(tidyverse, R2jags)


# Gloeckner (2012) ------------------------------------------------------------

paper <- read_rds("data/cpt/cpt_gloeckner12.rds.bz2")

## create data object for JAGS
# create Array 

preArrayData <- subset(paper, select = - c(index, problem,id))

preArrayData <- preArrayData %>%
  select(-paper, paper)

split_data <- split(preArrayData, preArrayData$subject)
nsub  # number of subjects
nprob # max number of problems
num_cols <- ncol(preArrayData) - 1  #number of columns without subject

# initiate Array
data_array <- array(NA, dim = c(nprob, num_cols, nsub))

# fill 3D-Array
for (i in 1:nsub) {
  subject_data <- as.matrix(split_data[[i]][, -1])  # Data without subject ID
  data_array[1:nrow(subject_data), 1:ncol(subject_data), i] <- subject_data  # Fill the data
}

# check dimension
dim(data_array)


# Converge array in list again to have NAs in the end of each element vector
# Due to control variable in model code the NAs won't be taken into account

nsolvedprob #double check number of solved problems per participant
nsub  #double check number of subjects

## load and prepare data ---------------------------------------------------

data_list <- list(
  nprob = nsolvedprob,          
  nsub =nsub,           
  HA = apply(t(data_array[, 1, ]),2, as.numeric), 
  LA = apply(t(data_array[, 3, ]),2, as.numeric), 
  HB = apply(t(data_array[, 5, ]),2, as.numeric),
  LB = apply(t(data_array[, 7, ]),2, as.numeric),
  sprobHA = apply(t(data_array[, 2, ]),2, as.numeric),
  sprobLA = apply(t(data_array[, 4, ]),2, as.numeric),
  sprobHB = apply(t(data_array[, 6, ]),2, as.numeric),
  sprobLB = apply(t(data_array[, 8, ]),2, as.numeric),  
  choice = apply(t(data_array[, 9, ]) ,2, as.numeric) 
)


## fitting ------------------------------------------------------------

params <- c("mu.alpha", "alpha", "mu.gamma", "gamma", "mu.delta", "delta", "mu.rho", "rho") # fitted parameters that should be shown in the results

# function for creating parameter values to initialize the MCMC chains with
# when rnorm(1, .4, .1), then .4 is the initial value for the parameter on the desired scale

params_init <- function(){
  list("mu.probit.alpha" = qnorm(rnorm(1, .4, .1)) , # hyper parameters (mu. prefix refers to  group level)
       "mu.probit.gamma" = qnorm(rnorm(1, .5, .1)) ,
       "mu.probit.delta" = qnorm(rnorm(1, .7, .1)) , 
       "mu.probit.rho" = qnorm(.01) ,  
       "probit.alpha" = qnorm(rnorm(nsub, .4, .1)) , # individual level parameters
       "probit.gamma" = qnorm(rnorm(nsub, .5, .1)) ,
       "probit.delta" = qnorm(rnorm(nsub, .7, .1)) , 
       "probit.rho" = qnorm(rep(.01, nsub) ))
}


# run JAGS model

mfit <- jags.parallel(
  
  data = data_list , # data list
  inits = params_init , # creates list of initial values for each chain
  parameters.to.save = params , 
  model.file = "code/models/CPT_hierarchical_array.txt" , # model code, see file
  n.chains = 6 , # number of MCMC chains
  n.iter = 2000 , # number of iterations (should be set much higher once it's clear that the model works)
  n.burnin = 1000 , # first 1000 samples of each chain are discarded
  n.thin = 1 , # with 1, every sample is stored, with 2, every 2nd sample is stored, ... to reduce autocorrelatons, use higher values. however, higher values require more iterations
  n.cluster = 6 , # compute MCMC chains in parallel (on different cores of the computer)
  DIC = TRUE # store all posterior samples
)

# sanity check 2
data_list$sprobHA+data_list$sprobLA
data_list$sprobHB+data_list$sprobLB


## store results -----------------------------------------------------------------

# traceplot(mfit) # can be cumbersome with many subjecyts


# summary of posterior distributions

estimates <- mfit$BUGSoutput$summary %>% round(4) %>%  as_tibble(rownames = "parameter")
posteriors <- mfit$BUGSoutput$sims.matrix %>% as_tibble()

papername <- unique(paper$paper)
write_rds(estimates, paste("data/cpt/cpt_", tolower(papername), "_estimates.rds.bz2", sep =""), compress = "bz2")
write_rds(posteriors, paste("data/cpt/cpt_", tolower(papername), "_posteriors.rds.bz2", sep =""), compress = "bz2")


