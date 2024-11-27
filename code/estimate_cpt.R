
# Data set: XY ------------------------------------------------------------


## load and prepare data ---------------------------------------------------

# load data


# create list with all objects

dat <- list(
  HA = problems$HA ,
  LA = problems$LA , 
  HB = problems$HB ,
  LB = problems$LB , 
  sprobHA = round(sprobHA, 2) , 
  sprobLA = round(sprobLA, 2) ,
  sprobHB = round(sprobHB, 2) , 
  sprobLB = round(sprobLB, 2) ,
  choice = choices ,
  nprob = nprob ,
  nsub = nsub
)



## fitting ------------------------------------------------------------

params <- c("mu.alpha", "alpha", "mu.gamma", "gamma", "mu.delta", "delta", "mu.rho", "rho") # fitted parameters that should be shown in the results

# function for creating parameter values to initialize the MCMC chains with
# when rnorm(1, .4, .1), then .4 is the initial value for the parameter on the desired scale

params_init <- function(){
  list("mu.probit.alpha" = qnorm(rnorm(1, .4, .1)) , # hyper parameters (mu. prefix refers to  group level)
       "mu.probit.gamma" = qnorm(rnorm(1, .5, .1)) ,
       "mu.probit.delta" = qnorm(rnorm(1, .7, .1)) , 
       "mu.probit.rho" = qnorm(.1) ,  
       "probit.alpha" = qnorm(rnorm(nsub, .4, .1)) , # individual level parameters
       "probit.gamma" = qnorm(rnorm(nsub, .5, .1)) ,
       "probit.delta" = qnorm(rnorm(nsub, .7, .1)) , 
       "probit.rho" = qnorm(rep(.01, nsub) ))
}


# run JAGS model

mfit <- jags.parallel(
  
  data = dat , # data list
  inits = params_init , # creates list of initial values for each chain
  parameters.to.save = params , 
  model.file = "code/models/CPT_hierarchical.txt" , # model code, see file
  n.chains = 5 , # number of MCMC chains
  n.iter = 2000 , # number of iterations (should be set much higher once it's clear that the model works)
  n.burnin = 1000 , # first 1000 samples of each chain are discarded
  n.thin = 1 , # with 1, every sample is stored, with 2, every 2nd sample is stored, ... to reduce autocorrelatons, use higher values. however, higher values require more iterations
  n.cluster = 6 , # compute MCMC chains in parallel (on different cores of the computer)
  DIC = TRUE # store all posterior samples
)


## store results -----------------------------------------------------------------

# traceplot(mfit) # can be cumbersome with many subjecyts


# summary of posterior distributions

fits <- mfit$BUGSoutput$summary %>% round(4) %>%  as_tibble(rownames = "parameter")
fits



# Data set: WZ --------------------------------------------------------------


