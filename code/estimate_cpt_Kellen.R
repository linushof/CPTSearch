# load packages
pacman::p_load(tidyverse, R2jags)

#load and preprocess data with the pre-processing script
paper <- read_rds(glue("data/PreprocessedPaperData/cpt_kellen16.rds.bz2"))

papername <- unique(paper$paper)

#Get vector of the number of solved problems per participant
nsolvedprob <- paper %>% 
  group_by(subject) %>% 
  summarise(nsolvedprob = n()) 

nsolvedprob <- nsolvedprob$nsolvedprob

#Create a list of subject and the respective switch rate
switch_rate_data <- paper %>%
  select(paper, subject, cat_switch, avg_switchrate) %>% 
  distinct() 

switch_rate_data$subject <- seq_len(nrow(switch_rate_data))

# prepare data for CPT --------------------------------------------------------------------


# compute number of problems and subjects in each data set 

nproblem <-paper %>% 
  group_by(paper) %>% 
  distinct(problem) %>% 
  summarise(nproblem = n()) 

nsubjects <- paper %>% 
  group_by(paper) %>% 
  distinct(subject) %>% 
  summarise(nsubject = n()) 

overview <- left_join(nproblem, nsubjects, by=join_by(paper)) %>% 
  mutate(nchoice = nproblem * nsubject)


# create data objects for JAGS

## data that is the same for all subjects (but see ISSUES above)

problems <- paper %>% distinct(problem, HA, LA, HB, LB)  
problems <- head(problems,42)

## data that differs between subjects 
## create data frames that store the trial/problem-specific choices and sampled probabilities for each subject
## in model code referred to as [j=subject=rows, i=problem=columns]

choices <- paper %>% select(subject, problem, choice) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = choice) %>% select(-subject)
sprobHA <- paper %>% select(subject, problem, sprobHA) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobHA) %>% select(-subject)
sprobLA <- paper %>% select(subject, problem, sprobLA) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobLA) %>% select(-subject)
sprobHB <- paper %>% select(subject, problem, sprobHB) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobHB) %>% select(-subject)
sprobLB <- paper %>% select(subject, problem, sprobLB) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobLB) %>% select(-subject)

nprob <- nrow(problems)
nsub <- nrow(choices)

# sanity checks (should be true)
sprobHA+sprobLA == 1
sprobHB+sprobLB == 1
nrow(choices) == nrow(sprobHA) 
nrow(sprobHA) ==  nrow(sprobLA)
nrow(sprobLA) ==  nrow(sprobHB)
nrow(sprobHB) ==  nrow(sprobLB)

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



# fitting ------------------------------------------------------------

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
  
  data = dat , # data list
  inits = params_init , # creates list of initial values for each chain
  parameters.to.save = params , 
  model.file = "code/models/CPT_hierarchical_Kellen.txt" , # model code, see file
  n.chains = 5 , # number of MCMC chains
  n.iter = 2000 , # number of iterations (should be set much higher once it's clear that the model works)
  n.burnin = 1000 , # first 1000 samples of each chain are discarded
  n.thin = 1 , # with 1, every sample is stored, with 2, every 2nd sample is stored, ... to reduce autocorrelatons, use higher values. however, higher values require more iterations
  n.cluster = 6 , # compute MCMC chains in parallel (on different cores of the computer)
  DIC = TRUE # store all posterior samples
)


# results -----------------------------------------------------------------

# traceplot(mfit) # can be cumbersome with many subjecyts


# summary of posterior distributions

fits <- mfit$BUGSoutput$summary %>% round(4) %>%  as_tibble(rownames = "parameter")


# plot weighting function

## weighting function for group level
mu.fits <- fits %>% 
  filter(parameter %in% c("mu.alpha", "mu.gamma", "mu.delta", "mu.rho")) %>% 
  mutate(subject = 0, 
         parameter = c("alpha", "delta", "gamma", "rho")) 

mu.weights <- mu.fits %>% 
  select(subject, parameter, mean) %>% 
  pivot_wider(names_from = parameter, values_from = mean) %>% 
  select(-c(alpha, rho)) %>%
  expand_grid(p = seq(0, 1, .01)) %>% # create vector of sampled relative frequencies
  mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) # comput

nsub
## weighting function for individual level
ind.fits <- fits %>% 
  filter(! parameter %in% c("mu.alpha", "mu.gamma", "mu.delta", "mu.rho", "deviance")) %>% 
  mutate(subject = rep(1:nsub, 4),
         parameter = c(rep("alpha", nsub), rep("delta", nsub),  rep("gamma", nsub), rep("rho", nsub))
  )

ind.weights <- ind.fits %>%
  select(subject, parameter, mean) %>% 
  pivot_wider(names_from = parameter, values_from = mean) %>% 
  select(-c(alpha, rho)) %>%
  expand_grid(p = seq(0, 1, .01)) %>% # create vector of sampled relative frequencies
  mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) # comput

# plot
mu.weights %>% 
  ggplot(aes(p, w)) +
  scale_x_continuous(breaks = seq(0, 1, length.out = 3)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 3)) +
  labs(x = "p",
       y = "w(p)") +
  geom_line(data=ind.weights, aes(group=subject), color = "gray") + 
  geom_abline(intercept=0, slope =1, linewidth = 1, "black", linetype = "dashed") +
  geom_line(linewidth = 1, color = "black") +
  theme_classic()

