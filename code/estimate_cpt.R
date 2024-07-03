pacman::p_load(tidyverse, R2jags)

source("code/preprocessing.R")

# prepare data for CPT --------------------------------------------------------------------

dat_cpt <- choices %>% 
  
  filter(type == "free", 
         dom == "Gain", 
         probA3 == 0 & probA4 == 0 & probA5 == 0, 
         probB3 == 0 & probB4 == 0 & probB5 == 0) %>% 
  
  mutate(
    
    HA = if_else(outA1 > outA2, outA1, outA2) , 
    LA = if_else(outA1 < outA2, outA1, outA2) ,
    HB = if_else(outB1 > outB2, outB1, outB2) , 
    LB = if_else(outB1 < outB2, outB1, outB2) ,
    sprobHA = if_else(outA1 > outA2, sprobA1, sprobA2) ,
    sprobLA = if_else(outA1 < outA2, sprobA1, sprobA2) , 
    sprobHB = if_else(outB1 > outB2, sprobB1, sprobB2) ,
    sprobLB = if_else(outB1 < outB2, sprobB1, sprobB2)
    
  ) %>%  
  
  select(index:subject , 
         HA, sprobHA, LA, sprobLA , 
         HB, sprobHB, LB, sprobLB , 
         choice) 

# select data set/paper

paper <- dat_cpt %>% filter(paper == "Kellen16")

# hierarchical ------------------------------------------------------------

# prepare data 

problems <- paper %>% distinct(problem, HA, LA,HB, LB)
choices <- paper %>% select(subject, problem, choice) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = choice) %>% select(-subject)
sprobHA <- paper %>% select(subject, problem, sprobHA) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobHA) %>% select(-subject)
sprobLA <- paper %>% select(subject, problem, sprobLA) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobLA) %>% select(-subject)
sprobHB <- paper %>% select(subject, problem, sprobHB) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobHB) %>% select(-subject)
sprobLB <- paper %>% select(subject, problem, sprobLB) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobLB) %>% select(-subject)

nprob <- nrow(problems)
nsub <- nrow(choices)

# sanity checks
sprobHA+sprobLA == 1
sprobHB+sprobLB == 1
nrow(choices) == nrow(sprobHA) 
nrow(sprobHA) ==  nrow(sprobLA)
nrow(sprobLA) ==  nrow(sprobHB)
nrow(sprobHB) ==  nrow(sprobLB)

# create data list for JAGS

dat <- list(
  HA = problems$HA ,
  LA = problems$LA , 
  HB = problems$HB ,
  LB = problems$LB , 
  sprobHA = sprobHA , 
  sprobLA = sprobLA ,
  sprobHB = sprobHB , 
  sprobLB = sprobLB ,
  choice = choices ,
  nprob = nprob ,
  nsub = nsub
)

params <- c("mu.alpha", "alpha", "mu.gamma", "gamma", "mu.delta", "delta", "mu.rho", "rho") # parameters

params_init <- function(){
  list("mu.probit.alpha" = rnorm(1, 0, .1) ,
       "mu.probit.gamma" = rnorm(1, 0, .1) ,
       "mu.probit.delta" = rnorm(1, 0, .1) , 
       "mu.probit.rho" = -5,  
       "probit.alpha" = rnorm(nsub, 0, .1) ,
       "probit.gamma" = rnorm(nsub, 0, .1) ,
       "probit.delta" = rnorm(nsub, 0, .1) , 
       "probit.rho" = rep(-5, nsub) )
}


## sample from posterior distributions using MCMC

mfit <- jags.parallel(
  
  data = dat ,
  inits = params_init ,
  parameters.to.save = params ,
  model.file = "code/CPT_hierarchical.txt" ,
  n.chains = 4 ,
  n.iter = 10000 ,
  n.burnin = 5000 ,
  n.thin = 10 ,
  n.cluster = 4 , # compute MCMC chains in parallel
  DIC = FALSE ,
  jags.seed = 56121
  
)

traceplot(mfit)

fits <- mfit$BUGSoutput$summary %>% round(4) %>%  as_tibble(rownames = "parameter")


mu.fits <- fits %>% 
  filter(parameter %in% c("mu.alpha", "mu.gamma", "mu.delta", "mu.rho")) %>% 
  mutate(subject = 0, 
         parameter = c("alpha", "delta", "gamma", "rho")) %>% 

ind.fits <- fits %>% 
  filter(! parameter %in% c("mu.alpha", "mu.gamma", "mu.delta", "mu.rho")) %>% 
  mutate(subject = rep(1:nsub, 4),
         parameter = c(rep("alpha", nsub), rep("delta", nsub),  rep("gamma", nsub), rep("rho", nsub))
  )

mu.weights <- mu.fits %>% 
  select(subject, parameter, mean) %>% 
  pivot_wider(names_from = parameter, values_from = mean) %>% 
  select(-c(alpha, rho)) %>%
  expand_grid(p = seq(0, 1, .01)) %>% # create vector of sampled relative frequencies
  mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) # comput


ind.weights <- ind.fits %>%
  select(subject, parameter, mean) %>% 
  pivot_wider(names_from = parameter, values_from = mean) %>% 
  select(-c(alpha, rho)) %>%
  expand_grid(p = seq(0, 1, .01)) %>% # create vector of sampled relative frequencies
  mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) # comput


mu.weights %>% 
  ggplot(aes(p, w)) +
  scale_x_continuous(breaks = seq(0, 1, length.out = 3)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 3)) +
  labs(x = "p",
       y = "w(p)") +
  geom_line(data=ind.weights, aes(group=subject), linewidth = .5, color = "gray") + 
  geom_line(linewidth = 1, color = "#ff02ff") +
  theme_minimal(base_size = 20)





# pooling --------------------------------------------------------

params <- c("alpha", "gamma", "delta", "rho") # free parameters
params_init <- function(){
  list("alpha" = rbeta(1, 20,20) ,
       "gamma" = rbeta(1, 20,20) ,
       "delta" = rbeta(1, 20,20) , 
       "rho" = rbeta(1, 1, 20))
}


  ## sample from posterior distributions using MCMC

mfit <- jags.parallel(
  
  data = dat ,
  inits = params_init ,
  parameters.to.save = params ,
  model.file = "code/CPT_pooled.txt" ,
  n.chains = 4 ,
  n.iter = 10000 ,
  n.burnin = 5000 ,
  n.thin = 1 ,
  n.cluster = 4 , # compute MCMC chains in parallel
  DIC = FALSE ,
  jags.seed = 61721
  
  )

mfit
traceplot(mfit)

fits <- mfit$BUGSoutput$summary %>% round(4) %>%  as_tibble(rownames = "parameter")

# show results

weights <- fits %>%
  select(parameter, mean) %>%
  pivot_wider(names_from = parameter, values_from = mean) %>% 
  select(-c(alpha, rho)) %>%
  expand_grid(p = seq(0, 1, .01)) %>% # create vector of sampled relative frequencies
  mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) # comput

weights %>% 
  ggplot(aes(p, w)) +
  scale_x_continuous(breaks = seq(0, 1, length.out = 3)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 3)) +
  labs(x = "p",
       y = "w(p)") +
  geom_line(linewidth = 2, color = "#ff02ff") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed") + 
  theme_minimal(base_size = 20)






