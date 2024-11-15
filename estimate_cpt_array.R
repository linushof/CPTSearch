# load packages
pacman::p_load(tidyverse, R2jags)
library(dplyr)
library(abind)
library(ggplot2)
library(dplyr)
library(patchwork)
library(purrr)

# load and preprocess data with the pre-processing script
source("preprocessing.R") 

# prepare data for CPT --------------------------------------------------------------------


dat_cpt <- choices %>% 
  
  # get problems where sampling is autonomous, outcomes are gains only, and with up to 2 outcomes per option
  
  filter(type == "free", 
         dom == "Gain", 
         probA3 == 0 & probA4 == 0 & probA5 == 0, 
         probB3 == 0 & probB4 == 0 & probB5 == 0) %>% 
  
  # identify higher and lower outcome for each option (as CPT uses rank-dependent transformations)
  
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
  
  # select only relevant variables
  
  select(index:subject , 
         HA, sprobHA, LA, sprobLA , 
         HB, sprobHB, LB, sprobLB , 
         choice) 


# compute number of problems and subjects in each data set 

nproblem <- dat_cpt %>% 
  group_by(paper) %>% 
  distinct(problem) %>% 
  summarise(nproblem = n()) 

nsubjects <- dat_cpt %>% 
  group_by(paper) %>% 
  distinct(subject) %>% 
  summarise(nsubject = n()) 

# overview of max problems and subjects per paper 
overview <- left_join(nproblem, nsubjects, by=join_by(paper)) %>% 
  mutate(nchoice = nproblem * nsubject)

# check for paper picks
table(dat_cpt$paper)

# function to create individual plot pair of sprob for each paper
create_plot <- function(df, paper_name) {
  p1 <- ggplot(df, aes(x = sprobHA)) +
    geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = paste0("sprobHA - ", paper_name), x = "sprobHA", y = "Frequency") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8),
      axis.text.x = element_text(size = 5)
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25))
  
  p2 <- ggplot(df, aes(x = sprobHB)) +
    geom_histogram(bins = 20, fill = "orange", color = "black", alpha = 0.7) +
    labs(title = paste0("sprobHB - ", paper_name), x = "sprobHB", y = "Frequency") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8),
      axis.text.x = element_text(size = 5)
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25))
  
  # side by side combination of the two plots for the same paper
  p1 | p2
}

# create the plots of sprob
plots <- dat_cpt %>%
  group_split(paper) %>%
  map(~ create_plot(.x, as.character(.x$paper[1])))  # Pass title explicitly for each group

# arrange all plot pairs in a 4x4 grid
overview_plot <- wrap_plots(plots, ncol = 4)

# display the final grid of plots
print(overview_plot)

# pick paper
papername <- "Gloeckner12"
paper <- dat_cpt %>% filter(paper == papername)

# check for NAs in sprob and remove the respective rows
NAcolumns <- paper  %>% filter(is.na(sprobHA) == TRUE | is.na(sprobLA) == TRUE |is.na(sprobHB) == TRUE|is.na(sprobLB) == TRUE)

paper <- setdiff(paper,NAcolumns)

# replace NAs with 0
paper[is.na(paper)] <-0

# get vector with solved problems per participant
nsolvedprob <- paper %>% 
  group_by(subject) %>% 
  summarise(nsolvedprob = n()) 

# merge nsolvedprob in original working table
paper <- merge(paper,nsolvedprob, by = "subject")

# vector of solved problems per participant
nsolvedprob <- nsolvedprob$nsolvedprob

# check duplicates (add only for Frey, Hau, Erev, Wullff)
paper_unique <- paper %>%
  group_by(subject) %>%
  mutate(unique_subject_problem = row_number()) %>%
  ungroup() %>%
  select(-problem) %>% # replace old problem column with unique ID
  rename(problem = unique_subject_problem)

# ensure right order of the columns ((add only for Frey, Hau, Erev, Wullff))
paper <- paper_unique
paper <- paper %>%
  select(1:4, problem, everything()) 

## data that is the same for all subjects (but see ISSUES above)

problems <- paper %>% distinct(problem) 

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
sprobHA+sprobLA
sprobHB+sprobLB == 1
nrow(choices) == nrow(sprobHA) 
nrow(sprobHA) ==  nrow(sprobLA)
nrow(sprobLA) ==  nrow(sprobHB)
nrow(sprobHB) ==  nrow(sprobLB)


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
  
  data = data_list , # data list
  inits = params_init , # creates list of initial values for each chain
  parameters.to.save = params , 
  model.file = "CPT_hierarchical_array.txt" , # model code, see file
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
  mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) 

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
  mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) 

# plot
mu.weights %>% 
  ggplot(aes(p, w)) +
  scale_x_continuous(breaks = seq(0, 1, length.out = 3)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 3)) +
  labs(title = papername, x = "p",
       y = "w(p)") +
  geom_line(data=ind.weights, aes(group=subject), color = "gray") + 
  geom_abline(intercept=0, slope =1, linewidth = 1, "black", linetype = "dashed") +
  geom_line(linewidth = 1, color = "black") +
  theme_classic()+theme(
    plot.title = element_text(
      size = 20,    
      face = "bold",  
      hjust = 0.5
    ))
    

