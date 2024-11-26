# load packages
pacman::p_load(tidyverse)

#load data
choices <- read_rds("data/trial_summaries.rds.bz2")

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

overview <- left_join(nproblem, nsubjects, by=join_by(paper)) %>% 
  mutate(nchoice = nproblem * nsubject)


# use Kellen (2016) data for demonstration (comprehensive data set with no/few issues, see also comment below)
# possible ISSUES with other data sets that need to be fixed : Missings, not all subjects solve the same set of problems, duplicate subject-problem pairings in different data sets (id) within a paper (e.g. Erev)
# none of these issues in Kellen16, however, for other data sets, code below needs to be adapted. 
# Ideally make code flexible to deal with all data sets
# moreover: fitting strategy might be adjusted later to group by switching behavior

paper <- dat_cpt %>% filter(paper == "Kellen16")


# create data objects for JAGS

## data that is the same for all subjects (but see ISSUES above)

problems <- paper %>% distinct(problem, HA, LA,HB, LB)  

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
