# Preparation -------------------------------------------------------------


# load packages
pacman::p_load(tidyverse , 
               R2jags ,
               abind , 
               patchwork
)

# load data
choices <- read_rds("data/trial_summaries.rds.bz2")


# Data Selection --------------------------------------------------------------------


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


# Data set: Gloeckner (2012) ----------------------------------------------------------------

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

write_rds(paper, paste("data/cpt/cpt_", tolower(papername), ".rds.bz2", sep =""), compress = "bz2")

# Dataset: XY -------------------------------------------------------------

