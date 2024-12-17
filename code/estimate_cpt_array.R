# Preparation -------------------------------------------------------------

# load pkgs
pacman::p_load(tidyverse, R2jags)
library(viridis)
library(glue)

# Gloeckner (2012) ------------------------------------------------------------

#valid papers
valid_papers <- c("gloeckner12", "hertwig04", "ungemach09", "rakow08")


parameters_df <- data.frame(
  paper = character(),
  subject = numeric(),
  alpha = numeric(),
  gamma = numeric(),
  rho = numeric(),
  delta = numeric()
)

for (p in valid_papers) {
  paper <- read_rds(glue("data/PreprocessedPaperData/cpt_{p}.rds.bz2"))
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

  ## create data object for JAGS
  # create Array 

  preArrayData <- subset(paper, select = - c(index, problem,id))

  preArrayData <- preArrayData %>%
    select(-paper, paper)

  split_data <- split(preArrayData, preArrayData$subject)

  #Get max number of problems 
  problems <- paper %>% distinct(problem) 
  nprob <- nrow(problems) 

  #Get number of subjects
  nsub <-length(unique(paper$subject)) 

  #number of columns without subject
  num_cols <- ncol(preArrayData) - 1 

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

  # traceplot(mfit) # can be cumbersome with many subjects

  # summary of posterior distributions

  estimates <- mfit$BUGSoutput$summary %>% round(4) %>%  as_tibble(rownames = "parameter")
  posteriors <- mfit$BUGSoutput$sims.matrix %>% as_tibble()

  papername <- unique(paper$paper)
  write_rds(estimates, paste("data/PostEst/cpt_", tolower(papername), "_estimates.rds.bz2", sep =""), compress = "bz2")
  write_rds(posteriors, paste("data/PostEst/cpt_", tolower(papername), "_posteriors.rds.bz2", sep =""), compress = "bz2")


  # Create summary dataframe with estimates per subject
  parameter_data <- estimates %>%
   mutate(subject = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", parameter)),  # Extrahiert Subjekt-ID
          param_name = gsub("\\[.*", "", parameter)) %>%                  # Extrahiert Parametername
   select(subject, param_name, mean) %>%                                 # Relevante Spalten
    pivot_wider(names_from = param_name, values_from = mean)              # Parameter als Spalten


  #plots
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
  # select(-c(alpha, rho)) %>%
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

  ind.weights <- left_join(ind.weights, switch_rate_data, by = "subject")

  #Weighted Probability plots
  median_switch_rate <- median(ind.weights$avg_switchrate, na.rm = TRUE)
  prob_weight_plot <- mu.weights %>% 
    ggplot(aes(p, w)) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25), 
      labels = scales::number_format(accuracy = 0.01)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = scales::number_format(accuracy = 0.01) 
    ) +
    labs(
      title = papername,
      x = "p",
      y = "w(p)"
    ) +
    geom_line(data = ind.weights, aes(group = subject, color = avg_switchrate), alpha = 0.8) +
    scale_color_viridis_c(
      option = "plasma",  
      name = "Average Switching Rate",  
      breaks = c(0, 0.25, 0.5, 0.75, 1),  
      labels = scales::number_format(accuracy = 0.01),
      limits = c(0, 1) 
    ) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, color = "black", linetype = "dashed") +
    geom_line(linewidth = 1.2, color = "black") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "bottom",  
      legend.direction = "horizontal",  
      legend.box = "horizontal",  
      legend.title.align = 0.5, 
      legend.text = element_text(size = 10, color = "black"),  
      legend.title = element_text(size = 12, face = "plain", hjust = 0.5, vjust =0.8),  
      legend.key.width = unit(2, "cm"),  
      legend.spacing.y = unit(2.5, "cm"),  
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA) ,
    )
  # Calculation median switching rate
  median_switch_rate <- median(ind.weights$avg_switchrate, na.rm = TRUE)


  # Plot with median switching rate
  prob_weight_plot <- mu.weights %>% 
    ggplot(aes(p, w)) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25), 
      labels = scales::number_format(accuracy = 0.01)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = scales::number_format(accuracy = 0.01) 
    ) +
    labs(
      title = papername,
      x = "p",
      y = "w(p)",
      caption = paste("Median Switching Rate:", round(median_switch_rate,2))
    ) +
    geom_line(data = ind.weights, aes(group = subject, color = avg_switchrate), alpha = 0.8) +
    scale_color_viridis_c(
      option = "plasma",  
      name = "Average Switching Rate",  
      breaks = c(0, 0.25, 0.5, 0.75, 1),  
      labels = scales::number_format(accuracy = 0.01),
      limits = c(0, 1) 
    ) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, color = "black", linetype = "dashed") +
    geom_line(linewidth = 1.2, color = "black"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "bottom",  
      legend.direction = "horizontal",  
      legend.box = "horizontal",  
      legend.title.align = 0.5, 
      legend.text = element_text(size = 10, color = "black"),  
      legend.title = element_text(size = 12, face = "plain", hjust = 0.5, vjust = 0.8),  
      legend.key.width = unit(2, "cm"),  
      legend.spacing.y = unit(2.5, "cm"),  
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      plot.caption = element_text(size = 12, hjust = 0.5, color = "black")  
    )


  prob_weight_plot

  ggsave(filename = paste0("plots/ProbWeighting/ProbWeighting_", papername, ".png"), plot = prob_weight_plot)

  #regression data preparations

  data_regression <- fits %>%
    mutate(
      subject = as.integer(str_extract(parameter, "\\d+")), 
      parameter = str_remove(parameter, "\\[\\d+\\]")        
    ) %>%
    select(subject, parameter, mean) %>%                     
    pivot_wider(
      names_from = parameter,                                
      values_from = mean                                     
    )

  data_regression <- left_join(data_regression, switch_rate_data, by = "subject")

  #mean_values <- data_regression %>%
   # summarise(
    #  mean_alpha = mean(alpha, na.rm = TRUE),
     # mean_delta = mean(delta, na.rm = TRUE),
    #  mean_gamma = mean(gamma, na.rm = TRUE),
    #  mean_rho = mean(rho, na.rm = TRUE)
    #)
  
  #mean_values <- mean_values %>%
   # mutate(paper = p) %>%
    #select(paper, everything())  
  
  
 # parameters_df <- bind_rows(parameters_df, mean_values)
  
  
  parameters_df_paper<- data_regression %>%
    select(subject, alpha, gamma, rho, delta) %>%
    mutate(paper = p)  
  
  # Saving parameters for each paper in dataframe
  parameters_df <- rbind( parameters_df, parameters_df_paper)
  
}
#Distribution Gamma
#Logistic Regression
logit_model <- glm(cat_switch ~ gamma + delta, 
                   data = data_regression, 
                   family = binomial)
summary(logit_model)

exp(coef(logit_model))
pR2(logit_model)

#Linear Regression
lin_model_gamma <- lm(gamma ~ avg_switchrate, 
                   data = data_regression)
lin_model_delta <- lm(delta ~ avg_switchrate, 
                data = data_regression)
lin_model_alpha <- lm(alpha ~ avg_switchrate, 
                data_regression)


#Summary
summary(lin_model_gamma)
summary(lin_model_delta)
summary(lin_model_alpha)

# Plot alpha values across papers
ggplot(parameters_df, aes(x = paper, y = alpha)) +
  geom_point(size = 4, color = "skyblue") +  
  labs(
    title = "alpha values across papers",
    x = "Paper",
    y = "Alpha"
  ) +
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12)  
  )

# Plot gamma values across papers
ggplot(parameters_df, aes(x = paper, y = gamma)) +
  geom_point(size = 4, color = "skyblue") +  
  labs(
    title = "gamma values across papers",
    x = "Paper",
    y = "Gamma"
  ) +
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12)  
  )

# Plot rho values across papers
ggplot(parameters_df, aes(x = paper, y = rho)) +
  geom_point(size = 4, color = "skyblue") +  
  labs(
    title = "rho values across papers",
    x = "Paper",
    y = "Rho"
  ) +
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12)  
  )

# Plot delta values across papers
ggplot(parameters_df, aes(x = paper, y = delta)) +
  geom_point(size = 4, color = "skyblue") +  
  labs(
    title = "delta values across papers",
    x = "Paper",
    y = "Delta"
  ) +
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12)  
  )




