
# preparation -------------------------------------------------------------


# non array ---------------------------------------------------------------

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


## weighting function for individual level
ind.fits <- fits %>% 
  filter(! parameter %in% c("mu.alpha", "mu.gamma", "mu.delta", "mu.rho")) %>% 
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


# array -------------------------------------------------------------------



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


