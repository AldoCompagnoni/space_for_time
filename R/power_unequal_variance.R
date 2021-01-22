# unequal variance
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggthemes)
library(cowplot)
library(MASS)


# 1. power analysis -----------------------------------
x_rep_grid <- expand.grid( x_ran = 2,
                           rep_x = seq(20,50,by=2),
                           uneq  = c(T, F) )

# bootstrap range
bootstrap_beta <- function(ii){
  
  # set the bounds
  x_bound <- x_rep_grid$x_ran[ii]
  
  x <- seq( -x_bound, x_bound,
            length.out = x_rep_grid$rep_x[ii] )
  
  boostrap_reg <- function(boot_i, unequal ){
    
    if( unequal ){
    
      y1 <- rnorm(x_rep_grid$rep_x[ii]/2, 0 + x * 0.05, 0.2)
      y2 <- rnorm(x_rep_grid$rep_x[ii]/2, 0 + x * 0.05, 0.4)
      y  <- c( y1, y2 )
    
    } else{
      
      y <- rnorm(x_rep_grid$rep_x[ii], 0 + x * 0.05, 0.2)
      
    }
    
    lm( y ~ x ) %>% coef %>% .[2]
    
  }
  
  beta_v <- sapply(1:1000, boostrap_reg, x_rep_grid$uneq)
  
  data.frame( beta_mean = mean(beta_v),
              beta_sd   = sd(beta_v,na.rm=T),
              beta_min  = min(beta_v),
              beta_max  = max(beta_v) )
  
}

# recovered betas
beta_sim <- lapply(1:nrow(x_rep_grid), bootstrap_beta)

# data frame for plotting
beta_df  <- beta_sim %>% 
  bind_rows %>% 
  bind_cols( x_rep_grid ) %>% 
  mutate( beta_range = beta_max - beta_min ) %>% 
  mutate( rep_x = as.factor(rep_x) )

# plot  
ggplot(beta_df) + 
  geom_point( aes( x = as.factor(uneq),
                  y = beta_mean,
                  color = rep_x) ) +
  theme_minimal() + 
  labs( y     = expression('sd('*beta*')'),
        x     = 'Range of climate anomalies',
        color = 'Replicates' ) 
  

beta_df %>% boxplot(beta_mean ~ data=.)