
design_df <- expand.grid( rep_n = seq(5,30,by=5), 
                          x_r   = seq(1,3,by=0.5),
                          sd_s  = seq(0.05, 0.15, by=0.025) )

# the power of increasing range of x predictor
power_range <- function( ii ){
  
  print( ii )
  
  x <- seq( (-design_df$x_r[ii]),
            design_df$x_r[ii], 
            length.out = design_df$rep_n[ii] )  
  
  sim_lm <- function( sim_i, x, rep_n, sd_s ){
    
    y <- rnorm( rep_n, x * 0.05, design_df$sd_s[ii] )
    lm( y ~ x ) %>% summary %>% .$coefficients %>% .[2,4]
    
  }

  p_vals <- sapply(1:1000, sim_lm, x, design_df$rep_n[ii] )

  design_df[ii,] %>% 
    mutate( power = sum(p_vals < 0.05) / 1000 )
  
}

# 
power_l   <- lapply(1:nrow(design_df), power_range)
power_df  <- power_l %>% bind_rows


power_df %>% 
  subset( x_r == 1 ) %>% 
  ggplot() +
  geom_jitter( aes(rep_n, power,
                   color = as.factor(sd_s),
                   size = 2),
               width = 0.1 ) 


power_df %>% 
  subset( )
  ggplot() +
  geom_jitter( aes(x_r, power,
                  color = rep_n,
                  pch = as.factor(sd_s),
                  size = 2),
               width = 0.1 ) +
  scale_color_viridis_c()

write.csv(power_df, 'PROVA_1.4.2020.csv', row.names=F)