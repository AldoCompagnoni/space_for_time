# power analysis: is it better to have more data points, 
# or sample more "climate space"? 
# answer: the latter
# 1. bootstrap range of x by sample size
# 2. power analysis 
# 3. cross-validation
# 4. cross-validation + correlation
# 5. hypothetical case: 5 years, 100 sites
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggthemes)
library(cowplot)
library(MASS)

# read in compadre
load('C:/CODE/plant_review_analyses/Data/COMPADRE_v.5.0.0.RData')

# average lambda by population
id <- which( compadre$metadata$MatrixComposite == 'Individual' )

# keep these SppAuthor/MatrixPopulation
keep_pop <- compadre$metadata[id,] %>% 
              subset( !is.na(MatrixStartYear) ) %>% 
              count(SpeciesAuthor, MatrixPopulation, MatrixStartYear) %>%  
              subset( n > 5 )  


keep_ids <- intersect( which(compadre$metadata$SpeciesAuthor %in% keep_pop$SpeciesAuthor),
                       which(compadre$metadata$MatrixPopulation %in% keep_pop$MatrixPopulation) )


calc_lam <- function(ii){
  
  mat <- compadre$mat[[keep_ids[ii]]]$matA
  
  data.frame( SpeciesAuthor    = compadre$metadata$SpeciesAuthor[keep_ids[ii]],
              MatrixPopulation = compadre$metadata$MatrixPopulation[keep_ids[ii]],
              log_lambda       = log(Re(eigen(mat)$value[1])) )
  
}

# store lambdas in one data frame
lam_df <- lapply(1:length(keep_ids), calc_lam) %>% bind_rows
  
# sd of lambdas
sd_df <- lam_df %>% 
          group_by( SpeciesAuthor,MatrixPopulation ) %>% 
          summarise( lam_sd = sd(log_lambda) )

# note: almost exactly the variance found in the plant review dataset!
sd_df$lam_sd %>% hist
sd_df$lam_sd %>% mean
sd_df$lam_sd %>% median


# 1. bootstrap range of x by sample size ---------------

# univariate normal
bootstrap_univ <- function(ii, n_rep, x_val){

  set.seed(ii)
  
  if( x_val == 'random'){
    x <- rnorm(n_rep)   
  } else {
    x <- seq(-1,1,length.out=n_rep)
  }
  
  y <- rnorm(n_rep, 0 + x * 0.05, 0.2)
  
  data.frame( 
              # beta = lm( y ~ x ) %>% coef %>% .[2],
              low  = as.numeric( range(x) )[1],
              high = as.numeric( range(x) )[2] )
  
}


# multivariate normal --------------------------

yr_rep     <- 5
pop_rep    <- 5

# two sites, uncorrelated anomalies
range_uncorr <- list()

bootstrap_cor <- function(ii,cor_val, yr_rep, pop_rep){
  
  # set up sigma
  mat          <- matrix(cor_val, pop_rep, pop_rep)
  diag(mat)    <- rep(1, pop_rep)
  
  MASS::mvrnorm(yr_rep, rep(0,pop_rep), Sigma=mat) %>% 
    as.numeric %>% 
    range
}

# multivariate, strongly correlated
range_cor    <- lapply(1:10000, bootstrap_cor, 0.9, 5, 2) %>% 
                    do.call( rbind, .)

# multivariate normal
range_uncorr <- lapply(1:10000, bootstrap_cor, 0, 5, 20) %>% 
                    do.call( rbind, .)
                    

# univariate normal 
rang_univ    <- lapply(1:10000, bootstrap_univ, rep_n, 'random') %>% 
                    do.call( rbind, .)

# compare ranges
apply(rang_univ,  2, mean)
apply(range_cor, 2, mean)
apply(range_uncorr, 2, mean)


# make response surface
surf_df <- expand.grid( pop_rep = c(2:10),
                        cor     = c(0.95, 0.9, 0.85, 0.8) )

range_by_rep <- list()

for(ii in 1:nrow(surf_df) ){
  
  range_uncorr <- lapply(1:10000, bootstrap_cor, 
                         surf_df$cor[ii], 5, surf_df$pop_rep[ii]) %>% 
                      do.call( rbind, .)
  
  df_out <- data.frame( cor     = surf_df$cor[ii],
                        pop_rep = surf_df$pop_rep[ii],
                        range   = apply(range_uncorr, 2, mean) %>% abs %>% mean )
  
  range_by_rep[[ii]]   <- df_out
  
}

range_df <- bind_rows( range_by_rep )

# change in range of x values with spatial replication
pn <- ggplot(range_df) +
        geom_tile( aes( x    = as.factor(pop_rep),
                        y    = cor,
                        fill = range ) ) +
        scale_fill_viridis() +
        labs( x    = "Number of populations",
              y    = "Climatic correlation",
              fill = "Anomalies range")

range_df



# 2. power analysis -----------------------------------
x_rep_grid <- expand.grid( x_ran = seq(1,3,length.out=20),
                           rep_x = seq(5,30,by=5) )

# bootstrap range
bootstrap_beta <- function(ii){
  
  # set the bounds
  x_bound <- x_rep_grid$x_ran[ii]
  
  x <- seq( -x_bound, x_bound,
            length.out = x_rep_grid$rep_x[ii] )
  
  boostrap_reg <- function(boot_i){
    y <- rnorm(x_rep_grid$rep_x[ii], 0 + x * 0.05, 0.2)
    lm( y ~ x ) %>% coef %>% .[2]
  }
  
  beta_v <- sapply(1:10000, boostrap_reg)
  
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
              subset( x_ran < 2.6 ) %>% 
              mutate( rep_x = as.factor(rep_x) )

# # write out results
# write.csv(beta_df, 'C:/CODE/plant_review_analyses/results/biases/power_analysis.csv',
#           row.names=F)

# read results
# beta_df <- read.csv('C:/CODE/plant_review_analyses/results/biases/power_analysis.csv')


p1 <- ggplot(beta_df) + 
  geom_point( aes( x = x_ran,
                  y = beta_sd,
                  color = rep_x) ) +
  scale_color_colorblind() +
  theme_minimal() + 
  labs( y     = expression('sd('*beta*')'),
        x     = 'Range of climate anomalies',
        color = 'Replicates' ) 
  

p2 <- ggplot(beta_df) + 
  geom_point( aes( x = rep_x,
                  y = beta_sd,
                  color = x_ran) ) +
  scale_color_viridis() +
  theme_minimal() + 
  labs( y     = expression('sd('*beta*')'),
        color = 'clim. range',
        x     = 'Replicates' ) 

p1 <- ggplot(beta_df) + 
  geom_point( aes( x = x_ran,
                   y = beta_sd,
                   color = rep_x) ) +
  scale_color_colorblind() +
  theme_minimal() + 
  labs( y     = expression('sd('*beta*')'),
        x     = 'Range of climate anomalies',
        color = 'Replicates' ) 


p3 <- ggplot(beta_df) + 
  geom_point( aes( x = rep_x,
                   y = beta_range,
                   color = x_ran) ) +
  scale_color_viridis() +
  theme_minimal() + 
  labs( y     = expression(beta*' range'),
        color = 'clim. range',
        x     = 'Replicates' ) 

p4 <- ggplot(beta_df) + 
  geom_point( aes( x = x_ran,
                   y = beta_range,
                   color = x_ran) ) +
  scale_color_viridis() +
  theme_minimal() + 
  labs( y     = expression(beta*' range'),
        color = 'Replicates',
        x     = 'Range of climate anomalies' ) 


pa_plot <- plot_grid( p1, p2, 
                      labels = 'AUTO',
                      label_size = 18,
                      align='h', 
                      nrow = 1, ncol = 2 )

ggsave('power_analysis.png', 
       pa_plot,
       width = 6.3, 
       height = 3.15)


# single plot with DISCRETE scale
beta_df %>% 
  mutate( rep_x = as.factor(rep_x) ) %>% 
  ggplot() + 
  geom_point( aes( x = x_ran,
                   y = beta_sd,
                   color = rep_x),
              size = 3 ) +
  scale_color_colorblind() +
  theme_minimal() + 
  labs( y     = expression('sd('*beta*')'),
        x     = 'Range of climate anomalies',
        color = 'Replicates' ) +
  ggsave('power_analysis_reps.tiff', 
         width = 6.3, 
         height = 6.3)


#3. cross-validation ------------------------------------------

x_rep_grid <- expand.grid( x_ran = seq(1,3,length.out=20),
                           rep_x = seq(5,30,by=5) )

# bootstrap range
crossval_sim <- function(ii, sd_p, b_p){
  
  # set the bounds
  x_bound <- x_rep_grid$x_ran[ii]
  
  x <- seq( -x_bound, x_bound,
            length.out = x_rep_grid$rep_x[ii] )
  
  
  sim <- function( ss, ii ){
    
    full_df <- data.frame( yr = c(1:x_rep_grid$rep_x[ii]),
                           x  = x,
                           y  = rnorm(x_rep_grid$rep_x[ii], 
                                      0 + x * b_p, 
                                      sd_p) 
                          )
    
    # crossvalidation
    crossval <- function( zz, ii ){
      
      train_df <- subset(full_df, !(yr %in% zz) )
      test_df  <- subset(full_df, yr == zz)
      
      pred_y   <- predict( lm( y ~ x, data = train_df ),
                           newdata = test_df )
      pred_0   <- predict( lm( y ~ 1, data = train_df ),
                           newdata = test_df )
      
      data.frame( pred_y = (pred_y - test_df$y)^2,
                  pred_0 = (pred_0 - test_df$y)^2 )
      
    }
    
    rmse <- lapply(1:x_rep_grid$rep_x[ii], crossval, ii) %>% 
              bind_rows %>% 
              apply(2, function(x) x^2) %>% 
              apply(2, mean) %>% 
              sqrt
      
    data.frame( rep    = x_rep_grid$rep_x[ii],
                ran    = x_rep_grid$x_ran[ii],
                pred_y = rmse[1],
                pred_0 = rmse[2] )
                
  }
  
  lapply(1:100, sim, ii) %>% bind_rows
  
}

cv_df1 <- lapply(1, crossval_sim, 0.2, 0)

cv_df1[[1]] %>% 
  mutate( best = pred_y < pred_0 ) %>% 
  mutate( best = as.numeric(best) ) %>% 
  group_by( rep, ran ) %>% 
  summarise( prop = sum(best)/100 )



cv_df1 <- lapply(1:50, crossval_sim)
cv_df2 <- lapply(51:100, crossval_sim)
cv_df3 <- lapply(101:120, crossval_sim)
cv_df  <- bind_rows(cv_df1, cv_df2, cv_df3) %>% 
            mutate( best = pred_y < pred_0 ) %>% 
            mutate( best = as.numeric(best) ) %>% 
            group_by( rep, ran ) %>% 
            summarise( prop = sum(best)/30 )  
  
  
write.csv( cv_df, 'results/crossval_1.31.2020.csv',
           row.names = F)



# 4. cross-validation + correlation -------------------------

x_corr_rep <- expand.grid( corr   = c(-0.5, 0, 0.5, 0.95),
                           site_n = c(1, seq(5,25,by=5) ),
                           yr     = c(3, 5) ) %>% 
                bind_rows( data.frame( corr   = c(.95,.95),
                                       site_n = c(1,1),
                                       yr     = c(20,25) ) )

# Produce sigma based
make_sigma <- function(site_n, corr){
  
  if( site_n == 1 ){
    sig <- matrix( c(1,corr,corr,1), 2, 2)
  }
  
  if( site_n > 1 ){
    reps  <- diag(site_n)
    reps[upper.tri(reps)] <- 0.98
    reps[lower.tri(reps)] <- 0.98
    
    off_d <- matrix( corr, site_n, site_n )
    
    sig <- rbind( cbind( reps,  off_d ),
                  cbind( off_d, reps  ) )
    
  }
  
  sig
  
}


# bootstrap range
bootstrap_beta <- function(ii){
  
  # var-covar of 
  sig <- make_sigma( x_corr_rep$site_n[ii],
                     x_corr_rep$corr[ii] )
  
  boostrap_reg <- function(boot_i){
    
    # set the bounds
    x <- mvrnorm( x_corr_rep$yr[ii], 
                  rep(0, x_corr_rep$site_n[ii] * 2),
                  Sigma = sig ) %>% 
          as.numeric
    
    y   <- rnorm( length(x) , 0 + as.numeric(x) * 0.043, 0.16) #0.043
    
    mod <- lm( y ~ x )
    
    data.frame( beta = mod %>% coef %>% .[2],
                sig  = mod %>% summary %>% .$coefficient %>% .[2,4]
               )
    
  }
  
  beta_df <- lapply(1:1000, boostrap_reg) %>% bind_rows
  
  data.frame( beta_mean = mean(beta_df$beta),
              beta_sd   = sd(beta_df$beta,na.rm=T),
              beta_min  = min(beta_df$beta),
              beta_max  = max(beta_df$beta),
              type_m    = sum(beta_df$beta > (0.043*2) ) / 1000 ,
              p_wrong_s = sum(beta_df$beta < 0)   / 1000 ,
              sig       = sum(beta_df$sig < 0.05) / 1000 )
  
}

# recovered betas
beta_sim <- lapply(1:nrow(x_corr_rep), bootstrap_beta)

# save.image( 'results/power_analysis.Rdata' )

# data frame for plotting
beta_df  <- beta_sim %>% 
  bind_rows %>% 
  bind_cols( x_corr_rep ) %>% 
  mutate( beta_range = beta_max - beta_min ) %>% 
  mutate( site_n = site_n * 2 ) %>% 
  mutate( rep_x = as.factor(site_n) )

# Statistical power 
ggplot(beta_df) + 
  geom_jitter( aes( x = rep_x,
                   y = sig,
                   # y = beta_sd,
                   color = as.factor(corr),
                   pch = as.factor(yr) ),
              size= 4.5,
              width = 0.1 ) +
  geom_hline( aes(yintercept=0.8),
              lty = 2 ) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, 1 ) +
  labs( y     = 'Power',
        color = 'Correlation',
        x     = 'Reps (spatial)',
        pch   = 'Reps (year)') +
  theme( text = element_text( size = 20 ) ) +
  ggsave( 'results/sig_by_corr_rep_effect_0.043.tiff',
          width=6.3, height=6.3, compression='lzw')

# Wrong sign
ggplot(beta_df) + 
  geom_jitter( aes( x = rep_x,
                    y = p_wrong_s,
                    color = as.factor(corr),
                    pch = as.factor(yr) ),
               size= 4.5,
               width = 0.1 ) +
  geom_hline( aes(yintercept=0.8),
              lty = 2 ) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, 1 ) +
  labs( y     = 'Power',
        color = 'Correlation',
        x     = 'Reps (spatial)',
        pch   = 'Reps (year)') +
  theme( text = element_text( size = 20 ) ) +
  ggsave( 'results/wrong_s_effect_0.043.tiff',
          width=6.3, height=6.3, compression='lzw')


# Standard deviation of beta
ggplot(beta_df) + 
  geom_jitter( aes( x = rep_x,
                    y = beta_sd,
                    color = as.factor(corr),
                    pch = as.factor(yr) ),
              size= 4.5,
              width = 0.1) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, max(beta_df$beta_sd) ) +
  labs( y     = expression('sd('*beta*')' ),
        color = 'Correlation',
        x     = 'Reps (spatial)',
        pch   = 'Reps (year)') +
  theme( text = element_text( size = 20 ) ) +
  ggsave( 'results/sd_beta_by_corr_rep_effect_0.043.tiff',
          width=6.3, height=5, compression='lzw')



p1 <- ggplot(beta_df) + 
  geom_point( aes( x = x_ran,
                   y = beta_sd,
                   color = rep_x) ) +
  scale_color_colorblind() +
  theme_minimal() + 
  labs( y     = expression('sd('*beta*')'),
        x     = 'Range of climate anomalies',
        color = 'Replicates' ) 


# 5. hypothetical cases: 3 years, 100 sites ----------------------

set.seed (11)

x <- rnorm( 3 )

rep_reg <- function(ii, n_rep){
  
  set.seed( ii + 1000 )
  
  # random site specific slope
  x1 <- rnorm(1, 0.05, 0.05)
  y  <- rnorm( n_rep, x * x1, 0.2 )
  data.frame( x = x,
              y = y )
}

df <- lapply(1:100, rep_reg, 3) %>% bind_rows

ggplot(df) +
  geom_point(aes(x,y)) +
  geom_smooth(aes(x,y), method='lm') +
  geom_jitter( aes(x,y), 
               width = 0.05 ) +
  theme_minimal() + 
  labs( y = expression('log('*lambda*')'),
        x = 'Climate anomaly' ) #+
  theme( text = element_text( size = 20 ) ) +
  ggsave( 'results/100_sites_random_slope.tiff',
          width=6.3, height=5, compression='lzw')

lm( y ~ x, data = df ) %>% coef
lm( y ~ x, data = df ) %>% summary


# power analysis for the random beta model
power_ran_beta <- function( ii, n_rep, n_yrs ){
  
  # generate the three climatic anomalies
  set.seed( ii )
  x <- rnorm( n_yrs )

  # generate the N replicate values  
  rep_reg <- function(ii, n_yrs){
    
    set.seed( ii )
    
    # random site specific slope
    x1 <- rnorm(1, 0.05, 0.05)
    y  <- rnorm( n_yrs, x * x1, 0.2 )
    data.frame( x = x,
                y = y )
  }
  
  # data frame
  df <- lapply(1:n_rep, rep_reg, n_yrs) %>% bind_rows

  # put out p value and coefficient
  data.frame( 
    p = lm( y ~ x, data = df ) %>% 
          summary %>% 
          .$coefficients %>% 
          .[2,4],
    c = lm( y ~ x, data = df ) %>% 
          summary %>% 
          .$coefficients %>% 
          .[2,1]
    )  
  
}

# design matrix 
design_df <- expand.grid( n_rep = c(25,50,75,100),
                          n_yrs = c(3:5) )

# simulate the design matrix
sim_design <- function( des_i, power_ran_beta ){
  
  raw_df <- lapply(1:1000, 
                power_ran_beta, 
                design_df$n_rep[des_i], 
                design_df$n_yrs[des_i]) %>% 
            bind_rows
  
  data.frame( p     = sum(raw_df$p < 0.05) / 1000,
              n_rep = design_df$n_rep[des_i],
              n_yrs = design_df$n_yrs[des_i] )
  
}

# Power analysis with "random beta"
power_ran_b_df <- lapply( 1:12, sim_design, power_ran_beta) %>% 
                    bind_rows %>% 
                    mutate( n_yrs = as.factor(n_yrs) )


ggplot( power_ran_b_df ) +
  geom_point( aes(n_rep,p,
                  group = n_yrs,
                  color = n_yrs),
              ) +
  labs( y = "Power",
        x = "Spatial replicates",
        color = 'Rep (yrs)') +
  theme_minimal() +
  ggthemes::scale_color_colorblind() +
  ggsave( 'results/power_random_slope.tiff',
          width=3.15, height=2.5, compression='lzw')

#  power analyses
p_25  <- lapply(1:10000, power_ran_beta, 25, 5) %>% bind_rows
p_50  <- lapply(1:10000, power_ran_beta, 50, 5) %>% bind_rows
p_75  <- lapply(1:10000, power_ran_beta, 75, 5) %>% bind_rows
p_100 <- lapply(1:10000, power_ran_beta, 100, 5) %>% bind_rows


# power analysis 
p_v <- c( sum(p_25$p  < 0.05) / 1000,  sum(p_50$p  < 0.05) / 1000,
          sum(p_75$p  < 0.05) / 1000,  sum(p_100$p < 0.05) / 1000,
          sum(p_125$p  < 0.05) / 1000, sum(p_150$p < 0.05) / 1000,
          sum(p_175$p  < 0.05) / 1000, sum(p_200$p < 0.05) / 1000
          )

# coefficients 
c_v <- c( mean(p_25$c),  mean(p_50$c),  mean(p_75$c),  mean(p_100$c), 
          mean(p_125$c), mean(p_150$c), mean(p_175$c), mean(p_200$c) )

# coefficients 
c_sd_v <- c( sd(p_25$c),  sd(p_50$c),  sd(p_75$c),  sd(p_100$c), 
             sd(p_125$c), sd(p_150$c), sd(p_175$c), sd(p_200$c) )


# 6. SAME SLOPE: 3 years, 100 sites; 100 years, 1 site ----------------------


seed_i  <- 13

set.seed( seed_i )
x_100   <- rnorm( 99 )
x_3     <- rnorm( 3 )

df_tmp  <- data.frame( x = x_100,
                       y = rnorm(length(x_100), x_100 * 0.05, 0.16 ),
                       type = 'Temporal replication',
                       stringsAsFactors = F )


sim_lam_rep <- function( ii ){
    
  data.frame( x = x_3,
              y = rnorm(3, x_3 * 0.05, 0.16 ),
              type = 'Mostly spatial replication',
              stringsAsFactors = F )
  
} 

df_spat <- lapply(1:33, sim_lam_rep) %>% bind_rows

# 
bind_rows( df_tmp, df_spat ) %>% 
  ggplot(  ) +
  geom_point( aes(x,y) ) +
  geom_smooth( aes(x,y),
               method = 'lm',
               color = 'black') +
  facet_wrap( ~ type ) +
  labs( x = 'Climatic anomaly',
        y = expression('Log('*lambda*')') ) +
  theme_minimal() +
  ggsave( 'results/spatial_versus_temporal_rep.tiff',
          width = 6.3, height = 4, compression = 'lzw' )


# 7: Correlated errors --------------------------------------------------

x_corr_rep <- expand.grid( corr  = c(0,0.5,0.95),
                           site_n = c(1, seq(5,30,by=5) ),
                           yr     = c(3, 5) )

# Produce sigma based
make_sigma <- function(site_n, corr){
  
  if( site_n == 1 ){
    sig <- matrix( c(1,corr,corr,1), 2, 2)
  }
  
  if( site_n > 1 ){
    reps  <- diag(site_n)
    reps[upper.tri(reps)] <- 0.98
    reps[lower.tri(reps)] <- 0.98
    
    off_d <- matrix( corr, site_n, site_n )
    
    sig <- rbind( cbind( reps,  off_d ),
                  cbind( off_d, reps  ) )
    
  }
  
  sig
  
}


# bootstrap range
bootstrap_beta <- function(ii){
  
  # var-covar of climate covariates
  sig_x <- make_sigma( x_corr_rep$site_n[ii],
                       0.98 )
  
  # var-covar of climate covariates
  sig_y <- matrix( x_corr_rep$corr[ii] * (0.16^2),
                   x_corr_rep$site_n[ii]*2, 
                   x_corr_rep$site_n[ii]*2 )
  diag(sig_y) <- 0.16^2
  
  # boostrap
  boostrap_reg <- function(boot_i){
    
    # set the bounds
    x     <- mvrnorm( x_corr_rep$yr[ii], 
                  rep(0, x_corr_rep$site_n[ii] * 2),
                  Sigma = sig_x ) %>% 
                as.numeric
    
    yhat  <- as.numeric(x) * 0.043
    y_eta <- mvrnorm( x_corr_rep$yr[ii], 
                      rep(0, x_corr_rep$site_n[ii] * 2),
                      Sigma = sig_y ) %>% 
                as.numeric
    
    y   <- yhat + y_eta
    
    mod <- lm( y ~ x )
    
    data.frame( beta = mod %>% coef %>% .[2],
                sig  = mod %>% summary %>% .$coefficient %>% .[2,4]
              )
    
  }
  
  beta_df <- lapply(1:1000, boostrap_reg) %>% bind_rows
  
  data.frame( beta_mean = mean(beta_df$beta),
              beta_sd   = sd(beta_df$beta,na.rm=T),
              beta_min  = min(beta_df$beta),
              beta_max  = max(beta_df$beta),
              type_m    = sum(beta_df$beta > (0.043*2) ) / 1000 ,
              p_wrong_s = sum(beta_df$beta < 0)   / 1000 ,
              sig       = sum(beta_df$sig < 0.05) / 1000 )
  
}

# recovered betas
beta_sim <- lapply(1:nrow(x_corr_rep), bootstrap_beta)

# data frame for plotting
beta_df  <- beta_sim %>% 
  bind_rows %>% 
  bind_cols( x_corr_rep ) %>% 
  mutate( beta_range = beta_max - beta_min ) %>% 
  mutate( site_n = site_n * 2 ) %>% 
  mutate( rep_x = as.factor(site_n) )


# Statistical power 
ggplot(beta_df) + 
  geom_jitter( aes( x = rep_x,
                    y = sig,
                    color = as.factor(corr),
                    pch = as.factor(yr) ),
               size= 4.5,
               width = 0.1 ) +
  geom_hline( aes(yintercept=0.8),
              lty = 2 ) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, 1 ) +
  labs( y     = 'Power',
        color = 'Correlation',
        x     = 'Reps (spatial)',
        pch   = 'Reps (year)') +
  theme( text = element_text( size = 20 ) ) 


ggplot(beta_df) + 
  geom_jitter( aes( x = rep_x,
                    y = p_wrong_s,
                    color = as.factor(corr),
                    pch = as.factor(yr) ),
               size= 4.5,
               width = 0.1 ) +
  geom_hline( aes(yintercept=0.8),
              lty = 2 ) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, 1 ) +
  labs( y     = 'Power',
        color = 'Correlation',
        x     = 'Reps (spatial)',
        pch   = 'Reps (year)') +
  theme( text = element_text( size = 20 ) ) 



ggplot(beta_df) + 
  geom_jitter( aes( x = rep_x,
                    y = beta_sd,
                    color = as.factor(corr),
                    pch = as.factor(yr) ),
               size= 4.5,
               width = 0.1 ) +
  geom_hline( aes(yintercept=0.8),
              lty = 2 ) +
  scale_color_colorblind() +
  theme_minimal() + 
  ylim( 0, 0.3 ) +
  labs( y     = 'Power',
        color = 'Correlation',
        x     = 'Reps (spatial)',
        pch   = 'Reps (year)') +
  theme( text = element_text( size = 20 ) ) 
