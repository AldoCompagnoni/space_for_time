# Try to quantify bias induced by multicollinearity
# Assumptions:
# 1. response variable is dependent on 

# year replication, and number of variables
yr_rep    <- 30
n_vars    <- 2


# simulate  regressions
simulate_reg <- function(ii, 
                         cor_val = -0.6, 
                         yr_rep  = 30, 
                         pop_rep = 2){
  
  # set up sigma and simulate predictors
  mat          <- matrix(cor_val, pop_rep, pop_rep)
  diag(mat)    <- rep(1, pop_rep)
  X            <- MASS::mvrnorm(yr_rep, rep(0,pop_rep), Sigma=mat) 
  
  # simulate dataset 
  y            <- rnorm( 30, 0 + 0.2*X[,1] + 0.4*X[,2], 0.2)
  
  # retrieve beta values in single and multiple regression
  data.frame(
    
    beta1_0 = lm( y ~ X[,1]) %>% coef %>% .[2] %>% as.numeric,
    beta1_i = lm( y ~ X) %>% coef %>% .[2] %>% as.numeric,
    
    beta2_0 = lm( y ~ X[,2]) %>% coef %>% .[2] %>% as.numeric,
    beta2_i = lm( y ~ X) %>% coef %>% .[3] %>% as.numeric

    )
   
}

# simulate retrieving 
betas_bias  <- lapply( 1:1000, bootstrap_cor ) %>% 
                  do.call( rbind, .)

betas_bias %>% 
  # dplyr::select( beta1_0, beta1_i) %>% 
  boxplot
