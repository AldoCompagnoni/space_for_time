# simulate how estimating beta is affected by bias on:
# 1. the y
# 2. the x
# 3. both x and y

# 1. the y
bias_y <- function(ii){

  set.seed(ii)
  x  <- rnorm(30)
  y  <- rnorm(length(x), x*0.5)
  y1 <- rnorm(length(x), y, 2)
  
  data.frame( beta     = lm( y ~ x) %>% coef %>% .[2],
              beta_att = lm( y1 ~ x) %>% coef %>% .[2] )

}

par(mar=c(3,3,0.1,0.1))
lapply(1:100, bias_y) %>% bind_rows %>% boxplot


# 2. the x
bias_x <- function(ii){
  
  set.seed(ii)
  x  <- rnorm(30)
  y  <- rnorm(length(x), x*0.5)
  x1 <- rnorm(length(x), x, 1)
  
  data.frame( beta     = lm( y ~ x) %>% coef %>% .[2],
              beta_att = lm( y ~ x1) %>% coef %>% .[2] )
  
}

par(mar=c(3,3,0.1,0.1))
lapply(1:100, bias_x) %>% bind_rows %>% boxplot


# 3. both x and y
bias_xy <- function(ii){
  
  set.seed(ii)
  x  <- rnorm(30)
  y  <- rnorm(length(x), x*0.5)
  x1 <- rnorm(length(x), x, 1)
  y1 <- rnorm(length(x), x, 3)
  
  data.frame( beta     = lm( y ~ x) %>% coef %>% .[2],
              beta_att = lm( y1 ~ x1) %>% coef %>% .[2] )
  
}

par(mar=c(3,3,0.1,0.1))
lapply(1:100, bias_xy) %>% bind_rows %>% boxplot
