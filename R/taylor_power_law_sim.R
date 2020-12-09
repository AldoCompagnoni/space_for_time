# simulation to see how taylor's power law affects errors in lambda
x1 <- rnorm( 1000, 100, 10 )
x2 <- rnorm( 1000, 1000, 300 )

df1 <- data.frame( t1 = x1,
                  t0 = lag(x1) )

df2 <- data.frame( t1 = x2,
                   t0 = lag(x2) )

log(df1$t1 / df1$t0) %>% sd(na.rm=T)
log(df2$t1 / df2$t0) %>% sd(na.rm=T)

sim_dilute <- function( ii, sampl ){
  
  x  <- rnorm(sampl)  
  x1 <- rnorm(length(x), x )  
  
}
