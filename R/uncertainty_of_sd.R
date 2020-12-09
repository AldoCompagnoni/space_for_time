# estimate SE of the standard deviation

# model of the mean
x <- rnorm(100)
summary(lm(x ~ 1))$coefficients[1,2]

# stan model.



# regressuion (se of the residual sigma)
x <- rnorm(100)
y <- rnorm(length(x), x, sd=0.43)
