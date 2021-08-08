
library(MASS)

set.seed(1337)
n <- 100

mu_S = rep(0, 5) # c(2,4,6,8,10)
sigma_S = diag(1,nrow = 5,ncol = 5)
S = mvrnorm(n, mu = mu_S, Sigma = sigma_S)

bets <- matrix(c(-2,6,3,0,0))
gams <- matrix(c(0,0,5,1,-10))

y <- rnorm(n, S %*% bets, exp(-3 + S %*% gams))




test_data <- cbind.data.frame(y,S)
names(test_data) <- c("y",paste0("s",1:5))
str(test_data)

m <- lmls(y ~ s1 + s2 + s3 + s4,
          scale = ~ s1 + s3 + s4 ,
          data = test_data,
          light =TRUE)

m


v0 <- 0.01
M <- 10000
stepsize <- .01



n <- 10000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
x4 <- runif(n)
x5 <- runif(n)
y <- rnorm(n, 0 + 1 * x1 + 2 * x3 - 6*x4, exp(-3 + 1 * x2 + 1 * x3))
test_data <- cbind.data.frame(y,x1,x2,x3,x4,x5)

m <- lmls(y ~ x1 + x2 + x3 + x4 + x5, ~ x2 + x3 + x4 + x5,
          data = test_data, light = FALSE)

res <- spike_slab(m)








