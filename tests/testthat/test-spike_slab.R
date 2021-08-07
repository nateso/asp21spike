
library(MASS)

set.seed(3487)
n <- 100
mu_S = c(2,4,6,8,10)
sigma_S = diag(1,nrow = 5,ncol = 5)
S = mvrnorm(n, mu = mu_S, Sigma = sigma_S)

bets <- matrix(c(-2,6,3,0,0))
gams <- matrix(c(0,0,5,1,-10))

y <- rnorm(n, S %*% bets, exp(S %*% gams))



test_data <- cbind.data.frame(y,S)
names(test_data) <- c("y",paste0("s",1:5))


m <- lmls(y ~ s1 + s2 + s3 + s4,
          scale = ~ s1 + s2 + s3 + s4 ,
          data = test_data,
          light =FALSE)



v0 <- 0.01
M <- 10000
stepsize <- .01






