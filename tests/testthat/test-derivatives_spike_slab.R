library(numDeriv)

set.seed(1337)

n <- 100
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- rnorm(n, 0 + 1 * x1 + 1 * x3, exp(-3 + 1 * x2 + 1 * x3))
m <- lmls(y ~ x1 + x3, ~ x2 + x3, light = FALSE)

# score -----------------------------------------------------------------------

beta <- gamma <- c(0, 0, 0)
m <- set_coef(m, "location", beta)
m <- set_coef(m, "scale", gamma)
Tau <- diag(3)
m$tau$location <- rep(1, 3)
m$tau$scale    <- rep(1, 3)

f <- function(x, predictor) {
  prop_log_prior <- mvtnorm::dmvnorm(x, mean = rep(0,nrow(Tau)), sigma = Tau, log = T)
  logLik(set_coef(m, predictor, x)) + prop_log_prior
}

test_that("score of beta is correct", {
  num_score <- numDeriv::grad(f, beta, predictor = "location")
  expect_exactly(score_spike(m, m$tau$location, "location"), num_score)
})

test_that("score of gamma is correct", {
  num_score <- numDeriv::grad(f, gamma, predictor = "scale")
  expect_exactly(score_spike(m, m$tau$scale, "scale"), num_score)
})

# fisher info -----------------------------------------------------------------

beta <- c(0, 1, 1)
gamma <- c(-3, 1, 1)
m <- set_coef(m, "location", beta)
m <- set_coef(m, "scale", gamma)

nsim <- 5000

reps <- replicate(nsim, {
  y <- rnorm(n, 0 + 1 * x1 + 1 * x3, exp(-3 + 1 * x2 + 1 * x3))
  m <- lmls(y ~ x1 + x3, ~ x2 + x3, light = FALSE)
  m <- set_coef(m, "location", beta)
  m <- set_coef(m, "scale", gamma)
  m$tau$location <- rep(1, 3)
  m$tau$scale    <- rep(1, 3)
  c(score_beta_spike(m, m$tau$location), score_gamma_spike(m, m$tau$scale))
})

num_info <- cov(t(reps))

test_that("fisher info of beta is correct", {
  expect_roughly(info_beta_spike(m, m$tau$location), num_info[1:3, 1:3])
})

test_that("fisher info of gamma is correct", {
  expect_roughly(info_gamma_spike(m, m$tau$scale), num_info[4:6, 4:6])
})

test_that("chol_info() works", {
  expect_equal(chol_info_spike(m, m$tau$location, "location"), 
               chol(info_beta_spike(m, m$tau$location)))
  expect_equal(chol_info_spike(m, m$tau$scale, "scale"), 
               chol(info_gamma_spike(m, m$tau$scale)))
})
