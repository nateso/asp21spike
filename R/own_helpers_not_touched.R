# stuff we do not need to touch

score <- function(m, predictor) {
  score_funs[[predictor]](m)
}

chol_info <- function(m, predictor) {
  chol_info_funs[[predictor]](m)
}

set_beta <- function(m, beta) {
  m$coefficients$location[] <- beta
  m$fitted.values$location <- drop(m$x %*% beta)
  m$residuals <- m$y - fitted(m, "location")
  m
}

set_gamma <- function(m, gamma) {
  m$coefficients$scale[] <- gamma
  m$fitted.values$scale <- exp(drop(m$z %*% gamma))
  m
}

set_coef_funs <- list(
  location = set_beta,
  scale = set_gamma
)

set_coef <- function(m, predictor, coef) {
  set_coef_funs[[predictor]](m, coef)
}

# old score functions 
score_beta_old <- function(m) {
  drop((resid(m) / fitted(m, "scale")^2) %*% m$x)
}

score_gamma_old <- function(m) {
  drop((resid(m, "pearson")^2 - 1) %*% m$z)
}

# old fisher information
info_beta_old <- function(m) {
  # this is simply the design matrix normalised with the fitted SCALE (!) values and
  # then t(matrix) %*% matrix (crossprod is slightly faster)
  
  # the the inverse of this is the variance of beta (weighted least squares)
  crossprod(m$x / fitted(m, "scale"))
}


#### The mulivariate Normal Distribution ####
rmvnorm <- function(n, mu = 0, chol_sig_inv) {
  dim <- nrow(chol_sig_inv)
  
  std_norm <- matrix(rnorm(dim * n), dim, n)
  scaled <- backsolve(chol_sig_inv, std_norm)
  shifted <- scaled + mu
  
  shifted
}

dmvnorm <- function(x, mu = 0, chol_sig_inv, log = FALSE) {
  std_norm <- drop(chol_sig_inv %*% (x - mu))
  correction <- sum(log(diag(chol_sig_inv)))
  
  log_prob <- dnorm(std_norm, log = TRUE)
  
  if (is.matrix(log_prob)) {
    log_prob <- colSums(log_prob) + correction
  } else {
    log_prob <- sum(log_prob) + correction
  }
  
  if (log) {
    log_prob
  } else {
    exp(log_prob)
  }
}


#### GIBBS SAMPLER ####

gibbs_update_beta <- function(curr_m) {
  fit <- lm.wfit(curr_m$x, curr_m$y, fitted(curr_m, "scale")^(-2))
  
  mu <- coef(fit)
  chol_sig_inv <- chol(info_beta(curr_m))
  next_beta <- rmvnorm(1, mu, chol_sig_inv)
  next_m <- set_beta(curr_m, next_beta)
  
  next_m
}

