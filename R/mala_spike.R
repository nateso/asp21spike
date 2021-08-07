#' @importFrom stats coef

mmala_propose_spike <- function(curr_m, predictor, stepsize) {
  curr_coef <- coef(curr_m, predictor)
  curr_score <- score_spike(curr_m, predictor)
  curr_chol_info <- chol_info_spike(curr_m, predictor)
  
  step <- backsolve(
    r = curr_chol_info,
    x = forwardsolve(
      l = curr_chol_info,
      x = curr_score,
      upper.tri = TRUE,
      transpose = TRUE
    )
  )
  
  mu <- curr_coef + stepsize^2 / 2 * step
  chol_sig_inv <- curr_chol_info / stepsize
  prop_coef <- drop(rmvnorm(1, mu, chol_sig_inv))
  forward <- dmvnorm(prop_coef, mu, chol_sig_inv, log = TRUE)
  list(prop_coef = prop_coef, forward = forward)
}

#' @importFrom stats coef

mmala_backward_spike <- function(curr_m, prop_m, predictor, stepsize) {
  curr_coef <- coef(curr_m, predictor)
  prop_coef <- coef(prop_m, predictor)
  prop_score <- score_spike(prop_m, predictor)
  prop_chol_info <- chol_info_spike(prop_m, predictor)
  
  step <- backsolve(
    r = prop_chol_info,
    x = forwardsolve(
      l = prop_chol_info,
      x = prop_score,
      upper.tri = TRUE,
      transpose = TRUE
    )
  )
  
  mu <- prop_coef + stepsize^2 / 2 * step
  chol_sig_inv <- prop_chol_info / stepsize
  dmvnorm(curr_coef, mu, chol_sig_inv, log = TRUE)
}

#' @importFrom stats logLik runif

mmala_update_spike <- function(curr_m, predictor, stepsize) {
  proposal <- mmala_propose_spike(curr_m, predictor, stepsize)
  prop_coef <- proposal$prop_coef
  forward <- proposal$forward
  
  prop_m <- set_coef(curr_m, predictor, prop_coef)
  backward <- mmala_backward_spike(curr_m, prop_m, predictor, stepsize)
  
  curr_coef <- coef(curr_m,predictor)
  Tau <- diag(curr_m$tau[[predictor]])
  
  prop_log_prior <- mvtnorm::dmvnorm(prop_coef, mean = rep(0,nrow(Tau)), sigma = Tau, log = T)
  curr_log_prior <- mvtnorm::dmvnorm(curr_coef, mean = rep(0,nrow(Tau)), sigma = Tau, log = T)
  
  alpha <- (logLik(prop_m) + prop_log_prior - logLik(curr_m) - curr_log_prior +
              backward - forward)
  
  if (log(runif(1)) <= alpha) {
    prop_m
  } else {
    curr_m
  }
}

# derivatives -----------------------------------------------------------------

#' @importFrom stats fitted resid

score_beta_spike <- function(m) {
  drop((resid(m) / fitted(m, "scale")^2) %*% m$x) - coef(m,'location')/m$tau$location)
}

#' @importFrom stats resid

score_gamma_spike <- function(m) {
  drop((resid(m, "pearson")^2 - 1) %*% m$z) - coef(m,'scale')/m$tau$scale)
}

#' @importFrom stats fitted

info_beta_spike <- function(m) {
  crossprod(m$x / fitted(m, "scale")) + diag(c(1/m$tau$location))
}

info_gamma_spike<- function(m) {
  2 * crossprod(m$z) + diag(c(1/m$tau$scale))
}

# wrappers ------------------------------------------------------
score_funs_spike <- list(
  location = score_beta_spike,
  scale = score_gamma_spike
)

score_spike <- function(m, predictor) {
  score_funs_spike[[predictor]](m)
}

chol_info_funs_spike <- list(
  location = function(m) chol(info_beta_spike(m)),
  scale = function(m) chol(info_gamma_spike(m))
)

chol_info_spike <- function(m, predictor) {
  chol_info_funs_spike[[predictor]](m)
}









