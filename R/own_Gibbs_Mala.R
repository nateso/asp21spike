# Gibbs and MALA sampler


#### MALA SAMPLER ####

mmala_propose <- function(curr_m, predictor, stepsize) {
  curr_coef <- coef(curr_m, predictor)
  curr_score <- score(curr_m, predictor)
  curr_chol_info <- chol_info(curr_m, predictor) 
  
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


mmala_backward <- function(curr_m, prop_m, predictor, stepsize) {
  curr_coef <- coef(curr_m, predictor)
  prop_coef <- coef(prop_m, predictor)
  prop_score <- score(prop_m, predictor)
  prop_chol_info <- chol_info(prop_m, predictor)
  
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

mmala_update <- function(curr_m, predictor, stepsize) {
  proposal <- mmala_propose(curr_m, predictor, stepsize)
  prop_coef <- proposal$prop_coef
  forward <- proposal$forward
  
  prop_m <- set_coef(curr_m, predictor, prop_coef)
  backward <- mmala_backward(curr_m, prop_m, predictor, stepsize)
  
  curr_coef <- coef(curr_m,predictor)
  Tau <- diag(curr_m$tau[[predictor]])

  prop_log_prior <- mvtnorm::dmvnorm(prop_coef, mean = rep(0,nrow(Tau)), sigma = Tau, log = T)
  curr_log_prior <- mvtnorm::dmvnorm(curr_coef, mean = rep(0,nrow(Tau)), sigma = Tau, log = T)
  
  alpha <- (logLik(prop_m) + prop_log_prior - logLik(curr_m) - curr_log_prior +
              backward - forward)
  thresh <- log(runif(1))

  if (thresh <= alpha) {
    prop_m
  } else {
    curr_m
  }
}

#### HELPER FUNCTIONS ####

#### Score functions ####

score_beta <- function(m){
  drop((resid(m) / fitted(m,"scale")^2) %*% m$x - coef(m,'location')/m$tau$location)
}

score_gamma <- function(m) {
  drop((resid(m, "pearson")^2 - 1) %*% m$z - coef(m,'scale')/m$tau$scale)
}


#### Fisher informations ####

info_beta <- function(m){
  crossprod(m$x / fitted(m,'scale')^2) + diag(c(1/m$tau$location))
}

info_gamma <- function(m) {
  2 * crossprod(m$z) + diag(c(1/m$tau$scale))
}


#### wrappers to call the infos ####
score_funs <- list(
  location = score_beta,
  scale = score_gamma
)

chol_info_funs <- list(
  location = function(m) chol(info_beta(m)),
  scale = function(m) m$chol_info_gamma
)











