#' @importFrom invgamma rinvgamma

sample_tau <- function(delta, a_tau, b_tau, v_0) {
  
  # 'INPUTS:
  # delta: status of inclusion
  # a_tau: hyper parameter for inverse gamma distribution
  # b_tau: hyper parameter for inverse gamma distribution
  # v_0: hyper parameter controlling the spike
  # OUTPUTS:
  # tau: sampled variances'
  
  n <- length(delta)
  tau <- delta * rinvgamma(n, a_tau, b_tau) + (1-delta) * rinvgamma(n, a_tau, b_tau * v_0)
  return(tau)
}

#' @importFrom invgamma dinvgamma
#' @importFrom stats rbinom

update_delta <- function(theta, tau, a_tau, b_tau, v_0) {
  
  # 'INPUTS:
  # theta: inclusion probability
  # a_tau: hyper parameter for inverse gamma distribution
  # b_tau: hyper paramerer for inverse gamma distribution
  # v_0: hyper parameter controling the spike
  # OUTPUTS:
  # delta: updated status of inclusion'
  n <- length(tau)
  
  num <- theta * dinvgamma(tau, a_tau, b_tau) 
  denum <- theta * dinvgamma(tau, a_tau, b_tau) + (1 - theta) * dinvgamma(tau, a_tau, b_tau * v_0)
  theta_new <- num / denum
  
  delta <- rbinom(n, 1, theta_new)
  return(delta)
}

#' @importFrom invgamma rinvgamma

update_tau <- function(delta, beta, a_tau, b_tau, v_0) {
  
  # 'INPUTS:
  # delta: status of inclusion
  # beta: coefficient of location/scale
  # a_tau: hyper parameter for inverse gamma distribution
  # b_tau: hyper paramerer for inverse gamma distribution
  # v_0: hyper parameter controling the spike
  # OUTPUTS:
  # tau: updated variances'
  n <- length(delta)
  
  tau <- delta * rinvgamma(n, a_tau + 0.5, b_tau + 0.5 * beta^2) +
    (1-delta) * rinvgamma(n, a_tau + 0.5, b_tau * v_0 + 0.5 * beta^2)

  return(tau)
}

#' @importFrom stats rbeta

update_theta <- function(delta, a_theta, b_theta) {
  
  # 'INPUTS:
  # delta: status of inclusion
  # a_theta: hyper parameter for beta distribution
  # b_theta: hyper parameter for beta distribution
  # OUTPUTS:
  # theta: updated inclusion probability'
   
  s <- sum(delta)
  k <- length(delta)
  theta = rbeta(1, a_theta + s, b_theta + k - s)
  return(theta)
}