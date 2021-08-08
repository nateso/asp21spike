#' @importFrom invgamma rinvgamma

sample_tau <- function(delta, a_tau, b_tau, v_0) {
  
  # 'INPUTS:
  # delta: status of inclusion
  # a_tau: hyper parameter for inverse gamma distribution
  # b_tau: hyper parameter for inverse gamma distribution
  # v_0: hyper parameter controlling the spike
  # OUTPUTS:
  # tau: sampled variances'
  
  k <- length(delta)
  tau <- rep(NA, k)
  for (jj in 1:k) {
    if (delta[jj] == 1) {
      tau[jj] <- rinvgamma(1, a_tau, b_tau)
    } else {
      tau[jj] <- rinvgamma(1, a_tau, b_tau * v_0)
    }
  }
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
  
  num <- theta * dinvgamma(tau, a_tau, b_tau) 
  denum <- theta * dinvgamma(tau, a_tau, b_tau) + (1 - theta) * dinvgamma(tau, a_tau, b_tau * v_0)
  theta_new <- num / denum
  
  delta <- rbinom(1, 1, theta_new)
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
  
  if (delta == 0) {
    tau <- rinvgamma(1, a_tau + 0.5, b_tau * v_0 + 0.5 * beta^2)
  } else {
    tau <- rinvgamma(1, a_tau + 0.5, b_tau + 0.5 * beta^2)
  }
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