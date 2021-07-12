## own helpers

sample_tau <- function(delta,a_tau,b_tau,v0 = 0.001){
  'INPUTS:
  delta: inclusion vector
  a_tau: hyper parameter for inverse gamma distribution
  b_tau: hyper paramerer for inverse gamma distribution
  v0: hyper parameter controling the spike
  
  OUTPUTS:
  out: vector with sampled variances'
  k <- length(delta)
  out <- rep(NA,k)
  for(jj in 1:k){
    if(delta[jj] == 1){
      out[jj] <- rinvgamma(1,a_tau,b_tau)
    }
    else{
      out[jj] <- rinvgamma(1,a_tau,b_tau*v0)
    }
  }
  return(out)
} 

update_delta <- function(theta,tau_j,a_tau,b_tau,v_0){
  
  num <- theta * dinvgamma(tau_j,a_tau,b_tau) 
  denom <- theta * dinvgamma(tau_j,a_tau,b_tau) + (1-theta) * dinvgamma(tau_j,a_tau,b_tau*v_0)
  theta_new_j = num/denom
  
  delta_j <- rbinom(1,1,theta_new_j)
  return(delta_j)
}

update_tau <- function(nsmp,a_tau,b_tau,v_0,beta_j,delta_j){
  if(delta_j == 0){
    tau_j <- rinvgamma(nsmp,a_tau + 0.5, b_tau * v_0 + 0.5 * beta_j^2)
  }
  else{
    tau_j <- rinvgamma(nsmp,a_tau + 0.5, b_tau + 0.5 * beta_j^2)
  }
  return(tau_j)
}

update_theta <- function(a_theta,b_theta,delta){
  s <- sum(delta)
  k <- length(delta)
  theta = rbeta(1,a_theta + s, b_theta + k - s)
  return(theta)
}
