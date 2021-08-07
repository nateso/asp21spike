


#' Title
#'
#' @param m 
#' @param a_theta_loc,b_theta_loc,a_tau_loc,b_tau_loc Hyperparameters for the location (see details). 
#' @param a_theta_scl,b_theta_scl,a_tau_scl,b_tau_scl Hyperparameters for the scale (see details).
#' @param v_0 Factor for the spike component.
#' @param nsim Number of iterations.
#' @inheritParams mcmc
#'
#' @export
spike_slab <- function(m,
                       a_theta_loc = 0.5,
                       b_theta_loc = 0.5,
                       a_tau_loc = 1,
                       b_tau_loc = 1,
                       a_theta_scl = a_theta_loc,
                       b_theta_scl = b_theta_loc,
                       a_tau_scl = a_tau_loc,
                       b_tau_scl = b_tau_loc,
                       nu0 = 1e-8,
                       nsim = 1000, 
                       stepsize = sqrt(3) * (m$df)^(-1/6)) {
  
  # initialise hyperparameters
  hyper <- list(a_theta = setNames(c(a_theta_loc, a_theta_scl), c("location", "scale")),
                b_theta = setNames(c(b_theta_loc, b_theta_scl), c("location", "scale")),
                a_tau = setNames(c(a_tau_loc, a_tau_scl), c("location", "scale")),
                b_tau = setNames(c(b_tau_loc, b_tau_scl), c("location", "scale")))
  
  # intitalise objects
  theta <- coefs <- tau <- delta <- list(location = matrix(NA, nrow = M, ncol = ncol(m$x)),
                                         scale = matrix(NA, nrow = M, ncol = ncol(m$z)))
  
  # initialise parameters
  for (k in 1:2) {
    n_params <- ncol(delta[[k]])
    theta[[k]][1, ] <- rbeta(n_params, hyper$a_theta[k], hyper$b_theta[k])
    for(l in 1:n_params){
      delta[[k]][1, l] <- rbinom(1, 1, theta[[k]][1, l])
      tau[[k]][1, l] <- sample_tau(delta[[k]][1, l], hyper$a_tau[k], hyper$b_tau[k], v_0)
    }
    m$tau[[k]] <- tau[[k]][1, ]
    
    param <- names(delta)[k]
    coefs[[k]][1, ] <- coef(mmala_update(m, param, stepsize = stepsize))[[param]]
  }
  
  # update parameters
  for (mm in 2:M) {
    for (kk in 1:2) {
      n_params <- ncol(coefs[[kk]])
      
      for (ll in 1:n_params) {
        ## update theta_lk
        # first make sure that we are using the most recent deltas
        current_delta <- delta[[kk]][mm, ]
        current_delta[is.na(current_delta)] <- delta[[kk]][mm-1, is.na(current_delta)]
        theta[[kk]][mm, ll] <- update_theta(hyper$a_theta[kk], 
                                           hyper$b_theta[kk], 
                                           current_delta)
        # this makes sure that we always include the intercept in the model
        # the intercept is always assigned to the slab. However, this does not 
        # imply a flat prior --> need to incorporate.
        has_intercept <- "(Intercept)" %in% names(m$coefficients[[kk]])
        if(has_intercept){
          theta[[kk]][mm, 1] <- 1
        }
        
        ## update delta_lk
        delta[[kk]][mm, ll] <- update_delta(theta = theta[[kk]][mm, ll], 
                                           tau_j = tau[[kk]][mm-1, ll], 
                                           a_tau = hyper$a_tau[kk], 
                                           b_tau = hyper$b_tau[kk], 
                                           v_0 = v_0)
        
        ## update tau_lk
        tau[[kk]][mm, ll] <- update_tau(a_tau = hyper$a_tau[kk], 
                                       b_tau = hyper$b_tau[kk], 
                                       v_0 = v_0, 
                                       beta_j = coefs[[kk]][mm-1, ll], 
                                       delta_j = delta[[kk]][mm, ll])
      }
      
      ## update taus in the model
      # as of now this first samples all other parameters first and then samples 
      # the coefficients all together. However, it makes more sense to sample the
      # coefficient values within the loop. 
      
      m$tau[[kk]] <- tau[[kk]][mm, ]
      
      # Lastly update coefficients
      param <- names(delta)[kk]
      m <- mmala_update_spike(m, param, stepsize = stepsize)
      coefs[[kk]][mm, ] <- coef(m)[[param]]
    }
  }
  
  # delete tau
  m$tau <- NULL
  
  # add results to model
  m$spike <- list("coefs" = coefs,
                  "tau" = tau,
                  "delta" = delta,
                  "theta" = theta)
  
  class(m) <- "lmls_spike"
  
  # return
  return(m)
}