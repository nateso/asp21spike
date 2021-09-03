#' Title
#'
#' @param a_theta_loc,b_theta_loc,a_tau_loc,b_tau_loc Hyperparameters for the location (see details). 
#' @param a_theta_scl,b_theta_scl,a_tau_scl,b_tau_scl Hyperparameters for the scale (see details).
#' @param v_0 Factor for the spike component.
#' @param burnin Number of samples used as burn in phase, will be set to zero if negative.
#' @param prog_bar Show a progress bar?
#' @param always_in_loc,always_in_scl Select the variables which are NOT subject to selection.
#' @param coef_init Single value. If not NULL (default), the coefficients will be initialised with the given value.
#' @inheritParams mcmc
#' 
#' @importFrom stats setNames
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
spike_slab_test <- function(m,
                       a_theta_loc = 0.5,
                       b_theta_loc = 0.5,
                       a_tau_loc = 1,
                       b_tau_loc = 1,
                       a_theta_scl = a_theta_loc,
                       b_theta_scl = b_theta_loc,
                       a_tau_scl = a_tau_loc,
                       b_tau_scl = b_tau_loc,
                       v_0 = 1e-2,
                       nsim = 1000,
                       burnin = round(nsim * 0.2), 
                       stepsize = sqrt(3) * (m$df)^(-1/6), 
                       prog_bar = TRUE,
                       always_in_loc = NULL,
                       always_in_scl = NULL,
                       coef_init = NULL,
                       seed = NULL){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(!is.numeric(burnin) | burnin < 0){
    burnin <- 0
  }
  
  M <- round(nsim + burnin)
  
  ## Normal prior for intercept----------------------------------------------
  # this enables us to implement normal priors for all variables which should
  # always be kept in the model, as well as for the intercept.
  always_in <- list(location = always_in_loc,
                    scale = always_in_scl)
  which_select <- which_always_in <- list(location = NULL,
                                          scale = NULL)
  n_select <- n_always_in <- c('location' = NA,
                               'scale' = NA)
  for(k in 1:2){
    has_intercept <- "(Intercept)" %in% names(m$coefficients[[k]])
    if(has_intercept == TRUE){
      always_in[[k]] <- c("(Intercept)", always_in[[k]])
    }
    if(length(always_in[[k]]) != 0){
      which_always_in[[k]] <- which(names(m$coefficients[[k]]) %in% always_in[[k]])
      which_select[[k]] <- which((names(m$coefficients[[k]]) %in% always_in[[k]]) == FALSE)
    }
    else{
      which_select[[k]] <- 1:length(names(m$coefficients[[k]]))
    }
    n_select[k] <- length(m$coefficients[[k]]) - length(which_always_in[[k]])
    n_always_in[k] <- length(which_always_in[[k]])
  }
  
  ## Initialise parameters------------------------------------------------------
  
  # initialise hyperparameter objects
  hyper <- list(a_theta = setNames(c(a_theta_loc, a_theta_scl), c("location", "scale")),
                b_theta = setNames(c(b_theta_loc, b_theta_scl), c("location", "scale")),
                a_tau = setNames(c(a_tau_loc, a_tau_scl), c("location", "scale")),
                b_tau = setNames(c(b_tau_loc, b_tau_scl), c("location", "scale")))
  
  # intitalise parameter objects
  coefs <- tau <- list(location = matrix(NA, nrow = M, ncol = ncol(m$x)),
                       scale = matrix(NA, nrow = M, ncol = ncol(m$z)))
  
  theta <- delta <- list(location = matrix(NA, nrow = M, ncol = n_select['location']),
                         scale = matrix(NA, nrow = M, ncol = n_select['scale'])) # only for variables subject to selection
  
  # initialise parameters
  for (k in 1:2){
    theta[[k]][1, ] <- rbeta(n_select[k], hyper$a_theta[k], hyper$b_theta[k])
    delta[[k]][1, ] <- rbinom(n_select[k], 1, theta[[k]][1, ])
    if(n_select[k] != 0){
      tau[[k]][1,which_select[[k]]] <- sapply(delta[[k]][1,],sample_tau,
                                              a_tau = hyper$a_tau[k],
                                              b_tau = hyper$b_tau[k],
                                              v_0 = v_0)
    }
    if(n_always_in[k] != 0){
      tau[[k]][1,which_always_in[[k]]] <- sample_tau_nosel(n_always_in[k],
                                                           hyper$a_tau[k],
                                                           hyper$b_tau[k]) # prior for hyper-variance if variable not subject to selection
    }
    param <- names(delta)[k]
    
    if(!is.null(coef_init)){
      coefs[[k]][1, ] <- coef_init # runif(ncol(coefs[[k]]),-100,100) # random initialization of coefficients
    } else {
      coefs[[k]][1, ] <- coef(mmala_update(m, param, stepsize = stepsize))[[param]] # ML initialization
    }
  }
  
  # update parameters-----------------------------------------------------------
  mod_tmp <- m
  if(prog_bar) {
    pb <- txtProgressBar(0, M, style = 3)
  }
  for (mm in 2:M) {
    for (kk in 1:2) {
      # Update coefficients
      param <- names(delta)[kk] # extract 'location' or 'scale'
      mod_tmp <- mmala_update_spike(curr_m = mod_tmp, 
                              curr_tau = tau[[kk]][mm - 1, ],
                              predictor = param, 
                              stepsize = stepsize)
      coefs[[kk]][mm, ] <- coef(mod_tmp)[[param]]

      # first update hyper-variances for coefficients, which are not subject to selection
      if(length(n_always_in[k]) != 0){
        tau[[kk]][mm,which_always_in[[kk]]] <- update_tau_nosel(beta = coefs[[kk]][mm,which_always_in[[kk]]],
                                                                a_tau = hyper$a_tau[kk],
                                                                b_tau = hyper$b_tau[kk])
      }
      # second, update parameters for all coefficients subject to selection
      for (ll in which_select[[kk]]){
        iter <- which(which_select[[kk]] == ll) # iter runs from 1 to ncol(delta), while ll iterates over all columns in tau and coef which are subject to selection
        
        # make sure that we are using the most recent deltas
        current_delta <- delta[[kk]][mm, ]
        current_delta[is.na(current_delta)] <- delta[[kk]][mm - 1, is.na(current_delta)]
        
        ## update tau_lk
        tau[[kk]][mm, ll] <- update_tau(a_tau = hyper$a_tau[kk], 
                                        b_tau = hyper$b_tau[kk], 
                                        v_0 = v_0, 
                                        beta = coefs[[kk]][mm, ll], 
                                        delta = current_delta[iter])
        ## update theta_lk
        theta[[kk]][mm, iter] <- update_theta(delta = current_delta,
                                              a_theta = hyper$a_theta[kk],
                                              b_theta = hyper$b_theta[kk])
        
        ## update delta_lk
        delta[[kk]][mm, iter] <- update_delta(theta = theta[[kk]][mm, iter], 
                                              tau = tau[[kk]][mm, ll], 
                                              a_tau = hyper$a_tau[kk], 
                                              b_tau = hyper$b_tau[kk], 
                                              v_0 = v_0)
      }
    }
    if(prog_bar) {
      setTxtProgressBar(pb, mm)
    }
  }
  if(prog_bar) {
    close(pb)
  }
  ## disregard burn-in samples---------------------------------------------------
  if (burnin > 0) {
    for (kk in 1:2) {
      coefs[[kk]] <- coefs[[kk]][-1:-burnin, ]
      tau[[kk]]   <- tau[[kk]][-1:-burnin, ]
      delta[[kk]] <- delta[[kk]][-1:-burnin, ]
      theta[[kk]] <- theta[[kk]][-1:-burnin, ]
    }
  }
  
  ## prepare data for export----------------------------------------------------
  for(k in 1:2){
    coefs[[k]] <- as.data.frame(coefs[[k]])
    tau[[k]] <- as.data.frame(tau[[k]])
    delta[[k]] <- as.data.frame(delta[[k]])
    theta[[k]] <- as.data.frame(theta[[k]])
    
    names(tau[[k]]) <- names(coefs[[k]]) <- names(m$coefficients[[k]])
    names(delta[[k]]) <- names(theta[[k]]) <- names(m$coefficients[[k]])[which_select[[k]]]
  }
  
  # add results to model
  m$spike <- list("coefs" = coefs,
                  "tau" = tau,
                  "delta" = delta,
                  "theta" = theta,
                  "hyper" = hyper,
                  "v_0" = v_0,
                  "nsim" = nsim,
                  "burnin" = burnin,
                  "stepsize" = stepsize)
  
  class(m) <- "lmls_spike"
  
  # return
  return(m)
}