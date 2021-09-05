# Set the current working directory (to save the results later)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
devtools::load_all()

# simulate some data -----------------------------------------------

set.seed(1509)

n <- list(1000,
          300)

bets <- list(matrix(c(2,6,3,3,3,0,0,0,0,0,0)),
             matrix(c(2,6,3,3,3,3,3,3,3,0,0)))

gams <-list(matrix(c(2,6,3,3,3,0,0,0,0,0,0)),
            matrix(c(2,6,3,3,3,3,3,3,3,0,0)))

snr <- list(10, 3)

all_combis <- expand.grid(n, bets, gams, snr)
all_combis <- all_combis[rep(1:nrow(all_combis), each = 200), ]
names(all_combis) <- c("n", "bets", "gams", "snr")

mu <- rep(1,11)
sigma <- rep(1,length(mu))

data_cuts <- round(seq(0, nrow(all_combis), length.out = 6))

save(all_combis, mu, sigma, data_cuts,
     file = "sim_data.RData")

pb <- txtProgressBar(0, nrow(all_combis), style = 3)
system.time({
  for(j in 2:length(data_cuts)){
    i_values <- seq((data_cuts[(j - 1)] + 1),
                    data_cuts[j])
    spsl <- m <- list()
    for(i in i_values){
      X <- matrix(NA,nrow = all_combis$n[[i]], ncol = length(mu))
      for(jj in 1:length(mu)){
        X[,jj] <- 0.1 * rnorm(all_combis$n[[i]],mu[jj],sigma[jj])
      }
      
      eta_b <- X %*% all_combis$bets[[i]]
      eta_g <- X %*% all_combis$gams[[i]]
      
      eps <- rt(all_combis$n[[i]], 5)
      
      y <- rnorm(all_combis$n[[i]], eta_b, exp(eta_g)) + (sd(eta_g)/all_combis$snr[[i]]) * eps
      ## summary(y)
      test_data <- cbind.data.frame(y,X)
      names(test_data) <- c("y",paste0("x",1:11))
      
      m[[i]] <- try(lmls(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
                         scale = ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
                         data = test_data,
                         light = FALSE,maxit = 1000), silent = TRUE)
      
      if(class(m[[i]]) == "lmls"){
        spsl[[i]] <- spike_slab(m[[i]],
                                v_0 = 0.1,
                                a_theta_loc = 1,
                                b_theta_loc = 1,
                                a_theta_scl = 1,
                                b_theta_scl = 1,
                                a_tau_loc = 5,
                                b_tau_loc = 25,
                                a_tau_scl = 5,
                                b_tau_scl = 25,
                                burnin = 10,
                                coef_init = 0,
                                #always_in_loc = 'x4',
                                nsim = 1000,
                                prog_bar = FALSE)
        
        ## cat("-----------------------------------------------------------------------------\n")
        ## print(all_combis[i,])
        ## summary(spsl[[i]])
      } else {
        spsl[[i]] <- NULL
        print("Fail")
      }
      setTxtProgressBar(pb, i)
    }
    save(m, spsl,
         file = paste0("sim_results_", j - 1, ".RData"))
  }
})
close(pb)


# calculate misslcassification rate or accuracy:

calc_acc <- function(mod,bets,gams){
  sel <- names(mod$spike$coefs$location) %in% names(mod$spike$delta$location)
  if("(Intercept)" %in% names(mod$spike$coefs$location)){
    sel <- sel[-1]
  }
  bets_sel <- bets[sel]
  sel <- names(mod$spike$coefs$scale) %in% names(mod$spike$delta$scale)
  if("(Intercept)" %in% names(mod$spike$coefs$scale)){
    sel <- sel[-1]
  }
  gams_sel <- gams[sel]
  actual_loc <- ifelse(bets_sel != 0, T,F)
  actual_scl <- ifelse(gams_sel != 0, T,F)
  
  TP <- sum(mod$spike$delta$location[,actual_loc]) + sum(mod$spike$delta$scale[,actual_scl])
  TN <- sum(mod$spike$delta$location[,actual_loc == F] == 0) + sum(mod$spik$delta$scale[,actual_scl == F] == 0)
  FN <- sum(mod$spike$delta$location[,actual_loc] == 0) + sum(mod$spike$delta$scale[,actual_scl] == 0)
  FP <- sum(mod$spike$delta$location[,actual_loc == F]) + sum(mod$spike$delta$scale[,actual_scl == F])
  
  out <- matrix(c(TP,FP,FN,TN),nrow = 2, ncol = 2,
                dimnames = list(c('actual_P','actual_N'),c('pred_P','pred_N')))
  n <- (length(actual_loc) + length(actual_scl)) * nrow(mod$spike$delta$location)
  out <- out / n
  return(list(confusion = out,
              accuracy = sum(diag(out)),
              missclass = 1 - sum(diag(out))))
}

for(i in 1:20){
  print(paste0(i, ": ", 
               round(calc_acc(spsl[[i]], all_combis$bets[[i]], all_combis$gams[[i]])$accuracy, 3)))
}

# View(all_combis)

all_combis
plot(test_data$x9, test_data$y)



