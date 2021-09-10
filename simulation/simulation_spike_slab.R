# Set the current working directory (to save the results later)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
devtools::load_all()

# simulate some data -----------------------------------------------

set.seed(1509)

n <- list(1000,
          300)

bets <- list(matrix(c(2,6,3,3,3,0,0,0,0,0,0)),
             matrix(c(2,6,3,3,3,3,3,3,3,0,0)))

gams <-list(matrix(c(0,0,0,0,0,0,0,0,0,0,0)),
            matrix(c(2,6,3,3,3,0,0,0,0,0,0)),
            matrix(c(2,6,3,3,3,3,3,3,3,0,0)))

snr <- list(10, 3)

all_combis <- expand.grid(n, bets, gams, snr)
all_combis <- all_combis[rep(1:nrow(all_combis), each = 200), ]
names(all_combis) <- c("n", "bets", "gams", "snr")

mu <- rep(0,11)
sigma <- rep(1,length(mu))

spsl <- m <- list()

data_cuts <- round(seq(0, nrow(all_combis), length.out = 6))

save(all_combis, mu, sigma, data_cuts,
     file = "../../simulation_study/sim_data.RData")

pb <- txtProgressBar(0, nrow(all_combis), style = 3)
system.time({
  for(j in 2:length(data_cuts)){
    i_values <- seq((data_cuts[(j - 1)] + 1),
                    data_cuts[j])
    for(i in i_values){
      X <- matrix(NA,nrow = all_combis$n[[i]], ncol = length(mu))
      for(jj in 1:length(mu)){
        X[,jj] <- 0.1 * rnorm(all_combis$n[[i]],mu[jj],sigma[jj])
      }
      
      eta_b <- X %*% all_combis$bets[[i]]
      eta_g <- X %*% all_combis$gams[[i]]
      
      eps <- rnorm(all_combis$n[[i]])
      
      y <- rnorm(all_combis$n[[i]], eta_b, exp(eta_g)) + (sd(eta_b)/all_combis$snr[[i]]) * eps
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
         file = paste0("../../simulation_study/sim_results_", j - 1, ".RData"))
  }
})
close(pb)
#     user   system  elapsed 
# 7679.278  166.338 7856.092 
# ca. 2h 10m

# calculate misslcassification rate or accuracy:

