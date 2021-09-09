# Set the current working directory (to save the results later)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
devtools::load_all()
getwd()
## Calculating accuracy ========================================================

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
  TP <- (sum(mod$spike$delta$location[,actual_loc]) +
         sum(mod$spike$delta$scale[,actual_scl]))
  TN <- (sum(mod$spike$delta$location[,actual_loc == F] == 0) +
         sum(mod$spik$delta$scale[,actual_scl == F] == 0))
  FN <- (sum(mod$spike$delta$location[,actual_loc] == 0) +
         sum(mod$spike$delta$scale[,actual_scl] == 0))
  FP <- (sum(mod$spike$delta$location[,actual_loc == F]) +
         sum(mod$spike$delta$scale[,actual_scl == F]))
  out <- matrix(c(TP, FP, FN, TN),
                nrow = 2, ncol = 2,
                dimnames = list(c('actual_P','actual_N'),
                                c('pred_P','pred_N')))
  n <- ((length(actual_loc) + length(actual_scl)) *
        nrow(mod$spike$delta$location))
  out <- out / n
  return(list(confusion = out,
              accuracy = sum(diag(out)),
              missclass = 1 - sum(diag(out))))
}

## Data ========================================================================

load("simulation_study/sim_data.RData")
load("simulation_study/sim_results.RData")


## Evaluation ==================================================================

eval <- data.frame("accuracy" = rep(NA, length(m)),
                   "error_loc" = NA,
                   "error_scl" = NA)

for(i in 1:length(m)){
  if(class(spsl[[i]]) == "lmls_spike"){
    # Accuracy
    eval$accuracy[i] <- calc_acc(spsl[[i]], 
                                 all_combis$bets[[i]],
                                 all_combis$gams[[i]])$accuracy
    # RMSE for location
    post_mean_loc <- apply(spsl[[i]]$spike$coefs$location, 2,
                           mean)
    eval$error_loc[i] <- sqrt(mean((post_mean_loc - c(0, all_combis$bets[[i]]))^2))
    # RMSE for scale
    post_mean_scl <- apply(spsl[[i]]$spike$coefs$scale, 2,
                           mean)
    eval$error_scl[i] <- sqrt(mean((post_mean_scl - c(0, all_combis$gams[[i]]))^2))
  }
}


hist(eval$accuracy)
hist(eval$error_loc)
hist(eval$error_scl)

# View(all_combis)

all_combis
plot(test_data$x9, test_data$y)



