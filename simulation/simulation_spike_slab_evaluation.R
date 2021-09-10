# Set the current working directory (to save the results later)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
devtools::load_all()
getwd()

library(ggplot2)

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
         if(any(actual_scl)){
           sum(mod$spike$delta$scale[,actual_scl])
           } else {
             0
           })
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

load("../../simulation_study/sim_data.RData")
load("../../simulation_study/sim_results.RData")


## Evaluation ==================================================================

eval <- data.frame("accuracy" = rep(NA, length(m)),
                   "error_loc" = NA,
                   "error_scl" = NA)

for(i in 1:length(m)){
  if(class(spsl[[i]]) == "lmls_spike"){
    # Accuracy
    eval$accuracy[i] <- calc_acc(mod = spsl[[i]], 
                                 bets = all_combis$bets[[i]],
                                 gams = all_combis$gams[[i]])$accuracy
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

bets_spars <- sapply(all_combis$bets, 
                     function(x){
                       ifelse(all(x == all_combis$bets[[1]]), 
                              "high", "low")
                     })

gams_spars <- sapply(all_combis$gams, 
                     function(x){
                       ifelse(all(x == all_combis$gams[[1]]), 
                              "high", "low")
                     })

boxplot(eval$accuracy ~ unlist(all_combis$snr) + unlist(all_combis$n))
boxplot(eval$error_loc ~ unlist(all_combis$snr) + unlist(all_combis$n))
boxplot(eval$error_scl ~ unlist(all_combis$snr) + unlist(all_combis$n))


boxplot(eval$accuracy ~ gams_spars + bets_spars)
boxplot(eval$error_loc ~ gams_spars + bets_spars)
boxplot(eval$error_scl ~ gams_spars + bets_spars)

boxplot(eval$accuracy ~ unlist(all_combis$snr) + unlist(all_combis$n) + bets_spars + gams_spars)

hist(eval$accuracy)
hist(eval$error_loc)
hist(eval$error_scl)

plot_data <- data.frame(eval,
                        "snr" = unlist(all_combis$snr),
                        "n" = unlist(all_combis$n),
                        "bets_spars" = bets_spars,
                        "gams_spars" = gams_spars)
str(plot_data)

plot_data$all_x     <- factor(apply(plot_data[, -1:-3], 1, paste0, collapse = " - "))
plot_data$snr_n     <- factor(apply(plot_data[, c("snr", "n")], 1, paste0, collapse = " - "))
plot_data$bets_gams <- factor(apply(plot_data[, c("bets_spars", "gams_spars")], 1, paste0, collapse = " - "))

ggplot(plot_data, 
       aes(y = accuracy, x = all_x)) +
  geom_violin(na.rm = TRUE, 
              draw_quantiles = 0.5, 
              scale = "width") +
  geom_boxplot(na.rm = TRUE,
               width = 0.2)

plot_data[which.max(plot_data$accuracy), ]
plot_data[which.min(plot_data$error_loc), ]
plot_data[which.max(plot_data$error_scl), ]

# View(all_combis)

all_combis
plot(test_data$x9, test_data$y)



