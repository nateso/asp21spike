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
         sum(mod$spike$delta$scale[,actual_scl == F] == 0))
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

gc()

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
bets_spars <- factor(bets_spars, levels = c("high", "low"))

gams_spars <- sapply(all_combis$gams, 
                     function(x){
                       ifelse(sum(x == 0) < 11, 
                              ifelse(sum(x == 0) < 6,
                                     "low",
                                     "high"), 
                              "max")
                     })
gams_spars <- factor(gams_spars, levels = c("max", "high", "low"))


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

all_x_levels <- apply(expand.grid(levels(gams_spars),
                                  levels(bets_spars),
                                  "n"   = sort(unique(plot_data$n)),
                                  "snr" = sort(unique(plot_data$snr)))[, 4:1],
                      1, paste0, collapse = " - ")
plot_data$all_x     <- factor(apply(plot_data[, -1:-3], 1, paste0, collapse = " - "), levels = all_x_levels)
plot_data$snr_n     <- factor(apply(plot_data[, c("snr", "n")], 1, paste0, collapse = " - "),
                              levels = c("3 - 300", "3 - 1000", "10 - 300", "10 - 1000"))
plot_data$bets_gams <- factor(apply(plot_data[, c("bets_spars", "gams_spars")], 1, paste0, collapse = " - "),
                              levels = paste0(rep(levels(bets_spars), each = length(levels(gams_spars))), " - ",
                                              rep(levels(gams_spars), length(levels(bets_spars)))))

plot_data_long <- data.frame("value" = c(plot_data[, 1], plot_data[, 2], plot_data[, 3]),
                             "what" = rep(c("Accuracy", "Error Location", "Error Scale"), each = nrow(plot_data)),
                             plot_data[rep(1:nrow(plot_data), 3), -1:-3])
str(plot_data_long)

col_n <- ifelse(length(unique(plot_data_long$all_x)) > 6, 4, 2)
col <- rep(letters[seq(1, length(unique(plot_data_long$all_x)) / col_n)], col_n)
plot_data_long$col <- col[as.numeric(plot_data_long$all_x)]

# ggplot(plot_data_long,
#               aes(all_x, value, fill = col)) +
#   geom_violin(na.rm = TRUE, 
#               draw_quantiles = 1:3 * 0.25, 
#               scale = "width") +
#   # theme_classic() +
#   theme(axis.text.x = element_text(angle = cust_angle, 
#                                    vjust = ifelse(cust_angle == 0, 0.5, 1), 
#                                    hjust = ifelse(cust_angle == 0, 0.5, 1)),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none") +
#   geom_vline(xintercept = 1:3 * 0.25 * length(unique(x)) + 0.5,
#              linetype = "longdash") + 
#   facet_grid(rows = vars(what),
#              scales = "free_y")


# default_plot <- function(data, x, y, x_lab, y_lab, cust_angle = 0, facet_rows = NULL){
#   # get_rect <- data.frame("xstart" = seq(0.5, length(levels(x)) - 0.5),
#   #                        "xend" = seq(1.5, length(levels(x)) + 0.5),
#   #                        "cols" = c("#33333333", "#00000000"),
#   #                        stringsAsFactors = FALSE)
#   # str(get_rect)
#   # str(get_rect)
#   col_n <- ifelse(length(unique(x)) > 6, 4, 2)
#   col <- rep(letters[seq(1, length(unique(x)) / col_n)], col_n)
#   data$col <- col[as.numeric(x)]
#   p <- ggplot(data, aes(x, y, fill = col)) +
#     geom_violin(na.rm = TRUE, 
#                 draw_quantiles = 1:3 * 0.25, 
#                 scale = "width") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = cust_angle, 
#                                      vjust = ifelse(cust_angle == 0, 0.5, 1), 
#                                      hjust = ifelse(cust_angle == 0, 0.5, 1)),
#           legend.position = "none") +
#     xlab(x_lab) +
#     ylab(y_lab) + 
#     geom_vline(xintercept = 0.5 * (length(unique(x)) + 1),
#                linetype = "longdash")
#   if(length(unique(x)) > 6){
#     p <- p + geom_vline(xintercept = c(0.25, 0.75) * (length(unique(x))) + 0.5,
#                         linetype = "longdash")
#   }
#   if(y_lab == "Accuracy"){ 
#     p <- p + ylim(c(min(y), 1))
#   }
#   if(!is.null(facet_rows)){
#     p <- p + facet_grid(rows = vars(facet_rows),
#                         scales = "free_y")
#   }
#   p
# }

pdf("../../simulation_study/violin_plots_all.pdf",
    width = 8,
    height = 8)
cust_angle <- 60
ggplot(plot_data_long,
       aes(all_x, value, fill = col)) +
  geom_violin(na.rm = TRUE, 
              draw_quantiles = 1:3 * 0.25, 
              scale = "width") +
  # theme_classic() +
  theme(axis.text.x = element_text(angle = cust_angle, 
                                   vjust = ifelse(cust_angle == 0, 0.5, 1), 
                                   hjust = ifelse(cust_angle == 0, 0.5, 1)),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("") + 
  geom_vline(xintercept = 1:3 * 0.25 * length(unique(plot_data_long$all_x)) + 0.5,
             linetype = "longdash") + 
  facet_grid(rows = vars(what),
             scales = "free_y")
dev.off()

pdf("../../simulation_study/violin_plots_sep.pdf",
    width = 4,
    height = 8)
cust_angle <- 30
ggplot(plot_data_long,
       aes(snr_n, value, fill = snr_n)) +
  geom_violin(na.rm = TRUE, 
              draw_quantiles = 1:3 * 0.25, 
              scale = "width") +
  # theme_classic() +
  theme(axis.text.x = element_text(angle = cust_angle, 
                                   vjust = ifelse(cust_angle == 0, 0.5, 1), 
                                   hjust = ifelse(cust_angle == 0, 0.5, 1)),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("") + 
  geom_vline(xintercept = 0.5 * length(unique(plot_data_long$snr_n)) + 0.5,
             linetype = "longdash") + 
  facet_grid(rows = vars(what),
             scales = "free_y")
## Sparsity -------------------------------------------------------------------
ggplot(plot_data_long,
       aes(bets_gams, value, fill = bets_gams)) +
  geom_violin(na.rm = TRUE, 
              draw_quantiles = 1:3 * 0.25, 
              scale = "width") +
  # theme_classic() +
  theme(axis.text.x = element_text(angle = cust_angle, 
                                   vjust = ifelse(cust_angle == 0, 0.5, 1), 
                                   hjust = ifelse(cust_angle == 0, 0.5, 1)),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("") + 
  geom_vline(xintercept = 0.5 * length(unique(plot_data_long$bets_gams)) + 0.5,
             linetype = "longdash") + 
  facet_grid(rows = vars(what),
             scales = "free_y")
dev.off()

plot_data[which.max(plot_data$accuracy), ]
plot_data[which.min(plot_data$error_loc), ]
plot_data[which.max(plot_data$error_scl), ]

# View(all_combis)

all_combis
plot(test_data$x9, test_data$y)


### ===========================================================================
### Evaluate for coefficients =================================================

calc_acc_coef <- function(mod, bets, gams, coef_val, coef_which, what){
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
  actual_loc <- bets_sel == coef_val
  actual_scl <- gams_sel == coef_val
  ## Accuracy -----------------------------------------------------------------
  if(what == "Accuracy"){
    if(coef_which == "Location"){
      if(coef_val == 0){
        delta <- 1 - mod$spike$delta$location
      } else {
        delta <- mod$spike$delta$location
      }
      TP <- sum(delta[, actual_loc])
      n <- ((sum(actual_loc)) *
              nrow(mod$spike$delta$location))
    } else if (coef_which == "Scale"){
      if(coef_val == 0){
        delta <- 1 - mod$spike$delta$scale
      } else {
        delta <- mod$spike$delta$scale
      }
      TP <- sum(delta[, actual_scl])
      n <- (sum(actual_scl) *
              nrow(mod$spike$delta$scale))
    }
    out <- TP / n
  } else if(what == "Error"){
    if(coef_which == "Location"){
      actual_coef <- c(FALSE, actual_loc)
    } else {
      actual_coef <- c(FALSE, actual_scl)
    }
    coefs <- apply(as.matrix(mod$spike$coefs[[tolower(coef_which)]][, actual_coef]),
                   2, mean)
    out <- sqrt(mean((coefs - coef_val)^2))
  }
  return(out)
}

spsl_sel <- (1:length(spsl))[sapply(spsl, function(x){class(x) == "lmls_spike"})]
eval_coef <- data.frame("coef" = rep(c(0, 2, 3, 6), 2 * 2),
                        "param" = rep(c("Location", "Scale"), each = 4 * 2),
                        "value" = NA,
                        "what" = rep(c("Accuracy", "Error"), each = 4),
                        stringsAsFactors = FALSE)
eval_coef <- eval_coef[rep(1:nrow(eval_coef), each = length(spsl_sel)), ]
eval_coef$spsl_sel <- spsl_sel

pb <- txtProgressBar(0, nrow(eval_coef), style = 3)
for(i in 1:nrow(eval_coef)){
  check <- any(all_combis[[ifelse(eval_coef$param[i] == "Location",
                                  "bets", "gams")]][[eval_coef$spsl_sel[i]]] %in% eval_coef$coef[i])
  if(check){
    eval_coef$value[i] <- calc_acc_coef(spsl[[eval_coef$spsl_sel[i]]], 
                                        all_combis$bets[[eval_coef$spsl_sel[i]]],
                                        all_combis$gams[[eval_coef$spsl_sel[i]]],
                                        coef_val = eval_coef$coef[i],
                                        coef_which = eval_coef$param[i], 
                                        what = eval_coef$what[i])
  }
  setTxtProgressBar(pb, i)
}
close(pb)

aggregate(eval_coef$value,
          by = as.list(eval_coef[, c(-3, -5)]), 
          FUN = mean, na.rm = TRUE)
summary(eval_coef)

eval_coef$coef  <- factor(eval_coef$coef, levels = c(0, 2, 3, 6))
eval_coef$param <- as.factor(eval_coef$param)
eval_coef$what  <- as.factor(eval_coef$what)

pdf("../../simulation_study/violin_plots_coef.pdf",
    width = 8,
    height = 5)
ggplot(eval_coef,
       aes(coef, value, fill = coef)) +
  geom_violin(na.rm = TRUE, 
              scale = "width",
              draw_quantiles = 1:3 * 0.25) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  facet_grid(rows = vars(what),
             cols = vars(param),
             scales = "free_y")
dev.off()


### ===========================================================================
### spikeSlabGAM comparison ===================================================

# Selecting results with worst accuracy
sel_worst <- (all_combis$snr == 3 &
                all_combis$n == 300 &
                bets_spars == "high" &
                gams_spars == "high")
# ... with best accuracy
sel_best  <-  (all_combis$snr == 10 &
               all_combis$n == 1000 &
               bets_spars == "low" &
               gams_spars == "low")
# ... with homoskedastic variance, for a fair comparison
sel_fair  <- (all_combis$snr == 10 &
                all_combis$n == 1000 &
                bets_spars == "low" &
                gams_spars == "max")
wb <- NA
wb[sel_worst] <- "worst"
wb[sel_best]  <- "best"
wb[sel_fair]  <- "fair"
wb <- factor(wb, levels = c("worst", "best", "fair"))

combis_ssg <- all_combis[sel_worst | sel_best | sel_fair, ]
combis_ssg$wb <- wb[sel_worst | sel_best | sel_fair]
combis_ssg$row <- (1:nrow(all_combis))[sel_best | sel_worst | sel_fair]

eval_ssg <- data.frame("row" = combis_ssg$row,
                       "accuracy" = NA)

unique(combis_ssg[, -6])

library(spikeSlabGAM)

set.seed(1509)
options(mc.cores = 4)

mod_ssg <- list()

pb <- txtProgressBar(0, nrow(combis_ssg), style = 3)
for(i in 1:nrow(combis_ssg)){
  mod_spsl <- spsl[[combis_ssg[i, ]$row]]
  if(class(mod_spsl) == "lmls_spike"){
    data_ssg <- data.frame("y" = mod_spsl$y,
                           mod_spsl$x[, -1])
    (system.time({
      sink("/dev/null")
      mod_ssg[[i]] <- spikeSlabGAM(y ~ lin(x1) + lin(x2) + lin(x3) + lin(x4) + lin(x5) + 
                                     lin(x6) + lin(x7) + lin(x8) + lin(x9) + lin(x10) + lin(x11),
                                   data = data_ssg)
      sink(NULL)
    }))
    sum_ssh <- summary(mod_ssg[[i]])
    ip <- sum_ssh$trmSummary
    ip <- ip[grep("lin", rownames(ip)), 1]
    true_coef <- c(combis_ssg$bets[[i]]) != 0
    eval_ssg$accuracy[i] <- mean(c(ip[true_coef], 1 - ip[!true_coef]))
  }
  setTxtProgressBar(pb, i)
}
close(pb)

boxplot(eval_ssg$accuracy ~ combis_ssg$wb)

plot_data_ssg <- data.frame(combis_ssg, 
                            "accuracy" = eval_ssg$accuracy,
                            "pack" = "spikeSlabGAM")
plot_data_ssg <- rbind.data.frame(plot_data_ssg,
                                  data.frame(plot_data_ssg[, -7:-8],
                                             "accuracy" = plot_data$accuracy[sel_worst | sel_best | sel_fair],
                                             "pack" = "asp21spike"))
str(plot_data_ssg, 1)

pdf("../../simulation_study/violin_plots_spikeslabgam.pdf",
    width = 8,
    height = 4)
ggplot(plot_data_ssg,
       aes(pack, accuracy, fill = pack)) +
  geom_violin(na.rm = TRUE, 
              scale = "width",
              draw_quantiles = 1:3 * 0.25) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  facet_grid(cols = vars(wb),
             scales = "free_y")
dev.off()

eval_coef_ssg <- combis_ssg[rep(1:nrow(combis_ssg), each = 4), c("wb", "row")]
eval_coef_ssg$coef <- c(0, 2, 3, 6)
str(eval_coef_ssg)
head(eval_coef_ssg)

