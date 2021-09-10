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

ggplot(plot_data_long,
              aes(all_x, value, fill = col)) +
  geom_violin(na.rm = TRUE, 
              draw_quantiles = 1:3 * 0.25, 
              scale = "width") +
  # theme_classic() +
  theme(axis.text.x = element_text(angle = cust_angle, 
                                   vjust = ifelse(cust_angle == 0, 0.5, 1), 
                                   hjust = ifelse(cust_angle == 0, 0.5, 1)),
        legend.position = "none") +
  xlab("SNR - n - sparsity Beta - sparsity Gamma") +
  ylab("") + 
  geom_vline(xintercept = 1:3 * 0.25 * length(unique(x)) + 0.5,
             linetype = "longdash") + 
  facet_grid(rows = vars(what),
             scales = "free_y")


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
        legend.position = "none") +
  xlab("SNR - n - sparsity Beta - sparsity Gamma") +
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
        legend.position = "none") +
  xlab("SNR - n") +
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
        legend.position = "none") +
  xlab("sparsity Beta - sparsity Gamma") +
  ylab("") + 
  geom_vline(xintercept = 1:3 * 0.25 * length(unique(plot_data_long$bets_gams)) + 0.5,
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



