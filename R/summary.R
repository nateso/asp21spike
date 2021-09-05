#' Summary for spike and slab priors
#' 
#' @param object Object of class lmls_spike
#' @param ... Currently not used
#' 
#' @method summary lmls_spike
#' 
#' @export
summary.lmls_spike <- function(object, ...) {
  x <- object
  coefs_loc <- x$spike$coefs$location
  means_loc <- colMeans(coefs_loc)
  delta_loc <- x$spike$delta$location
  probs_loc <- rep(NA, length(means_loc))
  sel <- names(coefs_loc) %in% names(delta_loc)
  probs_loc[sel] <- colMeans(delta_loc)
  sd_loc <- apply(coefs_loc, 2, sd)
  res_loc <- matrix(c(probs_loc, means_loc, sd_loc), 
                    ncol = 3, 
                    dimnames = list(colnames(coefs_loc),
                                    NULL))
  res_loc <- format(res_loc, digits = 3)
  res_loc <- gsub("NA", " -", res_loc)
  res_loc <- rbind(c("Inc. prob.", "Posterior mean", "Posterior sd"),
                   res_loc)
  res_loc <- format(res_loc, justify = "right")
  res_loc <- cbind(res_loc,
                   c("   ",
                     ifelse(!is.na(probs_loc) & probs_loc > 0.5,
                            ifelse(probs_loc > 0.75,
                                   ifelse(probs_loc > 0.9,
                                          "***",
                                          "** "),
                                   "*  "),
                            "   ")))
  colnames(res_loc) <- res_loc[1, ]
  res_loc <- res_loc[-1, ]
  
  # Scale
  coefs_scl <- x$spike$coefs$scale
  means_scl <- colMeans(coefs_scl)
  delta_scl <- x$spike$delta$scale
  probs_scl <- rep(NA, length(means_scl))
  sel <- names(coefs_scl) %in% names(delta_scl)
  probs_scl[sel] <- colMeans(delta_scl)
  sd_scl <- apply(coefs_scl, 2, sd)
  res_scl <- matrix(c(probs_scl, means_scl, sd_scl), 
                    ncol = 3, 
                    dimnames = list(colnames(coefs_scl),
                                    NULL))
  res_scl <- format(res_scl, digits = 3)
  res_scl <- gsub("NA", " -", res_scl)
  res_scl <- rbind(c("Inc. prob.", "Posterior mean", "Posterior sd"),
                   res_scl)
  res_scl <- format(res_scl, justify = "right")
  res_scl <- cbind(res_scl,
                   c("   ",
                     ifelse(!is.na(probs_scl) & probs_scl > 0.5,
                            ifelse(probs_scl > 0.75,
                                   ifelse(probs_scl > 0.9,
                                          "***",
                                          "** "),
                                   "*  "),
                            "   ")))
  colnames(res_scl) <- res_scl[1, ]
  res_scl <- res_scl[-1, ]
  
  cat("- Spike and slab prior variable selection -\n\n",
      x$nobs, " observations, ", x$spike$nsim, " iterations\n\n",
      "Location:\n", sep = "")
  prmatrix(res_loc, quote = FALSE)
  cat("\nScale:\n")
  prmatrix(res_scl, quote = FALSE)
  cat("\nIncl. prob.:  > 0.9: '***'  |  > 0.75: '**'  |  > 0.5: '*'\n")
  cat("Not selected: '-'\n")
}
