#' Print method for class 'lmls_spike'
#' 
#' @param x Object of class 'lmls_spike'
#' 
#' @method print lmls_spike
#' 
#' @export
print.lmls_spike <- function(x, ...) {
  probs_loc <- colMeans(x$spike$delta$location)
  names(probs_loc) <- colnames(x$spike$delta$location)
  probs_scl <- colMeans(x$spike$delta$scale)
  names(probs_scl) <- colnames(x$spike$delta$scale)
  cat("- Spike and slab prior variable selection -\n\n",
      x$nobs, " observations, ", x$spike$nsim, " iterations\n\n",
      "Inclusion probabilities: \nLocation:\n", sep = "")
  print(probs_loc)
  cat("\nScale:\n")
  print(probs_scl)
}
