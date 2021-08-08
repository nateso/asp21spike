
#' Plot Method for lmls_spike object
#' This is some description
#' 
#' @importFrom graphics abline barplot par hist lines
#' @param x object of class spike_lmls
#' @param parameter takes a string which coefficients should be plotted, location or scale (default location)
#' @param plot_posterior Boolean whether to plot the posterior distributions or not (default FALSE). If False it plots inclusion probabilities.
#' @export

plot.lmls_spike <- function(x, parameter = 'location', plot_posterior = FALSE){
  cov_names <- names(x$coefficients[[parameter]])
  if(plot_posterior == F){
    p <- colMeans(x$spike$delta[[parameter]])
    names(p) <- cov_names
    barplot(p,
            main = paste("inclusion probabilities for",parameter),
            ylab = 'average inclusion',
            ylim = c(0,1))
    abline(h = 0.5, col = "red", lty = 'dashed')
    abline(h = 0.25, col = 'red', lty = 'dashed')
    abline(h = 0.75, col = 'red', lty = 'dashed')
    }
  else{
    k <- ncol(x$spike$coefs[[parameter]])
    n <- nrow(x$spike$coefs[[parameter]])
    ncol <- ceiling(sqrt(k))
    nrow <- ceiling(k/ncol)
    
    op <- par(mfrow = c(nrow,ncol),
              no.readonly = TRUE)
    
    for(jj in 1:k){
      dat <- x$spike$coefs[[parameter]][,jj]
      post_mean <- mean(dat)
      hist(dat, 
           freq = F,
           breaks = sqrt(n),
           main = paste("Posterior distribution of",cov_names[[jj]]),
           xlab = cov_names[[jj]]
      )
      lines(x = density(x = dat), col = "black")
      
      abline(v = post_mean, col = 'red', lwd = 2)
      }
    par(op)
    }
}
  
  

 