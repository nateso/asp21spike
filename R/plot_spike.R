
#' Plot Method for lmls_spike object
#' 
#' @param x Object of class spike_lmls
#' @param parameter Takes a string which coefficients should be plotted, 'location' (default) or 'scale', can be abbreviated.
#' @param plot_type Can be one of 'inclusion' (default), 'posterior' or 'random walk', can be abbreviated.
#' @param ... Currently not used.
#' 
#' @importFrom graphics abline barplot par hist lines axis
#' @importFrom stats density sd
#' 
#' @export

plot.lmls_spike <- function(x,
                            parameter = "location", 
                            plot_type = "inclusion",
                            ...){
  parameter <- match.arg(parameter, c("location", "scale"))
  plot_type <- match.arg(plot_type, c("inclusion", "posterior", "random walk"))
  
  plot_posterior <- plot_type == "posterior"
  plot_random_walk <- plot_type == "random walk"
  
  if(all(plot_posterior == FALSE, plot_random_walk == FALSE)){
    cov_names <- names(x$spike$delta[[parameter]])
    p <- colMeans(x$spike$delta[[parameter]])
    names(p) <- cov_names
    barplot(p,
            main = paste("inclusion probabilities for",parameter),
            ylab = 'average inclusion',
            ylim = c(0,1),
            axes = FALSE)
    abline(h = c(0.5, 0.75, 0.9), col = "red", lty = "dashed")
    axis(2, at = c(0, 0.25, 0.5, 0.75, 0.9, 1))
  } else {
    cov_names <- names(x$spike$coefs[[parameter]])
    k <- ncol(x$spike$coefs[[parameter]])
    n <- nrow(x$spike$coefs[[parameter]])
    ncol <- ceiling(sqrt(k))
    nrow <- ceiling(k/ncol)
    
    op <- par(mfrow = c(nrow,ncol),
              no.readonly = TRUE)
    
    for(jj in 1:k){
      dat <- x$spike$coefs[[parameter]][,jj]
      if(plot_posterior == T){
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
      if(plot_random_walk == T){
        plot(x$spike$coefs[[parameter]][,jj], 
             type = "l",
             main = paste("Random Walk of",cov_names[[jj]]),
             xlab = "iteration",
             ylab = paste('coef of',cov_names[[jj]]))
      }
    }

    par(op)
  }
}
  
  

 