% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{mcmc}
\alias{mcmc}
\title{MCMC inference for location-scale regression}
\usage{
mcmc(m, nsim = 1000, stepsize = sqrt(3) * (m$df)^(-1/6))
}
\arguments{
\item{m}{A location-scale regression model from the \code{\link[=lmls]{lmls()}} function.}

\item{nsim}{The number of MCMC samples to draw.}

\item{stepsize}{The step size of the MMALA update.}
}
\description{
A Markov chain Monte Carlo (MCMC) sampler for location-scale regression
models from the \code{\link[=lmls]{lmls()}} function. The sampler uses a Gibbs update for the
location coefficients and the Riemann manifold Metropolis-adjusted Langevin
algorithm (MMALA) from Girolami and Calderhead (2011) with the Fisher-Rao
metric tensor for the scale coefficients. The priors for the regression
coefficients are assumed to be flat.
}
\references{
Girolami, M. and Calderhead, B. (2011), Riemann manifold Langevin and
Hamiltonian Monte Carlo methods. Journal of the Royal Statistical Society:
Series B (Statistical Methodology), 73: 123-214.
\url{https://doi.org/10.1111/j.1467-9868.2010.00765.x}
}
