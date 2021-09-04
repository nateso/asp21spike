# simulate some data -----------------------------------------------

set.seed(1509)

n <- list(200,
          1000)

bets <- list(matrix(c(rep(1,4), rep(2,4), rep(3,4), rep(0,4))),
             matrix(c(rep(1,4), rep(0, 12))))

gams <-list(matrix(c(rep(1,4), rep(2,4), rep(3,4), rep(0,4))),
            matrix(c(rep(1,4), rep(0, 12))))

snr <- list(2, 5)

all_combis <- expand.grid(n, bets, gams, snr)
names(all_combis) <- c("n", "bets", "gams", "snr")


for(i in 1:length(all_combis)){
  X <- matrix(NA,nrow = all_combis$n[[i]], ncol = 16)
  for(jj in 1:16){
    X[,jj] <- runif(n,-2,2)
  }
  
  eta_b <- X %*% all_combis$bets[[i]]
  eta_g <- X %*% all_combis$gams[[i]]
  
  eps <- rt(n, 5)
  
  y <- rnorm(n, eta_b, exp(eta_g)) + (sd(eta_g)/all_combis$snr[[i]]) * eps
  
  test_data <- cbind.data.frame(y,X)
  names(test_data) <- c("y",paste0("x",1:16))
  
  m <- lmls(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x12 + x13 + x14 + x15 + x16,
            scale = ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x12 + x13 + x14 + x15 + x16,
            data = test_data,
            light = FALSE,maxit = 1000)
  return(m)
}


# test spike_slab -------------------------------------------

spsl <- spike_slab(m,
                   v_0 = 0.05,
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
                   seed = 123, 
                   prog_bar = FALSE)
# summary(spsl)
# plot(spsl$spike$delta$location[, 4], type = "l")
# plot(spsl, "location", "inc")
# plot(spsl, "location", "post")
# plot(spsl, "location", "rand")
# barplot(colMeans(spsl$spike$delta$location))


# set.seed(243)
# spsl <- spike_slab(m,v_0 = 0.005,a_tau_loc = 3,b_tau_loc = 25,nsim = 5000)
# plot(spsl,'location')
# 
# set.seed(243)
# spsl <- spike_slab(m,v_0 = 0.005,a_tau_loc = 3,b_tau_loc = 25,nsim = 10000)
# plot(spsl,'location')
# 
# 
# set.seed(243)
# spsl <- spike_slab(m,v_0 = 0.005,a_tau_loc = 3,b_tau_loc = 25,nsim = 50000)
# plot(spsl,'location')

# set.seed(1312424)
# n <- 200
# snr <- 3
# sm1 <- runif(n)
# fsm1 <- dbeta(sm1, 7, 3)/2
# sm2 <- runif(n, 0, 1)
# f <- gl(3, n/3)
# ff <- as.numeric(f)/2
# fsm2f <- ff + ff * sm2 + ((f == 1) * -dbeta(sm2, 6, 4) +
#                                + (f == 2) * dbeta(sm2, 6, 9) + (f == 3) * dbeta(sm2,
#                                                                                 + 9, 6))/2
# lin <- matrix(rnorm(n * 3), n, 3)
# colnames(lin) <- paste("lin", 1:3, sep = "")
# noise1 <- sm1 + rnorm(n)
# noise2 <- runif(n)
# noise3 <- runif(n)
# noise4 <- sample(gl(4, n/4))
# eta <- fsm1 + fsm2f + lin %*% c(0.1, 0.2, 0.3)
# y <- eta + sd(eta)/snr * rt(n, df = 5)
# d <- data.frame(y, sm1, sm2, f, lin, noise1, noise2,
#                    noise3, noise4)
# f1 <- y ~ (sm1 + sm2 + f + lin1)^2 + lin2 + lin3 + noise1 + 
#   noise2 + noise3 + noise4
# 
# m <- lmls(location = f1,scale = f1,
#      data = d,
#      light = FALSE,maxit = 1000)
# 
# res <- spike_slab(m)
# colMeans(res$spike$delta$location)


# check helper functions for update
delta <- rep(c(0,1),2)
bets <- c(6,0.1,-6,-0.1)
a_tau <- 5 
b_tau <- 50
a_theta <- 1
b_theta <- 1
v_0 <- 0.005
theta <- 0.5

tau <- sample_tau(delta, a_tau, b_tau, v_0)

test_that('sample_tau() gives a vector of correct length',{
  expect_true(length(tau) == length(delta))
})

test_that('sample_tau() gives small values if delta = 0 and large if delta = 1',{
  expect_true(all(tau[delta == 0] < 0.1))
  expect_true(all(tau[delta == 1] > 0.1))
})

delta_update <- update_delta(theta,tau,a_tau,b_tau,v_0)
test_that("update_delta() yields one or zeros",{
  expect_true(all(delta_update %in% c(0,1)))
})

theta_update <- update_theta(delta,a_theta, b_theta)
test_that("update_theta() yields one value between 0 and 1",{
  expect_true(all(theta_update >= 0,theta_update <= 1))
  expect_true(length(theta_update) == 1)
})

tau_update <- update_tau(delta,bets,a_theta,b_theta,v_0)


