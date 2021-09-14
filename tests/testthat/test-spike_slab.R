# simulate some data -----------------------------------------------
set.seed(1509)
n <- 100

mu <- rep(1,5)
sigma <- rep(1,length(mu))

X <- matrix(NA,nrow = n, ncol = length(mu))
for(jj in 1:length(mu)){
  X[,jj] <- rnorm(n,mu[jj],sigma[jj])
}

bets <- matrix(c(2,6,3,0,0))
gams <- matrix(c(2,6,3,0,0))

y <- rnorm(n, X %*% bets, exp(-3 + X %*% gams))
y_noise <- y + rnorm(n,0,1)

y <- y_noise

test_data <- cbind.data.frame(y,X)
names(test_data) <- c("y",paste0("x",1:5))

m <- lmls(y ~ x1 + x2 + x3 + x4 + x5,
          scale = ~ x1 + x2 + x3 + x4 + x5,
          data = test_data,
          light = FALSE,maxit = 1000)

# test spike_slab -------------------------------------------

spsl <- spike_slab(m,
                   v_0 = 0.1,
                   a_theta_loc = 1,
                   b_theta_loc = 1,
                   a_theta_scl = 1,
                   b_theta_scl = 1,
                   a_tau_loc = 5,
                   b_tau_loc = 25,
                   a_tau_scl = 5,
                   b_tau_scl = 25,
                   burnin = 100,
                   coef_init = 0,
                   #always_in_loc = 'x4',
                   nsim = 1000,
                   seed = 123, 
                   prog_bar = FALSE)

summary(spsl)
plot(spsl$spike$delta$location[, 4], type = "l")
plot(spsl, "location", "inc")
plot(spsl, "location", "post")
plot(spsl, "location", "rand")

plot(test_data$x1, test_data$y)


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


