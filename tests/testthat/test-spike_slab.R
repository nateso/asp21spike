# simulate some data -----------------------------------------------

set.seed(1509)

n <- list(1000,
          300)

bets <- list(matrix(c(2,6,3,3,3,0,0,0,0,0,0)),
             matrix(c(2,6,3,3,3,3,3,3,3,0,0)))

gams <-list(matrix(c(2,6,3,3,3,0,0,0,0,0,0)),
            matrix(c(2,6,3,3,3,3,3,3,3,0,0)))

snr <- list(10, 3)

all_combis <- expand.grid(n, bets, gams, snr)
all_combis <- all_combis[rep(1:nrow(all_combis), each = 200), ]
names(all_combis) <- c("n", "bets", "gams", "snr")

spsl <-m <- list()

mu <- rep(1,11)
sigma <- rep(1,length(mu))

system.time({
for(i in 1:30){
  
  X <- matrix(NA,nrow = all_combis$n[[i]], ncol = length(mu))
  for(jj in 1:length(mu)){
    X[,jj] <- 0.1 * rnorm(all_combis$n[[i]],mu[jj],sigma[jj])
  }
  
  eta_b <- X %*% all_combis$bets[[i]]
  eta_g <- X %*% all_combis$gams[[i]]
  
  eps <- rt(all_combis$n[[i]], 5)
  
  y <- rnorm(all_combis$n[[i]], eta_b, exp(eta_g)) + (sd(eta_g)/all_combis$snr[[i]]) * eps
  summary(y)
  test_data <- cbind.data.frame(y,X)
  names(test_data) <- c("y",paste0("x",1:11))
  
  m[[i]] <- try(lmls(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
                     scale = ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11,
                     data = test_data,
                     light = FALSE,maxit = 1000), silent = TRUE)
  
  if(class(m[[i]]) == "lmls"){
    spsl[[i]] <- spike_slab(m[[i]],
                            v_0 = 0.1,
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
                            prog_bar = FALSE)
    
    print("-----------------------------------------------------------------------------")
    print(all_combis[i,])
    summary(spsl[[i]])
  } else {print("Fail")}
}})


# calculate misslcassification rate or accuracy:

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
  
  TP <- sum(mod$spike$delta$location[,actual_loc]) + sum(mod$spike$delta$scale[,actual_scl])
  TN <- sum(mod$spike$delta$location[,actual_loc == F] == 0) + sum(mod$spik$delta$scale[,actual_scl == F] == 0)
  FN <- sum(mod$spike$delta$location[,actual_loc] == 0) + sum(mod$spike$delta$scale[,actual_scl] == 0)
  FP <- sum(mod$spike$delta$location[,actual_loc == F]) + sum(mod$spike$delta$scale[,actual_scl == F])
  
  out <- matrix(c(TP,FP,FN,TN),nrow = 2, ncol = 2,
                dimnames = list(c('actual_P','actual_N'),c('pred_P','pred_N')))
  n <- (length(actual_loc) + length(actual_scl)) * nrow(mod$spike$delta$location)
  out <- out / n
  return(list(confusion = out,
              accuracy = sum(diag(out)),
              missclass = 1 - sum(diag(out))))
}

for(i in 1:30){
  print(i)
  print(calc_acc(spsl[[i]], all_combis$bets[[i]], all_combis$gams[[i]]))
}







View(all_combis)
for (i in 1:16) {
  summary()
}


all_combis
plot(test_data$x9, test_data$y)

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
                   burnin = 10,
                   coef_init = 0,
                   #always_in_loc = 'x4',
                   nsim = 1000,
                   seed = 123, 
                   prog_bar = FALSE)

summary(spsl)
plot(spsl$spike$delta$location[, 6], type = "l")
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


