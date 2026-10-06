# The model returned by TDS_MGWR() is the best sweep (by RMSE), not necessarily
# the last one (non-monotone backfitting since 1.4.1). Its standard errors must
# be computed from the same state as Betav, Shat, R_k, tS and AICc.

library(testthat)
library(mgwrsar)

test_that("tds_mgwr: sev is built from the returned R_k (best state)", {
  set.seed(2)  # n = 400: best sweep 8 of 11 with mgwrsar 1.4.2
  n <- 400L
  u <- runif(n); v <- runif(n)
  dat <- data.frame(X2 = rnorm(n), X3 = rnorm(n), X4 = rnorm(n))
  dat$Y <- 3 * (u + v) + 2 * u * sin(4 * v) * dat$X2 +
    4 * u * sin(6 * v) * dat$X3 + 2 * (u - 0.5) * dat$X4 + rnorm(n, sd = 1.5)
  m <- suppressMessages(TDS_MGWR(
    Y ~ X2 + X3 + X4, data = dat, coords = cbind(u, v), Model = "tds_mgwr",
    kernels = "gauss",
    control_tds = list(nns = 10, get_AIC = TRUE, ncore = 1, init_model = "OLS"),
    control = list(adaptive = TRUE, NN = n, Type = "GD", SE = TRUE)))
  sigma2 <- sum(m@residuals^2) / (n - m@tS)
  for (k in colnames(m@sev)) {
    se_k <- sqrt(sigma2 * rowSums(m@R_k[[k]]^2) / m@X[, k]^2)
    expect_equal(unname(m@sev[, k]), unname(se_k), tolerance = 1e-10, info = k)
  }
  expect_equal(m@tS, sum(diag(m@Shat)))
  expect_equal(m@AICc, mgwrsar:::aicc_f(m@residuals, m@tS, n))
})
