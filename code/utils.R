# MSE of posterior means (FLASH) ----------------------------------------

flash_mse <- function(fl, true_Y) {
  mean((flash_get_lf(fl) - true_flash_Y)^2)
}

# flash_pm_mse <- function(fl_samp, true_Y) {
#   n <- nrow(true_Y)
#   p <- ncol(true_Y)
#   nsamp <- length(fl_samp)
#
#   post_means <- matrix(0, nrow=n, ncol=p)
#   for (i in 1:nsamp) {
#     post_means <- post_means + fl_samp[[i]]
#   }
#   post_means <- post_means / nsamp
#   sum((post_means - true_Y)^2) / (n * p)
# }

# MSE for MASH ----------------------------------------------------------
mash_mse <- function(m, true_Y) {
  mean((get_pm(m) - t(true_Y))^2)
}


# 95% CI coverage for FLASH ---------------------------------------------
flash_ci <- function(fl_samp, true_Y) {
  n <- nrow(true_Y)
  p <- ncol(true_Y)
  nsamp <- length(fl_samp)

  flat_samp <- matrix(0, nrow=n*p, ncol=nsamp)
  for (i in 1:nsamp) {
    flat_samp[, i] <- as.vector(fl_samp[[i]])
  }
  CI <- t(apply(flat_samp, 1, function(x) {quantile(x, c(0.025, 0.975))}))
  mean((as.vector(true_Y) > CI[, 1]) & (as.vector(true_Y < CI[, 2])))
}

# 95% CI coverage for MASH ----------------------------------------------
mash_ci <- function(m, true_Y) {
  Y <- t(true_Y)
  mean((Y > get_pm(m) - 1.96 * get_psd(m))
      & (Y < get_pm(m) + 1.96 * get_psd(m)))
}


# LFSR for FLASH --------------------------------------------------------
flash_lfsr <- function(fl, fl_samp, true_Y, step=0.05) {
  n <- nrow(true_Y)
  p <- ncol(true_Y)
  nsamp <- length(fl_samp)

  pp <- matrix(0, nrow=n, ncol=p)
  pn <- matrix(0, nrow=n, ncol=p)
  for (i in 1:nsamp) {
    pp <- pp + (fl_samp[[i]] >= 0)
    pn <- pn + (fl_samp[[i]] <= 0)
  }
  lfsr <- 1 - pmin(pp, pn) / nsamp

  efsr_by_lfsr(flash_get_lf(fl), true_Y, lfsr, step)
}

# LFSR for MASH ---------------------------------------------------------
mash_lfsr <- function(m, true_Y, step=0.05) {
  efsr_by_lfsr(get_pm(m), t(true_Y), get_lfsr(m), step)
}

# empirical false sign rate vs. local false sign rate
efsr_by_lfsr <- function(pm, true_Y, lfsr, step) {
  pred_signs <- sign(pm)
  pred_zeros <- pred_signs == 0
  pred_signs[pred_zeros] <- sample(c(0, 1), length(pred_zeros), replace=T)

  gotitright <- (pred_signs == sign(true_Y))

  nsteps <- floor(.5 / step)
  efsr_by_lfsr <- rep(0, nsteps)
  for (k in 1:nsteps) {
    idx <- (lfsr >= (step * (k - 1)) & lfsr < (step * k))
    efsr_by_lfsr[k] <- ifelse(sum(idx) == 0, NA,
                              1 - sum(gotitright[idx]) / sum(idx))
  }
  efsr_by_lfsr
}
