# Evaluate methods based on MSE, CI coverage, and TPR vs. FPR -----------

flash_diagnostics <- function(fl, Y, true_Y, nsamp) {
  MSE <- flash_mse(fl, true_Y)

  # Sample from FLASH fit to estimate CI coverage and TPR vs. FPR
  fl_sampler <- flash_lf_sampler(Y, fl, ebnm_fn=ebnm_pn, fixed="loadings")
  fl_samp <- fl_sampler(nsamp)

  CI <- flash_ci(fl_samp, true_Y)
  ROC <- flash_roc(fl, fl_samp, true_Y)

  list(MSE = MSE, CI = CI, TP = ROC$TP, FP = ROC$FP,
       n_nulls = ROC$n_nulls, n_nonnulls = ROC$n_nonnulls)
}

mash_diagnostics <- function(m, true_Y) {
  MSE <- mash_mse(m, true_Y)
  CI <- mash_ci(m, true_Y)
  ROC <- mash_roc(m, true_Y)

  list(MSE = MSE, CI = CI, TP = ROC$TP, FP = ROC$FP,
       n_nulls = ROC$n_nulls, n_nonnulls = ROC$n_nonnulls)
}


# MSE of posterior means (FLASH) ----------------------------------------
flash_mse <- function(fl, true_Y) {
  mean((flash_get_lf(fl) - true_Y)^2)
}

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
  mean((as.vector(true_Y) >= CI[, 1]) & (as.vector(true_Y) <= CI[, 2]))
}

# 95% CI coverage for MASH ----------------------------------------------
mash_ci <- function(m, true_Y) {
  Y <- t(true_Y)
  mean((Y > get_pm(m) - 1.96 * get_psd(m))
      & (Y < get_pm(m) + 1.96 * get_psd(m)))
}


# LFSR for FLASH --------------------------------------------------------
flash_lfsr <- function(fl_samp) {
  nsamp <- length(fl_samp)
  n <- nrow(fl_samp[[1]])
  p <- ncol(fl_samp[[1]])

  pp <- matrix(0, nrow=n, ncol=p)
  pn <- matrix(0, nrow=n, ncol=p)
  for (i in 1:nsamp) {
    pp <- pp + (fl_samp[[i]] > 0)
    pn <- pn + (fl_samp[[i]] < 0)
  }
  1 - pmax(pp, pn) / nsamp
}


# Quantities for plotting ROC curves -----------------------------------
flash_roc <- function(fl, fl_samp, true_Y, step=0.01) {
  roc_data(flash_get_lf(fl), true_Y, flash_lfsr(fl_samp), step)
}

mash_roc <- function(m, true_Y, step=0.01) {
  roc_data(get_pm(m), t(true_Y), get_lfsr(m), step)
}

roc_data <- function(pm, true_Y, lfsr, step) {
  correct_sign <- pm * true_Y > 0
  is_null <- true_Y == 0
  n_nulls <- sum(is_null)
  n_nonnulls <- length(true_Y) - n_nulls

  ts <- seq(0, 1, by=step)
  tp <- rep(0, length(ts))
  fp <- rep(0, length(ts))

  for (t in 1:length(ts)) {
    signif <- lfsr <= ts[t]
    tp[t] <- sum(signif & correct_sign)
    fp[t] <- sum(signif & is_null)
  }

  list(ts = ts, TP = tp, FP = fp, n_nulls = n_nulls, n_nonnulls = n_nonnulls)
}


# empirical false sign rate vs. local false sign rate
# efsr_by_lfsr <- function(pm, true_Y, lfsr, step) {
#   pred_signs <- sign(pm)
#   pred_zeros <- pred_signs == 0
#   pred_signs[pred_zeros] <- sample(c(0, 1), length(pred_zeros), replace=T)
#
#   gotitright <- (pred_signs == sign(true_Y))
#
#   nsteps <- floor(.5 / step)
#   efsr_by_lfsr <- rep(0, nsteps)
#   for (k in 1:nsteps) {
#     idx <- (lfsr >= (step * (k - 1)) & lfsr < (step * k))
#     efsr_by_lfsr[k] <- ifelse(sum(idx) == 0, NA,
#                               1 - sum(gotitright[idx]) / sum(idx))
#   }
#   efsr_by_lfsr
# }
