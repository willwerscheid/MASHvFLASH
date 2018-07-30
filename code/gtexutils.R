get_missing <- function(data, seed, pct_missing = 0.05) {
  set.seed(seed)
  missing <- rbinom(length(data), 1, prob = pct_missing)
  return(matrix(missing, nrow=nrow(data), ncol=ncol(data)))
}

set_fl_data_with_missing <- function(data, missing) {
  data[missing] <- NA
  return(flash_set_data(data, S = 1))
}

set_m_data_with_missing <- function(data, missing) {
  data[missing] <- 0
  Shat <- matrix(1, nrow=nrow(data), ncol=ncol(data))
  Shat[missing] <- 1e6
  return(mash_set_data(data, Shat))
}

imputation_mse <- function(fitted, true_Y, missing) {
  sq_error <- (fitted - true_Y)^2
  sq_error <- sq_error[missing == 1]
  return(mean(sq_error))
}

flash_imputation_ci <- function(data, fl, true_Y, missing, nsamp = 200) {
  true_Y <- as.matrix(true_Y[missing == 1])

  fl_sampler <- flash_sampler(data, fl, fixed="loadings")
  fl_samp <- fl_sampler(nsamp)
  fl_samp <- lapply(fl_samp, function(samp){as.matrix(samp[missing == 1])})
  return(flash_ci(fl_samp, true_Y))
}

mash_imputation_ci <- function(m, true_Y, missing) {
  pm <- get_pm(m)[missing == 1]
  psd <- get_psd(m)[missing == 1]
  Y <- true_Y[missing == 1]
  mean((Y > pm - 1.96 * psd)
       & (Y < pm + 1.96 * psd))
}
