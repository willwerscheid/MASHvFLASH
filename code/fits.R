# Fit using FLASH -------------------------------------------------------
fit_flash <- function(Y, Kmax, method) {
  n <- nrow(Y)
  data <- flash_set_data(Y, S = 1)
  timing <- list()

  t0 <- Sys.time()
  if (method %in% c("OHF", "OHFplus")) {
    fl <- flash_add_fixed_l(data, diag(rep(1, n)))
    fl <- flash_backfit(data, fl, nullcheck = F, var_type = "zero")
    t1 <- Sys.time()
    timing$backfit <- t1 - t0
    fl <- flash_add_greedy(data, Kmax, fl, var_type = "zero")
    timing$greedy <- Sys.time() - t1
    if (method == "OHFplus") {
      t2 <- Sys.time()
      fl <- flash_backfit(data, fl, nullcheck = F, var_type = "zero")
      timing$backfit <- timing$backfit + (Sys.time() - t2)
    }
  } else {
    fl <- flash_add_greedy(data, Kmax, var_type = "zero")
    t1 <- Sys.time()
    timing$greedy <- t1 - t0
    if (method == "OHL") {
      fl <- flash_add_fixed_l(data, diag(rep(1, n)), fl)
    }
    fl <- flash_backfit(data, fl, nullcheck = F, var_type = "zero")
    timing$backfit <- Sys.time() - t1
  }

  timing$total <- Reduce(`+`, timing)

  list(fl = fl, timing = timing)
}


# Fit using MASH -------------------------------------------------------
fit_mash <- function(Y, ed=T) {
  data <- mash_set_data(t(Y))
  timing <- list()

  # time to create canonical matrices is negligible
  U = cov_canonical(data)

  if (ed) {
    t0 <- Sys.time()
    m.1by1 <- mash_1by1(data)
    strong <- get_significant_results(m.1by1, 0.05)
    U.pca <- cov_pca(data, 5, strong)
    U.ed <- cov_ed(data, U.pca, strong)
    U <- c(U, U.ed)
    timing$ed <- Sys.time() - t0
  }

  t0 <- Sys.time()
  m <- mash(data, U)
  timing$mash <- Sys.time() - t0

  timing$total <- Reduce(`+`, timing)

  list(m = m, timing = timing)
}
