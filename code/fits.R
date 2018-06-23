# Fit using FLASH -------------------------------------------------------
# Methods: 1. Vanilla, 2. Zero, 3. OHL, 4. OHF, 5. OHF+
fit_flash <- function(Y, Kmax, methods=1:5, ebnm_fn=ebnm_pn) {
  n <- nrow(Y)
  fits <- list()
  timing <- list()

  # Vanilla FLASH
  if (1 %in% methods) {
    timing$Vanilla <- list()
    data <- flash_set_data(Y)
    t0 <- Sys.time()
    fl <- flash_add_greedy(data, Kmax, var_type="by_column", ebnm_fn=ebnm_fn)
    t1 <- Sys.time()
    timing$Vanilla$greedy <- t1 - t0
    fits$Vanilla <- flash_backfit(data, fl, var_type="by_column", ebnm_fn=ebnm_fn)
    timing$Vanilla$backfit <- Sys.time() - t1
  }

  data <- flash_set_data(Y, S = 1)

  # Zero and OHL
  if (2 %in% methods || 3 %in% methods) {
    t0 <- Sys.time()
    fl <- flash_add_greedy(data, Kmax, var_type="zero", ebnm_fn=ebnm_fn)
    t1 <- Sys.time()
    if (2 %in% methods) {
      timing$Zero <- list()
      timing$Zero$greedy <- t1 - t0
      fits$Zero <- flash_backfit(data, fl, nullcheck=F, var_type="zero", ebnm_fn=ebnm_fn)
      timing$Zero$backfit <- Sys.time() - t1
    }
    if (3 %in% methods) {
      timing$OHL <- list()
      timing$OHL$greedy <- t1 - t0
      fl <- flash_add_fixed_l(data, diag(rep(1, n)), fl)
      t2 <- Sys.time()
      fits$OHL <- flash_backfit(data, fl, nullcheck=F, var_type="zero", ebnm_fn=ebnm_fn)
      timing$OHL$backfit <- Sys.time() - t2
    }
  }

  # OHF and OHF+
  if (4 %in% methods || 5 %in% methods) {
    t0 <- Sys.time()
    fl <- flash_add_fixed_l(data, diag(rep(1, n)))
    fl <- flash_backfit(data, fl, nullcheck=F, var_type="zero", ebnm_fn=ebnm_fn)
    t1 <- Sys.time()
    fl <- flash_add_greedy(data, Kmax, fl, var_type="zero", ebnm_fn=ebnm_fn)
    t2 <- Sys.time()
    if (4 %in% methods) {
      fits$OHF <- fl
      timing$OHF$backfit <- t1 - t0
      timing$OHF$greedy <- t2 - t1
    }
    if (5 %in% methods) {
      fits$OHFp <- flash_backfit(data, fl, nullcheck=F, var_type="zero", ebnm_fn=ebnm_fn)
      timing$OHFp$backfit <- (t1 - t0) + (Sys.time() - t2)
      timing$OHFp$greedy <- t2 - t1
    }
  }

  # timing$total <- Reduce(`+`, timing)

  list(fits = fits, timing = timing)
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

  # timing$total <- Reduce(`+`, timing)

  list(m = m, timing = timing)
}
