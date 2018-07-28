# Fit using FLASH -------------------------------------------------------
fit_flash_zero <- function(Y, Kmax, ebnm_fn="ebnm_pn",
                           init_fn="udv_svd", greedy=TRUE, backfit=TRUE,
                           warmstart=TRUE) {
  n <- nrow(Y)
  data <- flash_set_data(Y, S = 1)

  t0 <- Sys.time()
  if (greedy) {
    res <- flash_add_greedy(data, Kmax, var_type="zero",
                            ebnm_fn=ebnm_fn, init_fn=init_fn,
                            warmstart=warmstart)
    fl <- res$f
  } else {
    fl <- flash_add_factors_from_data(data, Kmax, init_fn=init_fn)
  }
  t1 <- Sys.time()
  if (backfit) {
    res <- flash_backfit(data, fl, var_type="zero",
                         ebnm_fn=ebnm_fn)
    fl <- res$f
  }
  t2 <- Sys.time()

  t.greedy <- t1 - t0
  t.backfit <- t2 - t1

  list(f = fl, t.greedy = t.greedy, t.backfit = t.backfit)
}


fit_flash_OHL <- function(Y, Kmax, ebnm_fn="ebnm_pn",
                          init_fn="udv_svd", greedy=TRUE, backfit=TRUE,
                          warmstart=TRUE) {
  n <- nrow(Y)
  data <- flash_set_data(Y, S = 1)
  canonical <- cbind(rep(1, n), diag(rep(1, n)))

  zero.res <- fit_flash_zero(Y, Kmax, ebnm_fn, init_fn,
                             greedy, backfit=FALSE, warmstart)

  fl <- flash_add_fixed_l(data, canonical, zero.res$f, init_fn=init_fn)

  t0 <- Sys.time()
  if (backfit) {
    res <- flash_backfit(data, fl,
                         var_type="zero", ebnm_fn=ebnm_fn,
                         nullcheck=FALSE)
    fl <- res$f
  } else {
    K <- flashr:::flash_get_k(fl)
    res <- flash_backfit(data, fl, kset=(K - ncol(canonical) + 1):K,
                         var_type="zero", ebnm_fn=ebnm_fn,
                         nullcheck=FALSE)
    fl <- res$f
  }
  t1 <- Sys.time()

  t.backfit <- Sys.time() - t0

  list(f = fl, t.greedy = zero.res$t.greedy, t.backfit = t.backfit)
}


fit_flash_OHF <- function(Y, Kmax, ebnm_fn="ebnm_pn",
                          init_fn="udv_svd", greedy=TRUE, backfit=TRUE,
                          warmstart=TRUE) {
  n <- nrow(Y)
  data <- flash_set_data(Y, S = 1)
  canonical <- cbind(rep(1, n), diag(rep(1, n)))

  t0 <- Sys.time()
  fl <- flash_add_fixed_l(data, canonical)
  res <- flash_backfit(data, fl, var_type="zero", ebnm_fn=ebnm_fn,
                       nullcheck=FALSE)
  fl <- res$f
  t1 <- Sys.time()

  if (greedy) {
    res <- flash_add_greedy(data, Kmax, fl, var_type="zero",
                            ebnm_fn=ebnm_fn, init_fn=init_fn,
                            warmstart=warmstart)
    fl <- res$f
  } else {
    fl <- flash_add_factors_from_data(data, Kmax, fl, init_fn=init_fn)
  }
  t2 <- Sys.time()
  K = flashr:::flash_get_k(fl)
  if (backfit && K > ncol(canonical)) {
    res <- flash_backfit(data, fl,
                         kset=(ncol(canonical) + 1):K,
                         var_type="zero", ebnm_fn=ebnm_fn,
                         nullcheck=FALSE)
    fl <- res$f
  }
  t3 <- Sys.time()

  t.greedy <- t2 - t1
  t.backfit <- (t1 - t0) + (t3 - t2)

  list(f = fl, t.greedy = t.greedy, t.backfit = t.backfit)
}


# Fit using MASH -------------------------------------------------------
fit_mash <- function(Y) {
  data <- mash_set_data(t(Y))
  timing <- list()

  # time to create canonical matrices is negligible
  U = cov_canonical(data)

  t0 <- Sys.time()
  m.1by1 <- mash_1by1(data)
  lvl <- 0.05
  strong <- get_significant_results(m.1by1, lvl)
  while (length(strong) < 5 && lvl < 0.5) {
    lvl <- lvl + 0.05
    strong <- get_significant_results(m.1by1, lvl)
  }
  if (length(strong) >= 5) {
    U.pca <- cov_pca(data, 5, strong)
    U.ed <- cov_ed(data, U.pca, strong)
    U <- c(U, U.ed)
    t.ed <- Sys.time() - t0
  } else {
    t.ed <- as.difftime(0, units="secs")
  }

  t0 <- Sys.time()
  m <- mash(data, U)
  t.mash <- Sys.time() - t0

  list(m = m, t.ed = t.ed, t.mash = t.mash)
}
