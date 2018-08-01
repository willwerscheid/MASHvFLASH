## @knitr gtex2
# Make sure to use branch "trackObj" when loading flashr.

# devtools::install_github("stephenslab/flashr", ref="trackObj")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

library(mashr)

source("./code/fits.R")
source("./code/utils.R")

fpath <- "./output/MASHvFLASHgtex2/"

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

strong <- t(gtex$strong.z)
random <- t(gtex$random.z)


# FLASH fits ------------------------------------------------------------

strong_data <- flash_set_data(strong, S = 1)
random_data <- flash_set_data(random, S = 1)

# methods to get DD: OHF; Zero(top5, top10, top20, full)
# LLs to use: canonical + OHF; Zero; canonical + top5/top10/top15; full

# Step 1. Learn data-driven loadings using strong tests.
t0 <- Sys.time()
fl_strong.OHF <- fit_flash_OHF(strong_data, Kmax=50, backfit=FALSE)
t_strong.OHF <- Sys.time() - t0

t0 <- Sys.time()
fl_strong.zero <- fit_flash_zero(strong_data, Kmax=50, backfit=FALSE)
t_strong.zero <- Sys.time() - t0

# Step 2. Fit the model to random tests to learn priors on factors.
LL <- list()
LL[[1]] <- fl_strong.OHF$f$EL

n <- nrow(strong)
canonical <- cbind(rep(1, n), diag(rep(1, n)))
LL[[2]] <- cbind(canonical, fl_strong.zero$f$EL[, 1:5])
LL[[3]] <- cbind(canonical, fl_strong.zero$f$EL[, 1:10])
LL[[4]] <- cbind(canonical, fl_strong.zero$f$EL[, 1:20])
LL[[5]] <- fl_strong.zero$f$EL

fl_random <- list()
t_random <- list()
for (i in 1:length(LL)) {
  t0 <- Sys.time()
  fl <- flash_add_fixed_l(random_data, LL[[i]])
  fl <- flash_backfit(random_data,
                      fl,
                      var_type="zero",
                      ebnm_fn = "ebnm_pn",
                      nullcheck = FALSE,
                      warmstart = TRUE,
                      verbose = TRUE)
  t_random[[i]] <- Sys.time() - t0
  fl_random[[i]] <- fl$f
}

# Step 3. Compute posteriors on strong tests, using priors from step 2.

fl_final <- list()
t_final <- list()
for (i in 1:length(LL)) {
  t0 <- Sys.time()
  fl <- flash_add_fixed_l(strong_data, LL[[i]])
  ebnm_param_f = lapply(fl_random[[i]]$gf, function(g) {list(g=g, fixg=TRUE)})
  fl <- flash_backfit(strong_data,
                      fl,
                      var_type="zero",
                      ebnm_fn = "ebnm_pn",
                      ebnm_param = list(f = ebnm_param_f, l = list()),
                      nullcheck = FALSE,
                      verbose = TRUE)
  t_final[[i]] <- Sys.time() - t0
  fl_final[[i]] <- fl$f
}

saveRDS(fl_final[[1]], paste0(fpath, "OHF.rds"))
saveRDS(fl_final[[2]], paste0(fpath, "top5.rds"))
saveRDS(fl_final[[3]], paste0(fpath, "top10.rds"))
saveRDS(fl_final[[4]], paste0(fpath, "top20.rds"))
saveRDS(fl_final[[5]], paste0(fpath, "zero.rds"))

t <- list()
t$strong <- list()
t$strong$OHF <- t_strong.OHF
t$strong$zero <- t_strong.zero
t$random <- t_random
t$final <- t_final


# Sample from posterior to get LFSR for FLASH fits ----------------------

nsamp <- 200
fl_lfsr <- list()
for (i in 1:length(fl_final)) {
  sampler <- flash_sampler(strong_data, fl_final[[i]], fixed="loadings")
  samp <- sampler(200)
  fl_lfsr[[i]] <- flash_lfsr(samp)
}

saveRDS(fl_lfsr, paste0(fpath, "fllfsr.rds"))


# MASH fit --------------------------------------------------------------

strong_data <- mash_set_data(t(strong), Shat = 1)
random_data <- mash_set_data(t(random), Shat = 1)

t$mash <- list()

# Step 1. Learn data-driven loadings using strong tests.
t0 <- Sys.time()
U.pca <- cov_pca(strong_data, 5)
U.ed <- cov_ed(strong_data, U.pca)
t$mash$ed <- Sys.time() - t0

# 2. Fit the model to random tests to learn mixture weights
t0 <- Sys.time()
U.c <- cov_canonical(random_data)
m_random <- mash(random_data, Ulist = c(U.ed,U.c))
t$mash$random <- Sys.time() - t0

# 3. Compute posterior summaries on the strong tests
t0 <- Sys.time()
m_final <- mash(strong_data, g=get_fitted_g(m_random), fixg=TRUE)
t$mash$final <- Sys.time() - t0

saveRDS(m_final, paste0(fpath, "m.rds"))
saveRDS(t, paste0(fpath, "t.rds"))
