# Make sure to use branch "trackObj" when loading flashr.

# devtools::install_github("stephenslab/flashr", ref="trackObj")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

library(mashr)

source("./code/fits.R")
source("./code/utils.R")
source("./code/gtexutils.R")

fpath <- "./output/MASHvFLASHgtex/"
# fpath <- "./output/MASHvFLASHrandom/"

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

data <- t(gtex$strong.z)
# data <- t(gtex$random.z)

seeds <- 1:4
pcts <- c(0.01, 0.05, 0.1, 0.25)

fl_fits <- c(function(fl_data) {fit_flash_zero(fl_data, Kmax = 50,
                                               ebnm_fn = "ebnm_pn",
                                               init_fn = "udv_si_svd",
                                               backfit = FALSE,
                                               warmstart = TRUE)},
             function(fl_data) {fit_flash_OHF(fl_data, Kmax = 50,
                                              ebnm_fn = "ebnm_pn",
                                              init_fn = "udv_si_svd",
                                              backfit = FALSE,
                                              warmstart = TRUE)},
             function(fl_data) {fit_flash_zero(fl_data, Kmax = 50,
                                               ebnm_fn = "ebnm_ash",
                                               init_fn = "udv_si_svd",
                                               backfit = FALSE,
                                               warmstart = TRUE)},
             function(fl_data) {fit_flash_OHF(fl_data, Kmax = 50,
                                              ebnm_fn = "ebnm_ash",
                                              init_fn = "udv_si_svd",
                                              backfit = FALSE,
                                              warmstart = TRUE)})
fl_fit_names <- c("pn.zero", "pn.OHF", "ash.zero", "ash.OHF")

fl_res <- list()
fl_diag <- list()
fl_t <- list()

for (i in 1:length(seeds)) {
  message(paste0("Running experiment #", i))
  missing <- get_missing(data, seed = seeds[i], pct_missing = pcts[i])
  fl_data <- set_fl_data_with_missing(data, missing)

  fl_res[[i]] <- list()
  fl_diag[[i]] <- list()
  fl_t[[i]] <- list()

  for (j in 1:length(fl_fits)) {
    message(paste("Fitting", fl_fit_names[j]))
    t0 <- Sys.time()
    fl_fit <- fl_fits[[j]](fl_data)
    fl_res[[i]][[fl_fit_names[j]]] <- fl_fit$f
    diagnostics <- list(mse = imputation_mse(flash_get_fitted_values(fl_fit$f),
                                             data, missing),
                        ci = flash_imputation_ci(fl_data, fl_fit$f,
                                                 data, missing))
    fl_diag[[i]][[fl_fit_names[j]]] <- diagnostics
    fl_t[[i]][[fl_fit_names[j]]] <- Sys.time() - t0
  }

  saveRDS(fl_res, paste0(fpath, "res.rds"))
  saveRDS(fl_diag, paste0(fpath, "fldiag.rds"))
  saveRDS(fl_t, paste0(fpath, "flt.rds"))
}

m_res <- list()
m_diag <- list()
m_t <- list()

for (i in 1:length(seeds)) {
  message(paste0("Running experiment #", i))
  missing <- get_missing(data, seed = seeds[i], pct_missing = pcts[i])
  m_data <- set_m_data_with_missing(t(data), t(missing))

  message("Fitting mash")
  t0 <- Sys.time()
  m_fit <- fit_mash(m_data)
  m_res[[i]] <- m_fit$m
  diagnostics <- list(mse = imputation_mse(get_pm(m_fit$m),
                                           t(data), t(missing)),
                      ci = mash_imputation_ci(m_fit$m, t(data),
                                              t(missing)))
  m_diag[[i]] <- diagnostics
  m_t[[i]] <- Sys.time() - t0
  saveRDS(m_res, paste0(fpath, "mres.rds"))
  saveRDS(m_diag, paste0(fpath, "mdiag.rds"))
  saveRDS(m_t, paste0(fpath, "mt.rds"))
}
