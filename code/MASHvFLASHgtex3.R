# devtools::install_github("stephenslab/flashr")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

library(mashr)
source("./code/utils.R")
source("./code/gtexanalysis.R")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

strong <- t(gtex$strong.z)
random <- t(gtex$random.z)

strong_data <- flash_set_data(strong, S = 1)
random_data <- flash_set_data(random, S = 1)

nn <- readRDS("./output/MASHvFLASHnn/fl.rds")
multi <- c(2, 5, 6, 8, 11:13, 17, 22:25, 31)
dd <- nn$EL[, multi]
dd <- dd / rep(apply(dd, 2, max), each=n) # Normalize data-driven loadings


## FLASH fit

n <- nrow(strong)
canonical <- cbind(rep(1, n), diag(rep(1, n)))
LL <- cbind(canonical, dd)

fl_random <- flash_add_fixed_loadings(random_data, LL)
system.time(
  fl_random <- flash_backfit(random_data,
                             fl_random,
                             var_type = "zero",
                             ebnm_fn = "ebnm_pn",
                             nullcheck = FALSE)
)
# 26 minutes (62 iterations)
saveRDS(fl_random$gf, "./output/MASHvFLASHgtex3/flgf.rds")

for (k in 2:45) {
  fl_random$gf[[k]]$a <- min(fl_random$gf[[k]]$a, 1)
}

fl_final <- flash_add_fixed_loadings(strong_data, LL)
ebnm_param_f = lapply(fl_random$gf, function(g) {list(g=g, fixg=TRUE)})
system.time(
  fl_final <- flash_backfit(strong_data,
                            fl_final,
                            var_type = "zero",
                            ebnm_fn = "ebnm_pn",
                            ebnm_param = list(f = ebnm_param_f, l = list()),
                            nullcheck = FALSE)
)
# 20 minutes (75 iterations)
saveRDS(fl_final, "./output/MASHvFLASHgtex3/fl.rds")

nsamp <- 200
system.time({
  sampler <- flash_sampler(strong_data, fl_final, fixed="loadings")
  samp <- sampler(200)
  fl_lfsr <- flash_lfsr(samp)
})
# 2 minutes
saveRDS(fl_lfsr, "./output/MASHvFLASHgtex3/fllfsr.rds")


# MASH fit

strong_data <- mash_set_data(t(strong), Shat = 1)
random_data <- mash_set_data(t(random), Shat = 1)

U.flash <- list()
for (i in 1:ncol(dd)) {
  U.flash[[i]] <- outer(dd[, i], dd[, i])
}
U.c <- cov_canonical(random_data)
system.time(m_random <- mash(random_data, Ulist = c(U.c, U.flash)))
# 4 minutes

system.time(
  m_final <- mash(strong_data, g=get_fitted_g(m_random), fixg=TRUE)
)
saveRDS(m_final, "./output/MASHvFLASHgtex3/m.rds")
# 19 seconds


m_lfsr <- t(get_lfsr(m_final))
m_pm <- t(get_pm(m_final))
fl_pm <- flash_get_fitted_values(fl_final)
mash.v.flash <- compare_methods(fl_lfsr, m_lfsr, fl_pm, m_pm)
for (n in mash.v.flash$diff_pms) {
  plot_test(n, fl_lfsr, fl_pm, "FLASH")
  plot_test(n, m_lfsr, m_pm, "MASH")
}
