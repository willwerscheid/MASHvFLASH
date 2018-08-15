# devtools::install_github("stephenslab/flashr")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

library(mashr)
source("./code/utils.R")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

strong <- t(gtex$strong.z)
random <- t(gtex$random.z)

strong_data <- flash_set_data(strong, S = 1)
random_data <- flash_set_data(random, S = 1)

# The multi-tissue effects from the previous analysis will serve as
#   the data-driven loadings:
nn <- readRDS("./output/MASHvFLASHnn/fl.rds")
multi <- c(2, 5, 6, 8, 11:13, 17, 22:25, 31)
n <- nrow(strong)
dd <- nn$EL[, multi]
dd <- dd / rep(apply(dd, 2, max), each=n) # normalize


## FLASH fit ------------------------------------------------------------

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

gf <- readRDS("./output/MASHvFLASHgtex3/flgf.rds")
fl_final <- flash_add_fixed_loadings(strong_data, LL)
ebnm_param_f = lapply(gf, function(g) {list(g=g, fixg=TRUE)})
system.time(
  fl_final <- flash_backfit(strong_data,
                            fl_final,
                            var_type = "zero",
                            ebnm_fn = "ebnm_pn",
                            ebnm_param = list(f = ebnm_param_f, l = list()),
                            nullcheck = FALSE)
)
# 21 minutes (75 iterations)
saveRDS(fl_final, "./output/MASHvFLASHgtex3/fl.rds")

nsamp <- 200
system.time({
  sampler <- flash_sampler(strong_data, fl_final, fixed="loadings")
  samp <- sampler(200)
  fl_lfsr <- flash_lfsr(samp)
})
# 2 minutes
saveRDS(fl_lfsr, "./output/MASHvFLASHgtex3/fllfsr.rds")


# MASH fit --------------------------------------------------------------

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
