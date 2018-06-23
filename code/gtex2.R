devtools::load_all("/Users/willwerscheid/GitHub/flashr2/")
library(mashr)
source("./code/fits.R")
source("./code/utils.R")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
strong <- gtex$strong.z
random <- gtex$random.z

# MASH
mdata.strong <- mash_set_data(strong, Shat = 1)
mdata.random <- mash_set_data(random, Shat = 1)

U.pca <- cov_pca(mdata.strong, 5)
U.ed <- cov_ed(mdata.strong, U.pca)
U.c <- cov_canonical(mdata.random)
m <- mash(mdata.random, Ulist = c(U.ed,U.c), outputlevel = 1)
m2 <- mash(mdata.strong, g=get_fitted_g(m), fixg=TRUE)
saveRDS(m2, "./output/gtex2mfit.rds")
m.lfsr <- get_lfsr(m2)

fldata.strong <- flash_set_data(t(strong), S = 1)
fldata.random <- flash_set_data(t(random), S = 1)
# 1. Learn "data-driven" loadings using "strong" tests (use ebnm_pn
#   for speed).
# No backfit...
fl.strong <- flash(fldata.strong, Kmax=50, var_type="zero", verbose=T)
saveRDS(fl.strong, "./output/gtexstrongfit.rds")

# 2. Fit the model to the random tests to learn the priors on the factors.
LL = cbind(flash_get_l(fl.strong), diag(rep(1, 44)))
fl.random <- flash_add_fixed_l(fl.random, LL = LL)
fl.random <- flash_backfit(fldata.random, fl.random, var_type="zero",
                           ebnm_fn = ebnm_ash, nullcheck=F, verbose=T)
saveRDS(fl.random, "./output/gtexrandomfit.rds")

# Compute posterior summaries on the strong tests, using g_f from step 2.
fl <- flash_add_fixed_l(fldata.strong, LL = LL)
fl <- flash_backfit(fldata.strong, fl, var_type="zero", ebnm_fn = ebnm_ash,
                    gf=flash_get_gf(fl.random), nullcheck=F, verbose=T)
saveRDS(fl, "./output/gtex2flfit.rds")

fl.sampler <- flash_lf_sampler(fldata.strong, fl, ebnm_fn=ebnm_ash,
                               fixed="loadings")
fl.samp <- fl.sampler(200)
fl.lfsr <- flash_lfsr(fl.samp)
saveRDS(fl.lfsr, "./output/gtex2lfsr.rds")

fl.pm <- flash_get_lf(fl)
m.pm <- t(get_pm(m2))
png("./output/gtex2compare.png")
plot(as.vector(fl.pm), as.vector(m.pm), xlab="FLASH PM", ylab="MASH PM",
     main="Posterior means on GTEx data", pch='.')
abline(0, 1, lty=2)
dev.off()
cor(as.vector(fl.pm), as.vector(m.pm)) # 0.98

# look at 5007: unique in thyroid; mash shrinks others to zero, flash doesn't

confusion_matrix <- function(t) {
  mash_signif <- m.lfsr <= t
  flash_signif <- fl.lfsr <= t
  round(table(mash_signif, flash_signif)
        / length(mash_signif), digits=3)
}
confusion_matrix(.05)
confusion_matrix(.01)
confusion_matrix(.001)
