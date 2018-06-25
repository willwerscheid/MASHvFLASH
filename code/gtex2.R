devtools::load_all("/Users/willwerscheid/GitHub/flashr2/")
library(mashr)
library(corrplot)
source("./code/fits.R")
source("./code/utils.R")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
strong <- gtex$strong.z
random <- gtex$random.z

# MASH ------------------------------------------------------------------
mdata.strong <- mash_set_data(strong, Shat = 1)
mdata.random <- mash_set_data(random, Shat = 1)

# 1. Learn "data-driven" loadings using "strong" tests
U.pca <- cov_pca(mdata.strong, 5)
U.ed <- cov_ed(mdata.strong, U.pca)

# 2. Fit the model to the random tests to learn mixture weights
U.c <- cov_canonical(mdata.random)
m <- mash(mdata.random, Ulist = c(U.ed,U.c))
m.pi <- get_estimated_pi(m)
m.pi[m.pi > .00005] * 20000 # not sparse!
corrplot(m$fitted_g$Ulist[["ED_tPCA"]])
corrplot(m$fitted_g$Ulist[["ED_PCA_2"]])
corrplot(m$fitted_g$Ulist[["simple_het_2"]])
random.lfsr <- get_lfsr(m)
sum (random.lfsr > .05) / length(random.lfsr) # 0.76

# 3. Compute posterior summaries on the strong tests
m2 <- mash(mdata.strong, g=get_fitted_g(m), fixg=TRUE)
saveRDS(m2, "./output/gtex2mfit.rds")
m.lfsr <- get_lfsr(m2)
sum(strong.lfsr > .05) / length(strong.lfsr) # 0.44
get_loglik(m2) -1353494


# FLASH -----------------------------------------------------------------
fldata.strong <- flash_set_data(t(strong), S = 1)
fldata.random <- flash_set_data(t(random), S = 1)
# 1. Learn "data-driven" loadings using "strong" tests (use ebnm_pn
#   for speed).
# No backfit...
fl.strong <- flash(fldata.strong, Kmax=50, var_type="zero", verbose=T)
saveRDS(fl.strong, "./output/gtexstrongfit.rds")
flash_get_pve(fl.strong) * 100 # 73% on first factor/loading

# 2. Fit the model to the random tests to learn the priors on the factors.
LL = cbind(flash_get_l(fl.strong), diag(rep(1, 44)))
fl.random <- flash_add_fixed_l(fldata.random, LL)
fl.random <- flash_backfit(fldata.random, fl.random, var_type="zero",
                           ebnm_fn = ebnm_ash, nullcheck=F, verbose=T)
saveRDS(fl.random, "./output/gtexrandomfit.rds")
flash_get_pve(fl.random) # now only 21%
flash_plot_pve(fl.random)

# 3. Compute posterior summaries on the strong tests, using g_f from step 2.
fl <- flash_add_fixed_l(fldata.strong, LL)
fl <- flash_backfit(fldata.strong, fl, var_type="zero", ebnm_fn = ebnm_ash,
                    gf=flash_get_gf(fl.random), fixgf=T, nullcheck=F, verbose=T)
# likelihood: -1361278
saveRDS(fl, "./output/gtex2flfit.rds")

fl.sampler <- flash_lf_sampler(fldata.strong, fl, ebnm_fn=ebnm_ash,
                               fixed="loadings")
fl.samp <- fl.sampler(200)
fl.lfsr <- flash_lfsr(fl.samp)
saveRDS(fl.lfsr, "./output/gtex2lfsr.rds")
sum(fl.lfsr > .05) / length(fl.lfsr) # 0.51

fl.pm <- flash_get_lf(fl)
m.pm <- t(get_pm(m2))
png("./output/gtex2compare.png")
plot(as.vector(fl.pm), as.vector(m.pm), xlab="FLASH PM", ylab="MASH PM",
     main="Posterior means on GTEx data", pch='.')
abline(0, 1, lty=2)
dev.off()
cor(as.vector(fl.pm), as.vector(m.pm)) # 0.98

fl.signif <- fl.lfsr <= .05
m.signif <- m.lfsr <= .05
flnotmash <- fl.signif & !t(m.signif)
which(colSums(flnotmash) >= 24) # 6129, 6897, 13684
mnotflash <- !fl.signif & t(m.signif)
which(colSums(mnotflash) >= 40) #2969, 4533, 10011, 11188

plot_comparison <- function(n) {
  clr = rep("peachpuff", 44)
  clr[fl.lfsr[, n] <= .05] <- "red4"
  clr[fl.lfsr[, n] <= .01] <- "red1"
  clr[m.lfsr[n, ] <= .05] <- "blue4"
  clr[m.lfsr[n, ] <= .01] <- "blue1"
  clr[fl.lfsr[, n] <= .05 & m.lfsr[n, ] <= .05] <- "purple"
  plot(strong[n, ], pch=19, col=clr)
  #plot((1 - pnorm(abs(strong[n, ])))*sign(strong[n, ]),
       #ylab="", pch=19, ylim=c(-1, 1))
  #points(fl.lfsr[, n], pch=2, ylab="")
  #points(m.lfsr[n, ], pch=3)
  #abline(.05, 0)
  #abline(-.05, 0)
}
plot_comparison(6129)
plot_comparison(6897)
plot_comparison(13684)
plot_comparison(2969)
plot_comparison(4533)
plot_comparison(10011)
plot_comparison(11188)

confusion_matrix <- function(t) {
  mash_signif <- m.lfsr <= t
  flash_signif <- fl.lfsr <= t
  round(table(mash_signif, flash_signif)
        / length(mash_signif), digits=3)
}
confusion_matrix(.05)
confusion_matrix(.01)
confusion_matrix(.001)
