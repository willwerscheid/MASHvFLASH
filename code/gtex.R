devtools::load_all("/Users/willwerscheid/GitHub/flashr2/")
library(mashr)

gtex <- readRDS("./data/MatrixEQTLSumStats.Portable.Z.rds")
data <- gtex$test.z
data <- t(data)
fl_data <- flash_set_data(data, S = 1)

source("./code/fits.R")
source("./code/sims.R")

gtex_mfit <- fit_mash(data)
saveRDS(gtex_mfit, "./output/gtexmfit.rds")

gtex_flfit <- fit_flash(data, Kmax = 40, add_onehots_first = FALSE)
saveRDS(gtex_flfit, "./output/gtexflfit.rds")

obj1 <- flash_get_objective(fl_data, gtex_flfit$fl) # -1259284


# Try OHF method of fitting FLASH object and compare likelihoods
gtex_flfit2 <- fit_flash(data, Kmax = 40, add_onehots_first = TRUE)
saveRDS(gtex_flfit2, "./output/gtexflfit2.rds")

obj2 <- flash_get_objective(fl_data, gtex_flfit2$fl) # -1635669

# Now do an additional backfit on the OHF fit
gtex_flfit3 <- list()
t0 <- Sys.time()
gtex_flfit3$fl <- flash_backfit(fl_data, gtex_flfit2$fl, var_type = "zero",
                             nullcheck = F, verbose = T)
t <- Sys.time() - t0
gtex_flfit3$timing <- gtex_flfit2$timing
gtex_flfit3$timing$backfit <- gtex_flfit3$timing$backfit + t
gtex_flfit3$timing$total <- gtex_flfit3$timing$total + t
saveRDS(gtex_flfit3, "./output/gtexflfit3.rds")

obj3 <- flash_get_objective(fl_data, gtex_flfit3$fl) # -1306242


# Use PM from each method as "true Y" and do diagnostics
# fl_pm <- flash_get_lf(gtex_flfit$fl)
# gtex_mres <- mash_diagnostics(gtex_mfit$m, fl_pm)
# saveRDS(gtex_mres, "./output/gtexmres.rds")
#
# m_pm <- t(get_pm(gtex_mfit$m))
# gtex_flres <- flash_diagnostics(gtex_flfit$fl, data, m_pm, nsamp = 200)
# saveRDS(gtex_flres, "./output/gtexflres.rds")


# Plot FLASH PM vs. MASH PM
fl_pm <- flash_get_lf(gtex_flfit$fl)
m_pm <- t(get_pm(gtex_mfit$m))
png("./output/gtexcompare.png")
plot(as.vector(fl_pm), as.vector(m_pm), xlab="FLASH PM", ylab="MASH PM",
     main="Posterior means on GTEx data", pch='.')
dev.off()
corr <- cor(as.vector(fl_pm), as.vector(m_pm))

# Use LFSR to get "significant" effects and get confusion matrices
m_lfsr <- t(get_lfsr(gtex_mfit$m))

fl_sampler <- flash_lf_sampler(data, gtex_flfit$fl, ebnm_fn=ebnm_pn, fixed="loadings")
fl_lfsr <- flash_lfsr(fl_sampler(200))

confusion_matrix <- function(t) {
  mash_signif <- m_lfsr <= t
  flash_signif <- fl_lfsr <= t
  round(table(mash_signif, flash_signif)
        / length(mash_signif), digits=3)
}
confusion_matrix(.05)
confusion_matrix(.01)
confusion_matrix(.001)
