## @knitr gtex

devtools::load_all("/Users/willwerscheid/GitHub/flashr2/")
library(mashr)

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

data <- gtex$random.z
data <- t(data)
fl_data <- flash_set_data(data, S = 1)

source("./code/fits.R")
source("./code/sims.R")
source("./code/utils.R")

gtex_mfit <- fit_mash(data)
saveRDS(gtex_mfit, "./output/gtexmfit.rds")

gtex_flfit <- fit_flash(data, Kmax = 40, methods=2:5)
saveRDS(gtex_flfit, "./output/gtexflfit.rds")

flash_get_objective(fl_data, gtex_flfit$fits$Zero) # -1277881
flash_get_objective(fl_data, gtex_flfit$fits$OHL) # -1277145
flash_get_objective(fl_data, gtex_flfit$fits$OHF) # -1315285
flash_get_objective(fl_data, gtex_flfit$fits$OHFp) # -1278991


# Use PM from each method as "true Y" and do diagnostics
# fl_pm <- flash_get_lf(gtex_flfit$fl)
# gtex_mres <- mash_diagnostics(gtex_mfit$m, fl_pm)
# saveRDS(gtex_mres, "./output/gtexmres.rds")
#
# m_pm <- t(get_pm(gtex_mfit$m))
# gtex_flres <- flash_diagnostics(gtex_flfit$fl, data, m_pm, nsamp = 200)
# saveRDS(gtex_flres, "./output/gtexflres.rds")


# Plot FLASH PM vs. MASH PM
fl_pm <- flash_get_lf(gtex_flfit$fits$OHL)
m_pm <- t(get_pm(gtex_mfit$m))
png("./output/gtexcompare.png")
plot(as.vector(fl_pm), as.vector(m_pm), xlab="FLASH PM", ylab="MASH PM",
     main="Posterior means on GTEx data", pch='.')
abline(0, 1, lty=2)
dev.off()
cor(as.vector(fl_pm), as.vector(m_pm)) # 0.952

# Use LFSR to get "significant" effects and get confusion matrices
m_lfsr <- t(get_lfsr(gtex_mfit$m))

fl_sampler <- flash_lf_sampler(data, gtex_flfit$fits$OHL, ebnm_fn=ebnm_pn, fixed="loadings")
fl_lfsr <- flash_lfsr(fl_sampler(200))
saveRDS(fl_lfsr, "./output/gtexfllfsr.rds")

confusion_matrix <- function(t) {
  mash_signif <- m_lfsr <= t
  flash_signif <- fl_lfsr <= t
  round(table(mash_signif, flash_signif)
        / length(mash_signif), digits=3)
}
confusion_matrix(.05)
confusion_matrix(.01)
confusion_matrix(.001)
