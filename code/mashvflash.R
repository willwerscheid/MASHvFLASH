run_sims <- function(sim_fn, nsims, plot_title, fpath) {
  #suppressMessages(
    #suppressWarnings(
      #capture.output(
        if (nsims == 1) {
          res = run_one_sim(sim_fn)
        } else {
          res = run_many_sims(sim_fn, nsims)
        }
      #)
    #)
  #)
  saveRDS(output_res_mat(res, plot_title), paste0(fpath, "res.rds"))
  if (!(plot_title == "Null simulation")) {
    png(paste0(fpath, "ROC.png"))
    plot_ROC(res, plot_title)
    dev.off()
  }
  png(paste0(fpath, "time.png"))
  plot_timing(res)
  dev.off()
}

run_many_sims <- function(sim_fn, nsims) {
  res <- list()
  combined_res <- list()

  for (i in 1:nsims) {
    res[[i]] <- run_one_sim(sim_fn)
  }
  list_elem <- names(res[[1]])
  for (elem in list_elem) {
    combined_res[[elem]] <- list()
    sub_elems <- names(res[[1]][[elem]])
    for (sub_elem in sub_elems) {
      tmp <- lapply(res, function(x) {x[[elem]][[sub_elem]]})
      combined_res[[elem]][[sub_elem]] <- Reduce(`+`, tmp)
      combined_res[[elem]][[sub_elem]] <- combined_res[[elem]][[sub_elem]] / nsims
    }
  }
  combined_res
}

run_one_sim <- function(sim_fn, Kmax = 10, nsamp=200) {
  data <- do.call(sim_fn, list())

  # If there are no strong signals, trying to run ED throws an error, so
  #   we need to do some error handling to fit the MASH object
  try(mfit <- fit_mash(data$Y))
  if (!exists("mfit")) {
    mfit <- fit_mash(data$Y, ed=F)
    mfit$timing$ed = as.difftime(0, units="secs")
  }

  flfit <- fit_flash(data$Y, Kmax, method = "FLASH")
  flfit1 <- fit_flash(data$Y, Kmax, method = "OHL")
  flfit2 <- fit_flash(data$Y, Kmax, method = "OHF")
  flfit3 <- fit_flash(data$Y, Kmax, method = "OHFplus")

  message("Running MASH diagnostics")
  mres <- mash_diagnostics(mfit$m, data$true_Y)
  message("Running FLASH diagnostics")
  flres <- flash_diagnostics(flfit$fl, data$Y, data$true_Y, nsamp)
  flres1 <- flash_diagnostics(flfit1$fl, data$Y, data$true_Y, nsamp)
  flres2 <- flash_diagnostics(flfit2$fl, data$Y, data$true_Y, nsamp)
  flres3 <- flash_diagnostics(flfit3$fl, data$Y, data$true_Y, nsamp)

  list(mash_timing = mfit$timing, mash_res = mres,
       flash_timing = flfit$timing, flash_res = flres,
       flash_OHL_timing = flfit1$timing, flash_OHL_res = flres1,
       flash_OHF_timing = flfit2$timing, flash_OHF_res = flres2,
       flash_OHFplus_timing = flfit3$timing, flash_OHFplus_res = flres3)
}

output_res_mat <- function(res, caption) {
  data.frame(MASH = c(res$mash_res$MSE, res$mash_res$CI),
             FLASH_vanilla = c(res$flash_res$MSE, res$flash_res$CI),
             FLASH_OHL = c(res$flash_OHL_res$MSE, res$flash_OHL_res$CI),
             FLASH_OHF = c(res$flash_OHF_res$MSE, res$flash_OHF_res$CI),
             FLASH_OHFp = c(res$flash_OHFplus_res$MSE,
                               res$flash_OHFplus_res$CI),
             row.names = c("MSE", "95% CI cov"))
}

plot_timing <- function(res) {
  data <- c(res$mash_timing$ed, res$mash_timing$mash,
            res$flash_timing$greedy, res$flash_timing$backfit,
            res$flash_OHL_timing$greedy, res$flash_OHL_timing$backfit,
            res$flash_OHF_timing$greedy, res$flash_OHF_timing$backfit,
            res$flash_OHFplus_timing$greedy,
            res$flash_OHFplus_timing$backfit)
  time_units <- units(data)
  data <- matrix(as.numeric(data), 2, 5)
  barplot(data, axes=T,
          main=paste("Average time to fit in", time_units),
          names.arg = c("MASH", "FLASH", "FL-OHL", "FL-OHF", "FL-OHF+"),
          legend.text = c("ED/Greedy", "MASH/Backfit"),
          ylim = c(0, max(colSums(data))*2))
  # (increasing ylim is easiest way to deal with legend getting in way)
}

plot_ROC <- function(res, main="ROC curve") {
  m_y <- res$mash_res$TP / res$mash_res$n_nonnulls
  m_x <- res$mash_res$FP / res$mash_res$n_nulls
  fl_y <- res$flash_res$TP / res$flash_res$n_nonnulls
  fl_x <- res$flash_res$FP / res$flash_res$n_nulls
  ohl_y <- res$flash_OHL_res$TP / res$flash_OHL_res$n_nonnulls
  ohl_x <- res$flash_OHL_res$FP / res$flash_OHL_res$n_nulls
  ohf_y <- res$flash_OHF_res$TP / res$flash_OHF_res$n_nonnulls
  ohf_x <- res$flash_OHF_res$FP / res$flash_OHF_res$n_nulls
  ohfp_y <- res$flash_OHFplus_res$TP / res$flash_OHFplus_res$n_nonnulls
  ohfp_x <- res$flash_OHFplus_res$FP / res$flash_OHFplus_res$n_nulls
  plot(m_x, m_y, xlim=c(0, 1), ylim=c(0, 1), type='l',
       xlab='FPR', ylab='TPR', main=main)
  lines(fl_x, fl_y, lty=2)
  lines(ohl_x, ohl_y, lty=3)
  lines(ohf_x, ohf_y, lty=4)
  lines(ohfp_x, ohfp_y, lty=5)
  legend("bottomright", c("MASH", "FLASH", "FLASH-OHL",
                          "FLASH-OHF", "FLASH-OHF+"), lty=1:5)
}

