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
  saveRDS(output_res_mat(res), paste0(fpath, "res.rds"))
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
  elems <- names(res[[1]])
  for (elem in elems) {
    combined_res[[elem]] <- list()
    subelems <- names(res[[1]][[elem]])
    for (subelem in subelems) {
      combined_res[[elem]][[subelem]] <- list()
      subsubelems <- names(res[[1]][[elem]][[subelem]])
      for (subsubelem in subsubelems) {
        tmp <- lapply(res, function(x) {x[[elem]][[subelem]][[subsubelem]]})
        combined_res[[elem]][[subelem]][[subsubelem]] <- Reduce(`+`, tmp) / nsims
      }
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

  flfits <- fit_flash(data$Y, Kmax)

  message("Running MASH diagnostics")
  mres <- mash_diagnostics(mfit$m, data$true_Y)

  message("Running FLASH diagnostics")
  flres <- list()
  methods <- names(flfits$fits)
  for (method in methods) {
    flres[[method]] <- flash_diagnostics(flfits$fits[[method]], data$Y,
                                         data$true_Y, nsamp)
  }

  # Give mash and flash the same structure so we can combine multiple sims
  mash_timing = list()
  mash_timing$mash <- mfit$timing
  mash_res = list()
  mash_res$mash <- mres
  list(mash_timing = mash_timing, mash_res = mash_res,
       flash_timing = flfits$timing, flash_res = flres)
}

output_res_mat <- function(res) {
  mat <- c(res$mash_res$mash$MSE, res$mash_res$mash$CI)
  methods <- names(res$flash_res)
  for (method in methods) {
    mat <- cbind(mat, c(res$flash_res[[method]]$MSE,
                        res$flash_res[[method]]$CI))
  }
  rownames(mat) = c("MSE", "95% CI cov")
  colnames(mat) = c("MASH", methods)
  mat
}

plot_timing <- function(res, units="secs") {
  data <- as.numeric(c(res$mash_timing$mash$ed,
                       res$mash_timing$mash$mash),
                     units=units)
  methods <- names(res$flash_timing)
  for (method in methods) {
    data <- cbind(data,
                  as.numeric(c(res$flash_timing[[method]]$greedy,
                               res$flash_timing[[method]]$backfit),
                             units=units))
  }
  barplot(data, axes=T,
          main=paste("Average time to fit in", units),
          names.arg = c("MASH", methods),
          legend.text = c("ED/Greedy", "MASH/Backfit"),
          ylim = c(0, max(colSums(data))*2))
  # (increasing ylim is easiest way to deal with legend getting in way)
}

plot_ROC <- function(res, main="ROC curve",
                     colors=c("gold", "pink", "red", "green", "blue")) {
  # Number of nulls and nonnulls are identical across methods
  n_nonnulls <- res$mash_res$mash$n_nonnulls
  n_nulls <- res$mash_res$mash$n_nulls

  m_y <- res$mash_res$mash$TP / n_nonnulls
  m_x <- res$mash_res$mash$FP / n_nulls
  plot(m_x, m_y, xlim=c(0, 1), ylim=c(0, 1), type='l',
       xlab='FPR', ylab='TPR', main=main)

  methods <- names(res$flash_res)
  idx <- 0
  for (method in methods) {
    idx <- idx + 1
    y <- res$flash_res[[method]]$TP / n_nonnulls
    x <- res$flash_res[[method]]$FP / n_nulls
    lines(x, y, col=colors[idx])
  }
  legend("bottomright", c("MASH", methods), lty=1,
         col=c("black", colors[1:idx]))
}

