# Functions for running simulations and combining results ---------------

run_sims <- function(sim_fn, nsims, plot_title, fpath, Kmax=50,
                     backfit=FALSE) {
  if (nsims == 1) {
    res = run_one_sim(sim_fn, Kmax, backfit)
  } else {
    res = run_many_sims(sim_fn, nsims, Kmax, backfit)
  }

  saveRDS(output_res_mat(res$res), paste0(fpath, "res.rds"))
  if (!(plot_title == "Null simulation")) {
    png(paste0(fpath, "ROC.png"))
    plot_ROC(res$res, plot_title)
    dev.off()
  }
  png(paste0(fpath, "time.png"))
  plot_timing(res$timing)
  dev.off()

  return(res)
}

run_one_sim <- function(sim_fn, Kmax, nsamp=200, backfit) {
  data <- do.call(sim_fn, list())

  # If there are no strong signals, trying to run ED throws an error, so
  #   we need to do some error handling to fit the MASH object
  mfit <- fit_mash(data$Y)

  flfits <- list()
  flfits$Zero <- fit_flash_zero(data$Y, Kmax, backfit=backfit)
  flfits$OHL <- fit_flash_OHL(data$Y, Kmax, backfit=backfit)
  flfits$OHF <- fit_flash_OHF(data$Y, Kmax, backfit=backfit)

  message("Running MASH diagnostics")
  mres <- mash_diagnostics(mfit$m, data$true_Y)

  message("Running FLASH diagnostics")
  methods <- names(flfits)
  flres <- list()
  for (method in methods) {
    flres[[method]] <- flash_diagnostics(flfits[[method]]$f, data$Y,
                                         data$true_Y, nsamp)
  }

  # Combine results:
  timing = lapply(flfits,
                  function(method) {
                    list(ed.or.greedy = method$t.greedy,
                         mash.or.backfit = method$t.backfit)
                  })
  timing$MASH <- list(ed.or.greedy = mfit$t.ed,
                      mash.or.backfit = mfit$t.mash)

  all_res <- flres
  all_res$MASH <- mres

  list(timing = timing, res = all_res)
}

run_many_sims <- function(sim_fn, nsims, Kmax, backfit) {
  res <- list()
  combined_res <- list()

  for (i in 1:nsims) {
    message(paste("  Simulating dataset", i))
    res[[i]] <- run_one_sim(sim_fn, Kmax, backfit=backfit)
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
  return(combined_res)
}


# Plotting and output functions -----------------------------------------

plot_timing <- function(timing, units="secs") {
  data <- sapply(timing, unlist)
  data <- data[, c(4, 1, 2, 3)]
  barplot(colSums(data), axes=T,
          main=paste("Average time to fit in", units),
          names.arg = colnames(data),
          #legend.text = c("ED/Greedy", "MASH/Backfit"),
          ylim = c(0, max(colSums(data))*2),
          cex.names = 0.8)
  # (increasing ylim is easiest way to deal with legend getting in way)
}

plot_ROC <- function(res, main="ROC curve") {
  # Number of nulls and nonnulls are identical across methods
  n_nonnulls <- res[[1]]$n_nonnulls
  n_nulls <- res[[1]]$n_nulls

  m_y <- res$MASH$TP / n_nonnulls
  m_x <- res$MASH$FP / n_nulls
  plot(m_x, m_y, xlim=c(0, 1), ylim=c(0, 1), type='l',
       xlab='FPR', ylab='TPR', main=main, col="orange")

  idx <- 1:3
  colors <- c("skyblue", "seagreen", "yellow2")
  for (i in idx) {
    y <- res[[i]]$TP / n_nonnulls
    x <- res[[i]]$FP / n_nulls
    lines(x, y, col=colors[i])
  }
  legend("bottomright", c("MASH", names(res)[idx]), lty=1,
         col=c("orange", colors[idx]), cex=0.8)
}

output_res_mat <- function(res) {
  mat <- rbind(sapply(res, function(method) {method$MSE}),
               sapply(res, function(method) {method$CI}))

  rownames(mat) = c("MSE", "95% CI cov")
  return(mat)
}
