## @knitr gtexanalysis

missing.tissues <- c(7, 8, 19, 20, 24, 25, 31, 34, 37)
gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", sep = '\t', comment.char = '')[-missing.tissues, 2]
OHF.colors <- c("tan4", "tan3")
zero.colors <- c("black", gray.colors(19, 0.2, 0.9),
                 gray.colors(17, 0.95, 1))

plot_test <- function(n, lfsr, pm, method_name) {
  plot(strong[, n], pch=1, col="black", xlab="", ylab="", cex=0.6,
       ylim=c(min(c(strong[, n], 0)), max(c(strong[, n], 0))),
       main=paste0("Test #", n, "; ", method_name))
  size = rep(0.6, 44)
  shape = rep(15, 44)
  signif <- lfsr[, n] <= .05
  shape[signif] <- 17
  size[signif] <- 1.35 - 15 * lfsr[signif, n]
  size <- pmin(size, 1.2)
  points(pm[, n], pch=shape, col=as.character(gtex.colors), cex=size)
  abline(0, 0)
}

plot_ohf_v_ohl_loadings <- function(n, ohf_fit, ohl_fit, ohl_name,
                                    legend_pos = "bottomright") {
  ohf <- abs(ohf_fit$EF[n, ] * apply(abs(ohf_fit$EL), 2, max))
  ohl <- -abs(ohl_fit$EF[n, ] * apply(abs(ohl_fit$EL), 2, max))
  data <- rbind(c(ohf, rep(0, length(ohl) - 45)),
                c(ohl[1:45], rep(0, length(ohf) - 45),
                  ohl[46:length(ohl)]))
  colors <- c("black",
              as.character(gtex.colors),
              OHF.colors,
              zero.colors[1:(length(ohl) - 45)])
  x <- barplot(data, beside=T, col=rep(colors, each=2),
               main=paste0("Test #", n, " loadings"),
               legend.text = c("OHF", ohl_name),
               args.legend = list(x = legend_pos, bty = "n", pch="+-",
                                  fill=NULL, border="white"))
  text(x[2*(46:47) - 1], min(data) / 10,
       labels=as.character(1:2), cex=0.4)
  text(x[2*(48:ncol(data))], max(data) / 10,
       labels=as.character(1:(length(ohl) - 45)), cex=0.4)
}

plot_ohl_v_zero_loadings <- function(n, ohl_fit, zero_fit, ohl_name,
                                    legend_pos = "topright") {
  ohl <- abs(ohl_fit$EF[n, ] * apply(abs(ohl_fit$EL), 2, max))
  # Combine equal effects and first data-driven loading
  ohl[1] <- ohl[1] + ohl[46]
  ohl <- ohl[-46]
  zero <- -abs(zero_fit$EF[n, ] * apply(abs(zero_fit$EL), 2, max))
  data <- rbind(c(ohl, rep(0, length(zero) - length(ohl) + 44)),
                c(zero[1], rep(0, 44), zero[2:length(zero)]))
  colors <- c("black", as.character(gtex.colors), zero.colors)
  x <- barplot(data, beside=T, col=rep(colors, each=2),
               main=paste0("Test #", n, " loadings"),
               legend.text = c(ohl_name, "Zero"),
               args.legend = list(x = legend_pos, bty = "n", pch="+-",
                                  fill=NULL, border="white"))
  text(x[2*(seq(46, ncol(data), by=2)) - 1], min(data) / 10,
       labels=as.character(seq(2, length(zero), by=2)), cex=0.4)
}

compare_methods <- function(lfsr1, lfsr2, pm1, pm2) {
  res <- list()
  res$first_not_second <- find_A_not_B(lfsr1, lfsr2)
  res$lg_first_not_second <- find_large_A_not_B(lfsr1, lfsr2)
  res$second_not_first <- find_A_not_B(lfsr2, lfsr1)
  res$lg_second_not_first <- find_large_A_not_B(lfsr2, lfsr1)
  res$diff_pms <- find_overall_pm_diff(pm1, pm2)
  return(res)
}

# Find tests where many conditions are significant according to
#   method A but not according to method B.
find_A_not_B <- function(lfsrA, lfsrB) {
  select_tests(colSums(lfsrA <= 0.05 & lfsrB > 0.05))
}

# Find tests where many conditions are highly significant according to
#   method A but are not significant according to method B.
find_large_A_not_B <- function(lfsrA, lfsrB) {
  select_tests(colSums(lfsrA <= 0.01 & lfsrB > 0.05))
}

find_overall_pm_diff <- function(pmA, pmB, n = 4) {
  pm_diff <- colSums((pmA - pmB)^2)
  return(order(pm_diff, decreasing = TRUE)[1:4])
}

# Get at least four (or min_n) "top" tests.
select_tests <- function(colsums, min_n = 4) {
  n <- 45
  n_tests <- 0
  while (n_tests < min_n && n > 0) {
    n <- n - 1
    n_tests <- sum(colsums >= n)
  }
  return(which(colsums >= n))
}
