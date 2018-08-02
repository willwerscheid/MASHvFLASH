## @knitr gtexanalysis

missing.tissues <- c(7, 8, 19, 20, 24, 25, 31, 34, 37)
gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", sep = '\t', comment.char = '')[-missing.tissues, 2]

plot_test <- function(n, lfsr, pm, method_name) {
  plot(strong[, n], pch=1, col="black", xlab="", ylab="", cex=0.6,
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

compare_methods <- function(lfsr1, lfsr2) {
  res <- list()
  res$first_not_second <- find_A_not_B(lfsr1, lfsr2)
  res$lg_first_not_second <- find_large_A_not_B(lfsr1, lfsr2)
  res$second_not_first <- find_A_not_B(lfsr2, lfsr1)
  res$lg_second_not_first <- find_large_A_not_B(lfsr2, lfsr1)
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
