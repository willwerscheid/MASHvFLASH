## @knitr sims2

devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
library(mashr)

source("./code/fits.R")
source("./code/sims.R")
source("./code/utils.R")

set.seed(1)

n <- 20
p <- 1200
nsims <- 10
nsamp <- 200 # for sampling lfsr (FLASH fits)
ncol <- 25 # number of columns that exhibit each variance type
t <- 0.05 # "significance" threshold

Sigma <- list()

# Independent (small)
Sigma[[1]] <- diag(2^2, n)
# Independent (large)
Sigma[[2]] <- diag(5^2, n)
# Independent (different sizes)
sizes <- seq(1, 5, length.out=n)
Sigma[[3]] <- diag(sizes^2)
# Identical (small)
Sigma[[4]] <- matrix(2^2, nrow=n, ncol=n)
# Identical (large)
Sigma[[5]] <- matrix(5^2, nrow=n, ncol=n)
# Rank-one
Sigma[[6]] <- outer(sizes, sizes)

zeros <- matrix(0, nrow=n, ncol=n)
for (j in 7:11) {
  Sigma[[j]] <- zeros
}
# Unique (small)
uniqsmidx <- 1
Sigma[[7]][uniqsmidx, uniqsmidx] <- 3^2
# Unique (large)
uniqlgidx <- 2
Sigma[[8]][uniqlgidx, uniqlgidx] <- 8^2
# Shared (3 conditions)
shar3idx <- 3:5
Sigma[[9]][shar3idx, shar3idx] <- matrix(3^2, nrow=3, ncol=3)
# Shared (10 conditions)
shar10idx <- 1:10
Sigma[[10]][shar10idx, shar10idx] <- matrix(2^2, nrow=10, ncol=10)
# Shared (5 conditions, different sizes)
shar5idx <- 6:10
sizes <- seq(2, 4, length.out=5)
Sigma[[11]][shar5idx, shar5idx] <- outer(sizes, sizes)

# Rank-5
A <- matrix(rnorm(n*5, 0, 2), nrow=5, ncol=n)
Sigma[[12]] <- t(A) %*% A
# Rank-10
A <- matrix(rnorm(n*10, 0, 2), nrow=10, ncol=n)
Sigma[[13]] <- t(A) %*% A
# Random
A <- matrix(rnorm(n*n, 0, 2), nrow=n, ncol=n)
Sigma[[14]] <- t(A) %*% A

# Large independent plus small identical
Sigma[[15]] <- Sigma[[2]] + Sigma[[4]]
# Small independent plus large unique
Sigma[[16]] <- Sigma[[1]] + Sigma[[8]]
# Small identical plus large unique
Sigma[[17]] <- Sigma[[4]] + Sigma[[8]]
ntypes <- 17

partnames <- c("IndSm", "IndLg", "IndDiff",
               "IdentSm", "IdentLg", "Rank1",
               "UniqSmNull", "UniqSmNonnull",
               "UniqLgNull", "UniqLgNonnull",
               "Shar3Null", "Shar3Nonnull",
               "Shar10Null", "Shar10Nonnull",
               "Shar5Null", "Shar5Nonnull",
               "Rank5", "Rank10", "Random",
               "IndIdent", "IndUniq", "IdentUniq",
               "Null")
partxidx <- list(1:n, 1:n, 1:n, 1:n, 1:n, 1:n,
                 setdiff(1:n, uniqsmidx), uniqsmidx,
                 setdiff(1:n, uniqlgidx), uniqlgidx,
                 setdiff(1:n, shar3idx), shar3idx,
                 setdiff(1:n, shar10idx), shar10idx,
                 setdiff(1:n, shar5idx), shar5idx,
                 1:n, 1:n, 1:n, 1:n, 1:n, 1:n, 1:n)
partyidx <- list(1:ncol, ncol + 1:ncol, 2*ncol + 1:ncol,
                 3*ncol + 1:ncol, 4*ncol + 1:ncol, 5*ncol + 1:ncol,
                 6*ncol + 1:ncol, 6*ncol + 1:ncol,
                 7*ncol + 1:ncol, 7*ncol + 1:ncol,
                 8*ncol + 1:ncol, 8*ncol + 1:ncol,
                 9*ncol + 1:ncol, 9*ncol + 1:ncol,
                 10*ncol + 1:ncol, 10*ncol + 1:ncol,
                 11*ncol + 1:ncol, 12*ncol + 1:ncol,
                 13*ncol + 1:ncol, 14*ncol + 1:ncol,
                 15*ncol + 1:ncol, 16*ncol + 1:ncol,
                 (17*ncol + 1):p)
nparts <- length(partnames)
partition_by_type <- function(X) {
  ret <- rep(0, nparts)
  for (i in 1:nparts) {
    ret[i] <- mean(X[partxidx[[i]], partyidx[[i]]])
  }
  names(ret) <- partnames
  ret
}

mses <- matrix(0, nrow=3, ncol=nparts)
ts <- c(seq(.005, .05, by=.005), seq(.06, .1, by=.01), seq(.2, 1.0, by=.1))

pr.pn <- matrix(0, nrow=length(ts), ncol=nparts)
pr.ash <- matrix(0, nrow=length(ts), ncol=nparts)
pr.m <- matrix(0, nrow=length(ts), ncol=nparts)

for (i in 1:nsims) {
  X <- matrix(0, nrow=n, ncol=p)
  for (j in 1:ntypes) {
    start_col = 1 + ncol*(j-1)
    end_col = ncol*j
    X[, start_col:end_col] <- t(MASS::mvrnorm(ncol, rep(0, n), Sigma[[j]]))
  }
  Y <- X + rnorm(n*p)

  fl.pn <- fit_flash(Y, Kmax=30, methods=3, ebnm_fn="ebnm_pn") # OHL
  fl.ash <- fit_flash(Y, Kmax=30, methods=3, ebnm_fn="ebnm_ash")
  m <- fit_mash(Y)

  base.se <- (Y - X)^2
  base.mse <- partition_by_type(base.se)

  fl.pn.se <- (flash_get_fitted_values(fl.pn$fits$OHL) - X)^2
  fl.pn.mse <- partition_by_type(fl.pn.se) / base.mse

  fl.ash.se <- (flash_get_fitted_values(fl.ash$fits$OHL) - X)^2
  fl.ash.mse <- partition_by_type(fl.ash.se) / base.mse

  m.se <- (t(get_pm(m$m)) - X)^2
  m.mse <- partition_by_type(m.se) / base.mse

  mses[1,] <- mses[1,] + fl.pn.mse
  mses[2,] <- mses[2,] + fl.ash.mse
  mses[3,] <- mses[3,] + m.mse


  fl.pn.sampler <- flash_sampler(Y, fl.pn$fits$OHL, fixed="loadings")
  fl.pn.samp <- fl.pn.sampler(nsamp)
  fl.pn.lfsr <- flash_lfsr(fl.pn.samp)

  fl.ash.sampler <- flash_sampler(Y, fl.ash$fits$OHL, fixed="loadings")
  fl.ash.samp <- fl.ash.sampler(nsamp)
  fl.ash.lfsr <- flash_lfsr(fl.ash.samp)

  m.lfsr <- t(get_lfsr(m$m))

  for (j in 1:length(ts)) {
    fl.pn.signif <- fl.pn.lfsr <= ts[j]
    fl.pn.pr <- partition_by_type(fl.pn.signif)

    fl.ash.signif <- fl.ash.lfsr <= ts[j]
    fl.ash.pr <- partition_by_type(fl.ash.signif)

    m.signif <- m.lfsr <= ts[j]
    m.pr <- partition_by_type(m.signif)

    pr.pn[j,] <- pr.pn[j,] + fl.pn.pr
    pr.ash[j,] <- pr.ash[j,] + fl.ash.pr
    pr.m[j,] <- pr.m[j,] + m.pr
  }

}

mses <- mses / nsims
pr.pn <- pr.pn / nsims
pr.ash <- pr.ash / nsims
pr.m <- pr.m / nsims

saveRDS(mses, "./output/sims2mse.rds")
saveRDS(pr.pn, "./output/sims2prpn.rds")
saveRDS(pr.ash, "./output/sims2prash.rds")
saveRDS(pr.m, "./output/sims2prm.rds")
