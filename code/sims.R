## SIMULATION FUNCTIONS -------------------------------------------------

# n is number of conditions, p is number of genes

# Noise is i.i.d. N(0, 1)
get_E <- function(n, p, sd = 1) {
  matrix(rnorm(n * p, 0, sd), n, p)
}

# Simulate from null model ----------------------------------------------

null_sim <- function(n, p, seed = NULL) {
  set.seed(seed)
  get_E(n, p)
}

# Simulate from MASH model ----------------------------------------------

# Sigma is list of covariance matrices
# pi[j] is probability that effect j has covariance Sigma[[j]]
# s is sparsity (percentage of null effects)
mash_sim <- function(n, p, Sigma, pi = NULL, s = 0.8, seed = NULL) {
  set.seed(NULL)
  if (is.null(pi)) {
    pi = rep(1, length(Sigma)) # default to uniform distribution
  }
  assertthat::are_equal(length(pi), length(Sigma))
  for (j in length(Sigma)) {
    assertthat::are_equal(dim(Sigma[j]), c(n, n))
  }

  pi <- pi / sum(pi) # normalize pi to sum to one
  which_sigma <- sample(1:length(pi), p, replace=TRUE, prob=pi)
  nonnull_fx <- sample(1:p, floor((1 - s)*p), replace=FALSE)

  X <- matrix(0, n, p)
  for (j in nonnull_fx) {
    X[, j] <- MASS::mvrnorm(1, rep(0, n), Sigma[[which_sigma[j]]])
  }
  Y <- X + get_E(n, p)
  list(Y = Y, nonnull = nonnull_fx)
}


# Simulate from FLASH model ---------------------------------------------

# fs is sparsity of factors (percentage of null effects)
# fvar is variance of effects (generated from normal distribution)
# ls is sparsity of loadings
# lvar is variance of loadings
# UVvar is variance of dense rank-one matrix included to mimic something
#   like unwanted variation (set it to 0 to ignore it)
flash_sim <- function(n, p, k, fs, fvar, ls, lvar, UVvar = 0, seed = NULL) {
  set.seed(seed)

  nonnull_ll <- matrix(sample(c(0, 1), n*k, TRUE, c(ls, 1 - ls)), n, k)
  LL <- nonnull_ll * matrix(rnorm(n*k, 0, sqrt(lvar)), n, k)

  nonnull_ff <- matrix(sample(c(0, 1), k*p, TRUE, c(fs, 1 - fs)), k, p)
  FF <- nonnull_ff * matrix(rnorm(k*p, 0, sqrt(fvar)), k, p)

  Y <- LL %*% FF + get_E(n, p)
  # add unwanted variation
  Y <- Y + outer(rnorm(n, 0, sqrt(UVvar)), rnorm(p, 0, sqrt(UVvar)))
  list(Y = Y, LL = LL, FF = FF)
}


## SIMULATIONS ----------------------------------------------------------

# I generate a total of seven datasets. One is null; three are from the
# MASH model; three are from the FLASH model.

n <- 25
p <- 1000
# 1. Everything is null
null_Y <- null_sim(n, p)

mashvar <- 100
Sigma <- list()
Sigma[[1]] <- diag(rep(mashvar, n))
# 2. Effects are independent across conditions
ind <- mash_sim(n, p, Sigma)

Sigma[[2]] <- matrix(mashvar, n, n)
# 3. Effects are either independent or shared
indsh <- mash_sim(n, p, Sigma)

for (j in 1:n) {
  Sigma[[2 + j]] <- matrix(0, n, n)
  Sigma[[2 + j]][j, j] <- mashvar
}
pi <- c(n, n, rep(1, n))
# 4. Effects are independent, shared, or unique to a single condition
mash <- mash_sim(n, p, Sigma)

# 5. Rank one model
fvar = 100
lvar = 100
rank1 <- flash_sim(n, p, 1, 0.8, fvar, 0.5, lvar)

# 6. Rank 5 model
rank5 <- flash_sim(n, p, 5, 0.8, fvar, 0.2, lvar)

UVvar = 25
# 7. Rank 3 model with unwanted variation
UV <- flash_sim(n, p, 3, 0.8, fvar, 0.3, lvar, UVvar)
