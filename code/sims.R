## SIMULATION FUNCTIONS -------------------------------------------------

# n is number of conditions, p is number of genes

# Noise is i.i.d. N(0, 1)
get_E <- function(n, p, sd = 1) {
  matrix(rnorm(n * p, 0, sd), n, p)
}


# Simulate from null model ----------------------------------------------

null_sim <- function(n, p, seed = NULL) {
  set.seed(seed)
  Y <- get_E(n, p)
  true_Y <- matrix(0, n, p)

  list(Y = Y, true_Y = true_Y)
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
  list(Y = Y, true_Y = X)
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
  LL <- nonnull_ll * matrix(rnorm(n*k, 0, sqrt(lvar)), nrow=n, ncol=k)

  nonnull_ff <- matrix(sample(c(0, 1), p, TRUE, c(fs, 1 - fs)),
                       nrow=k, ncol=p, byrow=TRUE)
  FF <- nonnull_ff * matrix(rnorm(k*p, 0, sqrt(fvar)), nrow=k, ncol=p)

  X <- LL %*% FF
  Y <- X + get_E(n, p)
  # add unwanted variation
  Y <- Y + outer(rnorm(n, 0, sqrt(UVvar)), rnorm(p, 0, sqrt(UVvar)))
  list(Y = Y, true_Y = X)
}


## SIMULATIONS ----------------------------------------------------------

# Functions to generate six types of datasets. One is null; three are
# from the MASH model; two are from the FLASH model.

sim_fns <- function(n, p, s,
                    indvar, shvar, uniqvar,
                    r1var, r5var) {

  # 1. Everything is null
  sim_null <- function(){ null_sim(n, p) }

  Sigma <- list()

  # 2. Effects are independent across conditions
  Sigma[[1]] <- diag(rep(indvar, n))
  Sigma_ind <- Sigma
  sim_ind <- function(){ mash_sim(n, p, Sigma_ind, s=s) }

  # 3. Effects are either independent or shared
  Sigma[[2]] <- matrix(shvar, n, n)
  Sigma_indsh <- Sigma
  sim_indsh <- function(){ mash_sim(n, p, Sigma_indsh, s=s) }

  # 4. Effects are independent, shared, or unique to a single condition
  for (j in 1:n) {
    Sigma[[2 + j]] <- matrix(0, n, n)
    Sigma[[2 + j]][j, j] <- uniqvar
  }
  pi <- c(n, n, rep(1, n))
  sim_mash <- function(){ mash_sim(n, p, Sigma, pi, s=s) }

  # 5. Rank one model
  sim_rank1 <- function(){ flash_sim(n, p, 1, s, r1var, 0.2, 1) }

  # 6. Rank 5 model
  sim_rank5 <- function(){ flash_sim(n, p, 5, s, r5var, 0.8, 1) }

  # 7. Rank 3 model with unwanted variation
  # sim_UV <- function(){ flash_sim(n, p, 3, s, r3var, 0.5, 1, UVvar) }

  c(sim_null, sim_ind, sim_indsh, sim_mash, sim_rank1, sim_rank5)
}

sim_names <- c("Null simulation",
               "All independent effects",
               "Independent and shared",
               "Independent, shared, and unique",
               "Rank 1 FLASH model",
               "Rank 5 FLASH model")
