# Make sure to use branch "trackObj" when loading flashr.

# devtools::install_github("stephenslab/flashr", ref="trackObj")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")
library(mashr)

source("./code/sims.R")
source("./code/fits.R")
source("./code/utils.R")

source("./code/mashvflash.R")

set.seed(1)

# s is sparsity (prop null genes); var parameters define sizes of effects
all_sims <- sim_fns(n=25, p=1000, s=0.8,
                    indvar=4, shvar=2, uniqvar=100,
                    r1var=4, r5var=1)

for (i in 1:length(all_sims)) {
  message(paste("  Beginning simulation #", i, sep=""))
  res <- run_sims(all_sims[[i]], nsims=10, plot_title=sim_names[[i]],
                  fpath = paste0("./output/MASHvFLASHsims/greedy/sim", i),
                  Kmax=25, backfit=FALSE)
}

for (i in 1:length(all_sims)) {
  message(paste("  Beginning simulation #", i, sep=""))
  res <- run_sims(all_sims[[i]], nsims=10, plot_title=sim_names[[i]],
                  fpath = paste0("./output/MASHvFLASHsims/backfit/sim", i),
                  Kmax=25, backfit=TRUE)
}
