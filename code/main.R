# library(flashr)
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
library(mashr)

source("./code/sims.R")
source("./code/fits.R")
source("./code/utils.R")
source("./code/mashvflash.R")

set.seed(1)

# s is sparsity (prop null genes); var parameters define size of effects
all_sims <- sim_fns(n=25, p=1000, s=0.8,
                    mashvar=1, fvar=1, lvar=1, UVvar=0.25)

for (i in 1:length(all_sims)) {
  message(paste("  Beginning simulation #", i, sep=""))
  res <- run_sims(all_sims[[i]], nsims=10, plot_title=sim_names[[i]],
                  fpath = paste0("./output/sim", i))
}
