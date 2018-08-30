devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

library(mashr)

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

strong <- gtex$strong.z
random <- gtex$random.z

# Step 1. Estimate correlation structure using MASH.

m_random <- mash_set_data(random, Shat = 1)
Vhat <- estimate_null_correlation(m_random)

# Step 2. Estimate data-driven loadings using FLASH.

# Step 2a. Fit Vhat.
n <- nrow(Vhat)
lambda.min <- min(eigen(Vhat, symmetric=TRUE, only.values=TRUE)$values)

data <- flash_set_data(Y, S = sqrt(lambda.min))

W.eigen <- eigen(Vhat - diag(rep(lambda.min, n)), symmetric=TRUE)
# The rank of W is at most n - 1, so we can drop the last eigenval/vec:
W.eigen$values <- W.eigen$values[-n]
W.eigen$vectors <- W.eigen$vectors[, -n, drop=FALSE]

fl <- flash_add_fixed_loadings(data,
                               LL=W.eigen$vectors,
                               init_fn="udv_svd",
                               backfit=FALSE)

ebnm_param_f <- lapply(as.list(W.eigen$values),
                       function(eigenval) {
                         list(g = list(a=1/eigenval, pi0=0), fixg = TRUE)
                       })
ebnm_param_l <- lapply(vector("list", n - 1),
                       function(k) {list()})
fl <- flash_backfit(data,
                    fl,
                    var_type="zero",
                    ebnm_fn="ebnm_pn",
                    ebnm_param=(list(f = ebnm_param_f, l = ebnm_param_l)),
                    nullcheck=FALSE)

# Step 2b. Add nonnegative factors.
ebnm_fn = list(f = "ebnm_pn", l = "ebnm_ash")
ebnm_param = list(f = list(warmstart = TRUE),
                  l = list(mixcompdist="+uniform"))
fl <- flash_add_greedy(data,
                       Kmax=50,
                       f_init=fl,
                       var_type="zero",
                       init_fn="udv_svd",
                       ebnm_fn=ebnm_fn,
                       ebnm_param=ebnm_param)
saveRDS(fl, "./output/MASHvFLASHVhat/2bGreedy.rds")

# Step 2c (optional). Backfit factors from step 2b.
fl <- flash_backfit(data,
                    fl,
                    kset=n:fl$nfactors,
                    var_type="zero",
                    ebnm_fn=ebnm_fn,
                    ebnm_param=ebnm_param,
                    nullcheck=FALSE)
saveRDS(fl, "./output/MASHvFLASHVhat/2cBackfit.rds")

# Step 2d (optional). Repeat steps 2b and 2c as desired.
fl <- flash_add_greedy(data,
                       Kmax=50,
                       f_init=fl,
                       var_type="zero",
                       init_fn="udv_svd",
                       ebnm_fn=ebnm_fn,
                       ebnm_param=ebnm_param)
fl <- flash_backfit(data,
                    fl,
                    kset=n:fl$nfactors,
                    var_type="zero",
                    ebnm_fn=ebnm_fn,
                    ebnm_param=ebnm_param,
                    nullcheck=FALSE)
saveRDS(fl, "./output/MASHvFLASHVhat/2dRepeat3.rds")
