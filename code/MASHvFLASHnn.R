devtools::load_all("/Users/willwerscheid/GitHub/flashr/")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
strong <- t(gtex$strong.z)
strong_data <- flash_set_data(strong, S = 1)

ebnm_param = list(f = list(), l = list(mixcompdist="+uniform"))

system.time(
  fl <- flash_add_greedy(strong_data,
                         100,
                         var_type="zero",
                         ebnm_fn="ebnm_ash",
                         ebnm_param=ebnm_param)
)

system.time(
  fl <- flash_backfit(strong_data,
                      f_init=fl,
                      var_type="zero",
                      ebnm_fn="ebnm_ash",
                      ebnm_param=ebnm_param,
                      tol=1)
)

# Repeat the following two steps until flash_add_greedy no longer adds
#   any factors:

system.time(
  fl <- flash_add_greedy(strong_data,
                         100,
                         f_init=fl,
                         var_type="zero",
                         ebnm_fn="ebnm_ash",
                         ebnm_param=ebnm_param)
)

system.time(
  fl <- flash_backfit(strong_data,
                      f_init=fl,
                      var_type="zero",
                      ebnm_fn="ebnm_ash",
                      ebnm_param=ebnm_param,
                      tol=1)
)

saveRDS(fl, "/Users/willwerscheid/GitHub/MASHvFLASH/output/MASHvFLASHnn/fl.rds")

# Tighten the tolerance and run a final backfit:

system.time(
  fl <- flash_backfit(strong_data,
                      f_init=fl,
                      var_type="zero",
                      ebnm_fn="ebnm_ash",
                      ebnm_param=ebnm_param,
                      tol=0.1)
)

saveRDS(fl, "/Users/willwerscheid/GitHub/MASHvFLASH/output/MASHvFLASHnn/fl.rds")
