devtools::load_all("~/GitHub/flashr/")
devtools::load_all("~/GitHub/ebnm/")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
strong <- gtex$strong.z
strong_data <- flash_set_data(strong, S = 1)

my_init_fn <- function(Y, K = 1) {
  ret = udv_si(Y, K)
  pos_sum = sum(ret$v[ret$v > 0])
  neg_sum = -sum(ret$v[ret$v < 0])
  if (neg_sum > pos_sum) {
    return(list(u = -ret$u, d = ret$d, v = -ret$v))
  } else
    return(ret)
}

# Using mixSQP maybe requires try/catch?

ebnm_fn = "ebnm_ash"
ebnm_param = list(l = list(mixcompdist = "normal",
                           optmethod = "mixSQP"),
                  f = list(mixcompdist = "+uniform",
                           optmethod = "mixSQP"))

# var_type = "zero" -----------------------------------------------------

t_g_zero <- system.time(
  fl_g_zero <- flashr:::flash_greedy_workhorse(strong_data,
                                               var_type = "zero",
                                               ebnm_fn = ebnm_fn,
                                               ebnm_param = ebnm_param,
                                               init_fn = "my_init_fn",
                                               stopping_rule = "factors",
                                               tol = 1e-3,
                                               verbose_output = "odF")
) # 34 factors in 5.5 min
saveRDS(fl_g_zero, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/fl_g_zero.rds")
saveRDS(t_g_zero, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/t_g_zero.rds")

t_b_zero <- system.time(
  fl_b_zero <- flashr::flash_backfit_workhorse(strong_data,
                                               f_init = fl_g_zero,
                                               var_type = "zero",
                                               ebnm_fn = ebnm_fn,
                                               ebnm_param = ebnm_param,
                                               stopping_rule = "factors",
                                               tol = 1e-3,
                                               verbose_output = "odF")
) # backfit in 15.0 minutes
saveRDS(fl_b_zero, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/fl_b_zero.rds")
saveRDS(t_b_zero, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/t_b_zero.rds")

t_g2_zero <- system.time(
  fl_g2_zero <- flashr:::flash_greedy_workhorse(strong_data,
                                                f_init = fl_b_zero,
                                                var_type = "zero",
                                                ebnm_fn = ebnm_fn,
                                                ebnm_param = ebnm_param,
                                                init_fn = "my_init_fn",
                                                stopping_rule = "factors",
                                                tol = 1e-3,
                                                verbose_output = "odF")
) # 6 more factors in 0.8 min
saveRDS(fl_g2_zero, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/fl_g2_zero.rds")
saveRDS(t_g2_zero, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/t_g2_zero.rds")

t_b2_zero <- system.time(
  fl_b2_zero <- flashr::flash_backfit_workhorse(strong_data,
                                                f_init = fl_g2_zero,
                                                var_type = "zero",
                                                ebnm_fn = ebnm_fn,
                                                ebnm_param = ebnm_param,
                                                stopping_rule = "factors",
                                                tol = 1e-3,
                                                verbose_output = "odF")
) # backfit in 20.7 minutes
saveRDS(fl_b2_zero, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/fl_b2_zero.rds")
saveRDS(t_b2_zero, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/t_b2_zero.rds")

t_g3_zero <- system.time(
  fl_g3_zero <- flashr:::flash_greedy_workhorse(strong_data,
                                                f_init = fl_b2_zero,
                                                var_type = "zero",
                                                ebnm_fn = ebnm_fn,
                                                ebnm_param = ebnm_param,
                                                init_fn = "my_init_fn",
                                                stopping_rule = "factors",
                                                tol = 1e-3,
                                                verbose_output = "odF")
)
# Nothing added this time!


# var_type = "const" ----------------------------------------------------

t_g_const <- system.time(
  fl_g_const <- flashr:::flash_greedy_workhorse(strong,
                                                var_type = "constant",
                                                ebnm_fn = ebnm_fn,
                                                ebnm_param = ebnm_param,
                                                init_fn = "my_init_fn",
                                                stopping_rule = "factors",
                                                tol = 1e-3,
                                                verbose_output = "odF")
) # 26 factors in 3.3 min
saveRDS(fl_g_const, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/fl_g_const.rds")
saveRDS(t_g_const, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/t_g_const.rds")

t_b_const <- system.time(
  fl_b_const <- flashr::flash_backfit_workhorse(strong,
                                                f_init = fl_g_const,
                                                var_type = "constant",
                                                ebnm_fn = ebnm_fn,
                                                ebnm_param = ebnm_param,
                                                stopping_rule = "factors",
                                                tol = 1e-3,
                                                verbose_output = "odF")
) # backfit in 13.1 min
saveRDS(fl_b_const, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/fl_b_const.rds")
saveRDS(t_b_const, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/t_b_const.rds")

t_g2_const <- system.time(
  fl_g2_const <- flashr:::flash_greedy_workhorse(strong,
                                                 f_init = fl_b_const,
                                                 var_type = "constant",
                                                 ebnm_fn = ebnm_fn,
                                                 ebnm_param = ebnm_param,
                                                 init_fn = "my_init_fn",
                                                 stopping_rule = "factors",
                                                 tol = 1e-3,
                                                 verbose_output = "odF")
) # 2 more factors in 0.4 min
saveRDS(fl_g2_const, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/fl_g2_const.rds")
saveRDS(t_g2_const, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/t_g2_const.rds")

t_b2_const <- system.time(
  fl_b2_const <- flashr::flash_backfit_workhorse(strong,
                                                 f_init = fl_g2_const,
                                                 var_type = "constant",
                                                 ebnm_fn = ebnm_fn,
                                                 ebnm_param = ebnm_param,
                                                 stopping_rule = "factors",
                                                 tol = 1e-3,
                                                 verbose_output = "odF")
) # backfit in 9.4 min
saveRDS(fl_b2_const, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/fl_b2_const.rds")
saveRDS(t_b2_const, "~/GitHub/MASHvFLASH/output/MASHvFLASHnn2/t_b2_const.rds")

t_g3_const <- system.time(
  fl_g3_const <- flashr:::flash_greedy_workhorse(strong,
                                                 f_init = fl_b2_const,
                                                 var_type = "constant",
                                                 ebnm_fn = ebnm_fn,
                                                 ebnm_param = ebnm_param,
                                                 init_fn = "my_init_fn",
                                                 stopping_rule = "factors",
                                                 tol = 1e-3,
                                                 verbose_output = "odF")
)
# Nothing new added.
