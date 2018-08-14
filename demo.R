library(BayesianITR)

## simulate data
K_u = 3
K_f = c(3,3)
W = c(60,70)

sim_data = SimulateData(I=200, W, trajectory_plot = TRUE)

tol_iter = 100
thin = 10
mh_sd = 0.3
prior_unit = 1

## U kernel
prior_mu_log_sd2_u0 = log(0.01)
prior_sigma_log_sd2_u0 = 0.09
prior_mu_log_sd2_u = c(log(0.01),log(0.01),log(0.01))
prior_mu_logit_rho_u = c(logit(0.1),logit(0.9),logit(0.5))

## f kernel
prior_sd1_f = 0.09
prior_u10_f = c(8,-8)
prior_sd10_f = 4

prior_sdr_f = 0.09
prior_ur0_f = c(10,20)
prior_sdr0_f = 4

prior_ub0_f = c(logit(0.01),logit(0.01))

prior_mu_log_sd2_f0 = log(0.01)
prior_sigma_log_sd2_f0 = 0.09

## run model fit
ptm = proc.time()
chain = model_fit(model = 'ITR', tol_iter, thin, mh_sd,
                  NULL, NULL, W,
                  prior_mu_log_sd2_u0, prior_sigma_log_sd2_u0,
                  0, prior_u10_f, prior_sd10_f,
                  0, prior_ur0_f, prior_sdr0_f,
                  prior_ub0_f,
                  prior_mu_log_sd2_f0, prior_sigma_log_sd2_f0,
                  prior_unit,
                  sim_data,
                  NULL,0)
runtime = proc.time() - ptm
chain$runtime = runtime;

## plot model fit
plot_model_fit(chain, start_iter = 5, tol_iter = 9, prior_unit = 1, train_data = sim_data)

## store the final results
save(chain, file = "chain.itr.sim_data.RData")
##################################################
## run model fit
ptm = proc.time()
K_u = 3
K_f = c(3,3)
chain = model_fit(model = 'DPITR', tol_iter, thin, mh_sd,
                  K_u, K_f, W,
                  prior_mu_log_sd2_u0, prior_sigma_log_sd2_u0,
                  prior_sd1_f, prior_u10_f, prior_sd10_f,
                  prior_sdr_f, prior_ur0_f, prior_sdr0_f,
                  prior_ub0_f,
                  prior_mu_log_sd2_f0, prior_sigma_log_sd2_f0,
                  prior_unit,
                  sim_data, NULL,0)
runtime = proc.time() - ptm
chain$runtime = runtime;

## plot model fit
plot_model_fit(chain, start_iter = 5, tol_iter = 9, prior_unit = 1, train_data = sim_data)

## store the final results
save(chain, file = "chain.dpitr.sim_data.RData")
##################################################
## run model fit
ptm = proc.time()
K_u = 3
K_f = c(3,3)
chain = model_fit(model = 'DPMITR', tol_iter, thin, mh_sd,
                  K_u, K_f, W,
                  prior_mu_log_sd2_u0, prior_sigma_log_sd2_u0,
                  prior_sd1_f, prior_u10_f, prior_sd10_f,
                  prior_sdr_f, prior_ur0_f, prior_sdr0_f,
                  prior_ub0_f,
                  prior_mu_log_sd2_f0, prior_sigma_log_sd2_f0,
                  prior_unit,
                  sim_data, NULL,0)
runtime = proc.time() - ptm
chain$runtime = runtime;

## plot model fit
plot_model_fit(chain, start_iter = 5, tol_iter = 9, prior_unit = 1, train_data = sim_data)
## store the final results
save(chain, file = "chain.dpmitr.sim_data.RData")
##################################################

## plot data
plot_sim_data(sim_data)
