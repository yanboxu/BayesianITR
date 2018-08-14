require(MASS)
require(reshape2)
require(statmod)
require(mnormt)
require(matrixcalc)
require(fGarch)
require(stats)
require(mvnfast)
require(sn)
require(MCMCpack)
require(Rcpp)
require(RcppArmadillo)

#' A Model fit Function
#'
#' This function allows you to obtain Bayesian etimation of baseline progressions (i.e. under no treatments) and treatment response curves (with treatments) from longitudinal observational data. Model can be trained only at indivdiual level (ITR), only at subgroup level (DPITR), and at both levels (DPMITR).
#' @param model Specify which model to fit; models include DPMITR, DPITR and ITR; defaults to ITR
#' @param K_u The max number of clusters for baseline progression; must be specified for model DPMITR and DPITR
#' @param K_f 1xD, the max number of clusters for each treatment type d; must be specified for model DPMITR and DPITR
#' @param W 1xD, window size for treatment type d
#'
#' @param tol_iter Total number of iterations for MCMC sampling
#' @param thin Number to thin in MCMC sampling
#' @param mh_sd Standard deviation for proposal normal distribution in metropolis-hastings
#'
#' @param prior_mu_log_sd2_u0 Mean of the sd in U's Gaussian process; defaults to log(1)
#' @param prior_sigma_log_sd2_u0 Sd of the sd in U's Gaussian process; defaults to 4
#'
#' @param prior_sd1_f Sd of the subgroup prior for peak-effect alpha_1; defaults to 0.01
#' @param prior_u10_f 1xD, mean of the base prior for peak-effect alpha_1, better specified with domain knowledge
#' @param prior_sd10_f Sd of the base prior for peak-effect alpha_1; defaults to 0.01
#'
#' @param prior_sdr_f Sd of the subgroup priors for time-to-peak-effect gamma; defaults to 0.09
#' @param prior_ur0_f 1xD, mean of the base prior for time-to-peak-effect gamma, better specified with domain knowledge
#' @param prior_sdr0_f Sd of the base prior for time-to-peak-effect gamma; defaults to 4
#'
#' @param prior_ub0_f 1xD, mean of the base prior for long-term-remaining-effect b, better specified with domain knowledge; if no long term effects, this can be set 0
#' @param prior_mu_log_sd2_f0 Mean of the sd in the noise Gaussian process
#' @param prior_sigma_log_sd2_f0 Sd of the sd in the noise Gaussian process; defaults to 4
#' @param prior_unit Unit of the prior distributions; defaults to 1
#' @param data Input data; please see SimulateData for detailed data format
#' @param model_fit Defaults to 'NULL'; otherwise specify the intermediate chain file as the current start
#' @param last_iter Defaults to 0; otherwise specify which iterations to continue with when start with an itermediate chain file
#' @references  Yanbo Xu, Yanxun Xu and Suchi Saria, 2016. A Bayesian Nonparametic Approach for Estimating Individualized Treatment-Response Curves. arXiv preprint arXiv:1608.05182. [tentatively accepted by Journal of Machine Learning Research]
#' @keywords heterogenous treatment effects
#' @keywords individualized treatment effects
#' @keywords treatment response curves
#' @author Yanbo Xu
#' @return A list of all the parameters defined in the model. \item{sigma2_u}{The individual-level estimated hyperparamters in the Gaussian process for baseline progression.}\item{rho_u}{The individual-level estimated hyperparamters in the Gaussian process for baseline progression.}\item{U}{The individual-level infered mean of the baseline progression.}\item{Beta}{The individual-level estimated Beta coefficients in the liner regression for the basline progression.}\item{a1}{The individual-level estimated peak effect alpha_1 defined in the treatment response curves.}\item{a2}{The individual-level estimated effect increasing slope alpha_2 defined in the treatment response curves.}\item{a3}{The individual-level estimated effect decreasing slope alpha_3 defined in the treatment response curves.}\item{r}{The individual-level estimated time to peak effect gamma defined in the treatment response curves.}\item{b}{The individual-level estimated remaining long-term effect b defined in the treatment response curves.}\item{rho_f}{The individual-level estimated hyperparameter in the Gaussian process prior of the treatment effect noises.}\item{sigma2_f}{The individual-level estimated hyperparameter in the Gaussian process prior of the treatment effect noises.}\item{f}{The indivdiual-level infered mean of the cumulative treatment effects}\item{Sigma_e}{The indivdiual-level estimated mean of the Gaussian noise}\item{ll}{The log-likelihood per iteration.} \item{C_u}{The estimated cluster assignments of baseline progressions; for model DPMITR and DPITR only.}\item{mu_log_sd2_u}{The subgroup-level estimated hyperparameter means in the Gaussian process prior for baseline progression; for model DPMITR and DPITR only.}\item{mu_log_rho_u}{The subgroup-level estimated hyperparameter means in the Gaussian process prior for baseline progression; for model DPMITR and DPITR only.}\item{Mu_b}{The subgroup-level estimated means for the prior of Beta; for model DPMITR and DPITR only.}\item{Sigma_b}{The subgroup-level estimated sd for the prior of Beta; for model DPMITR only.}\item{u1_f}{The subgroup-level estimated mean for the parameters defining treatment response curves; for model DPMITR and DPITR only.}\item{u2_f}{The subgroup-level estimated mean for the parameters defining treatment response curves; for model DPMITR and DPITR only.}\item{u3_f}{The subgroup-level estimated mean for the parameters defining treatment response curves; for model DPMITR and DPITR only.}\item{ur_f}{The subgroup-level estimated mean for the parameters defining treatment response curves; for model DPMITR and DPITR only.}\item{ub_f}{The subgroup-level estimated mean for the parameters defining treatment response curves; for model DPMITR and DPITR only.}\item{C_f}{The estimated cluster assignments for each treatment type d; for model DPMITR and DPITR only.}\item{M1}{The estimated hyperparameter for baseline progression clustering; for model DPMITR and DPITR only.}\item{M2}{The estimated hyperparameter for treatment effect noise clustering; for model DPMITR and DPITR only.}
#'
#' @examples
#' ## simulate data
#' K_u = 3
#' K_f = c(3,3)
#' W = c(60,70)
#'
#' sim_data = SimulateData(I=200, W, trajectory_plot = TRUE)
#'
#' tol_iter = 100
#' thin = 10
#' mh_sd = 0.3
#' prior_unit = 1
#'
#' ## U kernel
#' prior_mu_log_sd2_u0 = log(0.01)
#' prior_sigma_log_sd2_u0 = 0.09
#' prior_mu_log_sd2_u = c(log(0.01),log(0.01),log(0.01))
#' prior_mu_logit_rho_u = c(logit(0.1),logit(0.9),logit(0.5))
#'
#' ## f kernel
#' prior_sd1_f = 0.09
#' prior_u10_f = c(8,-8)
#' prior_sd10_f = 4
#'
#' prior_sdr_f = 0.01
#' prior_ur0_f = c(10,20)
#' prior_sdr0_f = 4
#'
#' prior_ub0_f = c(logit(0.01),logit(0.01))
#'
#' prior_mu_log_sd2_f0 = log(0.01)
#' prior_sigma_log_sd2_f0 = 0.09
#'
#' ## run model fit
#' ptm = proc.time()
#' chain = model_fit(model = 'DPMITR', tol_iter, thin, mh_sd,
#'                      K_u, K_f, W,
#'                      prior_mu_log_sd2_u0, prior_sigma_log_sd2_u0,
#'                      prior_sd1_f, prior_u10_f, prior_sd10_f,
#'                      prior_sdr_f, prior_ur0_f, prior_sdr0_f,
#'                      prior_ub0_f,
#'                      prior_mu_log_sd2_f0, prior_sigma_log_sd2_f0,
#'                      prior_unit,
#'                      data = sim_data,
#'                      model_fit = NULL,last_iter = 0)
#' runtime = proc.time() - ptm
#' chain$runtime = runtime;
#'
#' ## store the final results
#' save(chain, file = "chain.dpmitr.sim.RData");
#'
#' @useDynLib BayesianITR
#' @importFrom Rcpp sourceCpp
#' @export
model_fit <- function(model = 'ITR', tol_iter, thin, mh_sd,
                          K_u = NULL, K_f = NULL, W,
                          prior_mu_log_sd2_u0, prior_sigma_log_sd2_u0,
                          prior_sd1_f = 0.01, prior_u10_f, prior_sd10_f,
                          prior_sdr_f = 0.09, prior_ur0_f, prior_sdr0_f,
                          prior_ub0_f,
                          prior_mu_log_sd2_f0, prior_sigma_log_sd2_f0,
                          prior_unit = 1,
                          data, model_fit = NULL,last_iter = 1){
  ##### model includes DPMITR, DPITR, ITR ####
  if(!model %in% c('DPMITR','DPITR','ITR')){
    model = 'ITR'
  }
  print(paste('Fitting',model,'....'))
  
  #### Global variables ####
  tol_iter <- tol_iter
  mh_sd <- mh_sd
  K_u <- K_u
  K_f <- K_f
  prior_mu_log_sd2_u0 <- prior_mu_log_sd2_u0
  prior_sigma_log_sd2_u0 <- prior_sigma_log_sd2_u0
  prior_sd1_f <- prior_sd1_f
  prior_u10_f <- prior_u10_f
  prior_sd10_f <- prior_sd10_f
  prior_sdr_f <- prior_sdr_f
  prior_ur0_f <- prior_ur0_f
  prior_sdr0_f <- prior_sdr0_f
  prior_ub0_f <- prior_ub0_f
  prior_mu_log_sd2_f0 <- prior_mu_log_sd2_f0
  prior_sigma_log_sd2_f0 <- prior_sigma_log_sd2_f0
  prior_unit <- prior_unit
  model_fit <- model_fit
  last_iter <- last_iter
  
  ## Load data ##
  dimJ <- data$data_J; # maximum data points for each patient
  dimI <- data$data_I # number of patients
  dimL <- data$data_L # maximum treatment time points for each patient
  dimD <- data$data_D # teatment level
  dimP <- data$data_P; # dimension of fixed effect

  dataZ <- data$data_Z; # joint outcomes IxJ
  dataT <- data$data_T; # hospital time IxJ
  dataB <- data$data_B; # basic variables IxJxP
  dataTau <- data$data_Tau # treatment time IxL
  dataA <- data$data_A # treatments IxL
  #dataIU <- data.IU # baseline regimes
  J_I <- data$data_J_I; # data points per patient
  L_I <- data$data_L_I; # treatment points per patient

  dataIU = array(1, dim = c(dimI, dimJ));
  
  ######### Base Prior Initialization ##########
  prior_prob_sd2 <- 4; # Default sd if not secified for base priors
  
  # params in Normal(or NIW if DPM/DP) prior on Beta
  beta0 <- rep(0,dimP);
  S0 <- diag(1,dimP);
  
  # params in transformed Normal priors for U Gaussian process kernel params  
  mu_log_sd2_u0 <- log(1)
  if(!is.null(prior_mu_log_sd2_u0))
    mu_log_sd2_u0 <- prior_mu_log_sd2_u0
  sigma_log_sd2_u0 <- prior_prob_sd2
  if(!is.null(prior_sigma_log_sd2_u0))
    sigma_log_sd2_u0 <- prior_sigma_log_sd2_u0
  
  mu_logit_rho_u0 <- logit(0.5); # rho_u ~ 0.5
  sigma_logit_rho_u0 <- prior_prob_sd2;
  
  # params in transformed Normal priors for Treatment response curves params
  u10_f <- rep(0,dimD);
  sd10_f <- 1;
  if(!is.null(prior_u10_f))
    u10_f <- prior_u10_f
  if(!is.null(prior_sd10_f))
    sd10_f <- prior_sd10_f;

  u20_f <- rep(logit(0.5),dimD);
  sd20_f <- prior_prob_sd2;
  
  u30_f <- rep(logit(0.5),dimD);
  sd30_f <- prior_prob_sd2;

  ur0_f <- rep(5,dimD);
  sdr0_f <- 10;
  if(!is.null(prior_ur0_f))
    ur0_f <- prior_ur0_f
  if(!is.null(prior_sdr0_f))
    sdr0_f <- prior_sdr0_f;
  
  ub0_f <- rep(logit(0.2),dimD);
  sdb0_f <- prior_prob_sd2;
  if(!is.null(prior_ub0_f))
    ub0_f <- prior_ub0_f
  
  # params in transformed Normal priors for treatment response Gaussian process noises
  mu_log_sd2_f0 <- log(0.01);
  sigma_log_sd2_f0 <- prior_prob_sd2;
  if(!is.null(prior_mu_log_sd2_f0))
    mu_log_sd2_f0 <- prior_mu_log_sd2_f0;
  if(!is.null(prior_sigma_log_sd2_f0))
    sigma_log_sd2_f0 <- prior_sigma_log_sd2_f0;
  
  mu_logit_rho_f0 <- logit(0.9); # rho_u ~ 0.5
  sigma_logit_rho_f0 <- prior_prob_sd2;
  
  # IG params for Sigma_e[i]
  a_e <- rep(1,dimI);
  b_e <- rep(1,dimI);
  

  ######### Indivdiual-level Parameters Initialization ##########
  Beta <- array(0,dim=c(dimI,dimP));
  sigma2_u <- rep(0,dimI);
  rho_u <- rep(0,dimI);
  U <- array(-1,dim=c(dimI,dimJ));
  
  a1 <- array(0,dim=c(dimI,dimD));
  a2 <- array(0,dim=c(dimI,dimD));
  a3 <- array(0,dim=c(dimI,dimD));
  r <- array(0,dim=c(dimI,dimD));
  b <- array(0,dim=c(dimI,dimD));
  f <- array(0,dim=c(dimI,dimJ));
  
  for(i in 1:dimI){
    U[i,1:J_I[i]] = MASS::mvrnorm(1,rep(0,J_I[i]),diag(sigma2_u[i],J_I[i]));
    Beta[i,] = MASS::mvrnorm(1,beta0,S0);
    sigma2_u[i] = exp(stats::rnorm(1,mu_log_sd2_u0,sqrt(sigma_log_sd2_u0)));
    rho_u[i] = invlogit(stats::rnorm(1,mu_logit_rho_u0,sqrt(sigma_logit_rho_u0)));
    
    # a{1-3},r, parameters for g function
    for(d in 1:dimD){
      a1[i,d] = stats::rnorm(1,u10_f[d],sqrt(sd10_f));
      a2[i,d] = invlogit(stats::rnorm(1,u20_f[d],sqrt(sd20_f)));
      a3[i,d] = invlogit(stats::rnorm(1,u30_f[d],sqrt(sd30_f)));
      r[i,d] = stats::rnorm(1,ur0_f[d],sqrt(sdr0_f));
      a = a1[i,d] * (exp(a2[i,d]*r[i,d]/2) - 1)/(exp(a2[i,d]*r[i,d]/2) + 1);
      b[i,d] = invlogit(stats::rnorm(1,ub0_f[d],sqrt(sdb0_f))) * a;
    }
  }
  
  ######### Noise Initialization ##########
  sigma2_f <- exp(stats::rnorm(dimD, mu_log_sd2_f0, sqrt(sigma_log_sd2_f0)))
  rho_f <- invlogit(stats::rnorm(dimD, mu_logit_rho_f0, sqrt(sigma_logit_rho_f0)))
  
  Sigma_e_I <- rep(0,dimI);
  for(i in 1:dimI)
    Sigma_e_I[i] <- MCMCpack::rinvgamma(1,a_e[i],b_e[i])
  
  ############################ Read params from existing model_fit #######################
  if(!is.null(model_fit)){
    Sigma_e_I <- model_fit$Sigma_e[last_iter,];
    
    a1 <- model_fit$a1[last_iter,,];
    a2 <- model_fit$a2[last_iter,,];
    a3 <- model_fit$a3[last_iter,,];
    r <- model_fit$r[last_iter,,];
    b <- model_fit$b[last_iter,,];
    
    rho_f <- model_fit$rho_f[last_iter,];
    sigma2_f <- model_fit$sigma2_f[last_iter,];
    f <- model_fit$f[last_iter,,];
    
    Beta <- model_fit$Beta[last_iter,,];
    sigma2_u <- model_fit$sigma2_u[last_iter,];
    rho_u <- model_fit$rho_u[last_iter,];
    U <- model_fit$U[last_iter,,];

  }

  ########################### Prepare inputs for model sampling #####################
  dataB_3D = array(0, dim=c(dimP, dimI * dimJ));
  for(i in 1:dimI)
    for(j in 1:dimJ)
      for(p in 1:dimP)
        dataB_3D[p, (i-1)*dimJ + j] = dataB[i,j,p];

  data_list = list("dimI" = dimI, "dimJ" = dimJ, "dimL" = dimL, "dimD" = dimD,
                   "dimP" = dimP, "Z" = dataZ, "T" = dataT, "B" = dataB_3D,
                   "Tau" = dataTau, "A" = dataA, "IU" = dataIU,
                   "J_I" = J_I, "L_I" = L_I)
  data_list$A = data_list$A - 1;
  
  input_vars = list("tol_iter" = tol_iter, "thin" = thin, "mh_sd" = mh_sd,
                    "prior_unit" = prior_unit, "W" = W)

  hyper_param_list = list( "beta0" = beta0, "S0" = S0,
                          "a_e" = a_e, "b_e" = b_e, 
                          "mu_log_sd2_u0" = mu_log_sd2_u0, "sigma_log_sd2_u0" = sigma_log_sd2_u0,
                          "mu_logit_rho_u0" = mu_logit_rho_u0,
                          "sigma_logit_rho_u0" = sigma_logit_rho_u0,
                          "u10_f" = u10_f, "sd10_f" = sd10_f,
                          "u20_f" = u20_f, "sd20_f" = sd20_f,
                          "u30_f" = u30_f, "sd30_f" = sd30_f,
                          "ur0_f" = ur0_f, "sdr0_f" = sdr0_f,
                          "ub0_f" = ub0_f, "sdb0_f" = sdb0_f,
                          "mu_log_sd2_f0" = mu_log_sd2_f0, "sigma_log_sd2_f0" = sigma_log_sd2_f0,
                          "mu_logit_rho_f0" = mu_logit_rho_f0, "sigma_logit_rho_f0" = sigma_logit_rho_f0)

  param_list = list("Sigma_e_I" = Sigma_e_I, 
                    "Beta" = Beta, 
                    "sigma2_u" = sigma2_u, "rho_u" = rho_u, "U" = U, 
                    "a1" = a1, "a2" = a2, "a3" = a3, "r" = r, "b" = b,
                    "rho_f" = rho_f, "sigma2_f" = sigma2_f, "f" = f)
  
  ##### model is DPMITR or DPITR #########
  if(model == 'DPMITR'){
    if(is.null(K_u)){
      print(paste('Error: K_u needs to be specified for', model))
      return()
    }
    if(is.null(K_f)){
      print(paste('Error: K_f needs to be specified for', model))
      return()
    }
    
    K <- max(K_f)
    
    ########## Subgroup-level Prior Initialization #########
    prior_cluster_sd2 <- 0.09; # default sd within clusters
    prior_a1_cluster_sd2 <- 0.01; # sd for alpha_1 within clusters
    
    # NIW extra params for Mu_b, Sigma_b
    k0 <- 1;
    v0 <- dimP+1;
    
    # params in subgroup-level Normal priors for U Gaussian process kernel params 
    sigma_log_sd2_u <- prior_cluster_sd2;
    sigma_logit_rho_u <- prior_cluster_sd2;
    
    # params in subgroup-level Normal priors for Treatment response curves params
    sd1_f <- prior_a1_cluster_sd2;
    if(!is.null(prior_sd1_f))
      sd1_f <- prior_sd1_f;
    
    sd2_f <- prior_cluster_sd2;
    sd3_f <- prior_cluster_sd2;
    
    sdr_f <- prior_cluster_sd2;
    if(!is.null(prior_sdr_f))
      sdr_f <- prior_sdr_f
    
    sdb_f <- prior_cluster_sd2;
    
    # DPM parameter M for clustering
    M_c <- 2;
    M_d <- 2;
    
    ################## Subgroup-level Parameter initialization ##################
    # C_u[i], cluster indicator for U[i]
    C_u <- sample(1:K_u,dimI,replace = TRUE);
    
    # mu_b[k], Sigma_b[k], hyperparameters for kth cluster for Beta
    Mu_b <- array(0,dim=c(K_u,dimP));
    Sigma_b <- array(0,dim=c(K_u,dimP,dimP));
    for(k in 1:K_u)
      Sigma_b[k,,] <- diag(0.09,dimP); #xxxxxxxxxxxxxxxxxxx 0.2
    
    # u{1-3},r_f
    # C_f[i] cluster indicator for f
    u1_f <- array(0,dim=c(K,dimD));
    u2_f <- array(0,dim=c(K,dimD));
    u3_f <- array(0,dim=c(K,dimD));
    ur_f <- array(0,dim=c(K,dimD));
    ub_f <- array(0,dim=c(K,dimD));
    C_f <- array(0,dim=c(dimI,dimD))
    for(d in 1:dimD){
      u1_f[1:K_f[d],d] = stats::rnorm(K_f[d],u10_f[d],sqrt(sd10_f));
      u2_f[1:K_f[d],d] = stats::rnorm(K_f[d],u20_f[d],sqrt(sd20_f));
      u3_f[1:K_f[d],d] = stats::rnorm(K_f[d],u30_f[d],sqrt(sd30_f));
      ur_f[1:K_f[d],d] = stats::rnorm(K_f[d],ur0_f[d],sqrt(sdr0_f));
      ub_f[1:K_f[d],d] = stats::rnorm(K_f[d],ub0_f[d],sqrt(sdb0_f));
      C_f[,d] = sample(1:K_f[d],dimI,replace = TRUE);
    }
    
    # mu_log_sd2_u[k], mu_logit_rho_u[k], hyperparameters for kth cluster of sigma2_u, rho_u
    mu_log_sd2_u <- stats::rnorm(K_u,mu_log_sd2_u0,sqrt(sigma_log_sd2_u0));
    mu_logit_rho_u <- stats::rnorm(K_u,mu_logit_rho_u0,sqrt(sigma_logit_rho_u0));
    
    # M1, M2
    M1 <- stats::rgamma(1, M_c, M_d);
    M2 <- stats::rgamma(dimD, M_c, M_d);
    
    ######### Indivdiual-level Parameters Initialization ##########
    for(i in 1:dimI){
      # Beta
      Beta[i,] = MASS::mvrnorm(1,Mu_b[C_u[i],],Sigma_b[C_u[i],,]); 
      
      # sigma2_u, rho_u
      sigma2_u[i] = exp(stats::rnorm(1,mu_log_sd2_u[C_u[i]],sqrt(sigma_log_sd2_u)));
      rho_u[i] = invlogit(stats::rnorm(1,mu_logit_rho_u[C_u[i]],sqrt(sigma_logit_rho_u)));
      
      # U_i, within subject random effect
      U[i,1:J_I[i]] = MASS::mvrnorm(1,rep(0,J_I[i]),diag(sigma2_u[i],J_I[i]));
      
      # a{1-3},r, parameters for g function
      for(d in 1:dimD){
        a1[i,d] = stats::rnorm(1,u1_f[C_f[i,d],d],sqrt(sd1_f));
        a2[i,d] = invlogit(stats::rnorm(1,u2_f[C_f[i,d],d],sqrt(sd2_f)));
        a3[i,d] = invlogit(stats::rnorm(1,u3_f[C_f[i,d],d],sqrt(sd3_f)));
        r[i,d] = stats::rnorm(1,ur_f[C_f[i,d],d],sqrt(sdr_f));
        a = a1[i,d] * (exp(a2[i,d]*r[i,d]/2) - 1)/(exp(a2[i,d]*r[i,d]/2) + 1);
        b[i,d] = invlogit(stats::rnorm(1,ub_f[C_f[i,d],d],sqrt(sdb_f))) * a;
      }
    }
    
    C_f = C_f - 1
    C_u = C_u - 1
    ######### Read params from existing model_fit #########
    if(!is.null(model_fit)){
      if(model_fit$model == 'DPMITR'){
        C_u = model_fit$C_u[last_iter,];
        Mu_b = model_fit$Mu_b[last_iter,,];
        Sigma_b = model_fit$Sigma_b[last_iter,,,];
        
        mu_log_sd2_u = model_fit$mu_log_sd2_u[last_iter,];
        mu_logit_rho_u = model_fit$mu_logit_rho_u[last_iter,];
        
        C_f = model_fit$C_f[last_iter,,];
        u1_f = model_fit$u1_f[last_iter,,];
        u2_f = model_fit$u2_f[last_iter,,];
        u3_f = model_fit$u3_f[last_iter,,];
        ur_f = model_fit$ur_f[last_iter,,];
        ub_f = model_fit$ub_f[last_iter,,];
        
        M1 = model_fit$M1[last_iter];
        M2 = model_fit$M2[last_iter,];
      }
    }

    ######### Prepare inputs for model sampling #########  
    Sigma_b_3D = array(0, dim=c(dimP, dimP * K_u));
    for(k in 1:K_u)
      for(p in 1:dimP)
        for(pp in 1:dimP)
          Sigma_b_3D[p, (k-1)*dimP + pp] = Sigma_b[k,p,pp];
    
    input_vars = c(input_vars, list("K_u" = K_u, "K_f" = K_f))
    
    hyper_param_list = c(hyper_param_list, list("k0" = k0, "v0" = v0, "M_c" = M_c, "M_d" = M_d, "sigma_log_sd2_u" = sigma_log_sd2_u, "sigma_logit_rho_u" = sigma_logit_rho_u, "sd1_f" = sd1_f, "sd2_f" = sd2_f, "sd3_f" = sd3_f, "sdr_f" = sdr_f, "sdb_f" = sdb_f))
    
    param_list = c(param_list, list("C_u" = C_u, "Mu_b" = Mu_b, "Sigma_b_3D" = Sigma_b_3D, "mu_log_sd2_u" = mu_log_sd2_u, "mu_logit_rho_u" = mu_logit_rho_u, "C_f" = C_f, "u1_f" = u1_f, "u2_f" = u2_f, "u3_f" = u3_f, "ur_f" = ur_f, "ub_f" = ub_f, "M1" = M1, "M2" = M2))
  }

  if(model == 'DPITR'){
    if(is.null(K_u)){
      print(paste('Error: K_u needs to be specified for', model))
      return()
    }
    if(is.null(K_f)){
      print(paste('Error: K_f needs to be specified for', model))
      return()
    }
    
    K <- max(K_f)
    
    ################## Subgroup-level Parameter initialization ##################
    # DPM parameter M for clustering
    M_c <- 2;
    M_d <- 2;
    
    # C_u[i], cluster indicator for U[i]
    C_u <- sample(1:K_u,dimI,replace = TRUE);
    
    # mu_b[k], Sigma_b[k], hyperparameters for kth cluster for Beta
    Mu_b <- array(0,dim=c(K_u,dimP));
   
    # u{1-3},r_f
    # C_f[i] cluster indicator for f
    u1_f <- array(0,dim=c(K,dimD));
    u2_f <- array(0,dim=c(K,dimD));
    u3_f <- array(0,dim=c(K,dimD));
    ur_f <- array(0,dim=c(K,dimD));
    ub_f <- array(0,dim=c(K,dimD));
    C_f <- array(0,dim=c(dimI,dimD))
    for(d in 1:dimD){
      u1_f[1:K_f[d],d] = stats::rnorm(K_f[d],u10_f[d],sqrt(sd10_f));
      u2_f[1:K_f[d],d] = stats::rnorm(K_f[d],u20_f[d],sqrt(sd20_f));
      u3_f[1:K_f[d],d] = stats::rnorm(K_f[d],u30_f[d],sqrt(sd30_f));
      ur_f[1:K_f[d],d] = stats::rnorm(K_f[d],ur0_f[d],sqrt(sdr0_f));
      ub_f[1:K_f[d],d] = stats::rnorm(K_f[d],ub0_f[d],sqrt(sdb0_f));
      C_f[,d] = sample(1:K_f[d],dimI,replace = TRUE);
    }
    
    # mu_log_sd2_u[k], mu_logit_rho_u[k], hyperparameters for kth cluster of sigma2_u, rho_u
    mu_log_sd2_u <- stats::rnorm(K_u,mu_log_sd2_u0,sqrt(sigma_log_sd2_u0));
    mu_logit_rho_u <- stats::rnorm(K_u,mu_logit_rho_u0,sqrt(sigma_logit_rho_u0));
    
    # M1, M2
    M1 <- stats::rgamma(1, M_c, M_d);
    M2 <- stats::rgamma(dimD, M_c, M_d);
    
    ######### Indivdiual-level Parameters Initialization ##########
    for(i in 1:dimI){
      # Beta
      Beta[i,] = Mu_b[C_u[i],];
      
      # sigma2_u, rho_u
      sigma2_u[i] = exp(mu_log_sd2_u[C_u[i]]);
      rho_u[i] = invlogit(mu_logit_rho_u[C_u[i]]);
      
      # U_i, within subject random effect
      U[i,1:J_I[i]] = MASS::mvrnorm(1,rep(0,J_I[i]),diag(sigma2_u[i],J_I[i]));
      
      # a{1-3},r, parameters for g function
      for(d in 1:dimD){
        a1[i,d] = u1_f[C_f[i,d],d];
        a2[i,d] = invlogit(u2_f[C_f[i,d],d]);
        a3[i,d] = invlogit(u3_f[C_f[i,d],d]);
        r[i,d] = ur_f[C_f[i,d],d];
        a = a1[i,d] * (exp(a2[i,d]*r[i,d]/2) - 1)/(exp(a2[i,d]*r[i,d]/2) + 1);
        b[i,d] = invlogit(ub_f[C_f[i,d],d]) * a;
      }
    }
    
    C_f = C_f - 1
    C_u = C_u - 1
    ######### Read params from existing model_fit #########
    if(!is.null(model_fit)){
      if(model_fit$model != 'ITR'){
        C_u = model_fit$C_u[last_iter,];
        Mu_b = model_fit$Mu_b[last_iter,,];
        
        mu_log_sd2_u = model_fit$mu_log_sd2_u[last_iter,];
        mu_logit_rho_u = model_fit$mu_logit_rho_u[last_iter,];
        
        C_f = model_fit$C_f[last_iter,,];
        u1_f = model_fit$u1_f[last_iter,,];
        u2_f = model_fit$u2_f[last_iter,,];
        u3_f = model_fit$u3_f[last_iter,,];
        ur_f = model_fit$ur_f[last_iter,,];
        ub_f = model_fit$ub_f[last_iter,,];
        
        M1 = model_fit$M1[last_iter];
        M2 = model_fit$M2[last_iter,];
      }
    }
    
    ######### Prepare inputs for model sampling #########  
    
    input_vars = c(input_vars, list("K_u" = K_u, "K_f" = K_f))
    
    hyper_param_list = c(hyper_param_list, list("M_c" = M_c, "M_d" = M_d))
    
    param_list = c(param_list, list("C_u" = C_u, "Mu_b" = Mu_b, "mu_log_sd2_u" = mu_log_sd2_u, "mu_logit_rho_u" = mu_logit_rho_u, "C_f" = C_f, "u1_f" = u1_f, "u2_f" = u2_f, "u3_f" = u3_f, "ur_f" = ur_f, "ub_f" = ub_f, "M1" = M1, "M2" = M2))
  }
  
  
  ##### Sampling ####
  if(model == 'ITR'){
    modelFit = sampling_ind_C(input_vars,
                              data_list,
                              hyper_param_list,
                              param_list,
                              mvnfast::rmvn, mvnfast::dmvn, MCMCpack::rinvgamma)
    modelFit['model'] = 'ITR'
  }
  
  if(model == 'DPMITR'){
    modelFit = sampling_dpm_C(input_vars,
                              data_list,
                              hyper_param_list,
                              param_list,
                              mvnfast::rmvn, mvnfast::dmvn, MCMCpack::rinvgamma, MCMCpack::riwish)
    modelFit["model"] = "DPMITR"
  }
  
  if(model == 'DPITR'){
    modelFit = sampling_dp_C(input_vars,
                              data_list,
                              hyper_param_list,
                              param_list,
                              mvnfast::rmvn, mvnfast::dmvn, MCMCpack::rinvgamma, MCMCpack::riwish)
    modelFit["model"] = "DPITR"
  }
  
  ##### Sampling End ####
  return(modelFit);
}

#' logit function
#' 
#' @param p a probability
#' @return log(p/(1-p))
#' @examples 
#' logit(0.5)
#' @export
logit <- function(p){
  return(log(p/(1-p)));
}

#' inverse logit function
#' 
#' @param x a real number
#' @return exp(x)/(1+exp(x))
#' @examples 
#' invlogit(0.1)
#' @export
invlogit <- function(x){
  return(exp(x)/(1+exp(x)));
}
