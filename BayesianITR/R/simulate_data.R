require(MASS)
require(statmod)
require(mnormt)
require(MCMCpack)
require(fGarch)
require(mvtnorm)
require(sn)

#' A data simulation function
#' 
#' This functions generates heterogenous longitudinal data, in which baseline progressions are simulated from 3 clusters, two treatment types are generated (one type has positive effects, the other has negative effects), and the treatment responses for each type are respectively simulated from 3 clusters.
#' 
#' @param I Size of the cohort, i.e. number of trajectories
#' @param W 1xD, window size for treatment type d
#' @param trajectory_plot Defaults to TRUE; TRUE if need to plot upto 20 subsample trajectories from the simulated data. 
#' @author Yanbo Xu
#' @return A list of the simulated input data and the ground truth parameters that were used for data generation. \item{data_I}{Size of the cohort.}\item{data_D}{Number of treatment types.}\item{data_P}{Number of baseline variables in the linear regression.}\item{data_J_I}{1xI, the number of observations in each trajectory i.}\item{data_J}{The maximum number of observations that can be generated in each trajectory, i.e. max(data_J_I).}\item{data_L_I}{1xI, the number of treatments managed in each trajectory i.}\item{data_L}{The maximum number of treatments that can be generated in each trajectory, i.e. max(data_L_I).}\item{data_T_I}{1xI, the duration of each trajectory i.}\item{data_T}{IxJ, data_T[i,1:J_I[i]] contains the irregular time points at which observations were generated in trajectory i.}\item{data_Tau}{IxL, data_Tau[i,1:L_I[i]] contains the irregular time points at which treatments were managed in trajectory i.}\item{data_B}{IxJxP, the simulated baseline variables in the linear regression.}\item{data_A}{IxL, the sequence of treatment types managed in each trajectory i.}\item{data_Z}{IxJ, the sequene of the final observations for each trajectory i.}\item{data_avg_output}{Defaults to 0; the standardized mean of the data.}\item{data_sd_output}{Defaults to 1; the standardized sd of the data.}\item{param_a1}{IxD, the ground truth parameter for individual-level peak effects alpha_1.}\item{param_a2}{IxD, the ground truth parameter for individual-level effect increasing slope alpha_2.}\item{param_a3}{IxD, the ground truth parameter for individual-level effect decreasing slope alpha_3.}\item{param_r}{IxD, the ground truth parameter for individual-level time to peak effect gamma.}\item{param_b}{IxD, the ground truth parameter for individual-level remaining long-term effect b.}\item{param_C_f}{IxD, the true cluster assignment of the treatment effects in trajectory i.}\item{param_u1_f}{max(K_f)xD, the ground truth subgroup-level mean of the peak effects.}\item{param_u2_f}{max(K_f)xD, the ground truth subgroup-level mean of the effect increasing slope.}\item{param_u3_f}{max(K_f)xD, the ground truth subgroup-level mean of the effect decreasing slope.}\item{param_ur_f}{max(K_f)xD, the ground truth subgroup-level mean of the time to peak effect.}\item{param_ub_f}{max(K_f)xD, the ground truth subgroup-level mean of the remaining long-term effect.}\item{param_Beta}{IxP, the ground truth individual-level linear regression coefficients in the baseline progressions.}\item{param_sigma2_u}{1xI, the ground truth individual-level Gaussian process hyperparameters in the baseline progression}\item{param_rho_u}{1xI, the ground truth individual-level Gaussian process hyperparameters in the baseline progression}\item{param_U}{IxJ, the ground truth individual-level baseline progressions.}\item{param_C_u}{1xI, the true cluster assignment of the baseline progression in trajectory i.}\item{param_Mu_b}{K_uxP, the ground truth subgroup-level means of the linear regression coefficients in the baseline progression}\item{param_mu_log_sd2_u}{1xK_u, the ground truth subgroup-level means of the Gaussian process hyperparameters in the baseline progression}\item{param_mu_logit_rho_u}{1xK_u, the ground truth subgroup-level means of the Gaussian process hyperparameters in the baseline progression}\item{param_Sigma_e}{1xI, the ground truth Gaussian noises.}\item{param_sigma2_f}{1xD, the ground truth hyperparameters for the Gaussian process treatment effect noises.}\item{param_rho_f}{1xD, the ground truth hyperparameters for the Gaussian process treatment effect noises.}
#' @examples
#' SimulateData(I=100, W=c(60,70), trajectory_plot = TRUE)
#' @export
SimulateData <- function(I, W, trajectory_plot = TRUE){
  return(simulate_data(I, minT = 12, maxT = 24, unit = 60, minDur = 5, maxDur = 20, minADur = 60, maxADur = 80, W, prior_unit = 1, trajectory_plot))
}

#' A more detailed data simulation function
#' 
#' This functions generates heterogenous longitudinal data, in which baseline progressions are simulated from 3 clusters, two treatment types are generated (one type has positive effects, the other has negative effects), and the treatment responses for each type are respectively simulated from 3 clusters.
#' 
#' @param I Size of the cohort, i.e. number of trajectories
#' @param W 1xD, window size for treatment type d
#' @param trajectory_plot Defaults to TRUE; TRUE if need to plot upto 20 subsample trajectories from the simulated data. 
#' @param minT The minimum total duration of the trajectories generated in the data
#' @param maxT The maximum total duration of the trajectories generated in the data
#' @param unit The time unit of minT and maxT, e.g. unit = 60 means a unit of hour
#' @param minDur The minimum duration to generate next observation in a trajectory
#' @param maxDur The maximum duration to generate next observation in a trajectory
#' @param minADur The minimum duration to generate next treatment management in a trajectory
#' @param maxADur The maximum duration to generate next treatment management in a trajectory
#' @param prior_unit The unit of the prior distributions; defaults to 1
#' @author Yanbo Xu
#' @return A list of the simulated input data and the ground truth parameters that were used for data generation. \item{data_I}{Size of the cohort.}\item{data_D}{Number of treatment types.}\item{data_P}{Number of baseline variables in the linear regression.}\item{data_J_I}{1xI, the number of observations in each trajectory i.}\item{data_J}{The maximum number of observations that can be generated in each trajectory, i.e. max(data_J_I).}\item{data_L_I}{1xI, the number of treatments managed in each trajectory i.}\item{data_L}{The maximum number of treatments that can be generated in each trajectory, i.e. max(data_L_I).}\item{data_T_I}{1xI, the duration of each trajectory i.}\item{data_T}{IxJ, data_T[i,1:J_I[i]] contains the irregular time points at which observations were generated in trajectory i.}\item{data_Tau}{IxL, data_Tau[i,1:L_I[i]] contains the irregular time points at which treatments were managed in trajectory i.}\item{data_B}{IxJxP, the simulated baseline variables in the linear regression.}\item{data_A}{IxL, the sequence of treatment types managed in each trajectory i.}\item{data_Z}{IxJ, the sequene of the final observations for each trajectory i.}\item{data_avg_output}{Defaults to 0; the standardized mean of the data.}\item{data_sd_output}{Defaults to 1; the standardized sd of the data.}\item{param_a1}{IxD, the ground truth parameter for individual-level peak effects alpha_1.}\item{param_a2}{IxD, the ground truth parameter for individual-level effect increasing slope alpha_2.}\item{param_a3}{IxD, the ground truth parameter for individual-level effect decreasing slope alpha_3.}\item{param_r}{IxD, the ground truth parameter for individual-level time to peak effect gamma.}\item{param_b}{IxD, the ground truth parameter for individual-level remaining long-term effect b.}\item{param_C_f}{IxD, the true cluster assignment of the treatment effects in trajectory i.}\item{param_u1_f}{max(K_f)xD, the ground truth subgroup-level mean of the peak effects.}\item{param_u2_f}{max(K_f)xD, the ground truth subgroup-level mean of the effect increasing slope.}\item{param_u3_f}{max(K_f)xD, the ground truth subgroup-level mean of the effect decreasing slope.}\item{param_ur_f}{max(K_f)xD, the ground truth subgroup-level mean of the time to peak effect.}\item{param_ub_f}{max(K_f)xD, the ground truth subgroup-level mean of the remaining long-term effect.}\item{param_Beta}{IxP, the ground truth individual-level linear regression coefficients in the baseline progressions.}\item{param_sigma2_u}{1xI, the ground truth individual-level Gaussian process hyperparameters in the baseline progression}\item{param_rho_u}{1xI, the ground truth individual-level Gaussian process hyperparameters in the baseline progression}\item{param_U}{IxJ, the ground truth individual-level baseline progressions.}\item{param_C_u}{1xI, the true cluster assignment of the baseline progression in trajectory i.}\item{param_Mu_b}{K_uxP, the ground truth subgroup-level means of the linear regression coefficients in the baseline progression}\item{param_mu_log_sd2_u}{1xK_u, the ground truth subgroup-level means of the Gaussian process hyperparameters in the baseline progression}\item{param_mu_logit_rho_u}{1xK_u, the ground truth subgroup-level means of the Gaussian process hyperparameters in the baseline progression}\item{param_Sigma_e}{1xI, the ground truth Gaussian noises.}\item{param_sigma2_f}{1xD, the ground truth hyperparameters for the Gaussian process treatment effect noises.}\item{param_rho_f}{1xD, the ground truth hyperparameters for the Gaussian process treatment effect noises.}
#' @examples
#' simulate_data(I = 200, minT = 12, maxT = 24, unit = 60, minDur = 5, maxDur = 20, minADur = 60, maxADur = 80, W = c(60,70), 
#'               prior_unit = 1, trajectory_plot = TRUE)
#' @export
simulate_data <- function(I, minT, maxT, unit, minDur, maxDur, minADur, maxADur, W,
                          prior_unit = 1, trajectory_plot = TRUE){
  # basics
  D = 2;
  K_u = 3;
  K = 3;
  K_f = rep(K,D);
  P = 3;
  
  T_I = sample(minT:maxT,I, replace = TRUE) * unit; 
  L_I = rep(0, I);
  J_I = rep(0, I);
  
  J = maxT*unit/minDur;
  maxJ = J;
  L = maxT*unit/minADur;
  T = array(0, dim=c(I,J))
  Tau = array(0, dim=c(I,L))
  
  B = array(-1,dim =c(I,J,P));
  Z = array(0,dim=c(I,J));
  U = array(0,dim=c(I,J));
  f = array(0,dim=c(I,J));
  A = array(0,dim=c(I,L));
  
  C_u = sample(1:K_u,I,replace = TRUE);
  C_f = array(0,dim=c(I,D));
  for(d in 1:D){
    C_f[,d] = sample(1:(K_f[d]),I,replace = TRUE);
  }
  
  # initialize parameters
  Mu_b = array(0,dim=c(K_u,P))
  Beta = array(0,dim=c(I,P))
  
  mu_log_sd2_u = rep(0,K_u)
  mu_logit_rho_u = rep(0,K_u);
  sigma2_u = rep(0,I)
  rho_u = rep(0,I)
  
  u1_f = array(0,dim=c(K,D));
  u2_f = array(0,dim=c(K,D));
  u3_f = array(0,dim=c(K,D));
  ur_f = array(0,dim=c(K,D));
  ub_f = array(0,dim=c(K,D));
  
  a1 = array(0,dim=c(I,D));
  a2 = array(0,dim=c(I,D));
  a3 = array(0,dim=c(I,D));
  r = array(0,dim=c(I,D));
  b = array(0,dim=c(I,D));
  
  sigma2_f = rep(0,D)
  rho_f = rep(0,D)
  
  Sigma_e = rep(0.09,I)
  
  # specify paramters
  prior_sd2_cluster_sd2 = 0.01; 
  prior_rho_cluster_sd2 = 0.09;
  prior_a1_cluster_sd2 = 0.09;
  prior_r_cluster_sd2 = 0.09;
  
  Mu_b[1,] = c(5,5,3)
  Mu_b[2,] = c(30,-5,-3)
  Mu_b[3,] = c(10,-2,-1)
  
  mu_log_sd2_u = c(log(0.01),log(0.04),log(0.09))
  mu_logit_rho_u = c(logit(0.1),logit(0.9),logit(0.5))
  
  sigma2_f = c(0.01,0.01)
  rho_f = c(0.9,0.9)
  
  u1_f[1,1] = 10;
  u2_f[1,1] = logit(0.9);
  u3_f[1,1] = logit(0.4);
  ur_f[1,1] = 10;
  ub_f[1,1] = logit(0.001);
  
  u1_f[2,1] = 5;
  u2_f[2,1] = logit(0.9);
  u3_f[2,1] = logit(0.9);
  ur_f[2,1] = 5;
  ub_f[2,1] = logit(0.001);
  
  u1_f[3,1] = 8;
  u2_f[3,1] = logit(0.7);
  u3_f[3,1] = logit(0.7);
  ur_f[3,1] = 15;
  ub_f[3,1] = logit(0.001);
  
  u1_f[1,2] = -10;
  u2_f[1,2] = logit(0.9);
  u3_f[1,2] = logit(0.7);
  ur_f[1,2] = 20;
  ub_f[1,2] = logit(0.001);
  
  u1_f[2,2] = -6;
  u2_f[2,2] = logit(0.5);
  u3_f[2,2] = logit(0.5);
  ur_f[2,2] = 15;
  ub_f[2,2] = logit(0.001);
  
  u1_f[3,2] = -8;
  u2_f[3,2] = logit(0.4);
  u3_f[3,2] = logit(0.3);
  ur_f[3,2] = 25;
  ub_f[3,2] = logit(0.001);
  
  for(i in 1:I){
    Beta[i,] = mvrnorm(1,Mu_b[C_u[i],],diag(prior_sd2_cluster_sd2, P));
    sigma2_u[i] = exp(rnorm(1,mu_log_sd2_u[C_u[i]],sqrt(prior_sd2_cluster_sd2)));
    rho_u[i] = invlogit(rnorm(1,mu_logit_rho_u[C_u[i]],sqrt(prior_rho_cluster_sd2)));
    for(d in 1:D){
      a1[i,d] = rnorm(1,u1_f[C_f[i,d],d],sqrt(prior_a1_cluster_sd2));
      a2[i,d] = invlogit(rnorm(1,u2_f[C_f[i,d],d],sqrt(prior_rho_cluster_sd2)));
      a3[i,d] = invlogit(rnorm(1,u3_f[C_f[i,d],d],sqrt(prior_rho_cluster_sd2)));
      r[i,d] = rnorm(1,ur_f[C_f[i,d],d],sqrt(prior_r_cluster_sd2));
      a = a1[i,d] * (exp(a2[i,d]*r[i,d]/2) - 1)/(exp(a2[i,d]*r[i,d]/2) + 1);
      b[i,d] = invlogit(rnorm(1,ub_f[C_f[i,d],d],sqrt(prior_rho_cluster_sd2))) * a;
    }
  }
  
  avg_obs_after_A = 0;
  # generate data
  for(i in 1:I){
    curT = sample(minDur:maxDur,1);
    j = 0;
    while(curT <= T_I[i]){
      j = j+1;
      T[i,j] = curT;
      B[i,j,1] = 1;
      B[i,j,2] = curT/(maxJ*5);
      B[i,j,3] = B[i,j,2]^2;
      U[i,j] = B[i,j,] %*% Beta[i,];
      curT = curT + sample(minDur:maxDur,1);
    }
    J_I[i] = j;
    
    Kcov_u_i = array(0,dim=c(J_I[i],J_I[i]));
    for(j in 1:J_I[i]){
      for(jj in 1:J_I[i])
        Kcov_u_i[j,jj] = sigma2_u[i] * rho_u[i]^(abs(T[i,j] - T[i,jj])/60);
      Kcov_u_i[j,j] = Kcov_u_i[j,j];
    }
    U[i,1:J_I[i]] = U[i,1:J_I[i]] + mvrnorm(1, rep(0,J_I[i]), Kcov_u_i);
    Z[i,1:J_I[i]] = U[i,1:J_I[i]] + mvrnorm(1,rep(0,J_I[i]),diag(Sigma_e[i],J_I[i]));
    
    curTau = sample(minADur:maxADur,1);
    l = 0;
    while(curTau <= T[i,J_I[i]]){
      l = l + 1;
      Tau[i,l] = curTau;
      curTau = curTau + sample(minADur:maxADur,1);
    }
    L_I[i] = l;
  }
  
  for(i in 1:I){
    Kcov_f_i = array(0,dim=c(J_I[i],J_I[i]));
    for(l in 1:L_I[i]){
      d = sample(0:2, 1);
      A[i,l] = d;
      if(d > 0){
        Kcov_f_i = array(0,dim=c(J_I[i],J_I[i]));
        mu_f_i = rep(0,J_I[i]);
        for(j in 1:J_I[i]){
          deltaT1 = (T[i,j] -Tau[i,l])/prior_unit;
          if(deltaT1 <= 0) next;
          if(deltaT1 > W[d]) break;
          b0 = -a1[i,d]/(1+exp(a2[i,d]*r[i,d]/2));
          a0 = (a1[i,d]+2*b0-b[i,d])/(1+exp(-a3[i,d]*r[i,d]/2))
          if(deltaT1 <= r[i,d]){
            mu_f_i[j] = b0 + a1[i,d]/(1+exp(-a2[i,d]*(deltaT1 - r[i,d]/2)));
          }else{
            mu_f_i[j] = b[i,d] + a0/(1+exp(a3[i,d]*(deltaT1-3*r[i,d]/2)));
          }
          
          Kcov_f_i[j,j] = sigma2_f[d];
          if(j == J_I[i]) next;
          for(jj in (j+1):J_I[i]){
            deltaT2 = (T[i,jj] - Tau[i,l])/prior_unit;
            if(deltaT2 > W[d]) break;
            Kcov_f_i[j,jj] = sigma2_f[d] * (rho_f[d])^((T[i,jj] - T[i,j])/60)
            Kcov_f_i[jj,j] = Kcov_f_i[j,jj]
          }
        }
        idx_f = which(mu_f_i != 0);
        J_f = length(idx_f)
        avg_obs_after_A = avg_obs_after_A + J_f;
        if(J_f == 1)
          Fnoise = rnorm(1,0,sqrt(Kcov_f_i[idx_f,idx_f]))
        if(J_f > 1)
          Fnoise = mvrnorm(1,rep(0,J_f),Kcov_f_i[idx_f,idx_f])
        Z[i,idx_f] = Z[i,idx_f] + mu_f_i[idx_f] + Fnoise;
      }
    }
    
    idx_a = which(A[i,] > 0);
    L_I[i] = length(idx_a);
    if(L_I[i] > 0){
      A[i,1:L_I[i]] = A[i,idx_a]
      Tau[i,1:L_I[i]] = Tau[i,idx_a];
      if(L_I[i] < L){
        A[i,(L_I[i]+1):L] = 0
        Tau[i,(L_I[i]+1):L] = 0
      }
    }
  }
  
  tmpJ = max(J_I)
  tmpL = max(L_I)
  Z = Z[,-((tmpJ+1):J)]
  T = T[,-((tmpJ+1):J)]
  U = U[,-((tmpJ+1):J)]
  B = B[,-((tmpJ+1):J),]
  A = A[,-((tmpL+1):L)]
  Tau = Tau[,-((tmpL+1):L)]
  J = tmpJ
  L = tmpL
  
  if(trajectory_plot){
    if(I < 20){
      iids = 1:10
    }else{
      iids = sample(1:I,20)
    }
    for(i in iids){
      plot(x = T[i,1:J_I[i]]/60, y = Z[i,1:J_I[i]], type = "o", xlab = "hours", ylab = "outcome",
           ylim = c(min(Z[i,1:J_I[i]], U[i,1:J_I[i]]), max(Z[i,1:J_I[i]], U[i,1:J_I[i]])))
      lines(x = T[i,1:J_I[i]]/60, y = U[i,1:J_I[i]], col = "red", lty = 2)
      abline(v = Tau[i,1:L_I[i]]/60, col = A[i,1:L_I[i]]+3, lty = 2)
    }
  }
  
  output=list("data_Z" = Z, "data_T" = T, "data_B" = B, "data_Tau" = Tau, "data_A" = A, "data_I" = I, "data_J" = J, "data_L" = L, "data_D" = D, "data_P" = P, "data_J_I" = J_I, "data_L_I" = L_I, "data_avg_output" = 0, "data_sd_output" = 1, "param_a1" = a1, "param_a2" = a2, "param_a3" = a3, "param_r" = r, "param_b" = b, "param_C_u" = C_u, "param_C_f" = C_f, "param_Mu_b" = Mu_b, "param_Beta" = Beta, "param_mu_log_sd2_u" = mu_log_sd2_u, "param_mu_logit_rho_u" = mu_logit_rho_u, "param_sigma2_u" = sigma2_u, "param_rho_u" = rho_u, "param_u1_f" = u1_f, "param_u2_f" = u2_f, "param_u3_f" = u3_f, "param_ur_f" = ur_f, "param_ub_f" = ub_f, "param_sigma2_f" = sigma2_f, "param_rho_f" = rho_f, "param_Sigma_e" = Sigma_e, "param_U" = U, "param_maxJ" = maxJ)
  
  return(output)
}

logit <- function(p){
  return(log(p/(1-p)));
}

invlogit <- function(x){
  return(exp(x)/(1+exp(x)));
}