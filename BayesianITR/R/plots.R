require(tm)
require(R.utils)
require(ggplot2)
require(MCMCpack)
require(plyr)
require(SDMTools)
require(colorspace)

#' A Plot Function
#'
#' This function plots the heterogenous baseline progression and treatment response curves generated in the simulated data
#' 
#' @param sim_data The simulated data
#' @param W The treatment effective window
#' @author Yanbo Xu
#' @examples 
#' plot_sim_data(sim_data, W)
#' @export
plot_sim_data <- function(sim_data, W){
  treatment_ylim = rbind(c(0,10),c(-10,0))
  
  trainIndicator_D = array(0,dim=c(sim_data$data_I, sim_data$data_D))
  for(i in 1:sim_data$data_I){
    if(sim_data$data_L_I[i] > 0){
      treatment_idx = unique(sim_data$data_A[i,1:sim_data$data_L_I[i]])
      if(length(treatment_idx) > 0)
        trainIndicator_D[i, treatment_idx] = 1;
    }
  }
  
  # plot normal dynamics
  new_cluster = order(sim_data$param_Mu_b[,1])
  param.Mu_b = sim_data$param_Mu_b[new_cluster,]
  K_u = length(unique(sim_data$param_C_u))
  idx1 = which(sim_data$param_C_u == 1)
  idx2 = which(sim_data$param_C_u == 2)
  idx3 = which(sim_data$param_C_u == 3)
  param.C_u = sim_data$param_C_u
  param.C_u[idx1] = which(new_cluster == 1)
  param.C_u[idx2] = which(new_cluster == 2)
  param.C_u[idx3] = which(new_cluster == 3)
  
  end_T = max(sim_data$data_T);
  eval_T = 1:sim_data$data_J*(end_T/sim_data$data_J)
  eval_B = array(0,dim=c(sim_data$data_J,sim_data$data_P));
  eval_B[,1] = 1;
  eval_B[,2] = eval_T/(sim_data$param_maxJ*5);
  eval_B[,3] = eval_B[,2]^2;
  
  normal_dynamics = array(0,dim=c(sim_data$data_I,sim_data$data_J));
  for(i in 1:sim_data$data_I){
    normal_dynamics[i,1:sim_data$data_J_I[i]] = sim_data$param_U[i,1:sim_data$data_J_I[i]]
  }
  
  plot(1,type = "n",xlim = c(1,end_T/60), 
       ylim = c(min(normal_dynamics), max(normal_dynamics)),
       xlab = "hours", ylab = "Outcome",cex.axis = 2, main = "The heterogenous baseline progressions in data")
  for(i in 1:sim_data$data_I){
    lines(x = sim_data$data_T[i,1:sim_data$data_J_I[i]]/60, y = normal_dynamics[i,1:sim_data$data_J_I[i]],col=param.C_u[i] + 3,lwd = 0.5, lty = 2);
  }
  
  for(k in 1:K_u){
    ids = which(param.C_u == k)
    k_normal_dynamics = normal_dynamics[ids,]
    mean_normal_dynamics = eval_B %*% param.Mu_b[k,];
    mean_KK = array(0,dim = c(sim_data$data_J,sim_data$data_J))
    for(j in 1:sim_data$data_J)
      for(jj in 1:sim_data$data_J){
        mean_KK[j,jj] = exp(sim_data$param_mu_log_sd2_u[k]) * (invlogit(sim_data$param_mu_logit_rho_u[k]))^(abs(eval_T[j] - eval_T[jj])/60)
      }
    mean_normal_dynamics = MASS::mvrnorm(1,mean_normal_dynamics,mean_KK)
    lines(x = eval_T/60, y = mean_normal_dynamics,lwd = 2, col = k + 3)
  }
  
  # plot treatment level plots
  end_T = W;
  len_T = end_T + 1;
  eval_funcs = array(0, dim = c(sim_data$data_I,sim_data$data_D,max(len_T)))
  for(i in 1:sim_data$data_I){
    for(d in 1:sim_data$data_D){
      eval_T = 0:end_T[d];
      eval_funcs[i,d,eval_T+1] = g_function(eval_T,sim_data$param_a1[i,d],
                                            sim_data$param_a2[i,d],
                                            sim_data$param_a3[i,d],
                                            sim_data$param_r[i,d],
                                            sim_data$param_b[i,d])
    }  
  }
  
  param.u1_f = sim_data$param_u1_f
  param.u2_f = sim_data$param_u2_f
  param.u3_f = sim_data$param_u3_f
  param.ur_f = sim_data$param_ur_f
  param.ub_f = sim_data$param_ub_f
  param.C_f = sim_data$param_C_f
  
  for(d in 1:sim_data$data_D){
    new_cluster = order(param.u1_f[,d])
    param.u1_f[,d] = param.u1_f[new_cluster,d]
    param.u2_f[,d] = param.u2_f[new_cluster,d]
    param.u3_f[,d] = param.u3_f[new_cluster,d]
    param.ur_f[,d] = param.ur_f[new_cluster,d]
    param.ub_f[,d] = param.ub_f[new_cluster,d]
    
    idx1 = which(param.C_f[,d] == 1)
    idx2 = which(param.C_f[,d] == 2)
    idx3 = which(param.C_f[,d] == 3)
    param.C_f[idx1,d] = which(new_cluster == 1)
    param.C_f[idx2,d] = which(new_cluster == 2)
    param.C_f[idx3,d] = which(new_cluster == 3)
    
    eval_T = 0:end_T[d];
    idx_d = which(trainIndicator_D[,d] == 1)
    mean_eval_func = colMeans(eval_funcs[idx_d,d,eval_T+1]);
    
    y_min = min(eval_funcs[idx_d,d,eval_T+1]);
    y_max = max(eval_funcs[idx_d,d,eval_T+1]);
    
    ylab_str = paste("Response to treatment",d)
    plot(x = eval_T/60, y = mean_eval_func,lwd = 0.1,type="n",
         ylab = ylab_str, xlab = "hours",
         axes = TRUE, ylim = treatment_ylim[d,],cex.axis = 2, main = "The heterogenous treatment response curves in data")
    
    for(i in idx_d){
      lines(x = eval_T/60, y = eval_funcs[i,d,eval_T+1], col = param.C_f[i,d] + 3, lwd = 0.5, lty = 2)
    }
    ks = unique(param.C_f[idx_d,d]);
    for(k in ks){
      a1 = param.u1_f[k,d];
      a2 = invlogit(param.u2_f[k,d])
      a3 = invlogit(param.u3_f[k,d])
      r = param.ur_f[k,d]
      a = a1 * (exp(a2*r/2) - 1)/(exp(a2*r/2) + 1);
      b = invlogit(param.ub_f[k,d]) * a;
      k_eval_funcs = g_function(eval_T,a1,a2,a3,r,b)
      
      lines(x = eval_T/60, y = k_eval_funcs,lwd = 3,col = k+3)
    }
    box()
  } 
}

#' A Plot Function
#'
#' This function plots the estimated baseline progression and treatment response curves both at the individual level and subgroup level.
#' 
#' @param chain The model fit output from model_fit_dpmitr.
#' @param start_iter The thined iteration to start in chain for computing the MCMC estimation.
#' @param tol_iter The last thined iteration to end in chain for computing the MCMC estimation.
#' @param prior_unit The unit of the prior distribution; defaults to 1.
#' @param train_data The training data.
#' @param W The treatment effective window
#' @author Yanbo Xu
#' @examples 
#' plot_model_fit(chain, start_iter = 5, tol_iter = 9, prior_unit = 1, train_data = sim_data, W)
#' @export
plot_model_fit <-function(chain, start_iter, tol_iter, prior_unit, train_data, W){
  trainIndicator_D = array(0,dim=c(train_data$data_I, train_data$data_D))
  for(i in 1:train_data$data_I){
    if(train_data$data_L_I[i] > 0){
      treatment_idx = unique(train_data$data_A[i,1:train_data$data_L_I[i]])
      if(length(treatment_idx) > 0)
        trainIndicator_D[i, treatment_idx] = 1;
    }
  }
  
  model_fit = get_estimator_newT_C(chain,start_iter, tol_iter, prior_unit, train_data)
  
  end_T = max(train_data$data_T)
  eval_T = 1:train_data$data_J*(end_T/train_data$data_J)
  eval_B = array(0,dim=c(train_data$data_J,train_data$data_P));
  eval_B[,1] = 1;
  eval_B[,2] = eval_T/(24*60);
  eval_B[,3] = eval_B[,2]^2;
  
  # plot normal dynamics
  normal_dynamics = array(0,dim=c(train_data$data_I,train_data$data_J));
  for(i in 1:train_data$data_I){
    normal_dynamics[i,1:train_data$data_J_I[i]] = model_fit$Uz[i,1:train_data$data_J_I[i]]
  }
  
  plot(1,type = "n",xlim = c(1,end_T/60), 
       ylim = c(min(normal_dynamics), max(normal_dynamics)),
       xlab = "hours", ylab = "Outcome",cex.axis = 2, main = "The estimated baseline progressions")
  if(chain$model == 'ITR'){
    for(i in 1:train_data$data_I){
      lines(x = train_data$data_T[i,1:train_data$data_J_I[i]]/60, 
            y = normal_dynamics[i,1:train_data$data_J_I[i]],
            lwd = 0.5);
    }
  }else{
    new_cluster = order(model_fit$Mu_b[,1])
    model_fit$Mu_b = model_fit$Mu_b[new_cluster,]
    K_u = length(unique(model_fit$C_u))
    idx1 = which(model_fit$C_u == 1)
    idx2 = which(model_fit$C_u == 2)
    idx3 = which(model_fit$C_u == 3)
    model_fit$C_u[idx1] = which(new_cluster == 1)
    model_fit$C_u[idx2] = which(new_cluster == 2)
    model_fit$C_u[idx3] = which(new_cluster == 3)
    
    for(i in 1:train_data$data_I){
      lines(x = train_data$data_T[i,1:train_data$data_J_I[i]]/60, 
            y = normal_dynamics[i,1:train_data$data_J_I[i]],
            col=model_fit$C_u[i] + 3,lwd = 0.5, lty = 2);
    }
    
    for(k in 1:K_u){
      ids = which(model_fit$C_u == k)
      k_normal_dynamics = normal_dynamics[ids,]
      mean_normal_dynamics = eval_B %*% model_fit$Mu_b[k,];
      mean_KK = array(0,dim = c(train_data$data_J,train_data$data_J))
      for(j in 1:train_data$data_J)
        for(jj in 1:train_data$data_J){
          mean_KK[j,jj] = exp(model_fit$mu_log_sd2_u[k]) * (invlogit(model_fit$mu_logit_rho_u[k]))^(abs(eval_T[j] - eval_T[jj])/60)
        }
      mean_normal_dynamics = MASS::mvrnorm(1,mean_normal_dynamics,mean_KK)
      lines(x = eval_T/60, y = mean_normal_dynamics,lwd = 2, col = k + 3)
    }
  }
  
  # plot treatment level plots
  end_T = W;
  len_T = end_T + 1;
  eval_funcs = array(0, dim = c(train_data$data_I,train_data$data_D,max(len_T)))
  for(i in 1:train_data$data_I){
    for(d in 1:train_data$data_D){
      eval_T = 0:end_T[d];
      eval_funcs[i,d,eval_T+1] = g_function(eval_T,model_fit$a1[i,d],
                                            model_fit$a2[i,d],
                                            model_fit$a3[i,d],
                                            model_fit$r[i,d],
                                            model_fit$b[i,d])
    }  
  }
  
  for(d in 1:train_data$data_D){
    eval_T = 0:end_T[d];
    idx_d = which(trainIndicator_D[,d] == 1)
    mean_eval_func = colMeans(eval_funcs[idx_d,d,eval_T+1]);
    
    y_min = min(eval_funcs[idx_d,d,eval_T+1]);
    y_max = max(eval_funcs[idx_d,d,eval_T+1]);
    
    ylab_str = paste("Response to treatment",d)
    plot(x = eval_T/60, y = mean_eval_func,lwd = 0.1,type="l",
         ylab = ylab_str, xlab = "hours",
         axes = TRUE, ylim = c(y_min, y_max), cex.axis = 2, main = "The estimated treatment response curves")
    
    if(chain$model == 'ITR'){
      for(i in idx_d){
        lines(x = eval_T/60, y = eval_funcs[i,d,eval_T+1], lwd = 0.5)
      }
    }else{
      a1 = model_fit$u1_f[,d];
      a2 = invlogit(model_fit$u2_f[,d])
      a3 = invlogit(model_fit$u3_f[,d])
      r = model_fit$ur_f[,d]
      a = a1 * (exp(a2*r/2) - 1)/(exp(a2*r/2) + 1);
      new_cluster = order(a)
      model_fit$u1_f[,d] = model_fit$u1_f[new_cluster,d]
      model_fit$u2_f[,d] = model_fit$u2_f[new_cluster,d]
      model_fit$u3_f[,d] = model_fit$u3_f[new_cluster,d]
      model_fit$ur_f[,d] = model_fit$ur_f[new_cluster,d]
      model_fit$ub_f[,d] = model_fit$ub_f[new_cluster,d]
      
      idx1 = which(model_fit$C_f[,d] == 1)
      idx2 = which(model_fit$C_f[,d] == 2)
      idx3 = which(model_fit$C_f[,d] == 3)
      model_fit$C_f[idx1,d] = which(new_cluster == 1)
      model_fit$C_f[idx2,d] = which(new_cluster == 2)
      model_fit$C_f[idx3,d] = which(new_cluster == 3)
      
      for(i in idx_d){
        lines(x = eval_T/60, y = eval_funcs[i,d,eval_T+1], col = model_fit$C_f[i,d] + 3, lwd = 0.5, lty = 2)
      }
      ks = unique(model_fit$C_f[idx_d,d]);
      if(length(ks)){
        for(k in ks){
          a1 = model_fit$u1_f[k,d];
          a2 = invlogit(model_fit$u2_f[k,d])
          a3 = invlogit(model_fit$u3_f[k,d])
          r = model_fit$ur_f[k,d]
          a = a1 * (exp(a2*r/2) - 1)/(exp(a2*r/2) + 1);
          b = invlogit(model_fit$ub_f[k,d]) * a;
          k_eval_funcs = g_function(eval_T,a1,a2,a3,r,b)
          lines(x = eval_T/60, y = k_eval_funcs,lwd = 3,col = k+3)
        }
      }
    }
    box()
  } 
}

get_estimator_newT_C <- function(chain,start_iter, tol_iter, prior_unit, data){
  I = data$data_I
  J = data$data_J
  J_I = data$data_J_I
  L_I = data$data_L_I
  D = data$data_D
  T = data$data_T
  Tau = data$data_Tau;
  A = data$data_A;
  Z = data$data_Z;
  B = data$data_B;
  P = data$data_P;
  avg_output = data$data_avg_output;
  sd_output = data$data_sd_output;
  
  Uz = array(0,dim=c(I,J));
  U = array(0,dim=c(I,J));
  max_Uz = array(0,dim=c(I,J));
  min_Uz = array(0,dim=c(I,J));
  Z0 = array(0,dim=c(I,J));
  for(i in 1:I)
    for(j in 1:J_I[i]){
      #Z0[i,j] = mean(chain$Z0[(tol_iter/2):tol_iter,i,j]);
      tmp = rep(0,tol_iter - start_iter +  1)
      for(iter in start_iter:tol_iter)
        tmp[iter - start_iter + 1] = B[i,j,] %*% chain$Beta[, i, iter]
      Z0[i,j] = mean(tmp);
      U[i,j] = mean(chain$U[j,i,start_iter:tol_iter]);
      Uz[i,j] = mean(chain$U[j,i,start_iter:tol_iter] + tmp);
      max_Uz[i,j] = max(chain$U[j,i,start_iter:tol_iter] + tmp);
      min_Uz[i,j] = min(chain$U[j,i,start_iter:tol_iter] + tmp);
    }
  
  
  F = array(0,dim=c(I,J));
  max_F = array(0,dim=c(I,J));
  min_F = array(0,dim=c(I,J));
  for(i in 1:I)
    for(j in 1:J){
      F[i,j] = mean(chain$f[j,i,start_iter:tol_iter]);
      max_F[i,j] = max(chain$f[j,i,start_iter:tol_iter]);
      min_F[i,j] = min(chain$f[j,i,start_iter:tol_iter]);
    }
  
  rho_f = rowMeans(chain$rho_f[,start_iter:tol_iter])
  sigma2_f = rowMeans(chain$sigma2_f[,start_iter:tol_iter])
  
  Sigma_e = rowMeans(chain$Sigma_e[,start_iter:tol_iter]);
  
  rho_u = rowMeans(chain$rho_u[,start_iter:tol_iter]);
  sigma2_u = rowMeans(chain$sigma2_u[,start_iter:tol_iter]);
  
  a1 = array(0,dim = c(I,D))
  a2 = array(0,dim = c(I,D))
  a3 = array(0,dim = c(I,D))
  r = array(0,dim = c(I,D))
  b = array(0,dim = c(I,D))
  
  for(d in 1:D){
    for(i in 1:I){
      a1[i,d] = chain$a1[d,i,tol_iter];
      a2[i,d] = chain$a2[d,i,tol_iter];
      a3[i,d] = chain$a3[d,i,tol_iter];
      r[i,d] = chain$r[d,i,tol_iter];
      b[i,d] = chain$b[d,i,tol_iter];
      
    }
  }
  
  Beta = array(0,dim=c(I,P))
  for(p in 1:P)
    for(i in 1:I)
      Beta[i,p] = mean(chain$Beta[p,i,start_iter:tol_iter])
    
  fit_Z = Uz + F;
  max_Z = max_Uz + max_F;
  min_Z = min_Uz + min_F;
  
  ci_Z = array(0, dim = c(I,J))
  RMSE = 0;
  n = 0;
  for(i in 1:I){
    RMSE = RMSE + sum((fit_Z[i,1:J_I[i]]-Z[i,1:J_I[i]])^2)
    n = n + J_I[i];
  }
  RMSE = sqrt(RMSE/n)
  
  output = list("Uz" = Uz, "F" = F, "fit_Z" = fit_Z, "max_Z" = max_Z, "min_Z" = min_Z, 
                "RMSE" = RMSE * sd_output, "Beta" = Beta, "Z0" = Z0, "U" = U,
                "rho_f" = rho_f, "sigma2_f" = sigma2_f, "Sigma_e" = Sigma_e,
                "rho_u" = rho_u, "sigma2_u" = sigma2_u, "a1" = a1, "a2" = a2, "a3" = a3, "r" = r,"b" = b)
  
  if(chain$model != 'ITR'){
    C_u = (chain$C_u[,tol_iter])+1;
    
    C_f = array(0, dim = c(I,D));
    Kf = dim(chain$u1_f[,,tol_iter])[2];
    u1_f = array(0, dim = c(Kf,D));
    u2_f = array(0, dim = c(Kf,D));
    u3_f = array(0, dim = c(Kf,D));
    ub_f = array(0, dim = c(Kf,D));
    ur_f = array(0, dim = c(Kf,D));
    
    for(d in 1:D){
      for(i in 1:I)
        C_f[i,d] = (chain$C_f[d,i,tol_iter]) + 1;
      
      for(k in 1:Kf){
        u1_f[k,d] = chain$u1_f[d,k,tol_iter]
        u2_f[k,d] = chain$u2_f[d,k,tol_iter]
        u3_f[k,d] = chain$u3_f[d,k,tol_iter]
        ur_f[k,d] = chain$ur_f[d,k,tol_iter]
        ub_f[k,d] = chain$ub_f[d,k,tol_iter]
      }
    }
    mu_log_sd2_u = chain$mu_log_sd2_u[,tol_iter]
    mu_logit_rho_u = chain$mu_logit_rho_u[,tol_iter]
    
    Ku = length(mu_log_sd2_u)
    
    Mu_b = array(0,dim=c(Ku,P))
    
    for(p in 1:P){
      for(k in 1:Ku)
        Mu_b[k,p] = chain$Mu_b[p,k,tol_iter]
    }
    output = c(output, list("Mu_b" = Mu_b, "mu_log_sd2_u" = mu_log_sd2_u, "mu_logit_rho_u" = mu_logit_rho_u,
                            "C_u" = C_u, "C_f" = C_f,
                            "u1_f" = u1_f, "u2_f" = u2_f, "u3_f" = u3_f,"ur_f" = ur_f,"ub_f" = ub_f))
  }
  return(output)
} 

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

g_function <- function(t,a1,a2,a3,r,b){
  y = array(0,length(t))
  b0 = -a1/(1+exp(a2*r/2))
  a0 = (a1+2*b0-b) * (1+exp(-a3*r/2))
  for(i in 1:length(t)){
    if(t[i] < r){
      y[i] = a1/(1+exp(-a2*(t[i]-r/2))) + b0
    }else{
      y[i] = a0/(1+exp(a3*(t[i] - 3*r/2))) + b
    }
  }
  return(y);
}

logit <- function(p){
  return(log(p/(1-p)));
}

invlogit <- function(x){
  return(exp(x)/(1+exp(x)));
}