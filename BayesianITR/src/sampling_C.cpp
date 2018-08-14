#ifndef SAMPLING_UPDATE_H
#define SAMPLING_UPDATE_H

#include "sampling_utlities.h"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// ********* Declaration ***********//
// update functions
void Update_f(int, int, int, int, int, int,
              double, double, NumericVector,
              NumericMatrix, NumericMatrix, NumericMatrix,
              NumericMatrix, NumericMatrix,
              NumericMatrix, NumericMatrix,
              double, NumericVector, double,
              double, NumericVector, double,
              double, NumericVector, double,
              double, NumericVector, double,
              double, NumericVector, double,
              double, double,
              NumericVector, NumericVector,
              NumericVector, NumericVector,
              NumericMatrix &, NumericMatrix &, NumericMatrix &,
              NumericMatrix &, NumericMatrix &,
              NumericMatrix &,
              double &, double &,
              Function,
              Function); // Gibbs to update f, MH to update a1,a2,a3,r,b

void Update_C_f(int, int,
                NumericMatrix, NumericMatrix, NumericMatrix,
                NumericMatrix, NumericMatrix,
                double, double, double, double, double,
                NumericMatrix, NumericMatrix, NumericMatrix,
                NumericMatrix, NumericMatrix,
                IntegerVector, NumericVector,
                IntegerMatrix &, IntegerMatrix &,
                NumericMatrix &, NumericMatrix &, NumericMatrix &,
                NumericMatrix &, NumericMatrix &); // Gibbs to update C_f

void Update_Beta(int, int, int, int, int,
                 NumericMatrix, NumericMatrix,
                 double, double,
                 NumericMatrix, NumericMatrix,
                 IntegerVector, NumericMatrix, NumericMatrix,
                 NumericMatrix &,
                 Function); // Gibbs to update Beta

void Update_ITR_Beta(int, int, int, int, int,
                     NumericVector, NumericMatrix,
                     NumericMatrix, NumericMatrix,
                     double, double,
                     NumericMatrix, NumericMatrix,
                     NumericMatrix &,
                     Function); // Gibbs to update Beta in ITR

void Update_U(int, int, int, int, int,
              double, double, NumericVector,
              NumericMatrix, NumericMatrix, NumericMatrix, NumericMatrix,
              double, double, NumericMatrix,
              double, double, double,
              double, double, double,
              NumericMatrix, 
              NumericVector &, NumericVector &, NumericMatrix &,
              Function, Function); // Gibbs to update u, MH to update sigma2_u, rho_u

void Update_C_u(int, int, int, int,
                NumericMatrix, NumericMatrix, NumericMatrix,
                NumericVector, NumericVector,
                double, double,
                NumericVector, NumericVector,
                IntegerVector &, IntegerVector &,
                NumericMatrix &, NumericMatrix &,
                NumericVector &, NumericVector &,
                Function); // Gibbs to update C_u

void Update_e(int, int, int, int,
              NumericMatrix, NumericMatrix,
              NumericMatrix,
              NumericMatrix, NumericMatrix,
              NumericVector, NumericVector,
              NumericVector &,
              Function); // Gibbs to update e

void Update_f_param(int,
                    double, double, double, double,
                    NumericVector, NumericVector,
                    double, double,
                    NumericVector &, NumericVector &); // MH to update sigma2_f, rho_f

void Update_mu_u(int, int,
                 double, double, NumericVector, NumericMatrix,
                 IntegerVector,
                 double, double, double,
                 double, double, double,
                 NumericMatrix, NumericMatrix,
                 NumericVector, NumericVector,
                 NumericVector &, NumericVector &,
                 NumericMatrix &, NumericMatrix &,
                 Function,
                 Function); // Gibbs to update mu_logit_rho_u, mu_log_sd2_u

void Update_mu_f(int,
                 IntegerVector,
                 double, NumericVector, double,
                 double, NumericVector, double,
                 double, NumericVector, double,
                 double, NumericVector, double,
                 double, NumericVector, double,
                 IntegerMatrix,
                 NumericMatrix, NumericMatrix, NumericMatrix,
                 NumericMatrix, NumericMatrix,
                 NumericMatrix &, NumericMatrix &, NumericMatrix &,
                 NumericMatrix &, NumericMatrix &); // Gibbs to update u1_f, u2_f, u3_f, ur_f, ub_f

void Update_M1(int,
               double, double,
               int,
               IntegerVector,
               double &); // Gibbs to update M1

void Update_M2(int, int,
               double, double,
               IntegerVector,
               IntegerMatrix,
               NumericVector &); // Gibbs to update M2

//*********** Defination ****************//
// Gibbs to update f
// MH to update a1,a2,a3,r,b
// Prepare MH to update sigma2_f, rho_f
void Update_f(int i, int Ji, int Li, int dimD, int dimP, int dimJ,
              double mh_sd, double prior_unit, NumericVector W,
              NumericMatrix dataZ, NumericMatrix dataT, NumericMatrix dataB_3D,
              NumericMatrix dataTau, NumericMatrix dataA,
              NumericMatrix Beta, NumericMatrix U,
              double sd1_f, NumericVector u10_f, double sd10_f,
              double sd2_f, NumericVector u20_f, double sd20_f,
              double sd3_f, NumericVector u30_f, double sd30_f,
              double sdr_f, NumericVector ur0_f, double sdr0_f,
              double sdb_f, NumericVector ub0_f, double sdb0_f,
              double Sigma_e, double invSigma_e,
              NumericVector rho_f, NumericVector sigma2_f,
              NumericVector b_rho_f, NumericVector b_sd2_f,
              NumericMatrix & a1, NumericMatrix & a2, NumericMatrix & a3,
              NumericMatrix & r, NumericMatrix & b,
              NumericMatrix & f,
              double & Pa_f, double & Pb_f,
              Function rmvn,
              Function dmvn){
  NumericVector Yi(Ji);
  NumericMatrix KK_f(Ji, Ji);
  NumericVector mu_f(Ji);
  
  NumericMatrix b_KK_f(Ji, Ji); // later for updating sigma2_f, rho_f
  NumericVector b_mu_f(Ji); // for updating a1-a3,r,b
  NumericVector b_a1(dimD);
  NumericVector b_a2(dimD);
  NumericVector b_a3(dimD);
  NumericVector b_r(dimD);
  NumericVector b_b(dimD);
  
  double Pa = 0.0;
  double Pb = 0.0;
  for(int j = 0; j < Ji; j++){
    double BBij = 0.0;
    for(int p=0; p < dimP; p++){
      BBij += dataB_3D(p,i*dimJ + j) * Beta(i,p);
    }
    Yi(j) = dataZ(i,j) - U(i,j) - BBij;
  }
  
  for(int d = 0; d < dimD; d++){
    //Jacobian adjustment
    double g_r = a1(i,d) * (exp(a2(i,d)*r(i,d)/2.0) - 1.0)/(1.0+exp(a2(i,d)*r(i,d)/2.0));
    double zeta = std::abs(g_r/(b(i,d) * (g_r - b(i,d))));
    double b_g_r = 0;
    double b_zeta = 0;
    
    Pa += R::dnorm(a1(i,d),u10_f(d),sqrt(sd1_f + sd10_f),true);
    Pa += R::dnorm(logit(a2(i,d)),u20_f(d),sqrt(sd2_f + sd20_f),true);
    Pa += R::dnorm(logit(a3(i,d)),u30_f(d),sqrt(sd3_f + sd30_f),true);
    Pa += R::dnorm(r(i,d),ur0_f(d),sqrt(sdr_f + sdr0_f),true);
    Pa += R::dnorm(logit(b(i,d)/g_r),ub0_f(d),sqrt(sdb_f + sdb0_f),true);
    Pa += log(zeta);
    Pa = Pa - log(a2(i,d)) - log(1.0-a2(i,d)) - log(a3(i,d)) - log(1.0-a3(i,d));
    
    b_a1(d) = propose_norm(a1(i,d), mh_sd);
    b_a2(d) = propose_prob(a2(i,d),mh_sd/2.0);
    b_a3(d) = propose_prob(a3(i,d),mh_sd/2.0);
    b_r(d) = propose_norm(r(i,d), mh_sd);
    b_g_r = b_a1(d)*(exp(b_a2(d)*b_r(d)/2.0) - 1.0)/(1.0+exp(b_a2(d)*b_r(d)/2.0));
    b_b(d) = propose_prob(b(i,d)/g_r, mh_sd/2.0) * b_g_r;
    
    // Jacobian adjustment
    b_zeta = std::abs(b_g_r/(b_b(d) * (b_g_r - b_b(d))));
    
    Pb += R::dnorm(b_a1(d),u10_f(d),sqrt(sd1_f + sd10_f),true);
    Pb += R::dnorm(logit(b_a2(d)),u20_f(d),sqrt(sd2_f + sd20_f),true);
    Pb += R::dnorm(logit(b_a3(d)),u30_f(d),sqrt(sd3_f + sd30_f),true);
    Pb += R::dnorm(b_r(d),ur0_f(d),sqrt(sdr_f + sdr0_f),true);
    Pb += R::dnorm(logit(b_b(d)/b_g_r),ub0_f(d),sqrt(sdb_f + sdb0_f),true);
    Pb += log(b_zeta);
    Pb = Pb - log(b_a2(d)) - log(1.0-b_a2(d)) - log(b_a3(d)) - log(1.0-b_a3(d));
  }
  LogicalVector idx_f(Ji, false);
  
  for(int l = 0; l < Li; l++){
    int d = dataA(i,l);
    for(int j = 0; j < Ji; j++){
      double deltaT1 = (dataT(i,j) - dataTau(i,l))/prior_unit;
      if(deltaT1 <= 0)
        continue;
      if(deltaT1 > W(d))
        break;
      idx_f(j) = true;
      double b0 = -a1(i,d)/(1+exp(a2(i,d)*r(i,d)/2.0));
      double a0 = (a1(i,d)+2.0*b0-b(i,d))/(1.0+exp(-a3(i,d)*r(i,d)/2.0));
      if(deltaT1 <= r(i,d)){
        mu_f(j) = mu_f(j) + b0 + a1(i,d)/(1.0+exp(-a2(i,d)*(deltaT1 - r(i,d)/2.0)));
      }else{
        mu_f(j) = mu_f(j) + b(i,d) + a0/(1.0+exp(a3(i,d)*(deltaT1-3.0*r(i,d)/2.0)));
      }
      double b_b0 = -b_a1(d)/(1.0+exp(b_a2(d)*b_r(d)/2.0));
      double b_a0 = (b_a1(d)+2.0*b_b0-b_b(d))/(1.0+exp(-b_a3(d)*b_r(d)/2.0));
      if(deltaT1 <= b_r(d)){
        b_mu_f(j) = b_mu_f(j) + b_b0 + b_a1(d)/(1.0+exp(-b_a2(d)*(deltaT1 - b_r(d)/2.0)));
      }else{
        b_mu_f(j) = b_mu_f(j) + b_b(d) + b_a0/(1.0+exp(b_a3(d)*(deltaT1 - 3.0*b_r(d)/2.0)));
      }
      KK_f(j,j) += sigma2_f(d);
      b_KK_f(j,j) += b_sd2_f(d);
      
      for(int jj = j+1; jj < Ji; jj++){
        double deltaT2 = (dataT(i,jj) - dataTau(i,l))/prior_unit;
        if(deltaT2 > W(d))
          break;
        KK_f(j,jj) = KK_f(j,jj) + sigma2_f(d) * std::pow(rho_f(d),(dataT(i,jj) - dataT(i,j))/60.0);
        KK_f(jj,j) = KK_f(j,jj);
        
        b_KK_f(j,jj) = b_KK_f(j,jj) + b_sd2_f(d) * std::pow(b_rho_f(d),(dataT(i,jj) - dataT(i,j))/60.0);
        b_KK_f(jj,j) = b_KK_f(j,jj);
      }
    }
  }
  int cnt_idx_f = sum(idx_f);
  
  if(cnt_idx_f == 1){
    int uni_idx_f = 0;
    for(int j = 0; j < Ji; j++){
      if(idx_f(j)){
        uni_idx_f = j;
        break;
      }
    }
    
    double newYi = Yi(uni_idx_f);
    double newMu_f = mu_f(uni_idx_f);
    double newKK_f = KK_f(uni_idx_f, uni_idx_f);
    double b_newMu_f = b_mu_f(uni_idx_f);
    double b_newKK_f = b_KK_f(uni_idx_f, uni_idx_f);
    
    // Update f
    double invKK_f = 1/newKK_f;
    double hat_Sigma_f = 1/(invKK_f + invSigma_e);
    double hat_Mu_f = hat_Sigma_f * (invSigma_e * newYi + invKK_f * newMu_f);
    f(i,uni_idx_f) = rnorm(1, hat_Mu_f, sqrt(hat_Sigma_f))(0);
    
    // Update a1,a2,a3,r,b
    Pa += R::dnorm(newYi,newMu_f, sqrt(Sigma_e + newKK_f),true);
    Pb += R::dnorm(newYi,b_newMu_f, sqrt(Sigma_e + newKK_f),true);
    double Paccept = exp(Pb - Pa);
    if(runif(1,0,1)(0) < Paccept){
      a1(i,_) = b_a1;
      a2(i,_) = b_a2;
      a3(i,_) = b_a3;
      r(i,_) = b_r;
      b(i,_) = b_b;
      Pa_f += R::dnorm(newYi, b_newMu_f, sqrt(Sigma_e + newKK_f), true);
      Pb_f += R::dnorm(newYi, b_newMu_f, sqrt(Sigma_e + b_newKK_f), true);
    }else{
      Pa_f += R::dnorm(newYi, newMu_f, sqrt(Sigma_e + newKK_f), true);
      Pb_f += R::dnorm(newYi, newMu_f, sqrt(Sigma_e + b_newKK_f), true);
    }
  }
  if(cnt_idx_f > 1){
    arma::vec newYi = arma::zeros(cnt_idx_f);
    arma::vec newMu_f = arma::zeros(cnt_idx_f);
    arma::mat newKK_f = arma::zeros(cnt_idx_f,cnt_idx_f);
    arma::vec b_newMu_f = arma::zeros(cnt_idx_f);
    arma::mat b_newKK_f = arma::zeros(cnt_idx_f,cnt_idx_f);
    
    int idx_j = 0;
    int idx_jj = 0;
    for(int j = 0; j < Ji; j++){
      if(idx_f(j)){
        newYi(idx_j) = Yi(j);
        newMu_f(idx_j) = mu_f(j);
        b_newMu_f(idx_j) = b_mu_f(j);
        newKK_f(idx_j, idx_j) = KK_f(j,j);
        b_newKK_f(idx_j, idx_j) = b_KK_f(j,j);
        idx_jj = idx_j + 1;
        for(int jj = j + 1; jj < Ji; jj++){
          if(idx_f(jj)){
            newKK_f(idx_j, idx_jj) = KK_f(j,jj);
            newKK_f(idx_jj, idx_j) = newKK_f(idx_j, idx_jj);
            
            b_newKK_f(idx_j, idx_jj) = b_KK_f(j,jj);
            b_newKK_f(idx_jj, idx_j) = b_newKK_f(idx_j, idx_jj);
            idx_jj++;
          }
        }
        idx_j++;
      }
    }
    
    arma::mat eye_mat = arma::eye(cnt_idx_f, cnt_idx_f);
    arma::mat invKK_f = arma::inv(newKK_f);
    arma::mat hat_Sigma_f = arma::inv(invKK_f + eye_mat * invSigma_e);
    arma::vec hat_Mu_f = hat_Sigma_f * (eye_mat * invSigma_e * newYi + invKK_f * newMu_f);
    NumericVector tmp_f = rmvn(1, hat_Mu_f, hat_Sigma_f);
    idx_j = 0;
    for(int j = 0; j < Ji; j++){
      if(idx_f(j)){
        f(i,j) = tmp_f(idx_j);
        idx_j++;
      }
    }
    
    
    // Update a1,a2,a3,r,b
    Pa += as<double>(dmvn(newYi.t(), newMu_f.t(), eye_mat * Sigma_e + newKK_f,true));
    Pb += as<double>(dmvn(newYi.t(), b_newMu_f.t(), eye_mat * Sigma_e + newKK_f ,true));
    double Paccept = exp(Pb - Pa);
    
    if(runif(1,0,1)(0) < Paccept){
      a1(i,_) = b_a1;
      a2(i,_) = b_a2;
      a3(i,_) = b_a3;
      r(i,_) = b_r;
      b(i,_) = b_b;
      
      Pa_f += as<double>(dmvn(newYi.t(), b_newMu_f.t(), eye_mat * Sigma_e + newKK_f,true));
      Pa_f += as<double>(dmvn(newYi.t(), b_newMu_f.t(), eye_mat * Sigma_e + b_newKK_f,true));
    }else{
      Pa_f += as<double>(dmvn(newYi.t(), newMu_f.t(), eye_mat * Sigma_e + newKK_f,true));
      Pa_f += as<double>(dmvn(newYi.t(), newMu_f.t(), eye_mat * Sigma_e + b_newKK_f,true));
    }
  }
}

// Gibbs to update f
// pdate a1,a2,a3,r,b
// Prepare MH to update sigma2_f, rho_f, u1_f, ..., ub_f
void Update_DP_f(int i, int Ji, int Li, int dimD, int dimP, int dimJ, int K,
              double mh_sd, double prior_unit, NumericVector W,
              NumericMatrix dataZ, NumericMatrix dataT, NumericMatrix dataB_3D,
              NumericMatrix dataTau, NumericMatrix dataA,
              NumericMatrix Beta, NumericMatrix U,
              IntegerVector K_f, NumericVector M2,
              NumericVector u10_f, double sd10_f,
              NumericVector u20_f, double sd20_f,
              NumericVector u30_f, double sd30_f,
              NumericVector ur0_f, double sdr0_f,
              NumericVector ub0_f, double sdb0_f,
              double Sigma_e, double invSigma_e,
              NumericVector rho_f, NumericVector sigma2_f,
              NumericVector b_rho_f, NumericVector b_sd2_f,
              NumericMatrix u1_f, NumericMatrix u2_f, NumericMatrix u3_f,
              NumericMatrix ur_f, NumericMatrix ub_f,
              NumericMatrix b_u1_f, NumericMatrix b_u2_f, NumericMatrix b_u3_f,
              NumericMatrix b_ur_f, NumericMatrix b_ub_f,
              IntegerMatrix & N_f, IntegerMatrix & C_f,
              NumericMatrix & a1, NumericMatrix & a2, NumericMatrix & a3,
              NumericMatrix & r, NumericMatrix & b,
              NumericMatrix & f,
              double & Pa_uf, double & Pb_uf, double & Pa1_f, double & Pb1_f, double & Pa2_f, double & Pb2_f,
              Function rmvn,
              Function dmvn){
  NumericVector Yi(Ji);
  NumericMatrix KK_f(Ji, Ji);
  NumericVector mu_f(Ji);
  
  NumericMatrix b_KK_f(Ji, Ji); // later for updating sigma2_f, rho_f
  NumericVector b_mu_f(Ji); // later for updating u1-u3_f,ur_f,ub_f

  NumericMatrix mu_d_f(dimD, Ji); // later for updating C_f
  NumericMatrix new_mu_d_k_f(dimD, K*Ji); // later for upating C_f
  
  // later for updating u1_f,...,ub_f
  NumericVector b_a1(dimD);
  NumericVector b_a2(dimD);
  NumericVector b_a3(dimD);
  NumericVector b_r(dimD);
  NumericVector b_b(dimD);
  
  for(int d = 0; d < dimD; d++){
    // Jacobian adjustment
    double g_r = a1(i,d) * (exp(a2(i,d)*r(i,d)/2) - 1)/(1+exp(a2(i,d)*r(i,d)/2));
    double zeta = abs(g_r/(b(i,d) * (g_r - b(i,d))));
    
    Pa_uf += R::dnorm(u1_f(C_f(i,d),d),u10_f(d), sqrt(sd10_f),true);
    Pa_uf += R::dnorm(u2_f(C_f(i,d),d),u20_f(d), sqrt(sd20_f),true);
    Pa_uf += R::dnorm(u3_f(C_f(i,d),d),u30_f(d), sqrt(sd30_f),true);
    Pa_uf += R::dnorm(ur_f(C_f(i,d),d),ur0_f(d), sqrt(sdr0_f),true);
    Pa_uf += R::dnorm(ub_f(C_f(i,d),d),ub0_f(d), sqrt(sdb0_f),true);
    Pa_uf += log(zeta);
    Pa_uf -= log(a2(i,d)) - log(1-a2(i,d)) - log(a3(i,d)) - log(1-a3(i,d));
    
    b_a1(d) = b_u1_f(C_f(i,d),d);
    b_a2(d) = invlogit(b_u2_f(C_f(i,d),d));
    b_a3(d) = invlogit(b_u3_f(C_f(i,d),d));
    b_r(d) = b_ur_f(C_f(i,d),d);
    double b_g_r = b_a1(d)*(exp(b_a2(d)*b_r(d)/2) - 1)/(1+exp(b_a2(d)*b_r(d)/2));
    b_b(d) = invlogit(b_ub_f(C_f(i,d),d))*b_g_r;
    
    // Jacobian adjustment
    double b_zeta = abs(b_g_r/(b_b(d) * (b_g_r - b_b(d))));
    
    Pb_uf += R::dnorm(b_u1_f(C_f(i,d),d),u10_f(d), sqrt(sd10_f),true);
    Pb_uf += R::dnorm(b_u2_f(C_f(i,d),d),u20_f(d), sqrt(sd20_f),true);
    Pb_uf += R::dnorm(b_u3_f(C_f(i,d),d),u30_f(d), sqrt(sd30_f),true);
    Pb_uf += R::dnorm(b_ur_f(C_f(i,d),d),ur0_f(d), sqrt(sdr0_f),true);
    Pb_uf += R::dnorm(b_ub_f(C_f(i,d),d),ub0_f(d), sqrt(sdb0_f),true);
    Pb_uf += log(b_zeta);
    Pb_uf -= log(b_a2(d)) - log(1-b_a2(d)) - log(b_a3(d)) - log(1-b_a3(d));
  }

  for(int j = 0; j < Ji; j++){
    double BBij = 0.0;
    for(int p=0; p < dimP; p++){
      BBij += dataB_3D(p,i*dimJ + j) * Beta(i,p);
    }
    Yi(j) = dataZ(i,j) - U(i,j) - BBij;
  }
  
  LogicalVector idx_f(Ji, false);
  
  for(int l = 0; l < Li; l++){
    int d = dataA(i,l);
    for(int j = 0; j < Ji; j++){
      
      double deltaT1 = (dataT(i,j) - dataTau(i,l))/prior_unit;
      if(deltaT1 <= 0)
        continue;
      if(deltaT1 > W(d))
        break;
      idx_f(j) = true;
      
      double b0 = -a1(i,d)/(1+exp(a2(i,d)*r(i,d)/2.0));
      double a0 = (a1(i,d)+2.0*b0-b(i,d))/(1.0+exp(-a3(i,d)*r(i,d)/2.0));
      if(deltaT1 <= r(i,d)){
        mu_f(j) = mu_f(j) + b0 + a1(i,d)/(1.0+exp(-a2(i,d)*(deltaT1 - r(i,d)/2.0)));
        mu_d_f(d,j) = mu_d_f(d,j) + b0 + a1(i,d)/(1.0+exp(-a2(i,d)*(deltaT1 - r(i,d)/2.0)));
      }else{
        mu_f(j) = mu_f(j) + b(i,d) + a0/(1.0+exp(a3(i,d)*(deltaT1-3.0*r(i,d)/2.0)));
        mu_d_f(d,j) = mu_d_f(d,j) + b(i,d) + a0/(1.0+exp(a3(i,d)*(deltaT1-3.0*r(i,d)/2.0)));
      }
      
      double b_b0 = -b_a1(d)/(1+exp(b_a2(d)*b_r(d)/2));
      double b_a0 = (b_a1(d)+2*b_b0-b_b(d))/(1+exp(-b_a3(d)*b_r(d)/2));
      if(deltaT1 <= b_r(d)){
        b_mu_f(j) = b_mu_f(j) + b_b0 + b_a1(d)/(1+exp(-b_a2(d)*(deltaT1 - b_r(d)/2)));
      }else{
        b_mu_f(j) = b_mu_f(j) + b_b(d) + b_a0/(1+exp(b_a3(d)*(deltaT1 - 3*b_r(d)/2)));
      }
        
      for(int k = 0; k < K_f(d); k++){
        if(k == C_f(i,d)){
          new_mu_d_k_f(d,k*Ji + j) = mu_d_f(d,j);
          continue;
        }
        
        double new_a1 = u1_f(k,d);
        double new_a2 = invlogit(u2_f(k,d));
        double new_a3 = invlogit(u3_f(k,d));
        double new_r = ur_f(k,d);
        double new_g_r = new_a1 * (exp(new_a2*new_r/2) - 1)/(exp(new_a2*new_r/2) + 1);
        double new_b = invlogit(ub_f(k,d)) * new_g_r;
        
        double new_b0 = -new_a1/(1+exp(new_a2*new_r/2));
        double new_a0 = (new_a1+2*new_b0-new_b)/(1+exp(-new_a3*new_r/2));
          if(deltaT1 <= new_r){
            new_mu_d_k_f(d,k*Ji + j) = new_mu_d_k_f(d,k*Ji + j) + new_b0 + new_a1/(1+exp(-new_a2*(deltaT1 - new_r/2)));
          }else{
            new_mu_d_k_f(d,k*Ji + j) = new_mu_d_k_f(d,k*Ji + j) + new_b + new_a0/(1+exp(new_a3*(deltaT1-3*new_r/2)));
          }
      }
      
      KK_f(j,j) += sigma2_f(d);
      b_KK_f(j,j) += b_sd2_f(d);
      
      for(int jj = j+1; jj < Ji; jj++){
        double deltaT2 = (dataT(i,jj) - dataTau(i,l))/prior_unit;
        if(deltaT2 > W(d))
          break;
        KK_f(j,jj) = KK_f(j,jj) + sigma2_f(d) * std::pow(rho_f(d),(dataT(i,jj) - dataT(i,j))/60.0);
        KK_f(jj,j) = KK_f(j,jj);
        
        b_KK_f(j,jj) = b_KK_f(j,jj) + b_sd2_f(d) * std::pow(b_rho_f(d),(dataT(i,jj) - dataT(i,j))/60.0);
        b_KK_f(jj,j) = b_KK_f(j,jj);
      }
    }
  }
  int cnt_idx_f = sum(idx_f);
  
  if(cnt_idx_f == 1){
    int uni_idx_f = 0;
    for(int j = 0; j < Ji; j++){
      if(idx_f(j)){
        uni_idx_f = j;
        break;
      }
    }
    
    double newYi = Yi(uni_idx_f);
    double newMu_f = mu_f(uni_idx_f);
    double newKK_f = KK_f(uni_idx_f, uni_idx_f);
    double b_newMu_f = b_mu_f(uni_idx_f);
    double b_newKK_f = b_KK_f(uni_idx_f, uni_idx_f);
    
    // Update f
    double invKK_f = 1/newKK_f;
    double hat_Sigma_f = 1/(invKK_f + invSigma_e);
    double hat_Mu_f = hat_Sigma_f * (invSigma_e * newYi + invKK_f * newMu_f);
    f(i,uni_idx_f) = rnorm(1, hat_Mu_f, sqrt(hat_Sigma_f))(0);
    
    // Update C_f
    for(int d = 0; d < dimD; d++){
      N_f(C_f(i,d),d)--;
      NumericVector hatP_f(K_f(d));
      
      // stick-breaking
      NumericVector v(K_f(d));
      NumericVector pi(K_f(d));
      
      for(int k = 0; k < K_f(d); k++){
        if(k == K_f(d) - 1){
          v(k) = 1;
        }else{
          double tmp_sum = 0;
          for(int kk = k+1; kk < K_f(d); kk++)
            tmp_sum += N_f(d,kk);
          v(k) = rbeta(1, 1 + N_f(d,k), M2(d) + tmp_sum)(0);
        }
        
        if(k == 0){
          pi(k) = v(k);
        }else{
          double tmp_prod = 1;
          for(int kk = 0; kk < k; kk++)
            tmp_prod *= (1-v(kk));
          pi(k) = v(k) * tmp_prod;
        }
        
        hatP_f(k) = hatP_f(k) + log(pi(k));
        double new_mu = mu_f(uni_idx_f) - mu_d_f(d,uni_idx_f) + new_mu_d_k_f(d,k*Ji + uni_idx_f);
        hatP_f(k) = hatP_f(k) + R::dnorm(newYi, new_mu, sqrt(Sigma_e + newKK_f), true);
      }
      
      hatP_f = exp(hatP_f-max(hatP_f));
      if(sum(hatP_f) != 0){
        hatP_f = hatP_f/(sum(hatP_f));
        
        double tmp_rand = runif(1,0,1)(0);
        double tmp_sum = 0;
        int s;
        for(s = 0; s < K_f(d); s++){
          tmp_sum += hatP_f(s);
          if(tmp_sum > tmp_rand)
            break;
        }
        if(s < K_f(d))
          C_f(i,d) = s;
      }
      
      N_f(C_f(i,d),d)++;
      
      // Update a1,...,b
      a1(i,d) = u1_f(C_f(i,d),d);
      a2(i,d) = invlogit(u2_f(C_f(i,d),d));
      a3(i,d) = invlogit(u3_f(C_f(i,d),d));
      r(i,d) = ur_f(C_f(i,d),d);
      double a = a1(i,d) * (exp(a2(i,d)*r(i,d)/2) - 1)/(exp(a2(i,d)*r(i,d)/2) + 1);
      b(i,d) = invlogit(ub_f(C_f(i,d),d)) * a;
    }
    
    // Later for updating u1_f, .., ub_f
    Pa_uf += R::dnorm(newYi,newMu_f, sqrt(Sigma_e + newKK_f), true);
    Pb_uf += R::dnorm(newYi,b_newMu_f, sqrt(Sigma_e + newKK_f), true);
    
    // Later for updating sigma2_f and rho_f
    Pa1_f += R::dnorm(newYi, b_newMu_f, sqrt(Sigma_e + newKK_f),true);
    Pb1_f += R::dnorm(newYi, b_newMu_f, sqrt(Sigma_e + b_newKK_f),true);
    Pa2_f += R::dnorm(newYi, newMu_f, sqrt(Sigma_e + newKK_f),true);
    Pb2_f += R::dnorm(newYi, newMu_f, sqrt(Sigma_e + b_newKK_f),true);
  }
  
  if(cnt_idx_f > 1){
    arma::vec newYi = arma::zeros(cnt_idx_f);
    arma::vec newMu_f = arma::zeros(cnt_idx_f);
    arma::mat newKK_f = arma::zeros(cnt_idx_f,cnt_idx_f);
    arma::vec new_mu = arma::zeros(cnt_idx_f);
    arma::vec b_newMu_f = arma::zeros(cnt_idx_f);
    arma::mat b_newKK_f = arma::zeros(cnt_idx_f,cnt_idx_f);
    
    int idx_j = 0;
    int idx_jj = 0;
    for(int j = 0; j < Ji; j++){
      if(idx_f(j)){
        newYi(idx_j) = Yi(j);
        newMu_f(idx_j) = mu_f(j);
        newKK_f(idx_j, idx_j) = KK_f(j,j);
        b_newMu_f(idx_j) = b_mu_f(j);
        b_newKK_f(idx_j, idx_j) = b_KK_f(j,j);
        idx_jj = idx_j + 1;
        for(int jj = j + 1; jj < Ji; jj++){
          if(idx_f(jj)){
            newKK_f(idx_j, idx_jj) = KK_f(j,jj);
            newKK_f(idx_jj, idx_j) = newKK_f(idx_j, idx_jj);
            b_newKK_f(idx_j, idx_jj) = b_KK_f(j,jj);
            b_newKK_f(idx_jj, idx_j) = b_newKK_f(idx_j, idx_jj);
            idx_jj++;
          }
        }
        idx_j++;
      }
    }
    
    arma::mat eye_mat = arma::eye(cnt_idx_f, cnt_idx_f);
    arma::mat invKK_f = arma::inv(newKK_f);
    arma::mat hat_Sigma_f = arma::inv(invKK_f + eye_mat * invSigma_e);
    arma::vec hat_Mu_f = hat_Sigma_f * (eye_mat * invSigma_e * newYi + invKK_f * newMu_f);
    NumericVector tmp_f = rmvn(1, hat_Mu_f, hat_Sigma_f);
    idx_j = 0;
    for(int j = 0; j < Ji; j++){
      if(idx_f(j)){
        f(i,j) = tmp_f(idx_j);
        idx_j++;
      }
    }
    
    // Update C_f
    for(int d = 0; d < dimD; d++){
      N_f(C_f(i,d),d)--;
      NumericVector hatP_f(K_f(d));
      
      // stick-breaking
      NumericVector v(K_f(d));
      NumericVector pi(K_f(d));
      
      for(int k = 0; k < K_f(d); k++){
        if(k == K_f(d) - 1){
          v(k) = 1;
        }else{
          double tmp_sum = 0;
          for(int kk = k+1; kk < K_f(d); kk++)
            tmp_sum += N_f(d,kk);
          v(k) = rbeta(1, 1 + N_f(d,k), M2(d) + tmp_sum)(0);
        }
        
        if(k == 0){
          pi(k) = v(k);
        }else{
          double tmp_prod = 1;
          for(int kk = 0; kk < k; kk++)
            tmp_prod *= (1-v(kk));
          pi(k) = v(k) * tmp_prod;
        }
        
        hatP_f(k) = hatP_f(k) + log(pi(k));
        idx_j = 0;
        for(int j = 0; j < Ji; j++){
          if(idx_f(j)){
            new_mu(idx_j) = mu_f(j) - mu_d_f(d,j) + new_mu_d_k_f(d,k*Ji + j);
            idx_j++;
          }
        }
        
        hatP_f(k) += as<double>(dmvn(newYi.t(), new_mu.t(), eye_mat * Sigma_e + newKK_f, true));
      }
      
      hatP_f = exp(hatP_f-max(hatP_f));
      if(sum(hatP_f) != 0){
        hatP_f = hatP_f/(sum(hatP_f));
        
        double tmp_rand = runif(1,0,1)(0);
        double tmp_sum = 0;
        int s;
        for(s = 0; s < K_f(d); s++){
          tmp_sum += hatP_f(s);
          if(tmp_sum > tmp_rand)
            break;
        }
        if(s < K_f(d))
          C_f(i,d) = s;
      }
      
      N_f(C_f(i,d),d)++;
      
      // Update a1,...,b
      a1(i,d) = u1_f(C_f(i,d),d);
      a2(i,d) = invlogit(u2_f(C_f(i,d),d));
      a3(i,d) = invlogit(u3_f(C_f(i,d),d));
      r(i,d) = ur_f(C_f(i,d),d);
      double a = a1(i,d) * (exp(a2(i,d)*r(i,d)/2) - 1)/(exp(a2(i,d)*r(i,d)/2) + 1);
      b(i,d) = invlogit(ub_f(C_f(i,d),d)) * a;
    }
    
    // later for updating u1_f,...,ub_f
    Pa_uf += as<double>(dmvn(newYi.t(),newMu_f.t(), eye_mat * Sigma_e + newKK_f,true));
    Pb_uf += as<double>(dmvn(newYi.t(), b_newMu_f.t(), eye_mat * Sigma_e + newKK_f,true));
    
    // later for updating sigma2_f and rho_f
    Pa1_f += as<double>(dmvn(newYi.t(), b_newMu_f.t(), eye_mat * Sigma_e + newKK_f,true));
    Pb1_f += as<double>(dmvn(newYi.t(), b_newMu_f.t(), eye_mat * Sigma_e + b_newKK_f,true));
    Pa2_f += as<double>(dmvn(newYi.t(), newMu_f.t(), eye_mat * Sigma_e + newKK_f,true));
    Pb2_f += as<double>(dmvn(newYi.t(), newMu_f.t(), eye_mat * Sigma_e + b_newKK_f,true));
  }
}

// Gibbs to update C_f
void Update_C_f(int i, int dimD,
                NumericMatrix u1_f, NumericMatrix u2_f, NumericMatrix u3_f,
                NumericMatrix ur_f, NumericMatrix ub_f,
                double sd1_f, double sd2_f, double sd3_f, double sdr_f, double sdb_f,
                NumericMatrix a1, NumericMatrix a2, NumericMatrix a3,
                NumericMatrix r, NumericMatrix b,
                IntegerVector K_f, NumericVector M2,
                IntegerMatrix & N_f, IntegerMatrix & C_f,
                NumericMatrix & sum_a1, NumericMatrix & sum_a2, NumericMatrix & sum_a3,
                NumericMatrix & sum_r, NumericMatrix & sum_b){
  for(int d = 0; d < dimD; d++){
    N_f(d,C_f(i,d))--;
    NumericVector hatP_f(K_f(d));
    
    // Jacobian adjustment
    double g_r = a1(i,d) * (exp(a2(i,d)*r(i,d)/2) - 1)/(1+exp(a2(i,d)*r(i,d)/2));
    double zeta = std::abs(g_r/(b(i,d) * (g_r - b(i,d))));
    
    // stick-breaking
    NumericVector v(K_f(d));
    NumericVector pi(K_f(d));
    
    for(int k = 0; k < K_f(d); k++){
      if(k == K_f(d) - 1){
        v(k) = 1;
      }else{
        double tmp_sum = 0;
        for(int kk = k+1; kk < K_f(d); kk++)
          tmp_sum += N_f(d,kk);
        v(k) = rbeta(1, 1 + N_f(d,k), M2(d) + tmp_sum)(0);
      }
      
      if(k == 0){
        pi(k) = v(k);
      }else{
        double tmp_prod = 1;
        for(int kk = 0; kk < k; kk++)
          tmp_prod *= (1-v(kk));
        pi(k) = v(k) * tmp_prod;
      }
      
      hatP_f(k) = log(pi(k));
      hatP_f(k) += R::dnorm(a1(i,d),u1_f(k,d),sqrt(sd1_f),true);
      hatP_f(k) += R::dnorm(logit(a2(i,d)),u2_f(k,d),sqrt(sd2_f),true);
      hatP_f(k) += R::dnorm(logit(a3(i,d)),u3_f(k,d),sqrt(sd3_f),true);
      hatP_f(k) += R::dnorm(r(i,d),ur_f(k,d),sqrt(sdr_f),true);
      hatP_f(k) += R::dnorm(logit(b(i,d)/g_r),ub_f(k,d),sqrt(sdb_f),true);
      hatP_f(k) += log(zeta);
      hatP_f(k) = hatP_f(k) -log(a2(i,d)) - log(1-a2(i,d)) - log(a3(i,d)) - log(1-a3(i,d));
    }
    
    // Update C_f(i,d)
    hatP_f = exp(hatP_f-max(hatP_f));
    if(sum(hatP_f) != 0){
      hatP_f = hatP_f/(sum(hatP_f));
      
      double tmp_rand = runif(1,0,1)(0);
      double tmp_sum = 0;
      int s;
      for(s = 0; s < K_f(d); s++){
        tmp_sum += hatP_f(s);
        if(tmp_sum > tmp_rand)
          break;
      }
      if(s < K_f(d))
        C_f(i,d) = s;
    }
    
    N_f(d,C_f(i,d))++;
    sum_a1(d,C_f(i,d)) += a1(i,d);
    sum_a2(d,C_f(i,d)) += logit(a2(i,d));
    sum_a3(d,C_f(i,d)) += logit(a3(i,d));
    sum_r(d,C_f(i,d)) += r(i,d);
    sum_b(d,C_f(i,d)) += logit(b(i,d)/g_r);
  }
}

// Gibbs to update Beta
void Update_Beta(int i, int Ji, int dimD, int dimP, int dimJ,
                 NumericMatrix dataZ, NumericMatrix dataB_3D,
                 double Sigma_e, double invSigma_e,
                 NumericMatrix f, NumericMatrix U,
                 IntegerVector C_u, NumericMatrix Mu_b, NumericMatrix Sigma_b_3D,
                 NumericMatrix & Beta,
                 Function rmvn){
  arma::vec Yi = arma::zeros(Ji);
  for(int j = 0; j < Ji; j++)
    Yi(j) = dataZ(i,j) - U(i,j) - f(i,j);
  NumericVector origMu_b = Mu_b(C_u(i),_);
  arma::vec hatMu_b = as<arma::vec>(origMu_b);
  arma::mat hatSigma_b = arma::zeros(dimP,dimP);
  for(int p = 0; p < dimP; p++){
    hatSigma_b(p,p)= Sigma_b_3D(p, C_u(i) * dimP + p);
    for(int pp = p + 1; pp < dimP; pp++){
      hatSigma_b(p,pp) = Sigma_b_3D(p, C_u(i) * dimP + pp);
      hatSigma_b(pp,p) = hatSigma_b(p,pp);
    }
  }
  
  arma::mat invSigma_b = arma::inv(hatSigma_b);
  arma::mat tmpBB = arma::zeros(dimP, Ji);
  for(int p = 0; p < dimP; p++)
    for(int j = 0; j < Ji; j++)
      tmpBB(p,j) = dataB_3D(p, i * dimJ + j);
  arma::mat tmp_eye =  arma::eye(Ji,Ji) * invSigma_e;
  hatSigma_b = arma::inv(invSigma_b + tmpBB * tmp_eye * tmpBB.t());
  hatMu_b = hatSigma_b * (tmpBB * tmp_eye * Yi + invSigma_b * hatMu_b);
  Beta(i,_) = as<NumericVector>(rmvn(1,hatMu_b.t(),hatSigma_b));
}

// Gibbs to update Beta in ITR
void Update_ITR_Beta(int i, int Ji, int dimD, int dimP, int dimJ,
                 arma::vec beta0, arma::mat S0,
                 NumericMatrix dataZ, NumericMatrix dataB_3D,
                 double Sigma_e, double invSigma_e,
                 NumericMatrix f, NumericMatrix U,
                 NumericMatrix & Beta,
                 Function rmvn){
  arma::vec Yi = arma::zeros(Ji);
  for(int j = 0; j < Ji; j++)
    Yi(j) = dataZ(i,j) - U(i,j) - f(i,j);

  arma::vec hatMu_b = beta0;
  arma::mat hatSigma_b = S0;
  
  arma::mat invSigma_b = arma::inv(hatSigma_b);
  arma::mat tmpBB = arma::zeros(dimP, Ji);
  for(int p = 0; p < dimP; p++)
    for(int j = 0; j < Ji; j++)
      tmpBB(p,j) = dataB_3D(p, i * dimJ + j);
  arma::mat tmp_eye =  arma::eye(Ji,Ji) * invSigma_e;
  hatSigma_b = arma::inv(invSigma_b + tmpBB * tmp_eye * tmpBB.t());
  hatMu_b = hatSigma_b * (tmpBB * tmp_eye * Yi + invSigma_b * hatMu_b);
  Beta(i,_) = as<NumericVector>(rmvn(1,hatMu_b.t(),hatSigma_b));
}

// Gibbs to update u
// MH to update sigma2_u, rho_u
void Update_U(int i, int Ji, int dimD, int dimP, int dimJ,
              double mh_sd, double prior_unit, NumericVector W,
              NumericMatrix dataZ, NumericMatrix dataT, NumericMatrix dataB_3D,
              NumericMatrix dataIU,
              double Sigma_e, double invSigma_e, NumericMatrix f,
              double mu_log_sd2_u0, double sigma_log_sd2_u0, double sigma_log_sd2_u,
              double mu_logit_rho_u0, double sigma_logit_rho_u0, double sigma_logit_rho_u,
              NumericMatrix Beta,
              NumericVector & sigma2_u, NumericVector & rho_u,
              NumericMatrix & U,
              Function rmvn,
              Function dmvn){
  // double Pa = 0.0;
  // double Pb = 0.0;
  arma::vec Yi = arma::zeros(Ji);
  arma::vec IUi = arma::zeros(Ji);
  for(int j = 0; j < Ji; j++){
    double BBij = 0.0;
    for(int p=0; p < dimP; p++){
      BBij += dataB_3D(p,i*dimJ + j) * Beta(i,p);
    }
    Yi(j) = dataZ(i,j) - f(i,j) - BBij;
    IUi(j) = dataIU(i,j);
  }
  
  double b_sd2 = propose_positive(sigma2_u(i), mh_sd);
  double b_rho = propose_prob(rho_u(i), mh_sd/2.0);
  
  arma::mat KK = arma::zeros(Ji, Ji);
  arma::mat b_KK = arma::zeros(Ji, Ji);
  for(int j = 0; j < Ji; j++){
    KK(j,j) = sigma2_u(i);
    b_KK(j,j) = b_sd2;
    for(int jj = j+1; jj < Ji; jj++){
      if(IUi(jj) != IUi(j))
        break;
      KK(j,jj) = sigma2_u(i) * pow(rho_u(i), std::abs(dataT(i,j) - dataT(i,jj))/60.0);
      b_KK(j,jj) = b_sd2 * pow(b_rho, std::abs(dataT(i,j) - dataT(i,jj))/60.0);
      KK(jj,j) = KK(j,jj);
      b_KK(jj,j) = b_KK(j,jj);
    }
  }
  
  // Update U
  arma::mat eye_mat = arma::eye(Ji,Ji);
  arma::mat Sigma = arma::inv(arma::inv(KK) +  eye_mat * invSigma_e);
  arma::vec Mu = Sigma * (eye_mat * invSigma_e) * Yi;
  NumericVector tmpU = rmvn(1, Mu.t(), Sigma);
  for(int j = 0; j < Ji; j++)
    U(i,j) = tmpU(j);
  
  // Pa = as<double>(dmvn(Yi.t(),(arma::zeros(Ji)).t(), KK + eye_mat * Sigma_e,true)) +
  //   R::dnorm(log(sigma2_u(i)),mu_log_sd2_u0, sqrt(sigma_log_sd2_u + sigma_log_sd2_u0),true) - log(sigma2_u(i)) +
  //   R::dnorm(logit(rho_u(i)), mu_logit_rho_u0, sqrt(sigma_logit_rho_u + sigma_logit_rho_u0),true) - log(rho_u(i)) - log(1-rho_u(i)) -
  //   R::pnorm(sigma2_u(i), 0, mh_sd, true, true);
  // 
  // Pb = as<double>(dmvn(Yi.t(),(arma::zeros(Ji)).t(),b_KK + eye_mat * Sigma_e,true)) +
  //   R::dnorm(log(b_sd2),mu_log_sd2_u0, sqrt(sigma_log_sd2_u + sigma_log_sd2_u0),true) - log(b_sd2) +
  //   R::dnorm(logit(b_rho), mu_logit_rho_u0, sqrt(sigma_logit_rho_u + sigma_logit_rho_u0),true) - log(b_rho) - log(1-b_rho) -
  //   R::pnorm(b_sd2, 0, mh_sd, true, true);
  // 
  // // Update sigma2_u, rho_u
  // double Paccept = exp(Pb-Pa);
  // if(runif(1,0,1)(0) < Paccept){
  //   sigma2_u(i) = b_sd2;
  //   rho_u(i) = b_rho;
  // }
}

// Gibbs to update u
// MH to update sigma2_u, rho_u
void Update_DP_U(int i, int Ji, int dimD, int dimP, int dimJ,
              double mh_sd, double prior_unit, NumericVector W, 
              int K_u, int M1,
              NumericMatrix dataZ, NumericMatrix dataT, NumericMatrix dataB_3D,
              NumericMatrix dataIU,
              double Sigma_e, double invSigma_e, NumericMatrix f,
              double mu_log_sd2_u0, double sigma_log_sd2_u0, 
              double mu_logit_rho_u0, double sigma_logit_rho_u0,
              NumericMatrix Beta, NumericMatrix Mu_b,
              NumericVector mu_log_sd2_u, NumericVector mu_logit_rho_u,
              NumericVector b_mu_log_sd2_u, NumericVector b_mu_logit_rho_u,
              IntegerVector & C_u, IntegerVector & N_u,
              NumericVector & sigma2_u, NumericVector & rho_u,
              NumericMatrix & U,
              double & Pa_u, double & Pb_u,
              Function rmvn,
              Function dmvn){
  // Update Beta
  for(int p = 0; p < dimP; p++)
    Beta(i,p) = Mu_b(C_u(i),p);
  
  // Update sigma2_u and rho_u
  sigma2_u(i) = exp(mu_log_sd2_u(C_u(i)));
  rho_u(i) = invlogit(mu_logit_rho_u(C_u(i)));
  
  arma::vec Yi = arma::zeros(Ji);
  arma::vec IUi = arma::zeros(Ji);
  for(int j = 0; j < Ji; j++){
    double BBij = 0.0;
    for(int p=0; p < dimP; p++){
      BBij += dataB_3D(p,i*dimJ + j) * Beta(i,p);
    }
    Yi(j) = dataZ(i,j) - f(i,j) - BBij;
    IUi(j) = dataIU(i,j);
  }
  
  double b_sd2 = propose_positive(sigma2_u(i), mh_sd);
  double b_rho = propose_prob(rho_u(i), mh_sd/2.0);
  
  arma::mat KK = arma::zeros(Ji, Ji);
  arma::mat b_KK = arma::zeros(Ji, Ji);
  for(int j = 0; j < Ji; j++){
    KK(j,j) = sigma2_u(i);
    b_KK(j,j) = b_sd2;
    for(int jj = j+1; jj < Ji; jj++){
      if(IUi(jj) != IUi(j))
        break;
      KK(j,jj) = sigma2_u(i) * pow(rho_u(i), std::abs(dataT(i,j) - dataT(i,jj))/60.0);
      b_KK(j,jj) = b_sd2 * pow(b_rho, std::abs(dataT(i,j) - dataT(i,jj))/60.0);
      KK(jj,j) = KK(j,jj);
      b_KK(jj,j) = b_KK(j,jj);
    }
  }
  
  // Update U
  arma::mat eye_mat = arma::eye(Ji,Ji);
  arma::vec zero_vec = arma::zeros(Ji);
  arma::mat Sigma = arma::inv(arma::inv(KK) +  eye_mat * invSigma_e);
  arma::vec Mu = Sigma * (eye_mat * invSigma_e) * Yi;
  NumericVector tmpU = rmvn(1, Mu.t(), Sigma);
  for(int j = 0; j < Ji; j++)
    U(i,j) = tmpU(j);
  
  // Update C_u
  N_u(C_u(i))--;
  NumericVector hatP_u(K_u);
  
  // stick-breaking
  NumericVector v(K_u);
  NumericVector pi(K_u);
  
  for(int k = 0; k < K_u; k++){
    if(k == K_u - 1){
      v(k) = 1;
    }else{
      double tmp_sum = 0;
      for(int kk = k+1; kk < K_u; kk++)
        tmp_sum += N_u(kk);
      v(k) = rbeta(1, 1 + N_u(k), M1 + tmp_sum)(0);
    }
    
    if(k == 0){
      pi(k) = v(k);
    }else{
      double tmp_prod = 1;
      for(int kk = 0; kk < k; kk++)
        tmp_prod *= (1-v(kk));
      pi(k) = v(k) * tmp_prod;
    }
    
    double new_sigma2_u = exp(mu_log_sd2_u(k));
    double new_rho_u = invlogit(mu_logit_rho_u(k));
    arma::mat newKK = arma::zeros(Ji, Ji);
    for(int j = 0; j < Ji; j++){
      newKK(j,j) = new_sigma2_u;
      for(int jj = j+1; jj < Ji; jj++){
        newKK(j,jj) = new_sigma2_u * pow(new_rho_u, std::abs(dataT(i,j) - dataT(i,jj))/60.0);
        newKK(jj,j) = newKK(j,jj);
      }
    }
    
    arma::vec newYi = arma::zeros(Ji);
    for(int j = 0; j < Ji; j++){
      double BBij = 0.0;
      for(int p=0; p < dimP; p++){
        BBij += dataB_3D(p,i*dimJ + j) * Mu_b(k,p);
      }
      newYi(j) = dataZ(i,j) - f(i,j) - BBij;
    }
    
    hatP_u(k) = log(pi(k));
    hatP_u(k) += as<double>(dmvn(newYi.t(), zero_vec.t(), newKK + eye_mat * Sigma_e, true));
  }
  
  // Update C_u
  hatP_u = exp(hatP_u-max(hatP_u));
  if(sum(hatP_u) != 0){
    hatP_u = hatP_u/(sum(hatP_u));
    double tmp_rand = runif(1,0,1)(0);
    double tmp_sum = 0;
    int s;
    for(s = 0; s < K_u; s++){
      tmp_sum += hatP_u(s);
      if(tmp_sum > tmp_rand)
        break;
    }
    if(s < K_u)
      C_u(i) = s;
  }
  
  N_u(C_u(i))++;
  
  
  // Later for updating mu_log_sd2_u and mu_logit_rho_u
  Pa_u += as<double>(dmvn(Yi.t(), zero_vec.t(), KK + eye_mat * Sigma_e, true)) + 
    R::dnorm(mu_log_sd2_u(C_u(i)),mu_log_sd2_u0, sqrt(sigma_log_sd2_u0),true) - log(sigma2_u(i)) +
    R::dnorm(mu_logit_rho_u(C_u(i)), mu_logit_rho_u0, sqrt(sigma_logit_rho_u0),true) - log(rho_u(i)) - log(1-rho_u(i)) - 
    R::pnorm(sigma2_u(i), 0, mh_sd, true, true);
  
  Pb_u += as<double>(dmvn(Yi.t(), zero_vec.t(), b_KK + eye_mat * Sigma_e, true)) + 
    R::dnorm(b_mu_log_sd2_u(C_u(i)),mu_log_sd2_u0, sqrt(sigma_log_sd2_u0),true) - log(b_sd2) +
    R::dnorm(b_mu_logit_rho_u(C_u(i)), mu_logit_rho_u0, sqrt(sigma_logit_rho_u0),true) - log(b_rho) - log(1-b_rho) -
    R::pnorm(b_sd2, 0, mh_sd, true, true);
}

// Gibbs to update C_u
void Update_C_u(int i, int dimP,
                int K_u, int M1,
                NumericMatrix Beta, NumericMatrix Mu_b, NumericMatrix Sigma_b_3D,
                NumericVector mu_log_sd2_u, NumericVector mu_logit_rho_u,
                double sigma_log_sd2_u, double sigma_logit_rho_u,
                NumericVector sigma2_u, NumericVector rho_u,
                IntegerVector & C_u, IntegerVector & N_u,
                NumericMatrix & sum_beta, NumericMatrix & sum_S3D,
                NumericVector & sum_log_sd2_u, NumericVector & sum_logit_rho_u,
                Function dmvn){
  N_u(C_u(i))--;
  NumericVector hatP_u(K_u);
  
  // stick-breaking
  NumericVector v(K_u);
  NumericVector pi(K_u);
  
  for(int k = 0; k < K_u; k++){
    if(k == K_u - 1){
      v(k) = 1;
    }else{
      double tmp_sum = 0;
      for(int kk = k+1; kk < K_u; kk++)
        tmp_sum += N_u(kk);
      v(k) = rbeta(1, 1 + N_u(k), M1 + tmp_sum)(0);
    }
    
    if(k == 0){
      pi(k) = v(k);
    }else{
      double tmp_prod = 1;
      for(int kk = 0; kk < k; kk++)
        tmp_prod *= (1-v(kk));
      pi(k) = v(k) * tmp_prod;
    }
    
    arma::mat tmpS = arma::zeros(dimP,dimP);
    arma::vec tmpB = arma::zeros(dimP);
    arma::vec tmpMu = arma::zeros(dimP);
    for(int p = 0; p < dimP; p++){
      tmpB(p) = Beta(i,p);
      tmpMu(p) = Mu_b(k,p);
      tmpS(p,p) = Sigma_b_3D(p, k*dimP + p);
      for(int pp = p + 1; pp < dimP; pp++){
        tmpS(p,pp) = Sigma_b_3D(p, k * dimP + pp);
        tmpS(pp,p) = tmpS(p,pp);
      }
    }
    
    hatP_u(k) = log(pi(k));
    hatP_u(k) += as<double>(dmvn(tmpB.t(), tmpMu.t(), tmpS,true));
    hatP_u(k) += R::dnorm(log(sigma2_u(i)), mu_log_sd2_u(k), sqrt(sigma_log_sd2_u),true);
    hatP_u(k) += R::dnorm(logit(rho_u(i)), mu_logit_rho_u(k),sqrt(sigma_logit_rho_u),true);
    hatP_u(k) -= log(sigma2_u(i));
    hatP_u(k) = hatP_u(k) - log(rho_u(i)) - log(1-rho_u(i));
  }
  
  // Update C_u
  hatP_u = exp(hatP_u-max(hatP_u));
  if(sum(hatP_u) != 0){
    hatP_u = hatP_u/(sum(hatP_u));
    double tmp_rand = runif(1,0,1)(0);
    double tmp_sum = 0;
    int s;
    for(s = 0; s < K_u; s++){
      tmp_sum += hatP_u(s);
      if(tmp_sum > tmp_rand)
        break;
    }
    if(s < K_u)
      C_u(i) = s;
  }
  
  N_u(C_u(i))++;
  
  for(int p = 0; p < dimP; p++){
    sum_beta(p,C_u(i)) += Beta(i,p);
    sum_S3D(p, C_u(i) * dimP + p) += Beta(i,p) * Beta(i,p);
    for(int pp = p+1; pp < dimP; pp++){
      sum_S3D(p, C_u(i) * dimP + pp) += Beta(i,p) * Beta(i,pp);
      sum_S3D(pp, C_u(i) * dimP + p) += Beta(i,pp) * Beta(i,p);
    }
  }
  
  sum_log_sd2_u(C_u(i)) += log(sigma2_u(i));
  sum_logit_rho_u(C_u(i)) += logit(rho_u(i));
}

// Gibbs to update e
void Update_e(int i, int Ji, int dimP, int dimJ,
              NumericMatrix dataZ, NumericMatrix dataB_3D,
              NumericMatrix Beta,
              NumericMatrix U, NumericMatrix f,
              NumericVector a_e, NumericVector b_e,
              NumericVector & Sigma_e_I,
              Function rinvg){
  double sumYiYi = 0.0;
  for(int j = 0; j < Ji; j++){
    double BBij = 0.0;
    for(int p=0; p < dimP; p++){
      BBij += dataB_3D(p,i*dimJ + j) * Beta(i,p);
    }
    double Yi = dataZ(i,j) - U(i,j) - f(i,j) - BBij;
    sumYiYi += Yi * Yi;
  }
  Sigma_e_I(i) = as<double>(rinvg(1, a_e(i) + Ji/2.0, b_e(i) + sumYiYi/2.0));
}

// MH to update sigma2_f, rho_f
void Update_f_param(int dimD,
                    double mu_log_sd2_f0, double sigma_log_sd2_f0,
                    double mu_logit_rho_f0, double sigma_logit_rho_f0,
                    NumericVector b_rho_f, NumericVector b_sd2_f,
                    double Pa_f, double Pb_f,
                    NumericVector & rho_f, NumericVector & sigma2_f){
  for(int d = 0; d < dimD; d++){
    Pa_f = Pa_f + R::dnorm(log(sigma2_f(d)), mu_log_sd2_f0, sqrt(sigma_logit_rho_f0), true) - log(sigma2_f(d)) +
      R::dnorm(logit(rho_f(d)), mu_logit_rho_f0, sqrt(sigma_logit_rho_f0), true) - log(rho_f(d)) - log(1-rho_f(d));
    Pb_f = Pb_f + R::dnorm(log(b_sd2_f(d)), mu_log_sd2_f0, sqrt(sigma_logit_rho_f0), true) - log(b_sd2_f(d)) +
      R::dnorm(logit(b_rho_f(d)), mu_logit_rho_f0, sqrt(sigma_logit_rho_f0), true) - log(b_rho_f(d)) - log(1-b_rho_f(d));
  }
  
  double Paccept = exp(Pb_f - Pa_f);
  if(runif(1,0,1)(0) < Paccept){
    sigma2_f = b_sd2_f;
    rho_f = b_rho_f;
  }
}

// MH to update sigma2_f, rho_f
void Update_DP_f_param(int dimD, 
                    double mu_log_sd2_f0, double sigma_log_sd2_f0,
                    double mu_logit_rho_f0, double sigma_logit_rho_f0,
                    NumericVector b_rho_f, NumericVector b_sd2_f,
                    double Pa_f, double Pb_f,
                    double Pa_uf, double Pb_uf,
                    double Pa1_f, double Pb1_f, double Pa2_f, double Pb2_f,
                    NumericMatrix b_u1_f, NumericMatrix b_u2_f, NumericMatrix b_u3_f,
                    NumericMatrix b_ur_f, NumericMatrix b_ub_f,
                    NumericVector & rho_f, NumericVector & sigma2_f,
                    NumericMatrix & u1_f, NumericMatrix & u2_f, NumericMatrix & u3_f,
                    NumericMatrix & ur_f, NumericMatrix & ub_f){
  // Update u1_f, .., ub_f
  double Paccept = exp(Pb_uf - Pa_uf);
  if(runif(1,0,1)(0) < Paccept){
    u1_f = b_u1_f;
    u2_f = b_u2_f;
    u3_f = b_u3_f;
    ur_f = b_ur_f;
    ub_f = b_ub_f;
    Pa_f += Pa1_f;
    Pb_f += Pb1_f;
  }else{
    Pa_f += Pa2_f;
    Pb_f += Pb2_f;
  }
  
  // Update sigma2_f, rho_f
  for(int d = 0; d < dimD; d++){
    Pa_f = Pa_f + R::dnorm(log(sigma2_f(d)), mu_log_sd2_f0, sqrt(sigma_logit_rho_f0), true) - log(sigma2_f(d)) +
      R::dnorm(logit(rho_f(d)), mu_logit_rho_f0, sqrt(sigma_logit_rho_f0), true) - log(rho_f(d)) - log(1-rho_f(d));
    Pb_f = Pb_f + R::dnorm(log(b_sd2_f(d)), mu_log_sd2_f0, sqrt(sigma_logit_rho_f0), true) - log(b_sd2_f(d)) +
      R::dnorm(logit(b_rho_f(d)), mu_logit_rho_f0, sqrt(sigma_logit_rho_f0), true) - log(b_rho_f(d)) - log(1-b_rho_f(d));
  }
  
  Paccept = exp(Pb_f - Pa_f);
  if(runif(1,0,1)(0) < Paccept){
    sigma2_f = b_sd2_f;
    rho_f = b_rho_f;
  }
}

// Gibbs to update mu_logit_rho_u, mu_log_sd2_u
void Update_mu_u(int dimP, int K_u,
                 double k0, double v0, arma::vec beta0, arma::mat S0,
                 IntegerVector N_u,
                 double mu_log_sd2_u0, double sigma_log_sd2_u0, double sigma_log_sd2_u,
                 double mu_logit_rho_u0, double sigma_logit_rho_u0, double sigma_logit_rho_u,
                 NumericMatrix sum_beta, NumericMatrix sum_S3D,
                 NumericVector sum_log_sd2_u, NumericVector sum_logit_rho_u,
                 NumericVector & mu_log_sd2_u, NumericVector & mu_logit_rho_u,
                 NumericMatrix & Mu_b, NumericMatrix & Sigma_b_3D,
                 Function rmvn,
                 Function riwish){
  for(int k = 0; k < K_u; k++){
    // double hatSigma = 0;
    // double hatMu = 0;
    // // Update mu_logit_rho_u(k)
    //   hatSigma = 1/(1/sigma_logit_rho_u0 + N_u(k)/sigma_logit_rho_u);
    //   hatMu = (mu_logit_rho_u0/sigma_logit_rho_u0 + sum_logit_rho_u(k)/sigma_logit_rho_u) * hatSigma;
    //   mu_logit_rho_u(k) = R::rnorm(hatMu, sqrt(hatSigma));
    //
    // // Update mu_log_sd2_u(k)
    //   hatSigma = 1/(1/sigma_log_sd2_u0 + N_u(k)/sigma_log_sd2_u);
    //   hatMu = (mu_log_sd2_u0/sigma_log_sd2_u0 + sum_log_sd2_u(k)/sigma_log_sd2_u) * hatSigma;
    //   mu_log_sd2_u(k) = R::rnorm(hatMu,sqrt(hatSigma));
    
    // Update Mu_b(k), Sigma_b(k)
    double hat_k = k0 + N_u(k);
    //double hat_v = v0 + N_u(k);
    
    arma::vec sum_beta_k = arma::zeros(dimP);
    arma::mat sum_S_k = arma::zeros(dimP, dimP);
    for(int p = 0; p < dimP; p++){
      sum_beta_k(p) = sum_beta(p,k);
      sum_S_k(p,p) = sum_S3D(p, k*dimP + p);
      for(int pp = p+1; pp < dimP; pp++){
        sum_S_k(p,pp) = sum_S3D(p,k*dimP + pp);
        sum_S_k(pp,p) = sum_S_k(p,pp);
      }
    }
    
    arma::vec hat_beta = (k0 * beta0 + sum_beta_k)/hat_k;
    arma::mat hat_S = S0 + sum_S_k + k0 * beta0 * beta0.t() - hat_k * hat_beta * hat_beta.t();
    
    // NumericMatrix new_Sigma = as<NumericMatrix>(riwish(hat_v,hat_S));
    arma::mat new_Sigma = arma::eye(dimP, dimP) * 0.09;
    NumericVector new_Mu_b = as<NumericVector>(rmvn(1,hat_beta,new_Sigma/hat_k));
    for(int p = 0; p < dimP; p++){
      Mu_b(k,p) = new_Mu_b(p);
      Sigma_b_3D(p, k*dimP+p) = new_Sigma(p,p);
      for(int pp = p+1; pp < dimP; pp++){
        Sigma_b_3D(p, k*dimP+pp) = new_Sigma(p,pp);
        Sigma_b_3D(pp, k*dimP+p) = new_Sigma(pp,p);
      }
    }
  }
}

// Gibbs to update mu_logit_rho_u, mu_log_sd2_u
void Update_DP_mu_u(int dimP, int dimI, int dimJ, int K_u,
                 NumericVector J_I,
                 NumericMatrix dataZ, NumericMatrix dataB_3D,
                 arma::vec beta0, arma::mat S0,
                 IntegerVector N_u, IntegerVector C_u,
                 double mu_log_sd2_u0, double sigma_log_sd2_u0, 
                 double mu_logit_rho_u0, double sigma_logit_rho_u0, 
                 double Pa_u, double Pb_u,
                 NumericMatrix U, NumericMatrix f, NumericVector Sigma_e_I,
                 NumericVector b_mu_log_sd2_u, NumericVector b_mu_logit_rho_u,
                 NumericVector & mu_log_sd2_u, NumericVector & mu_logit_rho_u,
                 NumericMatrix & Mu_b, 
                 Function rmvn){
  
  // Update mu_log_sd2_u, mu_logit_rho_u
  double Paccept = exp(Pb_u-Pa_u);
  if(runif(1,0,1)(0) < Paccept){
    for(int k = 0; k < K_u; k++){
      mu_log_sd2_u(k) = b_mu_log_sd2_u(k);
      mu_logit_rho_u(k) = b_mu_logit_rho_u(k);
    }
  }
  
  // Update Mu_b
  for(int k = 0; k < K_u; k++){
    int sumJ = 0;
    for(int i = 0; i < dimI; i++){
      if(C_u(i) == k)
        sumJ += J_I(i);
    }
    
    if(sumJ == 0)
      continue;
    arma::vec bY = arma::zeros(sumJ);
    arma::mat bB = arma::zeros(dimP, sumJ);
    arma::mat tmp_eye = arma::eye(sumJ, sumJ);
    
    int idx_b = 0;
    for(int i = 0; i < dimI; i++){
      if(C_u(i) == k){
        for(int j = 0; j < J_I[i]; j++){
          bY(idx_b) = dataZ(i,j) - U(i,j) - f(i,j);
          for(int p = 0; p < dimP; p++)
            bB(p, idx_b) = dataB_3D(p, i*dimJ +j);
          tmp_eye(idx_b, idx_b) = (double)(1.0/Sigma_e_I(i));
          idx_b++;
        }
      }
    }
    
    arma::mat invSigma_b = arma::inv(S0);
    arma::mat hatSigma_b = arma::inv(invSigma_b + bB * tmp_eye * bB.t());
    arma::vec hatMu_b = hatSigma_b * (bB * tmp_eye * bY + invSigma_b * beta0);
    NumericVector new_Mu_b = as<NumericVector>(rmvn(1,hatMu_b,hatSigma_b));
    for(int p = 0; p < dimP; p++)
      Mu_b(k,p) = new_Mu_b(p);
  }
}

// Gibbs to update u1_f, u2_f, u3_f, ur_f, ub_f
void Update_mu_f(int dimD,
                 IntegerVector K_f,
                 double sd1_f, NumericVector u10_f, double sd10_f,
                 double sd2_f, NumericVector u20_f, double sd20_f,
                 double sd3_f, NumericVector u30_f, double sd30_f,
                 double sdr_f, NumericVector ur0_f, double sdr0_f,
                 double sdb_f, NumericVector ub0_f, double sdb0_f,
                 IntegerMatrix N_f,
                 NumericMatrix sum_a1, NumericMatrix sum_a2, NumericMatrix sum_a3,
                 NumericMatrix sum_r, NumericMatrix sum_b,
                 NumericMatrix & u1_f, NumericMatrix & u2_f, NumericMatrix & u3_f,
                 NumericMatrix & ur_f, NumericMatrix & ub_f){
  for(int d = 0 ; d < dimD; d++)
    for(int k = 0; k < K_f(d); k++){
      double hatMu = 0;
      double hatSigma = 0;
      hatSigma = 1.0/(1.0/sd10_f + N_f(d,k)/sd1_f);
      hatMu = (u10_f(d)/sd10_f + sum_a1(d,k)/sd1_f) * hatSigma;
      u1_f(k,d) = rnorm(1,hatMu,sqrt(hatSigma))(0);
      
      hatSigma = 1.0/(1.0/sd20_f + N_f(d,k)/sd2_f);
      hatMu = (u20_f(d)/sd20_f + sum_a2(d,k)/sd2_f) * hatSigma;
      u2_f(k,d) = rnorm(1,hatMu,sqrt(hatSigma))(0);
      
      hatSigma = 1.0/(1.0/sd30_f + N_f(d,k)/sd3_f);
      hatMu = (u30_f(d)/sd30_f + sum_a3(d,k)/sd3_f) * hatSigma;
      u3_f(k,d) = rnorm(1,hatMu,sqrt(hatSigma))(0);
      
      hatSigma = 1.0/(1.0/sdr0_f + N_f(d,k)/sdr_f);
      hatMu = (ur0_f(d)/sdr0_f + sum_r(d,k)/sdr_f) * hatSigma;
      ur_f(k,d) = rnorm(1,hatMu,sqrt(hatSigma))(0);
      
      hatSigma = 1.0/(1.0/sdb0_f + N_f(d,k)/sdb_f);
      hatMu = (ub0_f(d)/sdb0_f + sum_b(d,k)/sdb_f) * hatSigma;
      ub_f(k,d) = rnorm(1,hatMu,sqrt(hatSigma))(0);
    }
}

// Gibbs to update M1
void Update_M1(int dimI,
               double M_c, double M_d,
               int K_u,
               IntegerVector C_u,
               double & M1){
  double eta = R::rbeta(M1 + 1, dimI);
  IntegerVector k_ind(K_u);
  for(int i = 0; i < dimI; i++)
    k_ind(C_u(i)) = 1;
  int k = sum(k_ind);
  
  arma::vec hatP = arma::zeros(2);
  hatP(0) = (M_c + k -1)/(M_c + k - 1 + dimI*(M_d - log(eta)));
  hatP(1) = dimI*(M_d - log(eta))/(M_c + k - 1 + dimI*(M_d - log(eta)));
  
  arma::vec tmp_seq = arma::zeros(2);
  tmp_seq(1) = 1;
  int ind = RcppArmadillo::sample(tmp_seq,1,false,hatP)(0);
  if(ind == 0){
    M1 = rgamma(1, M_c+k, M_d-log(eta))(0);
  }else{
    M1 = rgamma(1, M_c+k-1, M_d-log(eta))(0);
  }
}

// Gibbs to update M2
void Update_M2(int dimD, int dimI,
               double M_c, double M_d,
               IntegerVector K_f,
               IntegerMatrix C_f,
               NumericVector & M2){
  for(int d = 0; d < dimD; d++){
    double eta = R::rbeta(M2(d) + 1, dimI);
    
    IntegerVector k_ind(K_f(d));
    for(int i = 0; i < dimI; i++){
      k_ind(C_f(i,d)) = 1;
    }
    int k = sum(k_ind);
    
    arma::vec hatP = arma::zeros(2);
    hatP(0) = (M_c + k -1)/(M_c + k - 1 + dimI*(M_d - log(eta)));
    hatP(1) = dimI*(M_d - log(eta))/(M_c + k - 1 + dimI*(M_d - log(eta)));
    
    arma::vec tmp_seq = arma::zeros(2);
    tmp_seq(1) = 1;
    int ind = RcppArmadillo::sample(tmp_seq,1,false,hatP)(0);
    if(ind == 0){
      M2(d) = rgamma(1, M_c+k, M_d-log(eta))(0);
    }else{
      M2(d) = rgamma(1, M_c+k-1, M_d-log(eta))(0);
    }
    
  }
}

#endif
