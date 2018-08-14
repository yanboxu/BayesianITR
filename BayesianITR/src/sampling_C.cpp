#include "sampling_update.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List sampling_dpm_C(const List input_vars,
                    const List data,
                    const List hyper_params,
                    List params,
                    Function rmvn, Function dmvn, Function rinvg, Function riwish){

  // *** std::string chainfile, int saveiter *** //
  List input_vars_list(input_vars);
  List data_list(data);
  List hyper_param_list(hyper_params);
  List param_list(params);
  
  const int dimI(data_list[0]);
  const int dimJ(data_list[1]);
  const int dimD(data_list[3]);
  const int dimP(data_list[4]);

  const NumericMatrix dataZ = as<NumericMatrix>(data_list[5]);
  const NumericMatrix dataT = as<NumericMatrix>(data_list[6]);
  const NumericMatrix dataB_3D = as<NumericMatrix>(data_list[7]);
  const NumericMatrix dataTau = as<NumericMatrix>(data_list[8]);
  const NumericMatrix dataA = as<NumericMatrix>(data_list[9]);
  const NumericMatrix dataIU = as<NumericMatrix>(data_list[10]);
  const NumericVector J_I(data_list[11]);
  const NumericVector L_I(data_list[12]);

  const int tol_iter(input_vars_list[0]);
  const int thin(input_vars_list[1]);
  const double mh_sd(input_vars_list[2]);
  const double prior_unit(input_vars_list[3]);
  const NumericVector W(input_vars_list[4]);
  const int K_u(input_vars_list[5]);
  const IntegerVector K_f(input_vars_list[6]);

  const NumericVector beta0(hyper_param_list[0]);
  const NumericMatrix S0 = as<NumericMatrix>(hyper_param_list[1]);
  const NumericVector a_e(hyper_param_list[2]);
  const NumericVector b_e(hyper_param_list[3]);

  const double mu_log_sd2_u0(hyper_param_list[4]);
  const double sigma_log_sd2_u0(hyper_param_list[5]);
  
  const double mu_logit_rho_u0(hyper_param_list[6]);
  const double sigma_logit_rho_u0(hyper_param_list[7]);
  
  const NumericVector u10_f(hyper_param_list[8]);
  const double sd10_f(hyper_param_list[9]);
  const NumericVector u20_f(hyper_param_list[10]);
  const double sd20_f(hyper_param_list[11]);
  const NumericVector u30_f(hyper_param_list[12]);
  const double sd30_f(hyper_param_list[13]);
  const NumericVector ur0_f(hyper_param_list[14]);
  const double sdr0_f(hyper_param_list[15]);
  const NumericVector ub0_f(hyper_param_list[16]);
  const double sdb0_f(hyper_param_list[17]);
  const double mu_log_sd2_f0(hyper_param_list[18]);
  const double sigma_log_sd2_f0(hyper_param_list[19]);
  const double mu_logit_rho_f0(hyper_param_list[20]);
  const double sigma_logit_rho_f0(hyper_param_list[21]);

  const double k0(hyper_param_list[22]);
  const double v0(hyper_param_list[23]);
  const double M_c(hyper_param_list[24]);
  const double M_d(hyper_param_list[25]);
  const double sigma_log_sd2_u(hyper_param_list[26]);
  const double sigma_logit_rho_u(hyper_param_list[27]);
  const double sd1_f(hyper_param_list[28]);
  const double sd2_f(hyper_param_list[29]);
  const double sd3_f(hyper_param_list[30]);
  const double sdr_f(hyper_param_list[31]);
  const double sdb_f(hyper_param_list[32]);
  
  NumericVector Sigma_e_I(param_list[0]);
  NumericMatrix Beta = as<NumericMatrix>(param_list[1]);
  NumericVector sigma2_u(param_list[2]);
  NumericVector rho_u(param_list[3]);
  NumericMatrix U = as<NumericMatrix>(param_list[4]);
  
  NumericMatrix a1 = as<NumericMatrix>(param_list[5]);
  NumericMatrix a2 = as<NumericMatrix>(param_list[6]);
  NumericMatrix a3 = as<NumericMatrix>(param_list[7]);
  NumericMatrix r = as<NumericMatrix>(param_list[8]);
  NumericMatrix b = as<NumericMatrix>(param_list[9]);
  NumericVector rho_f(param_list[10]);
  NumericVector sigma2_f(param_list[11]);
  NumericMatrix f = as<NumericMatrix>(param_list[12]);

  IntegerVector C_u(param_list[13]);
  NumericMatrix Mu_b = as<NumericMatrix>(param_list[14]);
  NumericMatrix Sigma_b_3D = as<NumericMatrix>(param_list[15]);
  NumericVector mu_log_sd2_u(param_list[16]);
  NumericVector mu_logit_rho_u(param_list[17]);
  IntegerMatrix C_f = as<IntegerMatrix>(param_list[18]);
  NumericMatrix u1_f = as<NumericMatrix>(param_list[19]);
  NumericMatrix u2_f = as<NumericMatrix>(param_list[20]);
  NumericMatrix u3_f = as<NumericMatrix>(param_list[21]);
  NumericMatrix ur_f = as<NumericMatrix>(param_list[22]);
  NumericMatrix ub_f = as<NumericMatrix>(param_list[23]);
  double M1(param_list[24]);
  NumericVector M2(param_list[25]);

  int _Thined_tol_iter = (int) (tol_iter/thin);
  int K = max(K_f);

  NumericMatrix __param_Sigma_e(dimI, _Thined_tol_iter);
  NumericMatrix __param_3DBeta(dimP, dimI * _Thined_tol_iter);
  NumericMatrix __param_3DMu_b(dimP, K_u * _Thined_tol_iter);
  NumericMatrix __param_4DSigma_b(dimP, dimP * K_u * _Thined_tol_iter);

  NumericMatrix __param_C_u(dimI, _Thined_tol_iter);
  NumericMatrix __param_sigma2_u(dimI, _Thined_tol_iter);
  NumericMatrix __param_rho_u(dimI, _Thined_tol_iter);
  NumericMatrix __param_mu_log_sd2_u(K_u, _Thined_tol_iter);
  NumericMatrix __param_mu_logit_rho_u(K_u, _Thined_tol_iter);
  NumericMatrix __param_3DU(dimJ, dimI * _Thined_tol_iter);

  NumericMatrix __param_3Da1(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Da2(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Da3(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Dr(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Db(dimD, dimI * _Thined_tol_iter);

  NumericMatrix __param_3Du1_f(dimD, K * _Thined_tol_iter);
  NumericMatrix __param_3Du2_f(dimD, K * _Thined_tol_iter);
  NumericMatrix __param_3Du3_f(dimD, K * _Thined_tol_iter);
  NumericMatrix __param_3Dur_f(dimD, K * _Thined_tol_iter);
  NumericMatrix __param_3Dub_f(dimD, K * _Thined_tol_iter);

  NumericMatrix __param_rho_f(dimD, _Thined_tol_iter);
  NumericMatrix __param_sigma2_f(dimD, _Thined_tol_iter);
  NumericMatrix __param_3DC_f(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Df(dimJ, dimI * _Thined_tol_iter);

  NumericVector __param_M1(_Thined_tol_iter);
  NumericMatrix __param_M2(dimD, _Thined_tol_iter);

  NumericVector __trainll(_Thined_tol_iter);

  IntegerVector N_u(K_u);
  IntegerMatrix N_f(dimD, K);
  for(int i = 0; i < dimI; i++){
    N_u(C_u(i))++;
    for(int d = 0; d < dimD; d++)
      N_f(d,C_f(i,d))++;
  }

  arma::vec a_beta0 = as<arma::vec>(beta0);
  arma::mat a_S0 = as<arma::mat>(S0);
  
  int __thined_iter = 0;
  for(int iter = 0; iter < tol_iter; iter++){
    NumericMatrix sum_beta(dimP, K_u);
    NumericMatrix sum_S3D(dimP, dimP * K_u);

    NumericVector sum_log_sd2_u(K_u);
    NumericVector sum_logit_rho_u(K_u);

    NumericMatrix sum_a1(dimD, K);
    NumericMatrix sum_a2(dimD, K);
    NumericMatrix sum_a3(dimD, K);
    NumericMatrix sum_r(dimD, K);
    NumericMatrix sum_b(dimD, K);

    double Pa_f = 0;
    double Pb_f = 0;
    NumericVector b_sd2_f(dimD);
    NumericVector b_rho_f(dimD);
    for(int d = 0; d < dimD; d++){
      b_sd2_f(d) = propose_positive(sigma2_f(d), mh_sd);
      b_rho_f(d) = propose_prob(rho_f(d), mh_sd/2);
      Pa_f -= R::pnorm(sigma2_f(d), 0, mh_sd, true, true);
      Pb_f -= R::pnorm(b_sd2_f(d), 0, mh_sd, true, true);
    }
    //Rcout << std::endl << "Iter " << iter;
    for(int i = 0; i < dimI; i++){
      // 1. Update f[i]
      if(L_I[i] > 0)
        Update_f(i, J_I(i), L_I(i), dimD, dimP, dimJ,
                 mh_sd, prior_unit, W,
                 dataZ, dataT, dataB_3D,dataTau, dataA,
                 Beta, U,
                 sd1_f, u10_f, sd10_f,
                 sd2_f, u20_f, sd20_f,
                 sd3_f, u30_f, sd30_f,
                 sdr_f, ur0_f, sdr0_f,
                 sdb_f, ub0_f, sdb0_f,
                 Sigma_e_I(i), (double) (1.0/Sigma_e_I(i)),
                 rho_f, sigma2_f,
                 b_rho_f, b_sd2_f,
                 a1, a2, a3, r, b, f, Pa_f, Pb_f,
                 rmvn, dmvn);

      // 2. Update C_f
      Update_C_f(i, dimD,
                 u1_f, u2_f, u3_f, ur_f, ub_f,
                 sd1_f, sd2_f, sd3_f, sdr_f, sdb_f,
                 a1, a2, a3, r, b,
                 K_f, M2,
                 N_f, C_f, sum_a1, sum_a2, sum_a3, sum_r, sum_b);

      // 3. Update Beta[i]
      Update_Beta(i, J_I(i), dimD, dimP, dimJ,
                  dataZ, dataB_3D,
                  Sigma_e_I(i), (double)(1.0/Sigma_e_I(i)),
                  f, U, C_u, Mu_b, Sigma_b_3D,
                  Beta,
                  rmvn);

      // 4. Update U[i], sigma2_u[i], rho_u[i]
      Update_U(i, J_I(i), dimD, dimP, dimJ,
               mh_sd, prior_unit, W,
               dataZ, dataT, dataB_3D, dataIU,
               Sigma_e_I(i), (double)(1.0/Sigma_e_I(i)), f,
               mu_log_sd2_u0, sigma_log_sd2_u0, sigma_log_sd2_u,
               mu_logit_rho_u0, sigma_logit_rho_u0, sigma_logit_rho_u,
               Beta, 
               sigma2_u, rho_u, U,
               rmvn, dmvn);

      // 5. Update C_u[i]
      Update_C_u(i, dimP,
                 K_u, M1,
                 Beta, Mu_b, Sigma_b_3D,
                 mu_log_sd2_u, mu_logit_rho_u,
                 sigma_log_sd2_u, sigma_logit_rho_u, sigma2_u, rho_u,
                 C_u, N_u, sum_beta, sum_S3D, sum_log_sd2_u, sum_logit_rho_u,
                 dmvn);

      // 6. Update Sigma_e
      Update_e(i, J_I(i), dimP, dimJ,
               dataZ, dataB_3D,
               Beta, U, f, a_e, b_e,
               Sigma_e_I,
               rinvg);

      if(iter % thin == 0){
        // Calculate the log-likelihood
        for(int j=0; j < J_I(i); j++){
          double Yij;
          double BBij = 0;
          for(int p=0; p < dimP; p++){
            BBij += dataB_3D(p,i*dimJ + j) * Beta(i,p);
          }
          Yij = dataZ(i,j) - BBij - U(i,j) - f(i,j);
          __trainll(__thined_iter) += R::dnorm(Yij,0,sqrt(Sigma_e_I(i)),true);
        }

        __param_Sigma_e(i, __thined_iter) = Sigma_e_I(i);
        __param_C_u(i, __thined_iter) = C_u(i);
        __param_sigma2_u(i, __thined_iter) = sigma2_u(i);
        __param_rho_u(i, __thined_iter) = rho_u(i);

        for(int p = 0; p < dimP; p++)
          __param_3DBeta(p, __thined_iter * dimI + i) = Beta(i,p);

        for(int j = 0; j < J_I(i); j++){
          __param_3DU(j, __thined_iter * dimI + i) = U(i,j);
          __param_3Df(j, __thined_iter * dimI + i) = f(i,j);
        }

        for(int d = 0; d < dimD; d++){
          __param_3Da1(d, __thined_iter * dimI + i) = a1(i,d);
          __param_3Da2(d, __thined_iter * dimI + i) = a2(i,d);
          __param_3Da3(d, __thined_iter * dimI + i) = a3(i,d);
          __param_3Dr(d, __thined_iter * dimI + i) = r(i,d);
          __param_3Db(d, __thined_iter * dimI + i) = b(i,d);
          __param_3DC_f(d, __thined_iter * dimI + i) = C_f(i,d);
        }
      }
    }
    // 7. Update sigma2_f, rho_f
    Update_f_param(dimD,
                   mu_log_sd2_f0, sigma_log_sd2_f0,
                   mu_logit_rho_f0, sigma_logit_rho_f0,
                   b_rho_f, b_sd2_f, Pa_f, Pb_f,
                   rho_f, sigma2_f);

    // 8. Update mu_log_sd2_u[k],mu_logit_rho_u[k]
    Update_mu_u(dimP, K_u,
                k0, v0, a_beta0, a_S0, N_u,
                mu_log_sd2_u0, sigma_log_sd2_u0, sigma_log_sd2_u,
                mu_logit_rho_u0, sigma_logit_rho_u0, sigma_logit_rho_u,
                sum_beta, sum_S3D,
                sum_log_sd2_u, sum_logit_rho_u,
                mu_log_sd2_u, mu_logit_rho_u, Mu_b, Sigma_b_3D,
                rmvn, riwish);

    // 9. Update u1-3,r_f[k,d]
    Update_mu_f(dimD, K_f,
                sd1_f, u10_f, sd10_f,
                sd2_f, u20_f, sd20_f,
                sd3_f, u30_f, sd30_f,
                sdr_f, ur0_f, sdr0_f,
                sdb_f, ub0_f, sdb0_f,
                N_f, sum_a1, sum_a2, sum_a3, sum_r, sum_b,
                u1_f, u2_f, u3_f, ur_f, ub_f);

    // 10. Update M1
    Update_M1(dimI,
              M_c, M_d, K_u, C_u,
              M1);
    //11. Update M2
    Update_M2(dimD, dimI,
              M_c, M_d, K_f, C_f,
              M2);

    if(iter % thin == 0){
      Rcout << "Thined iter: " << iter << std::endl;
      for(int k = 0; k < K_u; k++){
        for(int p = 0; p < dimP; p++){
          __param_3DMu_b(p, __thined_iter * K_u + k) = Mu_b(k,p);
          __param_4DSigma_b(_, __thined_iter * (K_u * dimP) + k * dimP + p) = Sigma_b_3D(_, k * dimP + p);
        }
        __param_mu_log_sd2_u(k, __thined_iter) = mu_log_sd2_u(k);
        __param_mu_logit_rho_u(k, __thined_iter) = mu_logit_rho_u(k);
      }
      for(int d = 0; d < dimD; d++){
        __param_rho_f(d, __thined_iter) = rho_f(d);
        __param_sigma2_f(d, __thined_iter) = sigma2_f(d);
        __param_M2(d, __thined_iter) = M2(d);
        for(int k = 0; k < K_f(d); k++){
          __param_3Du1_f(d, __thined_iter * K + k) = u1_f(k,d);
          __param_3Du2_f(d, __thined_iter * K + k) = u2_f(k,d);
          __param_3Du3_f(d, __thined_iter * K + k) = u3_f(k,d);
          __param_3Dur_f(d, __thined_iter * K + k) = ur_f(k,d);
          __param_3Dub_f(d, __thined_iter * K + k) = ub_f(k,d);
        }
      }
      __param_M1(__thined_iter) = M1;
      __thined_iter++;
    }
  }

  List output;
  __param_3DBeta.attr("dim") = IntegerVector::create(dimP,dimI,_Thined_tol_iter);
  __param_3DMu_b.attr("dim") = IntegerVector::create(dimP,K_u,_Thined_tol_iter);
  __param_4DSigma_b.attr("dim") = IntegerVector::create(dimP,dimP,K_u,_Thined_tol_iter);
  __param_3DU.attr("dim") = IntegerVector::create(dimJ,dimI,_Thined_tol_iter);

  __param_3Da1.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Da2.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Da3.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Dr.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Db.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);

  __param_3Du1_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  __param_3Du2_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  __param_3Du3_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  __param_3Dur_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  __param_3Dub_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);

  __param_3DC_f.attr("dim") = IntegerVector::create(dimD, dimI, _Thined_tol_iter);
  __param_3Df.attr("dim") = IntegerVector::create(dimJ, dimI, _Thined_tol_iter);

  output["Sigma_e"] = __param_Sigma_e;
  output["Beta"] = __param_3DBeta;
  output["Mu_b"] = __param_3DMu_b;
  output["Sigma_b"] = __param_4DSigma_b;

  output["C_u"] = __param_C_u;
  output["sigma2_u"] = __param_sigma2_u;
  output["rho_u"] = __param_rho_u;
  output["mu_log_sd2_u"] = __param_mu_log_sd2_u;
  output["mu_logit_rho_u"] = __param_mu_logit_rho_u;
  output["U"] = __param_3DU;

  output["a1"] = __param_3Da1;
  output["a2"] = __param_3Da2;
  output["a3"] = __param_3Da3;
  output["r"] = __param_3Dr;
  output["b"] = __param_3Db;

  output["u1_f"] = __param_3Du1_f;
  output["u2_f"] = __param_3Du2_f;
  output["u3_f"] = __param_3Du3_f;
  output["ur_f"] = __param_3Dur_f;
  output["ub_f"] = __param_3Dub_f;

  output["rho_f"] = __param_rho_f;
  output["sigma2_f"] = __param_sigma2_f;
  output["C_f"] = __param_3DC_f;
  output["f"] = __param_3Df;

  output["M1"] = __param_M1;
  output["M2"] = __param_M2;

  output["ll"] = __trainll;

  Rcout << "Sampling Done!" << std::endl;
  return(output);
}

// [[Rcpp::export]]
List sampling_ind_C(const List input_vars,
                    const List data,
                    const List hyper_params,
                    List params,
                    Function rmvn, Function dmvn, Function rinvg){
  
  // *** std::string chainfile, int saveiter *** //
  List input_vars_list(input_vars);
  List data_list(data);
  List hyper_param_list(hyper_params);
  List param_list(params);
  
  const int dimI(data_list[0]);
  const int dimJ(data_list[1]);
  const int dimD(data_list[3]);
  const int dimP(data_list[4]);
  
  const NumericMatrix dataZ = as<NumericMatrix>(data_list[5]);
  const NumericMatrix dataT = as<NumericMatrix>(data_list[6]);
  const NumericMatrix dataB_3D = as<NumericMatrix>(data_list[7]);
  const NumericMatrix dataTau = as<NumericMatrix>(data_list[8]);
  const NumericMatrix dataA = as<NumericMatrix>(data_list[9]);
  const NumericMatrix dataIU = as<NumericMatrix>(data_list[10]);
  const NumericVector J_I(data_list[11]);
  const NumericVector L_I(data_list[12]);
  
  const int tol_iter(input_vars_list[0]);
  const int thin(input_vars_list[1]);
  const double mh_sd(input_vars_list[2]);
  const double prior_unit(input_vars_list[3]);
  const NumericVector W(input_vars_list[4]);
  
  const NumericVector beta0(hyper_param_list[0]);
  const NumericMatrix S0 = as<NumericMatrix>(hyper_param_list[1]);
  const NumericVector a_e(hyper_param_list[2]);
  const NumericVector b_e(hyper_param_list[3]);
  
  const double mu_log_sd2_u0(hyper_param_list[4]);
  const double sigma_log_sd2_u0(hyper_param_list[5]);
  const double mu_logit_rho_u0(hyper_param_list[6]);
  const double sigma_logit_rho_u0(hyper_param_list[7]);
  
  const NumericVector u10_f(hyper_param_list[8]);
  const double sd10_f(hyper_param_list[9]);
  const NumericVector u20_f(hyper_param_list[10]);
  const double sd20_f(hyper_param_list[11]);
  const NumericVector u30_f(hyper_param_list[12]);
  const double sd30_f(hyper_param_list[13]);
  const NumericVector ur0_f(hyper_param_list[14]);
  const double sdr0_f(hyper_param_list[15]);
  const NumericVector ub0_f(hyper_param_list[16]);
  const double sdb0_f(hyper_param_list[17]);
  const double mu_log_sd2_f0(hyper_param_list[18]);
  const double sigma_log_sd2_f0(hyper_param_list[19]);
  const double mu_logit_rho_f0(hyper_param_list[20]);
  const double sigma_logit_rho_f0(hyper_param_list[21]);
  
  NumericVector Sigma_e_I(param_list[0]);
  NumericMatrix Beta = as<NumericMatrix>(param_list[1]);
  NumericVector sigma2_u(param_list[2]);
  NumericVector rho_u(param_list[3]);
  NumericMatrix U = as<NumericMatrix>(param_list[4]);
  
  NumericMatrix a1 = as<NumericMatrix>(param_list[5]);
  NumericMatrix a2 = as<NumericMatrix>(param_list[6]);
  NumericMatrix a3 = as<NumericMatrix>(param_list[7]);
  NumericMatrix r = as<NumericMatrix>(param_list[8]);
  NumericMatrix b = as<NumericMatrix>(param_list[9]);
  NumericVector rho_f(param_list[10]);
  NumericVector sigma2_f(param_list[11]);
  NumericMatrix f = as<NumericMatrix>(param_list[12]);
  
  int _Thined_tol_iter = (int) (tol_iter/thin);
  
  NumericMatrix __param_Sigma_e(dimI, _Thined_tol_iter);
  NumericMatrix __param_3DBeta(dimP, dimI * _Thined_tol_iter);
  NumericMatrix __param_sigma2_u(dimI, _Thined_tol_iter);
  NumericMatrix __param_rho_u(dimI, _Thined_tol_iter);
  NumericMatrix __param_3DU(dimJ, dimI * _Thined_tol_iter);
  
  NumericMatrix __param_3Da1(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Da2(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Da3(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Dr(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Db(dimD, dimI * _Thined_tol_iter);
  
  NumericMatrix __param_rho_f(dimD, _Thined_tol_iter);
  NumericMatrix __param_sigma2_f(dimD, _Thined_tol_iter);
  NumericMatrix __param_3Df(dimJ, dimI * _Thined_tol_iter);
  
  NumericVector __trainll(_Thined_tol_iter);
  
  arma::vec a_beta0 = as<arma::vec>(beta0);
  arma::mat a_S0 = as<arma::mat>(S0);
  
  int __thined_iter = 0;
  for(int iter = 0; iter < tol_iter; iter++){
    double Pa_f = 0;
    double Pb_f = 0;
    NumericVector b_sd2_f(dimD);
    NumericVector b_rho_f(dimD);
    for(int d = 0; d < dimD; d++){
      b_sd2_f(d) = propose_positive(sigma2_f(d), mh_sd);
      b_rho_f(d) = propose_prob(rho_f(d), mh_sd/2);
      Pa_f -= R::pnorm(sigma2_f(d), 0, mh_sd, true, true);
      Pb_f -= R::pnorm(b_sd2_f(d), 0, mh_sd, true, true);
    }
    //Rcout << std::endl << "Iter " << iter;
    for(int i = 0; i < dimI; i++){
      // 1. Update f[i]
      if(L_I[i] > 0)
        Update_f(i, J_I(i), L_I(i), dimD, dimP, dimJ,
                 mh_sd, prior_unit, W,
                 dataZ, dataT, dataB_3D,dataTau, dataA,
                 Beta, U,
                 0.0, u10_f, sd10_f,
                 0.0, u20_f, sd20_f,
                 0.0, u30_f, sd30_f,
                 0.0, ur0_f, sdr0_f,
                 0.0, ub0_f, sdb0_f,
                 Sigma_e_I(i), (double) (1.0/Sigma_e_I(i)),
                 rho_f, sigma2_f,
                 b_rho_f, b_sd2_f,
                 a1, a2, a3, r, b, f, Pa_f, Pb_f,
                 rmvn, dmvn);
      
      // 3. Update Beta[i]
      Update_ITR_Beta(i, J_I(i), dimD, dimP, dimJ,
                      a_beta0, a_S0,
                      dataZ, dataB_3D,
                      Sigma_e_I(i), (double)(1.0/Sigma_e_I(i)),
                      f, U,
                      Beta, rmvn);
      
      // 4. Update U[i], sigma2_u[i], rho_u[i]
      Update_U(i, J_I(i), dimD, dimP, dimJ,
               mh_sd, prior_unit, W,
               dataZ, dataT, dataB_3D, dataIU,
               Sigma_e_I(i), (double)(1.0/Sigma_e_I(i)), f,
               mu_log_sd2_u0, sigma_log_sd2_u0, 0.0,
               mu_logit_rho_u0, sigma_logit_rho_u0, 0.0,
               Beta,
               sigma2_u, rho_u, U,
               rmvn, dmvn);
      
      // 6. Update Sigma_e
      Update_e(i, J_I(i), dimP, dimJ,
               dataZ, dataB_3D,
               Beta, U, f, a_e, b_e,
               Sigma_e_I,
               rinvg);
      
      if(iter % thin == 0){
        // Calculate the log-likelihood
        for(int j=0; j < J_I(i); j++){
          double Yij;
          double BBij = 0;
          for(int p=0; p < dimP; p++){
            BBij += dataB_3D(p,i*dimJ + j) * Beta(i,p);
          }
          Yij = dataZ(i,j) - BBij - U(i,j) - f(i,j);
          __trainll(__thined_iter) += R::dnorm(Yij,0,sqrt(Sigma_e_I(i)),true);
        }
        
        __param_Sigma_e(i, __thined_iter) = Sigma_e_I(i);
        __param_sigma2_u(i, __thined_iter) = sigma2_u(i);
        __param_rho_u(i, __thined_iter) = rho_u(i);
        
        for(int p = 0; p < dimP; p++)
          __param_3DBeta(p, __thined_iter * dimI + i) = Beta(i,p);
        
        for(int j = 0; j < J_I(i); j++){
          __param_3DU(j, __thined_iter * dimI + i) = U(i,j);
          __param_3Df(j, __thined_iter * dimI + i) = f(i,j);
        }
        
        for(int d = 0; d < dimD; d++){
          __param_3Da1(d, __thined_iter * dimI + i) = a1(i,d);
          __param_3Da2(d, __thined_iter * dimI + i) = a2(i,d);
          __param_3Da3(d, __thined_iter * dimI + i) = a3(i,d);
          __param_3Dr(d, __thined_iter * dimI + i) = r(i,d);
          __param_3Db(d, __thined_iter * dimI + i) = b(i,d);
        }
      }
    }
    
    // 7. Update sigma2_f, rho_f
    Update_f_param(dimD,
                   mu_log_sd2_f0, sigma_log_sd2_f0,
                   mu_logit_rho_f0, sigma_logit_rho_f0,
                   b_rho_f, b_sd2_f, Pa_f, Pb_f,
                   rho_f, sigma2_f);
    
    if(iter % thin == 0){
      Rcout << "Thined iter: " << iter << std::endl;
      for(int d = 0; d < dimD; d++){
        __param_rho_f(d, __thined_iter) = rho_f(d);
        __param_sigma2_f(d, __thined_iter) = sigma2_f(d);
      }
      __thined_iter++;
    }
  }
  
  List output;
  __param_3DBeta.attr("dim") = IntegerVector::create(dimP,dimI,_Thined_tol_iter);
  __param_3DU.attr("dim") = IntegerVector::create(dimJ,dimI,_Thined_tol_iter);
  
  __param_3Da1.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Da2.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Da3.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Dr.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Db.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Df.attr("dim") = IntegerVector::create(dimJ, dimI, _Thined_tol_iter);
  
  output["Sigma_e"] = __param_Sigma_e;
  output["Beta"] = __param_3DBeta;
  
  output["sigma2_u"] = __param_sigma2_u;
  output["rho_u"] = __param_rho_u;
  output["U"] = __param_3DU;
  
  output["a1"] = __param_3Da1;
  output["a2"] = __param_3Da2;
  output["a3"] = __param_3Da3;
  output["r"] = __param_3Dr;
  output["b"] = __param_3Db;
  
  output["rho_f"] = __param_rho_f;
  output["sigma2_f"] = __param_sigma2_f;
  output["f"] = __param_3Df;
  
  output["ll"] = __trainll;
  
  Rcout << "Thined total iter:" << __thined_iter << std::endl;
  Rcout << "Sampling Done!" << std::endl;
  return(output);
}

// [[Rcpp::export]]
List sampling_dp_C(const List input_vars,
                    const List data,
                    const List hyper_params,
                    List params,
                    Function rmvn, Function dmvn, Function rinvg, Function riwish){
  
  // *** std::string chainfile, int saveiter *** //
  List input_vars_list(input_vars);
  List data_list(data);
  List hyper_param_list(hyper_params);
  List param_list(params);
  
  const int dimI(data_list[0]);
  const int dimJ(data_list[1]);
  const int dimD(data_list[3]);
  const int dimP(data_list[4]);
  
  const NumericMatrix dataZ = as<NumericMatrix>(data_list[5]);
  const NumericMatrix dataT = as<NumericMatrix>(data_list[6]);
  const NumericMatrix dataB_3D = as<NumericMatrix>(data_list[7]);
  const NumericMatrix dataTau = as<NumericMatrix>(data_list[8]);
  const NumericMatrix dataA = as<NumericMatrix>(data_list[9]);
  const NumericMatrix dataIU = as<NumericMatrix>(data_list[10]);
  const NumericVector J_I(data_list[11]);
  const NumericVector L_I(data_list[12]);
  
  const int tol_iter(input_vars_list[0]);
  const int thin(input_vars_list[1]);
  const double mh_sd(input_vars_list[2]);
  const double prior_unit(input_vars_list[3]);
  const NumericVector W(input_vars_list[4]);
  const int K_u(input_vars_list[5]);
  const IntegerVector K_f(input_vars_list[6]);
  
  const NumericVector beta0(hyper_param_list[0]);
  const NumericMatrix S0 = as<NumericMatrix>(hyper_param_list[1]);
  const NumericVector a_e(hyper_param_list[2]);
  const NumericVector b_e(hyper_param_list[3]);
  
  const double mu_log_sd2_u0(hyper_param_list[4]);
  const double sigma_log_sd2_u0(hyper_param_list[5]);
  
  const double mu_logit_rho_u0(hyper_param_list[6]);
  const double sigma_logit_rho_u0(hyper_param_list[7]);
  
  const NumericVector u10_f(hyper_param_list[8]);
  const double sd10_f(hyper_param_list[9]);
  const NumericVector u20_f(hyper_param_list[10]);
  const double sd20_f(hyper_param_list[11]);
  const NumericVector u30_f(hyper_param_list[12]);
  const double sd30_f(hyper_param_list[13]);
  const NumericVector ur0_f(hyper_param_list[14]);
  const double sdr0_f(hyper_param_list[15]);
  const NumericVector ub0_f(hyper_param_list[16]);
  const double sdb0_f(hyper_param_list[17]);
  const double mu_log_sd2_f0(hyper_param_list[18]);
  const double sigma_log_sd2_f0(hyper_param_list[19]);
  const double mu_logit_rho_f0(hyper_param_list[20]);
  const double sigma_logit_rho_f0(hyper_param_list[21]);
  const double M_c(hyper_param_list[22]);
  const double M_d(hyper_param_list[23]);
  
  NumericVector Sigma_e_I(param_list[0]);
  NumericMatrix Beta = as<NumericMatrix>(param_list[1]);
  NumericVector sigma2_u(param_list[2]);
  NumericVector rho_u(param_list[3]);
  NumericMatrix U = as<NumericMatrix>(param_list[4]);
  
  NumericMatrix a1 = as<NumericMatrix>(param_list[5]);
  NumericMatrix a2 = as<NumericMatrix>(param_list[6]);
  NumericMatrix a3 = as<NumericMatrix>(param_list[7]);
  NumericMatrix r = as<NumericMatrix>(param_list[8]);
  NumericMatrix b = as<NumericMatrix>(param_list[9]);
  NumericVector rho_f(param_list[10]);
  NumericVector sigma2_f(param_list[11]);
  NumericMatrix f = as<NumericMatrix>(param_list[12]);
  
  IntegerVector C_u(param_list[13]);
  NumericMatrix Mu_b = as<NumericMatrix>(param_list[14]);
  NumericVector mu_log_sd2_u(param_list[15]);
  NumericVector mu_logit_rho_u(param_list[16]);
  IntegerMatrix C_f = as<IntegerMatrix>(param_list[17]);
  NumericMatrix u1_f = as<NumericMatrix>(param_list[18]);
  NumericMatrix u2_f = as<NumericMatrix>(param_list[19]);
  NumericMatrix u3_f = as<NumericMatrix>(param_list[20]);
  NumericMatrix ur_f = as<NumericMatrix>(param_list[21]);
  NumericMatrix ub_f = as<NumericMatrix>(param_list[21]);
  double M1(param_list[23]);
  NumericVector M2(param_list[24]);
  
  int _Thined_tol_iter = (int) (tol_iter/thin);
  int K = max(K_f);
  
  NumericMatrix __param_Sigma_e(dimI, _Thined_tol_iter);
  NumericMatrix __param_3DBeta(dimP, dimI * _Thined_tol_iter);
  NumericMatrix __param_3DMu_b(dimP, K_u * _Thined_tol_iter);
  
  NumericMatrix __param_C_u(dimI, _Thined_tol_iter);
  NumericMatrix __param_sigma2_u(dimI, _Thined_tol_iter);
  NumericMatrix __param_rho_u(dimI, _Thined_tol_iter);
  NumericMatrix __param_mu_log_sd2_u(K_u, _Thined_tol_iter);
  NumericMatrix __param_mu_logit_rho_u(K_u, _Thined_tol_iter);
  NumericMatrix __param_3DU(dimJ, dimI * _Thined_tol_iter);
  
  NumericMatrix __param_3Da1(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Da2(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Da3(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Dr(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Db(dimD, dimI * _Thined_tol_iter);
  
  NumericMatrix __param_3Du1_f(dimD, K * _Thined_tol_iter);
  NumericMatrix __param_3Du2_f(dimD, K * _Thined_tol_iter);
  NumericMatrix __param_3Du3_f(dimD, K * _Thined_tol_iter);
  NumericMatrix __param_3Dur_f(dimD, K * _Thined_tol_iter);
  NumericMatrix __param_3Dub_f(dimD, K * _Thined_tol_iter);
  
  NumericMatrix __param_rho_f(dimD, _Thined_tol_iter);
  NumericMatrix __param_sigma2_f(dimD, _Thined_tol_iter);
  NumericMatrix __param_3DC_f(dimD, dimI * _Thined_tol_iter);
  NumericMatrix __param_3Df(dimJ, dimI * _Thined_tol_iter);
  
  NumericVector __param_M1(_Thined_tol_iter);
  NumericMatrix __param_M2(dimD, _Thined_tol_iter);
  
  NumericVector __trainll(_Thined_tol_iter);
  
  IntegerVector N_u(K_u);
  IntegerMatrix N_f(dimD, K);
  for(int i = 0; i < dimI; i++){
    N_u(C_u(i))++;
    for(int d = 0; d < dimD; d++)
      N_f(d,C_f(i,d))++;
  }
  
  arma::vec a_beta0 = as<arma::vec>(beta0);
  arma::mat a_S0 = as<arma::mat>(S0);
  
  int __thined_iter = 0;
  for(int iter = 0; iter < tol_iter; iter++){
    
    // Later for updating sigma2_f and rho_f
    double Pa_f = 0;
    double Pb_f = 0;
    NumericVector b_sd2_f(dimD);
    NumericVector b_rho_f(dimD);
    for(int d = 0; d < dimD; d++){
      b_sd2_f(d) = propose_positive(sigma2_f(d), mh_sd);
      b_rho_f(d) = propose_prob(rho_f(d), mh_sd/2);
      Pa_f -= R::pnorm(sigma2_f(d), 0, mh_sd, true, true);
      Pb_f -= R::pnorm(b_sd2_f(d), 0, mh_sd, true, true);
    }
    
    double Pa1_f = 0;
    double Pb1_f = 0;
    double Pa2_f = 0;
    double Pb2_f = 0;
    
    // Later for updating sigma2_u and rho_u
    double Pa_u = 0;
    double Pb_u = 0;
    NumericVector b_mu_log_sd2_u(K_u);
    NumericVector b_mu_logit_rho_u(K_u);
    for(int k = 0; k < K_u; k++){
      b_mu_log_sd2_u(k) = log(propose_positive(exp(mu_log_sd2_u(k)), mh_sd));
      b_mu_logit_rho_u(k) = logit(propose_prob(invlogit(mu_logit_rho_u(k)), mh_sd/2));
    }
    
    // Later for updating u1_f,..., ub_f
    double Pa_uf = 0;
    double Pb_uf = 0;
    NumericMatrix b_u1_f(K, dimD);
    NumericMatrix b_u2_f(K, dimD);
    NumericMatrix b_u3_f(K, dimD);
    NumericMatrix b_ur_f(K, dimD);
    NumericMatrix b_ub_f(K, dimD);
    for(int d = 0; d < dimD; d++){
      for(int k = 0; k < K_f(d); k++){
        b_u1_f(k,d) = propose_norm(u1_f(k,d), mh_sd);
        b_u2_f(k,d) = logit(propose_prob(invlogit(u2_f(k,d)), mh_sd));
        b_u3_f(k,d) = logit(propose_prob(invlogit(u3_f(k,d)), mh_sd));
        b_ur_f(k,d) = propose_norm(ur_f(k,d), mh_sd);
        b_ub_f(k,d) = logit(propose_prob(invlogit(ub_f(k,d)), mh_sd));
      }
    }
      
    //Rcout << std::endl << "Iter " << iter;
    for(int i = 0; i < dimI; i++){
      // 1. Update f, C_f, a1,.., b
      //Rcout << "1" << std::endl;
      if(L_I[i] > 0){
        Update_DP_f(i, J_I(i), L_I(i), dimD, dimP, dimJ, K,
                 mh_sd, prior_unit, W,
                 dataZ, dataT, dataB_3D,dataTau, dataA,
                 Beta, U, K_f, M2,
                 u10_f, sd10_f,
                 u20_f, sd20_f,
                 u30_f, sd30_f,
                 ur0_f, sdr0_f,
                 ub0_f, sdb0_f,
                 Sigma_e_I(i), (double) (1.0/Sigma_e_I(i)),
                 rho_f, sigma2_f,
                 b_rho_f, b_sd2_f,
                 u1_f, u2_f, u3_f, ur_f, ub_f,
                 b_u1_f, b_u2_f, b_u3_f, b_ur_f, b_ub_f,
                 N_f, C_f,
                 a1, a2, a3, r, b, f, Pa_uf, Pb_uf, Pa1_f, Pb1_f, Pa2_f, Pb2_f,
                 rmvn, dmvn);
      }
      
      //Rcout << "2" << std::endl;
      // 4. Update U[i], sigma2_u[i], rho_u[i]
      Update_DP_U(i, J_I(i), dimD, dimP, dimJ,
               mh_sd, prior_unit, W, K_u, M1,
               dataZ, dataT, dataB_3D, dataIU,
               Sigma_e_I(i), (double)(1.0/Sigma_e_I(i)), f,
               mu_log_sd2_u0, sigma_log_sd2_u0,
               mu_logit_rho_u0, sigma_logit_rho_u0,
               Beta, Mu_b,
               mu_log_sd2_u, mu_logit_rho_u,
               b_mu_log_sd2_u, b_mu_logit_rho_u,
               C_u, N_u,
               sigma2_u, rho_u, U,
               Pa_u, Pb_u,
               rmvn, dmvn);
      
      //Rcout << "3" << std::endl;
      // 6. Update Sigma_e
      Update_e(i, J_I(i), dimP, dimJ,
               dataZ, dataB_3D,
               Beta, U, f, a_e, b_e,
               Sigma_e_I,
               rinvg);
      
      if(iter % thin == 0){
        // Calculate the log-likelihood
        for(int j=0; j < J_I(i); j++){
          double Yij;
          double BBij = 0;
          for(int p=0; p < dimP; p++){
            BBij += dataB_3D(p,i*dimJ + j) * Beta(i,p);
          }
          Yij = dataZ(i,j) - BBij - U(i,j) - f(i,j);
          __trainll(__thined_iter) += R::dnorm(Yij,0,sqrt(Sigma_e_I(i)),true);
        }
        
        __param_Sigma_e(i, __thined_iter) = Sigma_e_I(i);
        __param_C_u(i, __thined_iter) = C_u(i);
        __param_sigma2_u(i, __thined_iter) = sigma2_u(i);
        __param_rho_u(i, __thined_iter) = rho_u(i);
        
        for(int p = 0; p < dimP; p++)
          __param_3DBeta(p, __thined_iter * dimI + i) = Beta(i,p);
        
        for(int j = 0; j < J_I(i); j++){
          __param_3DU(j, __thined_iter * dimI + i) = U(i,j);
          __param_3Df(j, __thined_iter * dimI + i) = f(i,j);
        }
        
        for(int d = 0; d < dimD; d++){
          __param_3Da1(d, __thined_iter * dimI + i) = a1(i,d);
          __param_3Da2(d, __thined_iter * dimI + i) = a2(i,d);
          __param_3Da3(d, __thined_iter * dimI + i) = a3(i,d);
          __param_3Dr(d, __thined_iter * dimI + i) = r(i,d);
          __param_3Db(d, __thined_iter * dimI + i) = b(i,d);
          __param_3DC_f(d, __thined_iter * dimI + i) = C_f(i,d);
        }
      }
    }
    //Rcout << "4" << std::endl;
    // 7. Update sigma2_f, rho_f
    Update_DP_f_param(dimD,
                   mu_log_sd2_f0, sigma_log_sd2_f0,
                   mu_logit_rho_f0, sigma_logit_rho_f0,
                   b_rho_f, b_sd2_f, Pa_f, Pb_f, Pa_uf, Pb_uf,
                   Pa1_f, Pb1_f, Pa2_f, Pb2_f,
                   b_u1_f, b_u2_f, b_u3_f, b_ur_f, b_ub_f,
                   rho_f, sigma2_f,
                   u1_f, u2_f, u3_f, ur_f, ub_f);
    
    //Rcout << "5" << std::endl;
    // 8. Update mu_log_sd2_u[k],mu_logit_rho_u[k]
    Update_DP_mu_u(dimP, dimI, dimJ, K_u, J_I,
                dataZ, dataB_3D,
                a_beta0, a_S0, N_u, C_u,
                mu_log_sd2_u0, sigma_log_sd2_u0, 
                mu_logit_rho_u0, sigma_logit_rho_u0, 
                Pa_u, Pb_u, U, f, Sigma_e_I,
                b_mu_log_sd2_u, b_mu_logit_rho_u,
                mu_log_sd2_u, mu_logit_rho_u, Mu_b,
                rmvn);
    
    //Rcout << "6" << std::endl;
    // 10. Update M1
    Update_M1(dimI,
              M_c, M_d, K_u, C_u,
              M1);
    
    //Rcout << "7" << std::endl;
    //11. Update M2
    Update_M2(dimD, dimI,
              M_c, M_d, K_f, C_f,
              M2);
    
    //Rcout << "8" << std::endl;
    if(iter % thin == 0){
      Rcout << "Thined iter: " << iter << std::endl;
      for(int k = 0; k < K_u; k++){
        for(int p = 0; p < dimP; p++){
          __param_3DMu_b(p, __thined_iter * K_u + k) = Mu_b(k,p);
        }
        __param_mu_log_sd2_u(k, __thined_iter) = mu_log_sd2_u(k);
        __param_mu_logit_rho_u(k, __thined_iter) = mu_logit_rho_u(k);
      }
      for(int d = 0; d < dimD; d++){
        __param_rho_f(d, __thined_iter) = rho_f(d);
        __param_sigma2_f(d, __thined_iter) = sigma2_f(d);
        __param_M2(d, __thined_iter) = M2(d);
        for(int k = 0; k < K_f(d); k++){
          __param_3Du1_f(d, __thined_iter * K + k) = u1_f(k,d);
          __param_3Du2_f(d, __thined_iter * K + k) = u2_f(k,d);
          __param_3Du3_f(d, __thined_iter * K + k) = u3_f(k,d);
          __param_3Dur_f(d, __thined_iter * K + k) = ur_f(k,d);
          __param_3Dub_f(d, __thined_iter * K + k) = ub_f(k,d);
        }
      }
      __param_M1(__thined_iter) = M1;
      __thined_iter++;
    }
  }
  
  List output;
  __param_3DBeta.attr("dim") = IntegerVector::create(dimP,dimI,_Thined_tol_iter);
  __param_3DMu_b.attr("dim") = IntegerVector::create(dimP,K_u,_Thined_tol_iter);
  __param_3DU.attr("dim") = IntegerVector::create(dimJ,dimI,_Thined_tol_iter);
  
  __param_3Da1.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Da2.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Da3.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Dr.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  __param_3Db.attr("dim") = IntegerVector::create(dimD,dimI,_Thined_tol_iter);
  
  __param_3Du1_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  __param_3Du2_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  __param_3Du3_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  __param_3Dur_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  __param_3Dub_f.attr("dim") = IntegerVector::create(dimD,K,_Thined_tol_iter);
  
  __param_3DC_f.attr("dim") = IntegerVector::create(dimD, dimI, _Thined_tol_iter);
  __param_3Df.attr("dim") = IntegerVector::create(dimJ, dimI, _Thined_tol_iter);
  
  output["Sigma_e"] = __param_Sigma_e;
  output["Beta"] = __param_3DBeta;
  output["Mu_b"] = __param_3DMu_b;
  
  output["C_u"] = __param_C_u;
  output["sigma2_u"] = __param_sigma2_u;
  output["rho_u"] = __param_rho_u;
  output["mu_log_sd2_u"] = __param_mu_log_sd2_u;
  output["mu_logit_rho_u"] = __param_mu_logit_rho_u;
  output["U"] = __param_3DU;
  
  output["a1"] = __param_3Da1;
  output["a2"] = __param_3Da2;
  output["a3"] = __param_3Da3;
  output["r"] = __param_3Dr;
  output["b"] = __param_3Db;
  
  output["u1_f"] = __param_3Du1_f;
  output["u2_f"] = __param_3Du2_f;
  output["u3_f"] = __param_3Du3_f;
  output["ur_f"] = __param_3Dur_f;
  output["ub_f"] = __param_3Dub_f;
  
  output["rho_f"] = __param_rho_f;
  output["sigma2_f"] = __param_sigma2_f;
  output["C_f"] = __param_3DC_f;
  output["f"] = __param_3Df;
  
  output["M1"] = __param_M1;
  output["M2"] = __param_M2;
  
  output["ll"] = __trainll;
  
  Rcout << "Sampling Done!" << std::endl;
  return(output);
}
