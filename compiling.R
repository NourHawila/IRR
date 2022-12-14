
stan_compile_full <- function(){
  write("functions {
      matrix make_corr_matrix_R(int k,real rho){
      matrix[k,k] m;
      for (i in 1:k) {
      m[i,i] = 1;
      }
      for (i in 1:(k-1)) {
      for (j in (i+1):k){
      m[i,j] = rho;
      m[j,i] = rho;
      }}
      return m;
      }

      matrix make_corr_matrix_T(int k,real rho){
      matrix[k,k] m = rep_matrix(0,k,k);
      for (i in 1:k) {
      m[i,i] = 1;
      }
      for (i in 1:(k-1)) {
      m[i,i+1] = rho;
      m[i+1,i] = rho;
      }
      return m;
      }

      }

      data {
      int N;
      int I;
      int J;
      int K;
      int p;
      int y[N];
      matrix[N,p] X;
      int nlevels[p];
      int subject[N];
      int rater[N];
      int time[N];

      int cov_T_str;
      int cov_R_str;
      real rho_T_eta;
      real rho_R_eta;

      real beta_a;
      real beta_b;
      real gamma_a;
      real gamma_b;
      vector [p] beta_mean;
      matrix [p,p] beta_sigma;
      }

      parameters {
      vector[p] beta_parm;
      vector[I] b_S;
      vector[J] b_R[I];
      vector[K] b_T[J,I];
      real<lower=0> tau_S;
      real<lower=0> tau_R;
      real<lower=0> tau_T;

      real <lower=-1, upper=1> rho_R;
      real <lower=-1, upper=1> rho_T_common;
      vector<lower=-1, upper=1> [J] rho_T;

      corr_matrix[J] Rho_R_unstr;
      corr_matrix[K] Rho_T_unstr[J];
      }


      transformed parameters {
      vector[J] Corr_T;
      matrix[K,K] Rho_T[J];
      real sigma_S = 1/sqrt(tau_S);
      real sigma_R = 1/sqrt(tau_R);
      real sigma_T = 1/sqrt(tau_T);

      real Corr_R = (sigma_S^2+rho_R*sigma_R^2) / (sigma_S^2+sigma_R^2+sigma_T^2);
      vector[J] Sigma_R = rep_vector(sigma_R,J);
      vector[K] Sigma_T = rep_vector(sigma_T,K);
      matrix[J,J] Rho_R = make_corr_matrix_R(J,rho_R);


      for(i in 1:J){
      Corr_T[i]= (sigma_S^2+sigma_R^2+rho_T[i]*sigma_T^2) / (sigma_S^2+sigma_R^2+sigma_T^2);
      }

      if(cov_T_str==1){
      for(i in 1:J) {
      Rho_T[i] = make_corr_matrix_T(K,rho_T[i]);
      Corr_T[i]= (sigma_S^2+sigma_R^2+rho_T[i]*sigma_T^2) / (sigma_S^2+sigma_R^2+sigma_T^2);
      }
      }

      if(cov_T_str==2){
      for(i in 1:J){
      Rho_T[i] = make_corr_matrix_T(K,rho_T_common);
      Corr_T[i]= (sigma_S^2+sigma_R^2+rho_T_common*sigma_T^2) / (sigma_S^2+sigma_R^2+sigma_T^2);
      }
      }
      }


      model {
      // probability model
      vector[N] prob;
      for (i in 1:N) {
        prob[i] = X[i,]*beta_parm + b_S[subject[i]] + b_R[subject[i], rater[i]] + b_T[rater[i], subject[i],time[i]];
        //prob[i] = inv_logit(prob[i]);
        prob[i] = Phi(prob[i]);
      }
      y ~ binomial(1, prob);


      //adaptive priors
      b_S ~ normal(0,sigma_S);
      if(cov_R_str==0){
      b_R ~ multi_normal(rep_vector(0, J),
                        quad_form_diag(Rho_R_unstr, Sigma_R));
      }
      if(cov_R_str!=0){
      b_R ~ multi_normal(rep_vector(0, J),
                        quad_form_diag(Rho_R, Sigma_R));
      }
      if(cov_T_str==0){
      for(i in 1:J){
      b_T[i] ~ multi_normal(rep_vector(0, K),
                        quad_form_diag(Rho_T_unstr[i], Sigma_T));
      }
      }
      if(cov_T_str!=0){
      for(i in 1:J){
      b_T[i] ~ multi_normal(rep_vector(0, K),
                        quad_form_diag(Rho_T[i], Sigma_T));
      }
      }


      //fixed priors
      beta_parm ~ multi_normal(beta_mean,beta_sigma);
      rho_R ~ beta(beta_a,beta_b);
      rho_T_common ~ beta(beta_a, beta_b);
      for(i in 1:J) rho_T[i] ~ beta(beta_a,beta_b);

      Rho_R_unstr ~ lkj_corr(rho_R_eta);
      for(i in  1:J) {
      Rho_T_unstr[i] ~ lkj_corr(rho_T_eta);
      }

      tau_S ~ gamma(gamma_a,gamma_b);
      tau_R ~ gamma(gamma_a,gamma_b);
      tau_T ~ gamma(gamma_a,gamma_b);
      }

      generated quantities {
      vector[N] log_lik;
      vector[N] y_hat;
      vector[N] mu_hat;
      for (n in 1:N) {
        //log_lik[n] = bernoulli_logit_lpmf(y[n]|X[n,]*beta_parm+ b_S[subject[n]] + b_R[subject[n], rater[n]] + b_T[rater[n], subject[n],time[n]]);
        log_lik[n] = bernoulli_lpmf(y[n]|Phi(X[n]*beta_parm+ b_S[subject[n]] + b_R[subject[n], rater[n]]+ b_T[rater[n], subject[n],time[n]]));
        y_hat[n] = bernoulli_rng(Phi(X[n]*beta_parm+ b_S[subject[n]] + b_R[subject[n], rater[n]]+ b_T[rater[n], subject[n],time[n]]));
        mu_hat[n] = Phi(X[n]*beta_parm+ b_S[subject[n]] + b_R[subject[n], rater[n]]+ b_T[rater[n], subject[n],time[n]]);
      }
      }",
        "jointmodel.stan")
  m1<-stan_model(file="jointmodel.stan")
  m1
}
compiled = stan_compile_full()
save(compiled,file="compiled_full.Rdata")


stan_compile_ind <- function(){
  write("data {
      int N;
      int I;
      int J;
      int K;
      int p;
      int y[N];
      matrix[N,p] X;
      int subject[N];
      int rater[N];
      int time[N];
      real gamma_a;
      real gamma_b;
      vector[p] beta_mean;
      matrix[p,p] beta_sigma;
      }

      parameters {
      vector[p] beta_parm;
      vector[I] b_S;
      vector[J] b_R;
      vector[K] b_T;

      real<lower=0> tau_S;
      real<lower=0> tau_R;
      real<lower=0> tau_T;

      }


      transformed parameters {
      real sigma_S = 1/sqrt(tau_S);
      real sigma_R = 1/sqrt(tau_R);
      real sigma_T = 1/sqrt(tau_T);

      real Corr_R = (sigma_S^2+sigma_T^2) / (sigma_S^2+sigma_R^2+sigma_T^2);
      real Corr_T = (sigma_S^2+sigma_R^2) / (sigma_S^2+sigma_R^2+sigma_T^2);
      }


      model {
      // probability model
      vector[N] prob;
      for (i in 1:N) {
        prob[i] = X[i,]*beta_parm+ b_S[subject[i]] + b_R[rater[i]] + b_T[time[i]];
        prob[i] = Phi(prob[i]);
      }
      y ~ bernoulli(prob);


      //adaptive priors
      b_S ~ normal(0,sigma_S);
      b_R ~ normal(0,sigma_R);
      b_T ~ normal(0,sigma_T);

      //fixed priors
      beta_parm ~ multi_normal(beta_mean,beta_sigma);

      tau_S ~ gamma(gamma_a,gamma_b);
      tau_R ~ gamma(gamma_a,gamma_b);
      tau_T ~ gamma(gamma_a,gamma_b);
      }

      generated quantities {
      vector[N] log_lik;
      vector[N] y_hat;
      vector[N] mu_hat;
      for (n in 1:N) {
        log_lik[n] = bernoulli_lpmf(y[n]|Phi(X[n]*beta_parm+ b_S[subject[n]] + b_R[rater[n]] + b_T[time[n]]));
        y_hat[n] = bernoulli_rng(Phi(X[n]*beta_parm+ b_S[subject[n]] + b_R[rater[n]] + b_T[time[n]]));
        mu_hat[n] = Phi(X[n]*beta_parm+ b_S[subject[n]] + b_R[rater[n]] + b_T[time[n]]);
      }
      }",
        "indmodel.stan")
  m1<-stan_model(file="indmodel.stan")
  m1
}
compiled = stan_compile_ind()
save(compiled,file="compiled_ind.Rdata")

stan_compile_partially <- function(){
  write("functions {
      matrix make_corr_matrix_T(int k,real rho){
      matrix[k,k] m = rep_matrix(0,k,k);
      for (i in 1:k) {
      m[i,i] = 1;
      }
      for (i in 1:(k-1)) {
      m[i,i+1] = rho;
      m[i+1,i] = rho;
      }
      return m;
      }
      }

      data {
      int N;
      int I;
      int J;
      int K;
      int p;
      int y[N];
      matrix[N, p] X;
      vector[N] Xk;
      int subject[N];
      int rater[N];
      int time[N];
      int cov_T_str;
      real rho_T_eta;
      real beta_a;
      real beta_b;
      real gamma_a;
      real gamma_b;
      vector[p] beta_mean;
      matrix[p,p] beta_sigma;
      real betak_mean;
      real betak_sigma;
      }

      parameters {
      vector[p] beta_parm;
      real beta_k;
      vector[I] b_S;
      vector[K] b_T[J];

      real<lower=0> tau_S;
      real<lower=0> tau_T;
      vector<lower=-1, upper=1> [J] rho_T;
      real <lower=-1, upper=1> rho_T_common;
      corr_matrix[K] Rho_T_unstr[J];
      }


      transformed parameters {
      real sigma_S = 1/sqrt(tau_S);
      real sigma_T = 1/sqrt(tau_T);
      vector[J] Corr_T;
      matrix[K,K] Rho_T[J];
      vector[K] Sigma_T = rep_vector(sigma_T,K);

      for(i in 1:J){
      Rho_T[i] = make_corr_matrix_T(K,rho_T[i]);
      Corr_T[i]= (sigma_S^2+rho_T[i]*sigma_T^2) / (sigma_S^2+sigma_T^2);
      }
      // separate  intrarater
      if(cov_T_str==1){
      for(i in 1:J) {
      Rho_T[i] = make_corr_matrix_T(K,rho_T[i]);
      Corr_T[i]= (sigma_S^2+rho_T[i]*sigma_T^2) / (sigma_S^2+sigma_T^2);
      }
      }

      // for common intrarater
      if(cov_T_str==2){
      for(i in 1:J){
      Rho_T[i] = make_corr_matrix_T(K,rho_T_common);
      Corr_T[i]= (sigma_S^2+rho_T_common*sigma_T^2) / (sigma_S^2+sigma_T^2);
      }
      }
      }


      model {
      // probability model
      vector[N] prob;
      for (i in 1:N) {
        prob[i] = X[i,]*beta_parm -(Xk[i]*beta_k+b_S[subject[i]] + b_T[rater[i], time[i]]);
        prob[i] = Phi(prob[i]);
      }
      y ~ bernoulli(prob);


      //adaptive priors
      b_S ~ normal(0,sigma_S);
      if(cov_T_str==0){
      for(i in 1:J){
      b_T[i] ~ multi_normal(rep_vector(0, K),
                        quad_form_diag(Rho_T_unstr[i], Sigma_T));
      }}

      if(cov_T_str!=0){
      for(i in  1:J){
      b_T[i] ~ multi_normal(rep_vector(0, K),
                        quad_form_diag(Rho_T[i], Sigma_T));
      }
      }


      //fixed priors
      beta_parm ~ multi_normal(beta_mean,beta_sigma);
      beta_k ~ normal(betak_mean,betak_sigma);
      for(i in  1:J) {
      Rho_T_unstr[i] ~ lkj_corr(rho_T_eta);
      }
      rho_T_common ~ beta(beta_a, beta_b);
      for(i in 1:K) rho_T[i] ~ beta(beta_a,beta_b);

      tau_S ~ gamma(gamma_a,gamma_b);
      tau_T ~ gamma(gamma_a,gamma_b);
      }
      generated quantities {
      vector[N] log_lik;
      vector[N] y_hat;
      vector[N] mu_hat;
      for (n in 1:N) {
        log_lik[n] = bernoulli_lpmf(y[n]|Phi(X[n]*beta_parm -(Xk[n]*beta_k+b_S[subject[n]] + b_T[rater[n], time[n]])));
        y_hat[n] = bernoulli_rng(Phi(X[n]*beta_parm -(Xk[n]*beta_k+b_S[subject[n]] + b_T[rater[n], time[n]])));
        mu_hat[n] = Phi(X[n]*beta_parm -(Xk[n]*beta_k+b_S[subject[n]] + b_T[rater[n], time[n]]));
      }
      }",
        "jointmodel.stan")
  m1<-stan_model(file="jointmodel.stan")
  m1
}
compiled = stan_compile_partially()
save(compiled,file="compiled_partial.Rdata")
