/*
     Project: Negative binomial models
              Power prior

     Date: June, 2024

          INPUT:

     PARAMETERS:
*/


data {
  int<lower=0>           N0;
  int<lower=0>           N1;
  array[N0] int<lower=0> Y0;
  array[N0] int<lower=0> Y1;

  real<lower=0, upper=1> WEIGHT;
  real                   THETA_EXT;
  real<lower=0>          SE_EXT;
  real                   THETA_WEAK;
  real<lower=0>          SE_WEAK;
  vector<lower=0>[2]     PRI_MU;
  vector<lower=0>[2]     PRI_K;
}

parameters {
  real                                    theta;
  real<lower=PRI_MU[1], upper=PRI_MU[2]>  mu0;
  real<lower=PRI_K[1],  upper=PRI_K[2]>   k;
}

model {
  mu0     ~ uniform(PRI_MU[1], PRI_MU[2]);
  k       ~ uniform(PRI_K[1],  PRI_K[2]);

  target += neg_binomial_lpmf(Y0 | k, k / mu0);
  target += neg_binomial_lpmf(Y1 | k, k / mu0 / exp(theta));
  target += normal_lpdf(theta | THETA_EXT, SE_EXT / sqrt(WEIGHT)) +
    normal_lpdf(theta | THETA_WEAK, SE_WEAK);
}
