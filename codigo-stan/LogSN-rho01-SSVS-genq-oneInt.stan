data {
  int<lower=1> K;                    // Number of groups
  int<lower=0> N_obs;                // Number of observed data
  int<lower=0> N_mis;                // Number of missing data
  int<lower=1> p;
  vector[N_obs] y_obs;                   // Observed dependent variable for all groups
  matrix[N_obs, p] X_obs;                // Matrix of observed explanatory variables for all groups
  matrix[N_mis, p] X_mis;             // Matrix of observed explanatory variables for all groups with missing y
  matrix[N_obs + N_mis, p] X_all;
  array[N_obs] int<lower=1, upper=K> group;
  array[N_mis] int<lower=1, upper=K> groupmis;
  array[N_obs+N_mis] int<lower=1, upper=K> group_all;
  real<lower=0> c;
  real<lower=0> tau;
  real epsilon;  // Threshold for significant coefficient
}


parameters {
  real<lower=0> sigma;               // sd of y
  vector<lower=0, upper=1>[K] rho;     // Correlation coefficient for each group
  real mu;                    // intercept for each region
  vector[p] b;                         // Regression coefficients for explanatory variables
  vector<lower=0>[K] z;                // Latent variable for each group
  // vector[N_mis] ymis;               // Non-observed dependent variable for all groups
  vector<lower=0, upper=1>[p] pr;        // probabilidad de inclusion
}


model {
  // Reference prior for rho
  target += sum(0.5 * log1p(rho.^2) - log1m(rho.^2));

  // https://sci-hub.se/https://www.jstor.org/stable/43601095
  // aproximacion de priori de Jeffrey
  // target += 0.5 * student_t_lupdf(rho | 1.0/2.0, 0, pi()/2.0);
  
  // Priors for regression coefficients b
  // b ~ normal(0, 100);  // Wide priors for regression coefficients
  
  // Flat prior for mu
  mu ~ normal(0, 10);
  
  // Prior for SSVS: seleccion of b and his inclusion probability pr.  
  pr ~ beta(0.5, 0.5);
  for (i in 1:p) {
    target += log_mix(pr[i],
                      normal_lpdf(b[i] | 0, tau * c),
                      normal_lpdf(b[i] | 0, tau));
  }

  // Prior for sigma
  // target += -sigma;   // log of 1/sigma
  sigma ~ inv_gamma(0.001, 0.001);

  // Prior for z (latent variable) for each group
  z ~ normal(0, sigma ./ sqrt(1 - 2/pi() * square(rho)) );  // Truncated normal prior for latent variables
  
  // Likelihood for observed y given z and b for each group
  
  // `sig` means the square root of variance parameter in the centered skew normal model
  
  vector[N_obs] sig = sigma * sqrt(1 - square(rho[group])) ./ sqrt(1 - 2/pi() * square(rho[group]));
  vector[N_obs] intercept = rho[group].* z[group] - sigma*rho[group] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[group]));
  target +=normal_id_glm_lupdf(y_obs | X_obs, mu + intercept, b, sig);
  
  // vector[N_mis] sigmis = sigma * sqrt(1 - square(rho[groupmis]))./ sqrt(1 - 2/pi() * square(rho[groupmis]));
  // vector[N_mis] interceptmis = rho[groupmis] .* z[groupmis]- sigma*rho[groupmis] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[groupmis]));
  // target +=normal_id_glm_lupdf(ymis | X_mis, mu + interceptmis, b, sigmis);
}

generated quantities {

   vector<lower=0, upper=1>[p] m_ind;
   for (j in 1:p) {
       if (abs(b[j]) > epsilon) {
           m_ind[j] = 1;
       } else {
           m_ind[j] = 0;
       }
   }

  int<lower=1, upper=N_obs+N_mis> j;
  // vector[N_mis] ymis;
  vector[N_obs+N_mis] y_all; // all y's
  real eta_ij;
  real<lower=0> sigma_i;
  //for(i in 1:N_mis){
  for(i in 1:(N_obs+N_mis)){
    j = group_all[i];
    //j = groupmis[i];
    //eta_ij = mu + X_mis[i, ]*b + rho[j] .* z[j]- sigma .* rho[j] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[j]));
    eta_ij = mu + X_all[i, ]*b + rho[j] .* z[j]- sigma .* rho[j] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[j]));
    sigma_i = sigma * sqrt(1 - square(rho[j]))./ sqrt(1 - 2/pi() * square(rho[j]));
    y_all[i] = normal_rng(eta_ij, sigma_i);
  }
}
