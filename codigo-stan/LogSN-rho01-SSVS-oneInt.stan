data {
  int<lower=1> K;                    // Number of groups
  int<lower=0> N_obs;                // Number of observed data
  int<lower=0> N_mis;                // Number of missing data
  int<lower=1> p;
  vector[N_obs] y_obs;                   // Observed dependent variable for all groups
  matrix[N_obs, p] X_obs;                // Matrix of observed explanatory variables for all groups
  matrix[N_mis, p] X_mis;             // Matrix of observed explanatory variables for all groups with missing y
  array[N_obs] int<lower=1, upper=K> group;
  array[N_mis] int<lower=1, upper=K> group_mis;
  real<lower=0> c;
  real<lower=0> tau;
  real epsilon;  // Threshold for significant coefficient
}

parameters {
  real<lower=0> sigma;               // sd of y
  vector<lower=0, upper=1>[K] rho;     // Correlation coefficient for each group
  real mu;                    // intercept for each region
  //real mu_0;
  //real<lower=0> sigma_mu;
  vector[p] b;                         // Regression coefficients for explanatory variables
  vector<lower=0>[K] z_std;                // Latent variable for each group
  vector[N_mis] y_mis;               // Non-observed dependent variable for all groups
  vector<lower=0, upper=1>[p] pr;        // probabilidad de inclusion
}

transformed parameters{
  vector<lower=0>[K] z = z_std .* sigma ./ sqrt(1 - 2/pi() * square(rho));
}


model {
  // Reference prior for rho
  target += sum( 0.5 * log1p(square(rho)) - log1m(square(rho)) );
   
  // Flat prior for mu
  mu ~ normal(0, 100);
  //mu ~ normal(mu_0, sigma_mu);
  //mu_0 ~ normal(0, 10);
  //sigma_mu ~ inv_gamma(0.001, 0.001);

  // Prior for sigma
  // target += -sigma;   // log of 1/sigma
  sigma ~ inv_gamma(0.001, 0.001); // almost Jeffrey prior for sigma
  
  // Prior for SSVS: seleccion of b and his inclusion probability pr.  
  pr ~ beta(0.5, 0.5);
  for (i in 1:p) {
    target += log_mix(pr[i],
                      normal_lpdf(b[i] | 0, tau * c),
                      normal_lpdf(b[i] | 0, tau));
  }

  // Latent z_std
  z_std ~ normal(0, 1);  // Truncated normal prior for latent variables
  
  // Likelihood for observed y given z and b for each group
  
  // `sig` means the square root of variance parameter in the centered skew normal model
  
  vector[N_obs] sig = sigma * sqrt(1 - square(rho[group])) ./ sqrt(1 - 2/pi() * square(rho[group]));
  vector[N_obs] intercept = rho[group].* z[group] - sigma*rho[group] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[group]));
  target +=normal_id_glm_lupdf(y_obs | X_obs, mu + intercept, b, sig);
  
  vector[N_mis] sigmis = sigma * sqrt(1 - square(rho[group_mis]))./ sqrt(1 - 2/pi() * square(rho[group_mis]));
  vector[N_mis] interceptmis = rho[group_mis] .* z[group_mis]- sigma*rho[group_mis] .* sqrt(2/pi()) ./ sqrt(1 - 2/pi() * square(rho[group_mis]));
  target +=normal_id_glm_lupdf(y_mis | X_mis, mu + interceptmis, b, sigmis);
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
}
