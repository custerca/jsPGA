functions{
  // Thermal response curve function
  matrix thermal_response_curve(vector temp, vector CTmax, vector Topt, vector sigma, int N, int K){
    matrix[N, K] TRC;
      for(k in 1:K){
        for(i in 1:N){
          if (temp[i] <= Topt[k]) 
              TRC[i,k] = exp(-((temp[i] - Topt[k])/(2*sigma[k]))^2);
          else if (CTmax[k] >= temp[i] && temp[i] > Topt[k])
              TRC[i,k] = 1 - ((temp[i] - Topt[k])/(Topt[k] - CTmax[k]))^2;
          else
              TRC[i,k] = 0.000001;
        }
      }
      return TRC;
  }
  
} 

// begin data section
data{
  // dimension parameters
  int<lower=0> N; // number of observations
  int<lower=0> J; // number of gears
  int<lower=0> K; // number of species
  int<lower=0> P; // number of abundance parameters
  vector<lower=1>[J] alpha; // Dirichlet prior parameter
  int M; // number of basis functions
  
  // data
  array[N, K] int<lower=0> Y; // total caught
  array[N*K] int<lower=0> y_vec; // vector of Y
  matrix<lower=0>[N, J] E; // effort
  matrix[N, P] X; // design matrix abundance
  vector[N] temp; // design matrix catchability
  matrix[N, M] Psi; // basis function matrix
  
  // priors for TRC...not actually priors anymore...fixed effect simulated through numerical integration
  vector<lower = 0>[K] CTmax_prior;
  vector<lower = 0>[K] Topt_prior;
  vector<lower = 0>[K] sigma; // scaling parameter for thermal response curves
    
}

// begin transformed data section
transformed data{
  vector[K] zeros = rep_vector(0, K);
}


// begin parameters section
parameters {
  row_vector[K] beta_0; // intercept
  matrix[P, K] beta; // abundance parameters
  array[K] simplex[J] theta; // catchability parameters
  real<lower=0> recip_phi; // inverse of negative binomial overdispersion parameter
  
  // optimization of code to estimate spatial basis coefficient's covariance structure
  // https://mc-stan.org/docs/2_19/stan-users-guide/multivariate-hierarchical-priors-section.html
  matrix[K,M] z_A; // standard normal
  cholesky_factor_corr[K] L_Sigma; // Cholesky decomposition of sigma
  vector<lower=0>[K] tau; // species scalings for covariance
}

// begin transformed parameters section
transformed parameters{

  matrix<lower=0,upper=1>[N, K] TRC; // thermal response curve
  matrix<lower=0>[N, K] Etilde; // scaled effort
  real<lower=0> phi; // negative binomial overdispersion parameter
  matrix[M, K] A; // spatial basis coefficients
  A = (diag_pre_multiply(tau,L_Sigma) * z_A)';
  
  phi = 1 / recip_phi; 
  
  TRC = thermal_response_curve(temp, CTmax_prior, Topt_prior, sigma, N, K);
  
  for(k in 1:K){
    Etilde[,k] = E * theta[k] * J;
  }
}

// begin model section
model {
  
  // define model parameters
  matrix[N, K] B0; // global intercept
  matrix[N, K] lambda; // mean intensity function
  

  // Global intercept
  B0 = rep_matrix(beta_0, N);

  // calculate relative abundance lambda - log scale
  lambda = B0 + (X * beta) + (Psi * A) + log(TRC);
  
  // priors 
  beta_0 ~ normal(0, 10); // beta 0 prior
  to_vector(beta) ~ normal(0, 10); // beta prior
  tau ~ cauchy(0,2.5); // prior for species scalings
  to_vector(z_A) ~ std_normal(); 
  L_Sigma ~ lkj_corr_cholesky(2); //https://distribution-explorer.github.io/multivariate_continuous/lkj.html
  recip_phi ~ cauchy(0.,5); // inverse of overdispersion parameter prior
  
  for(k in 1:K){
    theta[k] ~ dirichlet(alpha); // catchability prior
  }

  y_vec ~ neg_binomial_2_log(to_vector(log(Etilde) + lambda),phi); // modeling catch data
  
}

generated quantities{
  corr_matrix[K] Sigma_A; // Correlation matrix Sigma for spatial basis coefficients covariance structure
  cov_matrix[K] tau_Sigma; // Full covariance matrix T*A*T
  tau_Sigma = diag_pre_multiply(tau, L_Sigma) * diag_pre_multiply(tau, L_Sigma)';
  Sigma_A = multiply_lower_tri_self_transpose(L_Sigma); 
}
