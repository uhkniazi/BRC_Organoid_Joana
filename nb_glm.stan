data {
    int<lower=1> Ntotal ; // total number of observations
    int<lower=1> Ncols; // total number of levels for predictor
    matrix[Ntotal, Ncols] modmatrix; // model matrix
    int y[Ntotal]; // response variable
    real betaShape; // hyperparameters for the hierarchichal standard deviation parameter
    real betaRate;
  }
  transformed data {
    ## constants for parameter distributions
    
  }
  parameters {
    # parameters to track
    vector[Ncols] beta; // coefficients of the model
    real<lower=1, upper=1000> iSize; // size parameter for the nb distribution 
    real<lower=0.1> betaSigma; // standard deviation parameter for the joint prior for betas
  }
  transformed parameters {
    vector[Ntotal] iFitted; // the fitted value or mu
    iFitted = exp(modmatrix * beta);
  }
  model {
    // prior for standard deviations of coefficients
    betaSigma ~ gamma(betaShape, betaRate);
    // 2 betas
    beta[1] ~ normal( 0 , betaSigma ) ;
    beta[2] ~ normal( 0 , betaSigma ) ;
    y ~ neg_binomial_2(iFitted, iSize);  
  }
