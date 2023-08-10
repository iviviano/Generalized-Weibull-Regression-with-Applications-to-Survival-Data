functions {
  real exp_weib_log(vector Y, int N, int P, matrix X, vector trun, real kappa, real theta, real beta0, vector betas) {
    real lp = 0;
    real w;
    real a;
    real r = N - sum(trun);
    for (i in 1:N) {
      w = (Y[i] - beta0 - dot_product(X[i, ], betas)) / kappa;
      a = 1 - exp(-exp(w));
      if (trun[i]) {
        lp += log(1 - a ^ theta);
      }
      else
        lp += (theta - 1) * log(a) + w - exp(w);
    }
    lp += r * log(theta) - r * log(kappa);
    return lp;
  }
  
  vector exp_weib_rng(int N, int P, matrix X, real kappa, real theta, real beta0, vector betas) {
    vector[N] rns;
    for (i in 1:N) {
      real u = uniform_rng(0, 1);
      real sigma = exp(beta0 + dot_product(X[i, ], betas));
      real alpha = 1 / kappa;
      rns[i] = (-log(1 - u ^ (1 / theta))) ^ (1 / alpha) * sigma;
    }
    return rns;
  }
}

data {
  int<lower = 0, upper = 1> generate_samples;
  int<lower = 0, upper = 1> bayesian;
  int<lower = 0, upper = 1> QR;
  int<lower = 0, upper = 1> norm;

  int<lower = 0> N;
  int<lower = 0> P;
  vector<lower = 0>[N] T;
  vector<lower = 0, upper = 1>[N] trun;
  matrix[N, (P == 0 ? 1 : P)] X;
  
  int<lower = 0> N_new;
  matrix[N_new, (P == 0 ? 1 : P)] X_new;
}

transformed data {
  vector[N] Y = log(T);
  
  matrix[N, P == 0 ? 1 : P] Q = qr_thin_Q(X);
  matrix[P == 0 ? 1 : P, P == 0 ? 1 : P] R_inv = inverse(qr_thin_R(X));
//  matrix[N, P] Q = qr_thin_Q(X) * sqrt(N - 1);
//  matrix[P, P] R_inv = inverse(qr_thin_R(X) / sqrt(N - 1));
}

parameters {
  real beta0;
  vector[P == 0 ? 1 : P] thetas;
  real<lower = 0> kappa;
  real<lower = 0> theta;
  
  vector[generate_samples ? N_new : 0] Y_new;
}

transformed parameters {
  vector[P == 0 ? 1 : P] betas = QR ? R_inv * thetas : thetas;
  
}

model {
  if (bayesian) {
    if (norm) {
      thetas ~ cauchy(0, 2.5);
    }
    else {
      thetas ~ normal(0, 1);
    }
    kappa ~ exponential(1);
    theta ~ exponential(1); 
  }
  else {
    thetas ~ uniform(-3, 3);
    kappa ~ uniform(0, 3);
    theta ~ uniform(1, 5);
  }
  
  Y ~ exp_weib(N, P, QR ? Q : X, trun, kappa, theta, beta0, thetas);

  if (generate_samples)
    Y_new ~ exp_weib(N_new, P, X_new, rep_vector(0, N_new), kappa, theta, beta0, betas);
}

generated quantities {
  vector[N_new] T_new;
  if (generate_samples)
    T_new = exp_weib_rng(N_new, P, X_new, kappa, theta, beta0, betas);
}
