functions {
  real weib_lpdf(vector Y, int N, int P, matrix X, vector trun, real b, real beta0, vector betas) {
    real lp = 0;
    real w;
    for (i in 1:N) {
      w = (Y[i] - beta0 - dot_product(X[i, ], betas)) / b;
      if (trun[i])
        lp -= exp(w);
      else
        lp += w - exp(w);
    }
    lp -= (N - sum(trun)) * log(b);
    return lp;
  }

  vector weib_rng(int N, matrix X, real b, real beta0, vector betas) {
    vector[N] rns;
    for (i in 1:N) {
      real u = uniform_rng(0, 1);
      real alpha = 1 / b;
      real sigma = exp(beta0 + dot_product(X[i, ], betas));
      rns[i] = (-log(1 - u)) ^ (1 / alpha) * sigma;
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
  matrix[N, P == 0 ? 1 : P] X;

  int<lower = 0> N_new;
  matrix[N_new, P == 0 ? 1 : P] X_new;
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
  real<lower = 0> b;

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
    b ~ exponential(1);
  }
  else {
    thetas ~ uniform(-3, 3);
    b ~ uniform(0, 3);
  }

  Y ~ weib(N, P, QR ? Q : X, trun, b, beta0, thetas);

  if (generate_samples)
    Y_new ~ weib(N_new, P, X_new, rep_vector(0, N_new), b, beta0, betas);
}

generated quantities {
  vector[N_new] T_new;
  if (generate_samples)
    // T_new = weib_rng(N_new, X_new, b, beta0, betas);
    //T_new = exp(Y_new);
    for (i in 1:N_new)
      T_new[i] = weibull_rng(1 / b, exp(beta0 + dot_product(X[i, ], betas)));
}

// data {
//   int<lower = 0, upper = 1> generate_samples;
//   int<lower = 0, upper = 1> bayesian;
//   int<lower = 0, upper = 1> QR;
//   int<lower = 0, upper = 1> norm;
// 
//   int<lower = 0> N;
//   int<lower = 0> P;
//   vector<lower = 0>[N] T;
//   vector<lower = 0, upper = 1>[N] trun;
//   matrix[N, P == 0 ? 1 : P] X;
//   
//   int<lower = 0> N_new;
//   matrix[N_new, P == 0 ? 1 : P] X_new;
// }
// 
// transformed data {
//   vector[N] Y = log(T);
//   
//   matrix[N, P == 0 ? 1 : P] Q = qr_thin_Q(X);
//   matrix[P == 0 ? 1 : P, P == 0 ? 1 : P] R_inv = inverse(qr_thin_R(X));
// //  matrix[N, P] Q = qr_thin_Q(X) * sqrt(N - 1);
// //  matrix[P, P] R_inv = inverse(qr_thin_R(X) / sqrt(N - 1));
// }
// 
// parameters {
//   real beta0;
//   vector[P == 0 ? 1 : P] thetas;
//   real<lower = 0> alpha;
//   
//   vector[generate_samples ? N_new : 0] Y_new;
// }
// 
// transformed parameters {
//   vector[P == 0 ? 1 : P] betas = QR ? R_inv * thetas : thetas;
//   
// }
// 
// model {
//   if (bayesian) {
//     if (norm) {
//       thetas ~ cauchy(0, 2.5);
//     }
//     else {
//       thetas ~ normal(0, 1);
//     }
//     alpha ~ exponential(1);
//   }
//   else {
//     thetas ~ uniform(-3, 3);
//     alpha ~ uniform(0, 3);
//   }
//   
//   for (i in 1:N)
//     if (trun[i])
//       target += weibull_lccdf(T[i]| alpha, exp(beta0 + dot_product((QR ? Q : X)[i, ], thetas)));
//     else
//       target += weibull_lpdf(T[i]| alpha, exp(beta0 + dot_product((QR ? Q : X)[i, ], thetas)));
// 
//   if (generate_samples)
//     for (i in 1:N_new)
//       Y_new[i] ~ weibull(alpha, exp(beta0 + dot_product(X[i], betas)));
// }
// 
// generated quantities {
//   vector[N_new] T_new;
//   if (generate_samples)
//     for (i in 1:N_new)
//       T_new[i] = weibull_rng(alpha, exp(beta0 + dot_product(X[i, ], betas)));
// }
