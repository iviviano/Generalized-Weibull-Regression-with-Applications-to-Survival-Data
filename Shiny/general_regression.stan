functions {
  real gen_weib_lpdf(vector Y, int N, int P, matrix X, vector trun, real alpha, real lambda, real beta0, vector betas) {
    real lp = 0;
    real w;
    real power_cdf = 1 / lambda;
    real power_pdf =  1 / lambda - 1;
    for (i in 1:N) {
      w = (Y[i] - beta0 - dot_product(X[i, ], betas)) / alpha;
      if (trun[i])
        lp += power_cdf * log(1 - lambda * exp(w));
      else
        lp += w + power_pdf * log(1 - lambda * exp(w));
    }
    lp -= (N - sum(trun)) * log(alpha);
    return lp;
  }
  
  vector gen_weib_rng(int N, matrix X, real alpha, real lambda, real beta0, vector betas) {
    vector[N] rns;
    for (i in 1:N) {
      real u = uniform_rng(0, 1);
      real sigma = exp(beta0 + dot_product(X[i, ], betas));
      rns[i] = (((1 - (1 - u) ^ lambda) / lambda) ^ alpha) * sigma;
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
  matrix[N, P==0 ? 1 : P] X;
  
  int<lower = 0> N_new;
  matrix[N_new, P==0 ? 1 : P] X_new;
}

transformed data {
  vector[N] Y = log(T);
  
  matrix[N, P==0 ? 1 : P] Q = qr_thin_Q(X) * sqrt(N - 1);
  matrix[P==0 ? 1 : P, P==0 ? 1 : P] R_inv = inverse(qr_thin_R(X) / sqrt(N - 1));
}

parameters {
  real beta0;
  vector[P==0 ? 1 : P] thetas;
  real<lower = 0> alpha;
  real lambda;

  vector[generate_samples ? N_new : 0] Y_new;
}

transformed parameters {
  vector[P==0 ? 1 : P] betas = QR ? R_inv * thetas : thetas;
  
  // if (QR)
  //   betas = R_inv * thetas;
  // else
  //   betas = thetas;
}

model {
  if (bayesian) {
    if (norm) {
      thetas ~ cauchy(0, 2.5);
    }
    else {
      thetas ~ normal(0, 1);
    }
    alpha ~ exponential(1);
    lambda ~ normal(0, 1);
  }
  else {
    thetas ~ uniform(-3, 3);
    alpha ~ uniform(0, 2);
    lambda ~ uniform(-1, 1);
  }
  
  
  Y ~ gen_weib(N, P, QR ? Q : X, trun, alpha, lambda, beta0, thetas);

  if (generate_samples)
    Y_new ~ gen_weib(N_new, P, X_new, rep_vector(0, N_new), alpha, lambda, beta0, betas);
}

generated quantities {
  vector[N_new] T_new;
  if (generate_samples) {
    T_new = gen_weib_rng(N_new, X_new, alpha, lambda, beta0, betas);
  }
}

// functions {
//   real gen_weib_lpdf(vector Y, int N, int P, matrix X, vector trun, real alpha, real lambda, real beta0, vector betas) {
//     real lp = 0;
//     real w;
//     real power_cdf = 1 / lambda;
// 
//     real power_pdf =  1 / lambda - 1;
//     for (i in 1:N) {
//       w = (Y[i] - beta0 - dot_product(X[i, ], betas)) / alpha;
//       if (trun[i])
//         lp += power_cdf * log(1 - lambda * exp(w));
//       else
//         lp += w + power_pdf * log(1 - lambda * exp(w));
//     }
//     lp -= (N - sum(trun)) * log(alpha);
//     return lp;
//   }
// }
// 
// data {
//   int<lower = 0, upper = 1> test;
// 
//   int<lower = 0> N;
//   int<lower = 0> P;
//   vector<lower = 0>[N] T;
//   vector<lower = 0, upper = 1>[N] trun;
//   matrix[N, P] X;
// 
//   int<lower = 0> N_new;
//   matrix[N_new, P] X_new;
// }
// 
// transformed data {
//   vector[N] Y = log(T);
// 
//   matrix[N, P] Q = qr_thin_Q(X);
//   matrix[P, P] R_inv = inverse(qr_thin_R(X));
// }
// 
// parameters {
//   real beta0;
//   vector[P] thetas;
//   real<lower = 0> kappa;
//   real<lower = 0> theta;
// 
//   // vector<lower = 0>[N_new] T_new;
//   vector[N_new] Y_new;
// }
// 
// // transformed parameters {
// //   vector[N_new] Y_new = log(T_new);
// // }
// 
// model {
//   thetas ~ normal(0, 1);
//   kappa ~ exponential(1);
//   theta ~ exponential(1);
//   Y ~ gen_weib(N, P, Q, trun, kappa, theta, beta0, thetas);
// 
//   if (test)
//     Y_new ~ gen_weib(N_new, P, X_new, rep_vector(0, N_new), kappa, theta, beta0, R_inv * thetas);
// }
// 
// generated quantities {
//   vector[P] betas = R_inv * thetas;
// }
