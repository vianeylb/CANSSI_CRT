data {
  int<lower=0> N;
  int<lower=0> K;
  int y[N];
}

parameters {
  simplex[K] theta[K];
  // real lambda[K];
  positive_ordered[K] lambda;
}

model {
  // priors
  target+= gamma_lpdf(lambda[1] | 3, 1);
  target+= gamma_lpdf(lambda[2] | 10, 1);
  // forward algorithm
  {
  real acc[K];
  real gamma[N, K];
  for (k in 1:K)
    gamma[1, k] = poisson_lpmf(y[1] | lambda[k]);
  for (t in 2:N) {
    for (k in 1:K) {
      for (j in 1:K)
        acc[j] = gamma[t-1, j] + log(theta[j, k]) + poisson_lpmf(y[t] | lambda[k]);
      gamma[t, k] = log_sum_exp(acc);
    }
  }
  target += log_sum_exp(gamma[N]);
  }
}

// generated quantities {
//   int<lower=1,upper=K> z_star[N];
//   real log_p_z_star;
//   {
//     int back_ptr[N, K];
//     real best_logp[N, K];
//     for (k in 1:K)
//       best_logp[1, k] = poisson_lpmf(y[1] | lambda[k]);
//     for (t in 2:N) {
//       for (k in 1:K) {
//         best_logp[t, k] = negative_infinity();
//         for (j in 1:K) {
//           real logp;
//           logp = best_logp[t-1, j] + log(theta[j, k]) + poisson_lpmf(y[t] | lambda[k]);
//           if (logp > best_logp[t, k]) {
//             back_ptr[t, k] = j;
//             best_logp[t, k] = logp;
//           }
//         }
//       }
//     }
//     log_p_z_star = max(best_logp[N]);
//     for (k in 1:K)
//       if (best_logp[N, k] == log_p_z_star)
//         z_star[N] = k;
//     for (t in 1:(N - 1))
//       z_star[N - t] = back_ptr[N - t + 1, z_star[N - t + 1]];
//   }
// }