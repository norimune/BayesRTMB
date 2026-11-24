# Log-CDF and Log-CCDF Functions

Log cumulative distribution functions and log complementary cumulative
distribution functions for univariate distributions used in
\`rtmb_code()\`. These functions are useful for censored and
interval-censored likelihoods and are compatible with RTMB automatic
differentiation.

## Usage

``` r
normal_lcdf(x, mean, sd, sum = TRUE)

normal_lccdf(x, mean, sd, sum = TRUE)

lognormal_lcdf(x, meanlog, sdlog, sum = TRUE)

lognormal_lccdf(x, meanlog, sdlog, sum = TRUE)

exp_mod_normal_lcdf(x, mu, sigma, lambda, sum = TRUE)

exp_mod_normal_lccdf(x, mu, sigma, lambda, sum = TRUE)

exponential_lcdf(x, rate, sum = TRUE)

exponential_lccdf(x, rate, sum = TRUE)

half_normal_lccdf(x, sd, sum = TRUE)

half_normal_lcdf(x, sd, sum = TRUE)

beta_lcdf(x, a, b, sum = TRUE)

beta_lccdf(x, a, b, sum = TRUE)

gamma_lcdf(x, shape, rate, sum = TRUE)

gamma_lccdf(x, shape, rate, sum = TRUE)

inverse_gamma_lcdf(x, shape, scale, sum = TRUE)

inverse_gamma_lccdf(x, shape, scale, sum = TRUE)

cauchy_lcdf(x, location, scale, sum = TRUE)

cauchy_lccdf(x, location, scale, sum = TRUE)

student_t_lcdf(x, df, mu = 0, sigma = 1, sum = TRUE)

student_t_lccdf(x, df, mu = 0, sigma = 1, sum = TRUE)

laplace_lcdf(x, location = 0, scale = 1, sum = TRUE)

laplace_lccdf(x, location = 0, scale = 1, sum = TRUE)

logit_normal_lcdf(x, mu, sd, sum = TRUE)

logit_normal_lccdf(x, mu, sd, sum = TRUE)

weibull_lcdf(x, shape, scale, sum = TRUE)

weibull_lccdf(x, shape, scale, sum = TRUE)

uniform_lcdf(x, a, b, sum = TRUE)

uniform_lccdf(x, a, b, sum = TRUE)

normal_mixture_lcdf(x, pi_w, mean, sd, sum = TRUE)

normal_mixture_lccdf(x, pi_w, mean, sd, sum = TRUE)

mixture_lcdf(x, pi_w, lcdf_list, sum = TRUE)

mixture_lccdf(x, pi_w, lccdf_list, sum = TRUE)

bernoulli_lcdf(x, prob, sum = TRUE)

bernoulli_lccdf(x, prob, sum = TRUE)

binomial_lcdf(x, size, prob, sum = TRUE)

binomial_lccdf(x, size, prob, sum = TRUE)

bernoulli_logit_lcdf(x, eta, sum = TRUE)

bernoulli_logit_lccdf(x, eta, sum = TRUE)

binomial_logit_lcdf(x, size, eta, sum = TRUE)

binomial_logit_lccdf(x, size, eta, sum = TRUE)

poisson_lcdf(x, mean, sum = TRUE)

poisson_lccdf(x, mean, sum = TRUE)

neg_binomial_lcdf(x, size, prob, sum = TRUE)

neg_binomial_lccdf(x, size, prob, sum = TRUE)

neg_binomial_2_lcdf(x, mu, size, sum = TRUE)

neg_binomial_2_lccdf(x, mu, size, sum = TRUE)

ordered_logistic_lcdf(x, eta, cutpoints, sum = TRUE)

ordered_logistic_lccdf(x, eta, cutpoints, sum = TRUE)

sequential_logistic_lccdf(x, eta, cutpoints, sum = TRUE)

sequential_logistic_lcdf(x, eta, cutpoints, sum = TRUE)
```

## Arguments

- x:

  Quantile or vector of quantiles.

- mean:

  Mean parameter. For Poisson functions, the rate parameter.

- sd:

  Positive standard deviation.

- sum:

  Logical; return the sum when \`TRUE\`, otherwise return element-wise
  values.

- meanlog, sdlog:

  Mean and standard deviation on the log scale.

- mu:

  Location or mean parameter.

- sigma:

  Positive scale parameter.

- lambda:

  Positive exponential rate.

- rate:

  Positive rate parameter.

- a, b:

  Positive shape parameters, or lower and upper bounds for the uniform
  distribution.

- shape:

  Positive shape parameter.

- scale:

  Positive scale parameter.

- location:

  Location parameter.

- df:

  Positive degrees of freedom.

- pi_w:

  Mixture weights that sum to one.

- lcdf_list, lccdf_list:

  Lists of component-wise log-CDF or log-CCDF vectors, respectively.

- prob:

  Probability parameter.

- size:

  Number of trials or positive dispersion parameter.

- eta:

  Linear predictor.

- cutpoints:

  Ordered cutpoints.

## Value

A scalar sum of log probabilities when \`sum = TRUE\`, or a vector of
element-wise log probabilities when \`sum = FALSE\`.

## Details

Each function returns the sum of the element-wise log probabilities by
default. Set \`sum = FALSE\` to return one value per observation. The
naming follows Stan: \`\_lcdf\` computes \`log(P\[X \<= x\])\`, and
\`\_lccdf\` computes \`log(P\[X \> x\])\`. For discrete distributions
this distinction is important.

Supported continuous distributions are normal, lognormal, exponentially
modified normal, exponential, half-normal, beta, gamma, inverse-gamma,
Cauchy, Student-t, Laplace, logit-normal, Weibull, and uniform, as well
as normal and generic mixtures. Supported discrete distributions are
Bernoulli, binomial, Bernoulli-logit, binomial-logit, Poisson, negative
binomial, negative-binomial-2, ordered logistic, and sequential
logistic.

Multivariate and matrix distributions, nominal categorical
distributions, and diffusion models do not have a single generally
applicable scalar CDF and are therefore not included.

## See also

Other distributions:
[`distributions`](https://norimune.github.io/BayesRTMB/reference/distributions.md)
