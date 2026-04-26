# BayesRTMB: Bayesian Modeling with RTMB

The \`BayesRTMB\` package provides an intuitive, R-based Bayesian
modeling environment using the \`RTMB\` package as its computational
backend. It allows users to perform Maximum A Posteriori (MAP)
estimation, Markov Chain Monte Carlo (MCMC) sampling, and Automatic
Differentiation Variational Inference (ADVI) for hierarchical and
complex non-linear models through a unified and user-friendly interface.

## Details

The standard workflow of the package consists of three main steps:

1\. \*\*Model Definition\*\*: Use \[rtmb_code()\] to define the model
structure (likelihood, prior distributions, transformations, etc.) using
standard R syntax. 2. \*\*Model Instantiation\*\*: Pass the data and the
generated model code to \[rtmb_model()\] to create an \`RTMB_Model\` R6
object. 3. \*\*Estimation\*\*: Call the appropriate method on the model
object: \* \`\$optimize()\`: For MAP estimation and Laplace
approximation. \* \`\$sample()\`: For MCMC sampling using the No-U-Turn
Sampler (NUTS). \* \`\$variational()\`: For approximate Bayesian
inference using ADVI.

For more detailed examples and tutorials, please refer to the package
vignettes.

## See also

Useful links: \*
\[https://github.com/norimune/BayesRTMB\](https://github.com/norimune/BayesRTMB)
\* Report bugs at
\[https://github.com/norimune/BayesRTMB/issues\](https://github.com/norimune/BayesRTMB/issues)
