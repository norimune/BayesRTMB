\name{make_glmer_re_terms}
\alias{make_glmer_re_terms}
\title{Prepare GLMM Formula Components}
\usage{
make_glmer_re_terms(
  formula,
  data,
  family = "gaussian",
  resid_group = NULL,
  resid_time = NULL,
  within = NULL,
  factors = NULL
)
}
\arguments{
\item{formula}{lme4-style formula.}

\item{data}{Data frame.}

\item{family}{Character string of the distribution family.}

\item{resid_group}{Optional residual-correlation grouping variable.}

\item{resid_time}{Optional residual-correlation time variable.}

\item{within}{Optional list for wide-to-long conversion when the response is
written as \code{cbind(...)}.}

\item{factors}{Optional character vector of variables to treat as factors.}
}
\value{
A list containing \code{Y}, \code{X}, \code{trials}, \code{offset},
\code{N}, fixed-effect metadata, and random-effect terms.
}
\description{
Builds the response, fixed-effect model matrix, and lme4-style random-effect
design components used by \code{rtmb_glmer()}. This helper is the user-facing
entry point for reproducing the wrapper model with \code{rtmb_model()}.
}
