\name{constrained.estimate}
\alias{constrained.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate parameters with constrained approach
}
\description{
Estimate unknown parameters with constrained approach.
}
\usage{
constrained.estimate(m, p, q, n, betastart, bvarstart, psi.start, 
                     eta.start, gamma.start, Y, X, Y.c, ziMatrix, 
                     gn = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
Number of OTUs.
}
  \item{p}{
Number of covariates for count model (e.g., beta-binomial).
}
  \item{q}{
Number of covariates for zero model.
}
  \item{n}{
Number of samples.
}
  \item{betastart}{
Matrix of estimated betas, which are the effects/coefficients for the count model, with dimension p by m. It is used as initial values for the optimization procedure to estimate betas.
}
  \item{bvarstart}{
Matrix of variance of estimated betas with dimension p by m.
}
  \item{psi.start}{
Estimated vector of logit of overdispersion parameters with length m. And psi.start will be used as initial values for the optimization procedure to estimate psi.
}
  \item{eta.start}{
Matrix of estimated etas, which are the effects/coefficients for the zero model, with dimension q by m. It is used as initial values for the optimization procedure to estimate etas.
}
  \item{gamma.start}{
Estimation vector of the coefficients in the polynomial mean-overdispersion relationship in constrained approach.
}
  \item{Y}{
Count matrix with dimension n by m.
}
  \item{X}{
The design matrix (n by p, p is the number of covariates) for the count model (e.g., beta-binomial), and intercept is included.
}
  \item{Y.c}{
Vector of library size with length n.
}
  \item{ziMatrix}{
The design matrix (n by q) for the zero model, and intercept is included.
}
  \item{gn}{
We use a polynomial with degree of freedom gn to fit the mean-overdispersion relationship.
}
}
\value{
  \item{betahat}{Estimation matrix of beta (p by m).}
  \item{bvar}{Estimation matrix of the variance of estimated betahat (p by m).}
  \item{psi}{Estimation vector of the logit of the overdispersion parameters (with length m).}
  \item{eta}{Estimation matrix of eta (q by m).}
}
\author{
Tao Hu, Yihui Zhou
}
