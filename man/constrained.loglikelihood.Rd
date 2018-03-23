\name{constrained.loglikelihood}
\alias{constrained.loglikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Define the objective function in optimization procedure for estimating parameters with constrained approach
}
\description{
The objective function is the negative of log likelihood function.
}
\usage{
constrained.loglikelihood(para, X, Y.col, coeff, Y.c, ziMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{para}{
Vector of optimized parameters with length p+q, where p is the number of covariates for count model (e.g., beta-binomial), q is the number of covariates for zero model. The first p elements are betas which are the effects/coefficients for the count model. The last q elements are etas which are the effects/coefficients for the zero model.
}
  \item{X}{
The design matrix (n by p, p is the number of covariates) for the count model (e.g., beta-binomial), and intercept is included.
}
  \item{Y.col}{
Vector of counts corresponding to an OTU, with length n.
}
  \item{coeff}{
Vector of coefficients in the polynomial mean-overdispersion relationship in constrained approach.
}
  \item{Y.c}{
Vector of library size with length n.
}
  \item{ziMatrix}{
The design matrix (n by q) for the zero model, and intercept is included.
}
}
\author{
Tao Hu, Yihui Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
