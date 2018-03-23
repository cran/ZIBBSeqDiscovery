\name{regression}
\alias{regression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple linear regression
}
\description{
Use simple linear regression to find the initial value of betas in free approach, where beta are the effects/coefficients for the count model.
}
\usage{
regression(y, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The response vector.
}
  \item{x}{
The design matrix (n by p, p is the number of covariates) for the count model (e.g., beta-binomial), and intercept is included.
}
}
\value{
  \item{betahat}{Initial values for beta estimation matrix (p by m, p is the number of involved covariates and m is the number of OTUs).}
  \item{se}{Standard errors for the estimated betas.}
}
\author{
Tao Hu, Yihui Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
