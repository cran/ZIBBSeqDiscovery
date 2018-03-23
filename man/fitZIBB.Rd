\name{fitZIBB}
\alias{fitZIBB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The main function to fit ZIBB model
}
\description{
We use zero-inflated beta-binomial (ZIBB) to account for overdispersion and excessive zeros in the microbiome count data. The parameter estimation method is maximum likelihood. Two approaches are proposed to estimate the overdispersion parameters: free approach and constrained approach. For free approach, user does not need to provide initial values for unknown parameters because our program will come up with naive initial values automatically. For constrained approach, user should provide the estimations from free approach as the initial values. This function gives the estimations of the parameters, as well as the corresponding p values, which can be used to test the association between the composition of microbiome counts data and the interested covariates.
}
\usage{
fitZIBB(dataMatrix, X, ziMatrix, mode = "free", gn = 3, 
        betastart = matrix(NA, 0, 0), 
        psi.start = vector(mode = "numeric", length = 0), 
        eta.start = matrix(NA, 0, 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrix}{
The count matrix (m by n, m is the number of OTUs and n is the number of samples).
}
  \item{X}{
The design matrix (n by p, p is the number of covariates) for the count model (e.g., beta-binomial), and intercept is included. The second column is assumed to be the covariate of interest.
}
  \item{ziMatrix}{
The design matrix (n by q, q is the number of covariates) for the zero model, and intercept is included.
}
  \item{mode}{
Indicates which approach is used to estimate overdispersion parameters. mode can be set as either "free" or "constrained".
}
  \item{gn}{
In constrained approach, we use a polynomial with degree of freedom gn to fit the mean-overdispersion relationship. The default value for gn is 3. Note that gn is only valid when mode = "constrained".
}
  \item{betastart}{
Initial values for beta estimation matrix (p by m), where beta are the effects (or coefficients) for the count model. betastart is required only in constrained approach, and it should be assigned as the beta estimation matrix from free approach.
}
  \item{psi.start}{
Initial values for the logit of overdispersion parameters vector (with length m). psi.start is required only in constrained approach, and it should be assigned as the psi estimation vector from free approach.
}
  \item{eta.start}{
Initial values for eta estimation matrix (q by m), where eta are the effects (or coefficients) for the zero model. eta.start is required only in constrained approach, and it should be assigned as the eta estimation matrix from free approach.
}
}
\details{
In this package, we refer the covariate of interest as phenotype (only one phenotype is assumed currently, thouugh we can extend it to include multiple phenotypes), and the phenotype is assumed to corresonding to the second column of the degisn matrix \eqn{X} (note that the first column corresonds to the intercept). Assuming the parameters corresonding to the phenotype are \eqn{\{\beta_{1i}\}_{i=1,\ldots,m}}, this function tests the null hypothesis \eqn{H_0: \beta_{1i} = 0} for each OTU \eqn{i}.
}
\value{
  \item{betahat}{Estimation matrix of beta (p by m) for count model.}
  \item{bvar}{Estimation matrix of the variance of estimated betahat (p by m).}
  \item{p}{The vector (with length of m) of p values corresponding to the phenotype (aka, covariate of interest). In this package, we assume it corresponds to the second column of design matrix X (because intercept is included), though you can always include multiple phenotypes or other covariates. The \eqn{i}'th p value corresponds to the hypothesis test of \eqn{H_0: \beta_{1i} = 0} for OTU \eqn{i}}
  \item{psi}{Estimation vector of the logit of the overdispersion parameters (with length m).}
  \item{zeroCoef}{Estimation matrix of eta (q by m) for zero model.}
  \item{gamma}{Estimation vector of the coefficients in the mean-overdispersion relationship in constrained approach (with length gn+1). So gamma is only available when mode="constrained".}
}
\author{
Tao Hu, Yihui Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Load the data
## data.Y is a count matrix with 100 OTUs and 20 samples randomly selected  
##   from kostic data
data(data.Y)

## set random seed
set.seed(1)

## construct design matrix for count model
## data.X is a 20-by-2 matrix, phenotype is group, and the first 10 samples
##   come from group 1 and the rest samples come from group 2
data.X <- matrix(c(rep(1, 20), rep(0,10), rep(1, 10)), 20, 2)

## construct design matrix for zero model
## data.ziMatrix is a 20-by-2 matrix, the covariate is log of library size
data.ziMatrix <- matrix(1, 20, 2)
data.ziMatrix[, 2] <- log(colSums(data.Y))

## fit ZIBB with free approach
out.free <- fitZIBB(data.Y, data.X, data.ziMatrix, mode = "free")

## fit ZIBB with constrained approach
out.constrained <- fitZIBB(data.Y, data.X, data.ziMatrix, 
                           mode = "constrained", gn = 3, 
                           betastart = out.free$betahat, 
                           psi.start = out.free$psi, 
                           eta.start = out.free$zeroCoef)

## print OTUs which has p values smaller than 0.05
out.constrained$p[which(out.constrained$p < 0.05)]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fitZIBB}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
