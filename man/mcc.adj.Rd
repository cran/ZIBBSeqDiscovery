\name{mcc.adj}
\alias{mcc.adj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Using MCC method to replace NAs in the p values
}
\description{
When fitting the ZIBB model, some parameter estimations may fail due to numerical issues. In that case, a NA will be given as the corresponding p value. Here, a Moment Corrected Correlation (MCC) approach is employed to replace the NAs in the p values.
}
\usage{
mcc.adj(out.fitZIBB, dataMatrix, X, ziMatrix, K = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out.fitZIBB}{
  The output from function \code{\link{fitZIBB}}.
  }
  \item{dataMatrix}{
The count matrix (m by n, m is the number of OTUs and n is the number of samples).
}
  \item{X}{
The design matrix (n by p, p is the number of covariates) for the count model (e.g., beta-binomial), and intercept is included. The second column is assumed to be the covariate of interest.
}
  \item{ziMatrix}{
The design matrix (n by q, q is the number of covariates) for the zero model, and intercept is included.
}
  \item{K}{
  Divide covariate in ziMatrix (second colunm in default) into K stratum, under the requirement of MCC approach. The default value of K is 4.
  }
}
\value{
  The output has the exact same format as function fitZIBB, with corrected p values.
}
\author{
Tao Hu, Yihui Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{
Zhou, Y. H., & Wright, F. A. (2015). Hypothesis testing at the extremes: fast and robust association for high-throughput data. Biostatistics, 16(3), 611-625.
}

\seealso{
\code{\link{fitZIBB}}
}

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

## count how many NAs in the p values
sum(is.na(out.free$p))

## MCC adjustment
out.free.mcc <- mcc.adj(out.free, data.Y, data.X, data.ziMatrix, K=4)

## count how many NAs in the p values after MCC adjustment
sum(is.na(out.free.mcc$p))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
