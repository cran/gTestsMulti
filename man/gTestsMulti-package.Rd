\name{gTestsMulti}
\alias{gTestsMulti}
\title{New graph-based multi-sample tests}
\description{This package can be used to determine whether multiple samples are from the same distribution.}
\author{
Hoseung Song and Hao Chen

Maintainer: Hoseung Song (hosong@ucdavis.edu)
}
\references{
	Song, H. and Chen, H. (2022). New graph-based multi-sample tests for high-dimensional and non- Euclidean data. arXiv:2205.13787
}

\seealso{
  \code{\link{gtestsmulti}}
}
\examples{

## Mean difference in Gaussian distribution.
d = 50
mu = 0.2
sam = 50

set.seed(500)
X1 = matrix(rnorm(d*sam), sam)
X2 = matrix(rnorm(d*sam,mu), sam)
X3 = matrix(rnorm(d*sam,2*mu), sam)

data_list = list(X1, X2, X3)

# We use 'mstree' in 'ade4' package to construct the minimum spanning tree.
require(ade4)
x = rbind(X1, X2, X3)
E = mstree(dist(x))

a = gtestsmulti(E, data_list, perm = 1000)
# output results based on the permutation and the asymptotic results
# the test statistic values can be found in a$teststat
# p-values can be found in a$pval

}
