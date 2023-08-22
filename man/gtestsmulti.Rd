\name{gtestsmulti}
\alias{gtestsmulti}
\title{New graph-based multi-sample tests}
\description{This function provides graph-based multi-sample tests.}
\usage{
gtestsmulti(E, data_list, perm=0)
}
\arguments{
  \item{E}{The edge matrix for the similarity graph. Each row contains the node indices of an edge.}
  \item{data_list}{The list of multivariate matrices corresponding to the K different classes. The length of the list is K. Each element of the list is a matrix containing observations as the rows and features as the columns.}
  \item{perm}{The number of permutations performed to calculate the p-value of the test. The default value is 0, which means the permutation is not performed and only approximated p-value based on the asymptotic theory is provided. Doing permutation could be time consuming, so be cautious if you want to set this value to be larger than 10,000.}
}

\value{
  Returns a list \code{teststat} with each test statistic value and a list \code{pval} with p-values of the tests. See below for more details.
  \item{S}{The value of the test statistic \eqn{S}.}
  \item{S_A}{The value of the test statistic \eqn{S^{A}}.}
  \item{S_appr}{The approximated p-value of \eqn{S} based on asymptotic theory with a Bonferroni procedure.}
  \item{S_A_appr}{The approximated p-value of \eqn{S^{A}} based on asymptotic theory.}
  \item{S_perm}{The permutation p-value of \eqn{S} when argument `perm' is positive.}
  \item{S_A_perm}{The permutation p-value of \eqn{S^{A}} when argument `perm' is positive.}
}
\seealso{
  \code{\link{gTestsMulti-package}}
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

