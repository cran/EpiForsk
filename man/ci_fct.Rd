% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kija_fct_confint.R
\name{ci_fct}
\alias{ci_fct}
\title{solve optimization problem for  CI bounds}
\usage{
ci_fct(i, f, xtx_red, beta_hat, which_parm, level, n_grid, k)
}
\arguments{
\item{i}{An index for the point at which to solve for confidence limits.}

\item{f}{A function taking the parameter vector as its single argument, and
returning a numeric vector.}

\item{xtx_red}{Reduced form of matrix \emph{X^TX}.}

\item{beta_hat}{Vector of parameter estimates.}

\item{which_parm}{Vector indicating which parameters to include.}

\item{level}{The confidence level required.}

\item{n_grid}{Either \code{NULL} or an integer vector of length 1 or the number of
\code{TRUE}/indices in which_parm. Specifies the number of grid points in each
dimension of a grid with endpoints defined by len. If \code{NULL} or \code{0L}, will
instead sample k points uniformly on a sphere.}

\item{k}{If n_grid is \code{NULL} or \code{0L}, the number of points to sample
uniformly from a sphere.}
}
\value{
One row tibble with estimate and confidence limits.
}
\description{
solve optimization problem for each coordinate of f, to obtain the uniform
limit.
}
\examples{
1+1
}
