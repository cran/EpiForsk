% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kija_rate_omnibus_test.R
\name{RATEOmnibusTest}
\alias{RATEOmnibusTest}
\title{RATE based omnibus test of heterogeneity}
\usage{
RATEOmnibusTest(
  forest,
  level = 0.95,
  target = c("AUTOC", "QINI"),
  q = seq(0.1, 1, 0.1),
  R = 500,
  num.threads = 1,
  seed = NULL,
  honesty = TRUE,
  stabilize.splits = TRUE,
  ...
)
}
\arguments{
\item{forest}{An object of class \code{causal_forest}, as returned by
\link[grf]{causal_forest}(), with binary treatment.}

\item{level}{numeric, level of RATE confidence interval.}

\item{target}{character, see \link[grf]{rank_average_treatment_effect}.}

\item{q}{numeric, see \link[grf]{rank_average_treatment_effect}.}

\item{R}{integer, see \link[grf]{rank_average_treatment_effect}}

\item{num.threads}{passed to \link[grf]{causal_forest}. Number of threads used in
training. Default value is 1.}

\item{seed}{numeric, either length 1, in which case the same seed is used for
both new forests, or length 2, to train each forest with a different seed.
Default is \code{NULL}, in which case two seeds are randomly sampled.}

\item{honesty}{Boolean, \code{TRUE} if forest was trained using honesty. Otherwise \code{FALSE}.
Argument controls if honesty is used to train the new forests on the random
half-samples, so misspecification will lead to invalid results. Default is
\code{TRUE}, the default in \link[grf]{causal_forest}.}

\item{stabilize.splits}{Boolean, \code{TRUE} if forest was trained taking treatment into
account when determining the imbalance of a split. Otherwise \code{FALSE}.
Argument controls if treatment is taken into account when determining the
imbalance of a split during training of the new forests on the random
half-samples, so misspecification will lead to invalid results. Default is
\code{TRUE}, the default in \link[grf]{causal_forest}.}

\item{...}{additional arguments for \link[grf]{causal_forest}. By default, the
arguments used by forest will be used to train new forests on the random
half-samples. Arguments provided through \code{...} will override these. Note that
sample.weights and clusters are passed to both \link[grf]{causal_forest} and
\link[grf]{rank_average_treatment_effect.fit}.}
}
\value{
A list of class \code{rank_average_treatment_effect} with elements
\itemize{
\item estimate: the RATE estimate.
\item std.err: bootstrapped standard error of RATE.
\item target: the type of estimate.
\item TOC: a data.frame with the Targeting Operator Characteristic curve
estimated on grid q, along with bootstrapped SEs.
\item confint: a data.frame with the lower and upper bounds of the RATE
confidence interval.
\item pval: the p-value for the test that RATE is non-positive.
}
}
\description{
Provides the P-value for a formal test of heterogeneity based on the RATE
statistic by Yadlowsky et al.
}
\details{
RATE evaluates the ability of a provided prioritization rule to
prioritize treatment to subjects with a large benefit. In order to test for
heterogeneity, we want estimated CATE's to define the prioritization rule.
However, to obtain valid inference the prioritization scores must be
constructed independently of the evaluating forest training data. To
accomplice this, we split the data and train separate forests on each part.
Then we estimate double robust scores on the observations used to train
each forest, and obtain prioritization scores by predicting CATE's with
each forest on the samples not used for training.
}
\examples{
\donttest{
n <- 800
p <- 3
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5)
event_prob <- 1 / (1 + exp(2 * (pmax(2 * X[, 1], 0) * W - X[, 2])))
Y <- rbinom(n, 1, event_prob)
clusters <- sample(1:4, n, replace = TRUE)
cf <- grf::causal_forest(X, Y, W, clusters = clusters)
rate <- RATEOmnibusTest(cf, target = "QINI")
rate
}

}
\references{
Yadlowsky S, Fleming S, Shah N, Brunskill E, Wager S. Evaluating
Treatment Prioritization Rules via Rank-Weighted Average Treatment Effects.
2021. http://arxiv.org/abs/2111.07966.
}
\author{
KIJA
}
