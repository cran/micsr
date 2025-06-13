% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pscore.R
\name{pscore}
\alias{pscore}
\alias{summary.pscore}
\alias{print.pscore}
\alias{print.summary.pscore}
\alias{nobs.pscore}
\alias{nobs.summary.pscore}
\alias{rg}
\alias{rg.pscore}
\alias{rg.summary.pscore}
\alias{stdev}
\alias{mean.pscore}
\alias{mean.summary.pscore}
\alias{stdev.pscore}
\alias{stdev.summary.pscore}
\title{Propensity scores}
\usage{
pscore(formula, data, maxiter = 4, tol = 0.005, link = c("logit", "probit"))

\method{summary}{pscore}(object, ...)

\method{print}{pscore}(
  x,
  ...,
  digits = getOption("digits"),
  var_equal = c("none", "strata", "group", "both")
)

\method{print}{summary.pscore}(
  x,
  ...,
  digits = getOption("digits"),
  step = c("all", "strata", "covariates", "atet")
)

\method{nobs}{pscore}(object, ..., smpl = c("total", "cs"))

\method{nobs}{summary.pscore}(object, ..., smpl = c("total", "cs"))

rg(object, ...)

\method{rg}{pscore}(object, ..., smpl = c("total", "cs"))

\method{rg}{summary.pscore}(object, ..., smpl = c("total", "cs"))

stdev(object, ...)

\method{mean}{pscore}(x, ..., var_equal = c("none", "strat", "group", "both"))

\method{mean}{summary.pscore}(x, ...)

\method{stdev}{pscore}(object, ..., var_equal = c("none", "strata", "group", "both"))

\method{stdev}{summary.pscore}(object, ..., var_equal = c("none", "strata", "group", "both"))
}
\arguments{
\item{formula}{a Formula object; the left-hand side should contain
two variables (\code{x1 + x2}), where x1 is the group variable and
x2 the outcome. The group variable can be either a dummy for
treated individuals or a factor with levels \code{"treated"} and
\code{"control"}}

\item{data}{a data frame}

\item{maxiter}{the maximum number of iterations}

\item{tol}{stratas are cut in halves as long as the hypothesis of
equal means is rejected at the \code{tol} level,}

\item{link}{the link for the binomial glm estimation, either
\code{"logit"} or \code{"probit"}}

\item{...}{further arguments}

\item{x, object}{a \code{"pscore"} or a \code{"summary.pscore"} object}

\item{digits}{number of digits for the \code{print} methods}

\item{var_equal}{to compute the variance of the ATET, variances can
be computed at the class/group level (\code{var_equal = "none"}), at
the class level (\code{var_equal = "group"}), at the group level
(\code{var_equal = "strata"}) or globally (\code{var_equal = "both"})}

\item{step}{for the \code{print.summary} method, the step of the test to
be printed: one of \code{"all"} (the default), \code{strata},
\code{covariates} and \code{atet}}

\item{smpl}{the sample to use, either the whole sample (\code{smpl = "total"}) or the sample with common support (\code{smpl = "cs"})}
}
\value{
an object of class \code{"pscore"}, with the following elements:
\itemize{
\item \code{strata}: a tibble containing the stratas, the frequencies, the
means and the variances of the propensity scores for treated and
controled observations
\item \code{cov_balance}: a tibble containing the results of the balancing
tests for every covariate; the results for the class with the
lowest p-value is reported
\item \code{unchecked_cov}: a character vector containing the names of the
covariates for which the balancing test could be computed
\item \code{model}: a tibble containing the original data, with
supplementary columns: \code{.gp} for the groups, \code{.resp} for the
outcome and \code{.cls} for the stratas
\item \code{pscore}: the glm model fitted to compute the propensity scores
}
}
\description{
Propensity scores estimation, using an algorithm that checks the
balancing hypothesis using strata and enable the estimation of the
treatment effect using stratification methods
}
\examples{
data_tuscany <- twa |>
                subset(region == "Tuscany") |>
                transform(dist2 = dist ^ 2,
                livselfemp = I((city == "livorno") * (occup == "selfemp")),
                perm = ifelse(outcome == "perm", 1, 0))
formula_tuscany <- perm + group ~ city + sex + marital + age +
   loc + children + educ + pvoto + training +
   empstat + occup + sector + wage + hour + feduc + femp + fbluecol +
   dist + dist2 + livselfemp
pscore(formula_tuscany, data_tuscany)
}
\references{
\insertRef{DEHE:WAHB:02}{micsr}

\insertRef{BECK:ICHI:02}{micsr}
}
\keyword{models}
