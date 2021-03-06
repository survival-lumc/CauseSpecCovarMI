% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction-helpers.R
\name{cuminc_weib}
\alias{cuminc_weib}
\title{Compute true cumulative incidence for one of two competing events}
\usage{
cuminc_weib(alph_ev, lam_ev, alph_comp, lam_comp, t)
}
\arguments{
\item{alph_ev}{Numeric - shape of event of interest}

\item{lam_ev}{Numeric - rate of event of interest}

\item{alph_comp}{Numeric - shape of competing event}

\item{lam_comp}{Numeric - rate of competing event}

\item{t}{Scalar of vector of (positive) timepoints}
}
\value{
Scalar, value of cumulative incidence
}
\description{
Based on Weibull hazards from both competing events. Integration is
done via \code{stats::integrate}.
}
\examples{
time <- seq(0.1, 10, by = 0.1)

# Use baseline parametrisation from sim study
cuminc <- cuminc_weib(
alph_ev = 0.58,
lam_ev = 0.19,
alph_comp = 0.53,
lam_comp = 0.21,
t = time
)

plot(
x = time, 
y = cuminc, 
type = "l", 
xlab = "Time", 
ylab = "Cumulative incidence"
)

}
