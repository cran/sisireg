% File sisireg/man/psplot.Rd
\name{psplot}
\alias{psplot}
\title{Partial Sum Plot}
\description{Plots the Partial Sums with their quantiles for a given set of observations und the corresponding regression function.}
\usage{
psplot(dat, mu, text = 'Sample')
}
\arguments{
  \item{dat}{observations.}
  \item{mu}{regression function.}
  \item{text}{title of the chart.}
}
\value{No explicit return value: a plot is generated}
\examples{
psplot(sin(seq(-pi, pi, length.out = 200))+rnorm(200),
  sin(seq(-pi, pi, length.out = 200)), text='Test')
}
\author{
  Dr. Lars Metzner
}
\encoding{UTF-8}
\references{
  Dr. Lars Metzner (2021) \emph{Adäquates Maschinelles Lernen}.
  Independently Published (to appear).
}
\keyword{model}
\keyword{regression}
\keyword{nonparametric}
