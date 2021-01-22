% File sisireg/man/psplot3d.Rd
\name{psplot3d}
\alias{psplot3d}
\title{Partial Sum Plot for 2-dimensional coordinates}
\description{Plots the partial sum statistic for the 3-dimensional SSR model}
\usage{
psplot3d(koord, z, mu, text = "Sample")
}
\arguments{
  \item{koord}{data frame with coordinates.}
  \item{z}{vector of observations.}
  \item{mu}{vector of discrete regression function.}
  \item{text}{optional: title for the plot.}
}
\value{No explicit return value: a plot is generated}
\examples{
\donttest{
# generate data
set.seed(1234)
x <- rnorm(900)
y <- rnorm(900)
xy <- data.frame(x=x, y=y)
z <- rnorm(900) + atan2(x, y)
# Training
df_model <- ssr3d(xy, z, k = 4, fn = 8)
# plot partial sum statistic
psplot3d(xy, z, df_model$mu, 'ssr3d')
}
}
\author{
  Dr. Lars Metzner
}
\encoding{UTF-8}
\references{
  Dr. Lars Metzner (2021) \emph{Adäquates Maschinelles Lernen}.
  Independently Published.
}
\keyword{model}
\keyword{regression}
\keyword{nonparametric}
