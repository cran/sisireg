% File sisireg/man/psplotnd.Rd
\name{psplotnd}
\alias{psplotnd}
\title{Partial Sum Plot for the multidimensional coordinates}
\description{
plots the partial sum statistic for the general n-dimensional SSR-model
}
\usage{
psplotnd(koord, dat, mu, text = "Sample")
}
\arguments{
  \item{koord}{data frame with coordinates.}
  \item{dat}{data frame of observations.}
  \item{mu}{list of discrete regression function.}
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
z <-data.frame(z=rnorm(900) + atan2(x, y))
# Training
df_model <- ssrnd(xy, z, k = 4, fn = 8)
# plot partial sum statistic
psplotnd(xy, z, df_model$mu, 'ssr3d')
}
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
