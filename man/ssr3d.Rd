% File sisireg/man/ssr3d.Rd
\name{ssr3d}
\alias{ssr3d}
\title{3-dimensional SSR model}
\description{
Calculates the regression function for the 3-dimensional SSR-model.
}
\usage{
ssr3d(koord, dat, k = NULL, fn = NULL, iter = 1000)
}
\arguments{
  \item{koord}{data frame with 2-dimensional coordinates.}
  \item{dat}{vector with observations.}
  \item{k}{optional: maxumum run length for the model.}
  \item{fn}{optional: quantile for partial sums.}
  \item{iter}{optional: number of iterations for the numeric solver.}
}
\value{
  \item{df}{data frame with the relevant model data.}
}
\examples{
\donttest{
# generate data
set.seed(1234)
x <- rnorm(900)
y <- rnorm(900)
xy <- data.frame(x=x, y=y)
z <- rnorm(900) + atan2(x, y)
# Training
df_model <- ssr3d(xy, z)
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
