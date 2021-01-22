% File sisireg/man/ssrnd.Rd
\name{ssrnd}
\alias{ssrnd}
\title{Multi-dimensional SSR model}
\description{
Calculates the multi-dimensional SSR model
}
\usage{
ssrnd(koord, dat, k = NULL, fn = NULL, iter = 1000)
}
\arguments{
  \item{koord}{data frame with n-dimensional coordinates.}
  \item{dat}{data frame with observations.}
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
x <- rnorm(300)
y <- rnorm(300)
xy <- data.frame(x=x, y=y)
z <-data.frame(z=rnorm(300) + atan2(x, y))
# Training
df_model <- ssrnd(xy, z)
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
