% File sisireg/man/ssr3d_predict.Rd
\name{ssr3d_predict}
\alias{ssr3d_predict}
\title{3-dimensional SSR model prediction}
\description{
Calculates the prediction for a given 3-dimensional SSR model.
}
\usage{
ssr3d_predict(df_model, xy, ms = FALSE)
}
\arguments{
  \item{df_model}{data frame with model coordinates.}
  \item{xy}{data frame with coordinates for prediction.}
  \item{ms}{optional: boolean value to use the minimal surface algorithm.}
}
\value{
  \item{z}{array with predictions.}
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
# Prediction
xx <- c(c(0,1), c(-1,1), c(1,-1))
xx <- matrix(xx, ncol = 2)
yy <- ssr3d_predict(df_model, xx)
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
