% File sisireg/man/ssrnd_predict.Rd
\name{ssrnd_predict}
\alias{ssrnd_predict}
\title{Prediction function for the multi-dimensional SSR model}
\description{
Calculates the prediction for a given multi-dimensional SSR model
}
\usage{
ssrnd_predict(xx, mdl)
}
\arguments{
  \item{xx}{data frame with coordinates for prediction.}
  \item{mdl}{data frame with model coordinates.}
}
\value{
  \item{z}{list with predictions.}
}
\examples{
\donttest{
# generate data
set.seed(1234)
x_ <- y_ <- seq(-3, 3, length=20)
xy <- expand.grid(x=x_,y=y_)
x <- xy$x
y <- xy$y
z <- rnorm(400) + atan2(x, y)
X <- matrix(cbind(x, y), ncol = 2)
Y <- as.double(z)
# Training
df_model <- ssrnd(X, Y, k = 4, fn = 8)
# Prediction
xx <- c(c(0,1), c(-1,1), c(1,-1))
xx <- matrix(xx, ncol = 2)
yy <- ssrnd_predict(xx, df_model)
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
