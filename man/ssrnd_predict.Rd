% File sisireg/man/ssrnd_predict.Rd
\name{ssrnd_predict}
\alias{ssrnd_predict}
\title{Prediction function for the multi-dimensional SSR model}
\description{
Calculates the prediction for a given multi-dimensional SSR model
}
\usage{
ssrnd_predict(df_model, xx)
}
\arguments{
  \item{df_model}{data frame with model coordinates.}
  \item{xx}{data frame with coordinates for prediction.}
}
\value{
  \item{z}{list with predictions.}
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
# Prediction
xx <- c(c(0,1), c(-1,1), c(1,-1))
xx <- matrix(xx, ncol = 2)
yy <- ssrnd_predict(df_model, xx)
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
