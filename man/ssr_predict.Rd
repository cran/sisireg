% File sisireg/man/ssr_predict.Rd
\name{ssr_predict}
\alias{ssr_predict}
\title{SSR model Prediction}
\description{
  Calculates the prediction for a given SSR model.
}
\usage{
ssr_predict(df, xx)
}
\arguments{
  \item{df}{dataframe containing two series with x- und y-values.}
  \item{xx}{array containung locations for predictions.}
}
\value{
  \item{yy}{array containung the predicted values.}
}
\examples{
  set.seed(1234)
  df <- data.frame(x=runif(500, min=-1, max=1)*pi)
  df$y <- sin(df$x)*20 + rnorm(nrow(df), mean=0, sd=10)
  plot(df, xlim=c(-4, 4))
  dfl1 <- ssr(df)
  lines(dfl1)
  xx <- c(-4, -1, 0, 1, 4)
  yy <- ssr_predict(dfl1, xx)
  points(xx,yy, pch='+', col='red', cex=2)
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
