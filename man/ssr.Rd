% File sisireg/man/ssr.Rd
\name{ssr}
\alias{ssr}
\title{Onedimensional SSR-model calculation}
\description{
  Calculates L1- and L2-functions satisfiying the partial sum criterium.
}
\usage{
ssr(df, y1=NULL, yn=NULL, fn=0, iter=10000, minStat=FALSE, ne=TRUE, l1=TRUE, ps=TRUE)
}
\arguments{
  \item{df}{data frame with two-dimensional data.}
  \item{y1}{optional: fixed value left.}
  \item{yn}{optional: fixed value right.}
  \item{fn}{optional: partial-sum-quantile (standard: generic calculation from data).}
  \item{iter}{optional: maximum number of iterations.}
  \item{minStat}{optional: boolean value for the minimum statistic.}
  \item{ne}{optional: boolean value for non-equidistant observations.}
  \item{l1}{optional: boolean value for function type.}
  \item{ps}{optional: sign criterium (partial sum or run).}
}
\value{
  \item{mu}{SSR-function as array.}
}
\examples{
# generate equidistant data
set.seed(1234)
x <- seq(0, 2*pi, length.out = 200)
y <- 4*sin(x) + rnorm(200)
df <- data.frame(x=x, y=y)
# calculate regression functions
l1 <- ssr(df, ne=FALSE, ps=FALSE)
l2 <- ssr(df, ne=FALSE, l1=FALSE)
lmin <- ssr(df, ne=FALSE, minStat=TRUE, ps=FALSE)
# plot results
plot(x, y, main = 'Sign-Simplicity-Regression', xlab = 't', ylab = 'sin(t)+noise')
lines(x, l1, col = 'blue')
lines(x, l2, col = 'red')
lines(x, lmin, col = 'purple')
legend("topleft", inset=c(0.01,0.01), legend=c("L1 run-crit.", "L2 ps-crit.", "L1 min-stat."),
       col=c("blue", "red", "purple"), lty=1:1)

# generate nonequidistant data
df <- data.frame(x=runif(500, min=-1, max=1)*pi)
df$y <- sin(df$x)*20 + rnorm(nrow(df), mean=0, sd=10)
# calculate regression function
dfl1 <- ssr(df, fn = 5)
# plot results
plot(df)
lines(dfl1, col = 'red')
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
