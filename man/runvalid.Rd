% File sisireg/man/runvalid.Rd
\name{runvalid}
\alias{runvalid}
\title{Maximum Run Validity Check}
\description{Checks, if a given regression function is adequate with respect to the maximum run criterium.}
\usage{
runvalid(dat,mu,k=NULL)
}
\arguments{
  \item{dat}{obervations.}
  \item{mu}{discrete regression function.}
  \item{k}{optional: maximum run length.}
}
\value{
  \item{valid}{function is valid?}
}
\examples{
runvalid(sin(seq(-pi, pi, length.out = 200))+rnorm(200)/2, 
          sin(seq(-pi, pi, length.out = 200)))
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
