% File sisireg/man/psvalid.Rd
\name{psvalid}
\alias{psvalid}
\title{Partial Sum Validity Check}
\description{Checks, if a given regression function is adequate with respect to the partial sum criterium.}
\usage{
psvalid(dat,mu)
}
\arguments{
  \item{dat}{obervations.}
  \item{mu}{discrete regression function.}
}
\value{
  \item{valid}{function is valid?}
}
\examples{
psvalid(sin(seq(-pi, pi, length.out = 200))+rnorm(200), 
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
