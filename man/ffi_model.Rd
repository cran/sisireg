% File sisireg/man/ffi_model.Rd
\name{ffi_model}
\alias{ffi_model}
\title{
Factor-wise Influence Indicator (Model-FFI) for a given ssrmlp model
}
\description{
The Model-FFI depicts the overall influence of the input factors on the resulting trained ssrmlp model. For computation a unit matrix is used to accumulate the weights for each factor separately.
}
\usage{
ffi_model(W)
}
\arguments{
  \item{W}{a trained ssrmlp model}
}
\value{
  \item{ffi}{array of influence indicators}
}
\examples{
\donttest{
# generate data
set.seed(42)
x <- rnorm(300)
y <- rnorm(300)
z <- rnorm(300) + atan2(x, y)
# coordinates
X <- matrix(cbind(x,y), ncol = 2)
Y <- as.double(z)
# Training
W <- ssrmlp_train(X, Y)
ffi_model(W)
}

}
\author{
Dr. Lars Metzner
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
