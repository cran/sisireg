% File sisireg/man/ffi_prediction.Rd
\name{onnx_save}
\alias{onnx_save}
\title{
Saving a ssrmlp model in onnx format to file
}
\description{
Saving a ssrmlp model in onnx format to file (also see onnx.ai). This function uses the onnx python implementation, hence a python environment including modules onnx and numpy is required.
}
\usage{
onnx_save(ssrmlp_model, filename)
}
\arguments{
  \item{ssrmlp_model}{a trained ssrmlp model}
  \item{filename}{fully qualified file name}
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
ssrmlp_model <- ssrmlp_train(X, Y)
# prediction
p <- t(c(0.25, 0.25))
pred <- ssrmlp_predict(p, ssrmlp_model) 
# only if python is available
if (reticulate::py_module_available('onnx')) {
  tryCatch( 
    {
      # save in ONNX format
      onnx_save(ssrmlp_model, 'file.onnx')
      # and now test with onnxruntime.ai
    },
    error = function(cond) {
      message(conditionMessage(cond))
    },
    finally = {
      # cleanup
      file.remove('file.onnx')
      # to avoid NOTE in R CHECK
      tempfile <-  reticulate::import("tempfile")
      tmp <- tempfile$gettempdir()
      if (dir.exists(file.path(tmp, "__pycache__"))) {
        unlink(file.path(tmp, "__pycache__"), recursive = TRUE, force = TRUE)
      }
      tmp_py_files <- list.files(tmp, 
                                 pattern = "^__autograph_generated_file.*py$", full.names = TRUE)
      file.remove(tmp_py_files)
      print("done")
    }
  )
}

}
}
\author{
Dr. Lars Metzner
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
