% File sisireg/man/ffi_prediction.Rd
\name{onnx_load}
\alias{onnx_load}
\title{
Loading a ssrmlp model from ONNX file}
\description{
Loading a ssrmlp model from ONNX file (also see onnx.ai). This function uses the onnx python implementation, hence a python environment including modules onnx and numpy is required.
}
\usage{
onnx_load(filename)
}
\arguments{
  \item{filename}{fully qualified file name}
}
\value{
  \item{W}{parameters of ssrmlp model}
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
# only if python is available
if (reticulate::py_module_available('onnx')) {
  # save in ONNX format
  onnx_save(ssrmlp_model, 'file.onnx')
  # loading the file in ONNX format
  W <- onnx_load('file.onnx')
  # prediction with original implementation
  p <- t(c(0.25, 0.25))
  pred <- ssrmlp_predict(p, W) 
  # predict with TensorFlow
  onnx <- reticulate::import("onnx")
  backend <-  reticulate::import('onnx_tf.backend')
  model = onnx$load('file.onnx') 
  tf_rep = backend$prepare(model)       
  tf_y = tf_rep$run(p)                
  print(tf_y$Y-pred)
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
}
}

}
\author{
Dr. Lars Metzner
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
