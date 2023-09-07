#
# Handling of ONNX protocoll (onnx.ai) for ssrmlp models
#


# using python functionality
npR <- NULL
onnxR <- NULL
.onLoad <- function(libname, pkgname) {
  npR   <<- reticulate::import("numpy", convert = FALSE, delay_load = TRUE)
  onnxR <<- reticulate::import("onnx", delay_load = TRUE)
}


#
# save a given ssrmlp model to onnx file with file name filename
#
onnx_save <- function(ssrmlp_model, filename) {
  # initializers: model parameters 
  B  <- onnxR$numpy_helper$from_array(npR$array(as.matrix(1,1,1), dtype=npR$float64), name='B')
  W0 <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$W0, dtype=npR$float64), name='W0')
  W1 <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$W1, dtype=npR$float64), name='W1')
  W2 <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$W2, dtype=npR$float64), name='W2')
  minX <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$minX, dtype=npR$float64), name='minX')
  maxX <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$maxX, dtype=npR$float64), name='maxX')
  minY <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$minY, dtype=npR$float64), name='minY')
  maxY <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$maxY, dtype=npR$float64), name='maxY')
  
  # Input und Output
  X <- onnxR$helper$make_tensor_value_info('X', onnxR$TensorProto$FLOAT, list(1L, as.integer(ncol(ssrmlp_model$W0)-1)))
  Y <- onnxR$helper$make_tensor_value_info('Y', onnxR$TensorProto$FLOAT, list(1L, 1L))
  
  # create the graph
  nodes <- list()
  #X <- t(apply(X, 1, function(x) (x-W$minX)/(W$maxX-W$minX)))
  nodes <- append(nodes, onnxR$helper$make_node('Transpose', list('X'), list('Xt')))
  nodes <- append(nodes, onnxR$helper$make_node('Sub', list('Xt', 'minX'), list('X1')))
  nodes <- append(nodes, onnxR$helper$make_node('Sub', list('maxX', 'minX'), list('X2')))
  nodes <- append(nodes, onnxR$helper$make_node('Div', list('X1', 'X2'), list('X3')))
  #X <- cbind(X, rep(1, nrow(X)))
  nodes <- append(nodes, onnxR$helper$make_node('Concat', list('X3', 'B'), list('XB'), axis=1L))
  #O1 <- sigmoidR(W0 %*% t(X))
  nodes <- append(nodes, onnxR$helper$make_node('Transpose', list('XB'), list('XBt')))
  nodes <- append(nodes, onnxR$helper$make_node('MatMul', list('W0', 'XBt'), list('O1')))
  nodes <- append(nodes, onnxR$helper$make_node('Sigmoid', list('O1'), list('O1s')))
  #O2 <- sigmoidR(W1 %*% O1)
  nodes <- append(nodes, onnxR$helper$make_node('MatMul', list('W1', 'O1s'), list('O2')))
  nodes <- append(nodes, onnxR$helper$make_node('Sigmoid', list('O2'), list('O2s')))
  #y <- t(W2 %*% O2)
  nodes <- append(nodes, onnxR$helper$make_node('MatMul', list('W2', 'O2s'), list('Y0')))
  nodes <- append(nodes, onnxR$helper$make_node('Transpose', list('Y0'), list('Y0t')))
  #y = y * (W$maxY - W$minY)  + W$minY
  nodes <- append(nodes, onnxR$helper$make_node('Sub', list('maxY', 'minY'), list('Y1')))
  nodes <- append(nodes, onnxR$helper$make_node('Mul', list('Y0t', 'Y1'), list('Y2')))
  nodes <- append(nodes, onnxR$helper$make_node('Add', list('Y2', 'minY'), list('Y')))
  
  # all together in a graph
  graph <- onnxR$helper$make_graph(
    nodes,  # nodes
    'ssrmlp', # a name
    list(X),  # inputs
    list(Y),  # outputs
    list(B, W0, W1, W2, minX, maxX, minY, maxY))    # parameter
  
  model <- onnxR$helper$make_model(graph)
  onnxR$checker$check_model(model)
  onnxR$save_model(model, filename)
  
}


#
# load ssrmlp model parameters from onnx file with file name filename
#
onnx_load <- function(filename) {
  model <- onnxR$load(filename)
  
  param <- model$graph$initializer

  # Extraction of the parameters
  # list(B, W0, W1, W2, minX, maxX, minY, maxY)
  W <- list()
  for (i in 0:length(param)-1) {
    if (length(param[i]$dims) == 0) { # scalar values
      W[[param[i]$name]] <- onnxR$numpy_helper$to_array(param[i])[1]
    } else {
      W[[param[i]$name]] <- onnxR$numpy_helper$to_array(param[i])
    }
  }
  return(W)
}