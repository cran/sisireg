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
  W0 <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$W0, dtype=npR$float64), name='W0')
  W1 <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$W1, dtype=npR$float64), name='W1')
  W2 <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$W2, dtype=npR$float64), name='W2')
  minX <- onnxR$numpy_helper$from_array(npR$array(matrix(ssrmlp_model$minX,length(ssrmlp_model$minX),1), dtype=npR$float64), name='minX')
  maxX <- onnxR$numpy_helper$from_array(npR$array(matrix(ssrmlp_model$maxX,length(ssrmlp_model$maxX),1), dtype=npR$float64), name='maxX')
  minY <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$minY, dtype=npR$float64), name='minY')
  maxY <- onnxR$numpy_helper$from_array(npR$array(ssrmlp_model$maxY, dtype=npR$float64), name='maxY')
  
  # Input und Output (FLOAT=type 1, DOUBLE=type 11)
  X <- onnxR$helper$make_tensor_value_info('X', onnxR$TensorProto$DOUBLE, list(1L, as.integer(ncol(ssrmlp_model$W0)-1)))
  Y <- onnxR$helper$make_tensor_value_info('Y', onnxR$TensorProto$DOUBLE, list(1L, 1L))
  
  # create the nodes 
  nodes <- list()
  #X <- t(apply(X, 1, function(x) (x-W$minX)/(W$maxX-W$minX)))
  nodes <- append(nodes, onnxR$helper$make_node('Transpose', list('X'), list('Xt')))
  nodes <- append(nodes, onnxR$helper$make_node('Sub', list('Xt', 'minX'), list('X1')))
  nodes <- append(nodes, onnxR$helper$make_node('Sub', list('maxX', 'minX'), list('X2')))
  nodes <- append(nodes, onnxR$helper$make_node('Div', list('X1', 'X2'), list('X3')))
  #X <- cbind(X, rep(1, nrow(X)))
  nodes <- append(nodes, onnxR$helper$make_node('Constant', list(), list('B'), 
                                                value=onnxR$helper$make_tensor('const', onnxR$TensorProto$DOUBLE, list(1L, 1L), matrix(1,1))))
  nodes <- append(nodes, onnxR$helper$make_node('Concat', list('X3', 'B'), list('XB'), axis=-2L))
  #O1 <- sigmoidR(W0 %*% t(X))
  nodes <- append(nodes, onnxR$helper$make_node('MatMul', list('W0', 'XB'), list('O1')))
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
  
  # assembling the graph
  graph <- onnxR$helper$make_graph(
    nodes,  # nodes
    'ssrmlp', # a name
    list(X),  # inputs
    list(Y),  # outputs
    list(W0, W1, W2, minX, maxX, minY, maxY))    # parameter
  # completing the model
  model <- onnxR$helper$make_model(graph)
  
  # check the generated model
  inferred_model = onnxR$shape_inference$infer_shapes(model)
  print(inferred_model$graph$value_info)  # gives the dimension of each node
  onnxR$checker$check_model(model, full_check=TRUE)

  # and save the model
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