#' Create random binary matrices with prescribed row and column sums.
#'
#' The function \code{SampleFixed} produces \code{N} random binary matrices whose row and column sums match a given pair of integer vectors.
#' Each matrix is drawn uniformly from the set of all binary matrices possessing a fixed set of \code{rowsum}s and \code{colsum}s.
#' For this purpose, the user may choose one of four sampling algorithms by specifying the \code{method} parameter.
#' As most of these methods are based on the simulation of a Markov chain, they require an assumption on the number of \code{steps} that is neccessary to produce random samples. 
#' In contrast, an \code{exact} method may be used to produce random samples without having to rely on the number of steps.
#' 
#' @param rowsum Integer vector defining the row sums.
#' @param colsum Integer vector defining the column sums.
#' @param samples Number of samples to produce.
#' @param steps Number of steps. (Larger is better.)
#'   If \code{method} is "exact", no \code{steps} argument is required.
#' @param method Sampling algorithm. Must be one of the following
#'   "exact": Use a sampling algorithm that is proven to produce uniformly distributed samples. However, this algorithm may have an unpleasent running time.
#'   "ktc-switch": Simulation of the Markov chain defined by Kannan et al.: 'Simple Markov-chain algorithms for generating bipartite graphs and tournaments'. Random Structures andAlgorithms 14 (1997), 293–308.
#'   "edge-switch": Variant of the ktc-switch algorithm, based on an informed edge selection, at the cost of a larger memory consumption.
#'   "curveball": (Default.) This method implements the Curveball algorithm suggested by Strona et al.: 'A fast and unbiased procedure to randomize ecological binary matrices with fixed row and column totals.' Nature communications 5 (2014).
#' @return A list of binary matrices.
#' @examples
#' # First, we need to define the row and column sums of the binary matrices
#' # we want to produce, e.g.
#'
#' rowsum <- c(2, 1, 2)
#' colsum <- c(1, 1, 1, 2).
#'
#' # To use the 'exact' sampling method (which may be costly for large instances), we just
#' # pass the row and column sums in addition to the number of samples we require. Here, we want to create 100 independent matrices
#'
#' l <- SampleFixed(rowsum, colsum, 100, method="exact")
#' 
#' For large instances, the running time of the exact sampling method may become unpleasently high. In such cases, it is better to use 
#' one of the other sampling algorithms. In this case, we need to specify a number of 'steps'. 
#' It is hard to define a general rule on how to choose this parameter. In general, the larger this value is, the more independent will the binary matrices be. 
#' 
#' l <- SampleFixed(rowsum, colsum, 100, method="edge-switch", steps=100)
#'
SampleFixed <- function(rowsum,
                        colsum,
                        samples,
                        method = "exact",
                        steps = NULL) {
  
  supportedMethods <- c("exact", "curveball", "edge-switch", "ktv-switch")
  
  ## check 'samples'
  if (!is.numeric(samples))
    stop("The parameter 'samples' must be numeric.")
  if (samples < 0)
    stop("The parameter 'samples' must be postive.")
  
  ## check 'rowsum' and 'colsum'
  if (!is.numeric(rowsum))
    stop("'rowsum' must be a numeric vector.")
  if (!is.numeric(colsum))
    stop("'colsum' must be a numeric vector.")
  if (sum(rowsum < 0) > 0)
    stop("Every element in 'rowsum' must be non-negative.")
  if (sum(colsum < 0) > 0)
    stop("Every element in 'colsum' must be non-negative.")
  
  if (!CppIsRealizableFixed(rowsum, colsum))
    stop("Error! Empty sample space!")
  
  ## check 'method'
  if (!is.character(method))
    
    method <- tolower(method)
  if (!is.element(method, lapply(supportedMethods, tolower)))
    stop(paste("The parameter 'method' must be one of:", paste(supportedMethods)))
  
  ## check 'steps'
  if (is.null(steps)) {
    if(method != "exact")
      stop("A number of 'steps' is required.")
    else
      steps <- 0
  } else {
    if (!is.numeric(steps))
      stop("The parameter 'steps' must be an integer.")
    if (steps < 1)
      stop("The parameter 'steps' must not be negative.")
    if (method == "exact")
      warning("The parameter 'method' is set to 'exact'. Will ignore the parameter 'steps'.")
  }
  
  ## ok, we just sample binary matrices and return
  res <- CppSampleBinaryMatricesFixed(rowsum,
                                      colsum,
                                      samples,
                                      steps,
                                      method)
  return(res)
}

#' Create random binary matrices whose row and column sums are bounded by prescribed integers.
#'
#' The function \code{SampleInterval} produces \code{N} random binary matrices whose row and column sums lie in prescribed intervals.
#' Each matrix is drawn uniformly from the set of binary matrices that meet these conditions.
#' The sampling method can be specified by the \code{method} parameter. As the "simple" and "informed" sampling methods are based on the simulation of a Markov chain, 
#' they require an assumption on the number of \code{steps} that is neccessary to produce random samples. 
#' In contrast, an \code{exact} method may be used to produce random samples without having to rely on the number of steps.

#' 
#' @param rowsum.lower Integer vector defining the lower bounds on each row sum.
#' @param rowsum.upper Integer vector defining the upper bounds on each row sum.
#' @param colsum.lower Integer vector defining the lower bounds on each column sum.
#' @param colsum.upper Integer vector defining the upper bounds on each column sum.
#' @param samples Number of samples to produce.
#' @param steps Number of steps. (Larger is better.)
#'   If \code{method} is "exact", no \code{steps} argument is required.
#' @param method Sampling algorithm. Must be one of the following
#'   "exact": Use a sampling algorithm that is proven to produce uniformly distributed samples. However, this algorithm may have an unpleasent running time.
#'   "simple": Simulation of the first Markov chain defined in Rechner et al.: 'Uniform sampling of bipartite graphs with degrees in prescribed intervals'. Journal of Complex Networks (2017).
#'   "informed": Simulation of the second Markov chain defined in Rechner et al.: 'Uniform sampling of bipartite graphs with degrees in prescribed intervals'. Journal of Complex Networks (2017).
#' @return A list of binary matrices.
#' @examples
#' # First, we need to define the lower and upper bounds on the row and column sums of the binary matrices
#' # we want to produce, e.g.
#'
#' rowsum.lower <- c(1, 0, 1)
#' rowsum.upper <- c(2, 1, 1)
#' colsum.lower <- c(0, 1, 1, 0)
#' colsum.upper <- c(1, 1, 2, 3).
#'
#' # To use the 'exact' sampling method (which may be costly for large instances), we just
#' # pass the lower and upper bounds on the row and column sums in addition to the number of samples we require. Here, we want to create 100 independent matrices
#'
#' l <- SampleInterval(rowsum.lower, rowsum.upper, colsum.lower, colsum.upper, 100, method="exact")
#' 
#' For large instances, the running time of the exact sampling method may become unpleasently high. In such cases, it is better to use 
#' one of the other sampling algorithms. In this case, we need to specify a number of 'steps'. 
#' It is hard to define a general rule on how to choose this parameter. In general, the larger this value is, the more independent will the binary matrices be. 
#' 
#' l <- SampleInterval(rowsum.lower, rowsum.upper, colsum.lower, colsum.upper, 100, method="simple", steps=100)
#'
SampleInterval <- function(rowsum.lower,
                           rowsum.upper,
                           colsum.lower,
                           colsum.upper,
                           samples,
                           method = "exact",
                           steps = NULL) {
  
  supportedMethods <- c("exact", "simple", "informed")
  
  ## check 'samples'
  if (!is.numeric(samples))
    stop("The parameter 'samples' must be numeric.")
  if (samples < 0)
    stop("The parameter 'samples' must be postive.")
  
  ## check 'rowsum' and 'colsum'
  if (!is.numeric(rowsum.lower))
    stop("'rowsum_lower' must be a numeric vector.")
  if (!is.numeric(rowsum.upper))
    stop("'rowsum_upper' must be a numeric vector.")
  if (!is.numeric(colsum.lower))
    stop("'colsum_lower' must be a numeric vector.")
  if (!is.numeric(colsum.upper))
    stop("'colsum_upper' must be a numeric vector.")
  if (sum(rowsum.lower > rowsum.upper) > 0)
    stop("Empty sample space!")
  if (sum(colsum.lower > colsum.upper) > 0)
    stop("Empty sample space!")
  
  if (!CppIsRealizableInterval(rowsum.lower, rowsum.upper, colsum.lower, colsum.upper))
    stop("Error! Empty sample space!")
  
  ## check 'method'
  if (!is.character(method))
    stop(paste("The parameter 'method' must be one of:", paste(supportedMethods)))
  
  method <- tolower(method)
  if (!is.element(method, lapply(supportedMethods, tolower)))
    stop(paste("The parameter 'method' must be one of:", paste(supportedMethods)))
  
  ## check 'steps'
  if (is.null(steps)) {
    if(method != "exact")
      stop("A number of 'steps' is required!")
    else
      steps <- 0
  } else {
    if (!is.numeric(steps))
      stop("The parameter 'steps' must be an integer.")
    if (steps < 1)
      stop("The parameter 'steps' must not be negative.")
    if (method == "exact")
      warning("The parameter 'method' is set to 'exact'. Will ignore the parameter 'steps'.")
  }
  
  ## ok, we just sample binary matrices and return
  res <- CppSampleBinaryMatricesInterval(rowsum.lower,
                                         rowsum.upper,
                                         colsum.lower,
                                         colsum.upper,
                                         samples,
                                         steps,
                                         method)
  return(res)
}