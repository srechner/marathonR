% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SampleBinaryMatrices.R
\name{SampleFixed}
\alias{SampleFixed}
\title{Create random binary matrices with prescribed row and column sums.}
\usage{
SampleFixed(rowsum, colsum, samples, method = "exact", steps = NULL)
}
\arguments{
\item{rowsum}{Integer vector defining the row sums.}

\item{colsum}{Integer vector defining the column sums.}

\item{samples}{Number of samples to produce.}

\item{method}{Sampling algorithm. Must be one of the following
"exact": Use a sampling algorithm that is proven to produce uniformly distributed samples. However, this algorithm may have an unpleasent running time.
"ktc-switch": Simulation of the Markov chain defined by Kannan et al.: 'Simple Markov-chain algorithms for generating bipartite graphs and tournaments'. Random Structures andAlgorithms 14 (1997), 293–308.
"edge-switch": Variant of the ktc-switch algorithm, based on an informed edge selection, at the cost of a larger memory consumption.
"curveball": (Default.) This method implements the Curveball algorithm suggested by Strona et al.: 'A fast and unbiased procedure to randomize ecological binary matrices with fixed row and column totals.' Nature communications 5 (2014).}

\item{steps}{Number of steps. (Larger is better.)
If \code{method} is "exact", no \code{steps} argument is required.}
}
\value{
A list of binary matrices.
}
\description{
The function \code{SampleFixed} produces \code{N} random binary matrices whose row and column sums match a given pair of integer vectors.
Each matrix is drawn uniformly from the set of all binary matrices possessing a fixed set of \code{rowsum}s and \code{colsum}s.
For this purpose, the user may choose one of four sampling algorithms by specifying the \code{method} parameter.
As most of these methods are based on the simulation of a Markov chain, they require an assumption on the number of \code{steps} that is neccessary to produce random samples. 
In contrast, an \code{exact} method may be used to produce random samples without having to rely on the number of steps.
}
\examples{
# First, we need to define the row and column sums of the binary matrices
# we want to produce, e.g.

rowsum <- c(2, 1, 2)
colsum <- c(1, 1, 1, 2).

# To use the 'exact' sampling method (which may be costly for large instances), we just
# pass the row and column sums in addition to the number of samples we require. Here, we want to create 100 independent matrices

l <- SampleFixed(rowsum, colsum, 100, method="exact")

For large instances, the running time of the exact sampling method may become unpleasently high. In such cases, it is better to use 
one of the other sampling algorithms. In this case, we need to specify a number of 'steps'. 
It is hard to define a general rule on how to choose this parameter. In general, the larger this value is, the more independent will the binary matrices be. 

l <- SampleFixed(rowsum, colsum, 100, method="edge-switch", steps=100)

}
