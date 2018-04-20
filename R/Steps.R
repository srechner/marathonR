
stepsKTVSwitch <- function(rowsum, colsum) {
  return( 10 * length(rowsum) * length(colsum))
}

stepsEdgeSwitch <- function(rowsum, colsum) {
  return( 10 * length(rowsum) * length(colsum))
}

stepsCurveball <- function(rowsum, colsum) {
  return( 10 * length(rowsum) )
}

stepsSimple <- function(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper) {
  return (10 * length(rowsum) * length(colsum))
}

stepsInformed <- function(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper) {
  k <-(length(rowsum) + length(colsum)) / 2
  return (10 * k)
}
