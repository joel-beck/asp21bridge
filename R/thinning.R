### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### returns thinned numeric vector or matrix                                ####

thin_element <- function(samples, freq) {
  if (is.matrix(samples)) {
    samples[seq(from = 1, to = nrow(samples), by = freq), , drop = FALSE]
  } else if (is.vector(samples, mode = "numeric")) {
    samples[seq(from = 1, to = length(samples), by = freq)]
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector,",
      "matrix or list of vectors and matrices!"
    ))
  }
}


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### returns thinned object of the same type (vector, matrix or list)        ####

thinning <- function(samples, freq) {
  if (is.matrix(samples) || is.vector(samples, mode = "numeric")) {
    thin_element(samples = samples, freq = freq)
  } else if (is.list(samples)) {
    lapply(X = samples, thin_element, freq = freq)
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector,",
      "matrix or list of vectors and matrices!"
    ))
  }
}
