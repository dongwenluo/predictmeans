# ZOSull() gives a matrix for one value, and O'Sullivan predictions outside the basis range stop

    Code
      semipred(f, covariate = "x", at = list(x = 1.2))
    Condition
      Error:
      ! 'x' = 1.2 is outside the range of the O'Sullivan basis [0.0118805, 1.0088]: the smooth is not defined there. Predict within the data, or fit with s(x, bs = "os", xt = list(range = c(lo, hi))) covering the values.

