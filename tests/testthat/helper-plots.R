## Layer data of the first layer of a ggplot with the given geom (tests that
## do not depend on the order of the layers).
pm_layer <- function(g, geom, which = 1) {
  i <- which(vapply(g$layers, function(l) inherits(l$geom, geom), logical(1)))[which]
  ggplot2::layer_data(g, i)
}
