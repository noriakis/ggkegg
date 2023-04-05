#' @export
#'
geom_node_rect <- function(mapping = NULL, data = NULL, position = 'identity',
                            show.legend = NA, ...) {
  mapping <- c(mapping, aes(xmin = xmin,
                            ymin = ymin,
                            xmax = xmax,
                            ymax = ymax))
  class(mapping) <- "uneval"
  layer(
    data = data, mapping = mapping, stat = StatFilter, geom = GeomRect,
    position = position, show.legend = show.legend, inherit.aes = FALSE,
    params = list(na.rm = FALSE, ...)
  )
}

#' geom_node_rect_kegg
#' 
#' 
#' @param type type to be plotted (gene, map, compound ...)
#' @param rect_fill rectangular fill
#' @export
geom_node_rect_kegg <- function(type=NULL, rect_fill="grey") {
  ## [TODO] implement ggproto
  structure(list(type=type, rect_fill=rect_fill),
            class = "geom_node_rect")
}

#' ggplot_add.geom_node_rect
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_node_rect
#' @export
ggplot_add.geom_node_rect <- function(object, plot, object_name) {
  if (is.null(object$type)){
    type <- unique(plot$data$type)
    type <- type[type!="group"]
  } else {
    type <- object$type
  }
  if (!is.null(plot$data$undefined)) {
    plot <- plot + geom_node_rect(aes(filter=.data$undefined),
                fill="transparent", color="red")
    plot <- plot + geom_node_rect(aes(filter=.data$undefined & .data$type %in% type),
      fill=object$rect_fill, color="black")

  } else {
    plot <- plot + geom_node_rect(aes(filter=.data$type %in% type),
      fill=object$rect_fill, color="black")
  }
  plot
}