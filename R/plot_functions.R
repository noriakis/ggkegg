#' geom_node_rect
#' @param g ggraph object
#' @param with_group grouping rectangle will be drawn
#' @export
geom_node_rect <- function(with_group=TRUE) {
  structure(list(with_group=with_group),
            class = "geom_node_rect")
}

#' ggplot_add.geom_node_rect
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_node_rect
#' @export
ggplot_add.geom_node_rect <- function(object, plot, object_name) {
  plot <- plot +  geom_rect(aes(xmin=xmin, ymin=ymin,
                     xmax=xmax, ymax=ymax),
                 plot$data[!plot$data$undefined,],
                 color="black",
                 fill="grey")
  if (!is.null(plot$data$undefined)) {
    plot <- plot + geom_rect(aes(xmin=xmin, ymin=ymin,
                    xmax=xmax, ymax=ymax),
                data=plot$data[plot$data$undefined,],
                fill="transparent", color="red")
  }
  plot
}