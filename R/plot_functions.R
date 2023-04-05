#' geom_node_rect
#' @param g ggraph object
#' @param with_group grouping rectangle will be drawn
#' @export
geom_node_rect <- function(type=NULL, rect_fill="grey") {
  ## [TODO] implement ggproto
  structure(list(type=type,rect_fill=rect_fill),
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
    plot <- plot + geom_rect(aes(xmin=xmin, ymin=ymin,
                    xmax=xmax, ymax=ymax),
                data=plot$data[plot$data$undefined,],
                fill="transparent", color="red")
    plot <- plot +  geom_rect(aes(xmin=xmin, ymin=ymin,
                     xmax=xmax, ymax=ymax, fill=eval(parse(text=object$rect_fill))),
                 plot$data[!plot$data$undefined & plot$data$type %in% type,],
                 color="black")

  } else {
    plot <- plot +  geom_rect(aes(xmin=xmin, ymin=ymin,
                     xmax=xmax, ymax=ymax, fill=eval(parse(text=object$rect_fill))),
                data=plot$data[plot$data$type==type,],
                color="black")    
  }
  plot
}