#' stamp
#' 
#' place stamp on the specified node
#' 
#' @param name name of the nodes
#' @param color color of the stamp
#' @param which_column which node column to search
#' @param xval adjustment value for x-axis
#' @param yval adjustment value for y-axis
#' @export
#' @return ggplot2 object
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  stamp("hsa:6737")
stamp <- function(name, color="red", which_column="name", xval=2, yval=2) {
  structure(list(name=name, color=color, which_column=which_column, xval=xval, yval=yval),
            class="stamp")
}

#' ggplot_add.stamp
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param ... The other arguments
#' @export ggplot_add.geom_node_rect_kegg
#' @export
#' @return ggplot2 object
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  stamp("hsa:6737")
ggplot_add.stamp <- function(object, plot, ...) {
    plot <- plot + geom_node_rect(aes(xmin=.data$xmin-object$xval, xmax=.data$xmax+object$xval, 
                                      ymin=.data$ymin-object$yval, ymax=.data$ymax+object$yval, 
                                      filter=.data[[object$which_column]] %in% object$name),
                                  fill="transparent", color=object$color)
    plot
}
