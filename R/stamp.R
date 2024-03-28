#' stamp
#' 
#' place stamp on the specified node
#' 
#' @param name name of the nodes
#' @param color color of the stamp
#' @param which_column which node column to search
#' @export
#' @return ggplot2 object
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  stamp("hsa:6737")
stamp <- function(name, color="red", which_column="name") {
  structure(list(name=name, color=color, which_column=which_column),
            class="stamp")
}

#' ggplot_add.stamp
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_node_rect_kegg
#' @export
#' @return ggplot2 object
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  stamp("hsa:6737")
ggplot_add.stamp <- function(object, plot, object_name) {
    plot <- plot + geom_node_rect(aes(xmin=.data$xmin-2, xmax=.data$xmax+2, 
                                      ymin=.data$ymin-2, ymax=.data$ymax-2, 
                                      filter=.data[[object$which_column]] %in% object$name),
                                  fill="transparent", color=object$color)
    plot
}
