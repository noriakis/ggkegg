#' plot_module_text
#' plot the text representation of KEGG modules
#' @import tidygraph
#' @import patchwork
plot_module_text <- function(plot_list, show_name="name") {
  panel_list <- list()
  for (concat in seq_along(plot_list)) {
    plot_list[[concat]]$name <- plot_list[[concat]][[show_name]]
    g <- tbl_graph(nodes=plot_list[[concat]])
    panel_list[[concat]] <- ggraph(g, x=x, y=1) +
     geom_node_rect(aes(filter=.data$koflag), fill=plot_list[[concat]][plot_list[[concat]]$koflag,]$color,
      alpha=0.5, color="black")+
     geom_node_rect(aes(filter=!.data$koflag & !.data$conflag), fill="transparent", color="black")+
     geom_node_text(aes(label=name,filter=.data$koflag | .data$conflag))+
     theme_void()
  }
  wrap_plots(panel_list, ncol=1)
}


#' plot_module_steps
#' wrapper function for plotting module definition steps
#' @export
plot_module_steps <- function(all_steps, layout="kk") {
  allnodes <- unique(V(all_steps)$name)
  if (sum(startsWith(allnodes, "K"))==length(allnodes)) {
    stop("all nodes are KO.")
  }
  ggraph(all_steps, layout=layout) + 
    geom_node_point(size=4, aes(filter=!startsWith(name,"STEP") &
                                  !startsWith(name,"G") &
                                !startsWith(name,"CS"))) + 
    geom_node_point(size=2, shape=21, aes(filter=startsWith(name,"STEP"))) + 
    geom_node_point(size=2, shape=21, aes(filter=startsWith(name,"CS") |
                                            startsWith(name,"G"))) + 
    geom_edge_link(aes(filter=type %in% c("step_transition","rel")),
                       arrow = arrow(length = unit(3, 'mm')),
                       end_cap=circle(5, 'mm'),start_cap=circle(5,"mm"))+
    geom_edge_link(aes(filter=!type %in% c("step_transition","rel","instep"))) + 
    geom_edge_link(aes(label=type,
                           filter=!startsWith(type,"in") & 
                             !type %in% c("step_transition","rel")),
                       angle_calc="along",
                       label_dodge = unit(2, 'mm')) + 
    geom_node_text(aes(label=name,
                       filter=startsWith(name,"K")),
                   repel=TRUE, size=4, bg.colour="white")+
    theme_void()
}



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
            class = "geom_node_rect_kegg")
}

#' ggplot_add.geom_node_rect_kegg
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_node_rect_kegg
#' @export
ggplot_add.geom_node_rect_kegg <- function(object, plot, object_name) {
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