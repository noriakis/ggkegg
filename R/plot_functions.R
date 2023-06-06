#' multi_pathway_native
#' 
#' If you want to combine multiple KEGG pathways with their native coordinates,
#' supply this function a vector of pathway IDs and row number. This returns the
#' joined graph or list of graphs in which the coordinates are altered to panel
#' the pathways.
#' @param pathways pathway vector
#' @param row_num row number
#' @param return_list return list of graphs instead of joined graph
#' @export
#' @return graph adjusted for the position
#' 
multi_pathway_native <- function(pathways, row_num=2, return_list=FALSE) {
  plen <- length(pathways)
  if (plen %% 2) {col_num <- plen / row_num; addit <- 0} else {
    col_num <- as.integer(plen / row_num); addit <- plen %% row_num}
  
  tot_row <- 1
  tot_col <- 1
  miny <- 0
  gls <- list()
  for (pp in seq_len(pathways |> length())) {
    
    tot_row <- tot_row + 1

    g <- pathway(pathways[pp])
    g <- g |> mutate(x=(x/max(x)+tot_col-1),
                           y=y-miny)
    gls[[pp]] <- g

    edf <- g |> activate(nodes) |> data.frame()
    miny <- miny - min(edf$y)
    
    if (tot_row > row_num) {
      tot_row <- 1
      tot_col <- tot_col + 1
      miny <- 0
    }
  }
  if (return_list) {return(gls)}
  Reduce(graph_join, gls)
}

#' plot_module_text
#' plot the text representation of KEGG modules
#' @param plot_list the result of `module_text()`
#' @param show_name name column to be plotted
#' @importFrom tidygraph tbl_graph
#' @import patchwork
#' @return ggplot2 object
#' @export
plot_module_text <- function(plot_list, show_name="name") {
  panel_list <- list()
  for (concat in seq_along(plot_list)) {
    plot_list[[concat]]$name <- plot_list[[concat]][[show_name]]
    g <- tbl_graph(nodes=plot_list[[concat]])
    panel_list[[concat]] <- ggraph(g, x=x, y=1) +
     geom_node_rect(aes(filter=.data$koflag),
      fill=plot_list[[concat]][plot_list[[concat]]$koflag,]$color,
      alpha=0.5, color="black")+
     geom_node_rect(aes(filter=!.data$koflag & !.data$conflag),
      fill="transparent", color="black")+
     geom_node_text(aes(label=name,filter=.data$koflag | .data$conflag))+
     theme_void()
  }
  wrap_plots(panel_list, ncol=1)
}


#' plot_module_blocks
#' wrapper function for plotting module definition blocks
#' @param all_steps the result of `obtain_sequential_module_definition()`
#' @param layout ggraph layout parameter
#' @export
#' @return ggplot2 object
plot_module_blocks <- function(all_steps, layout="kk") {
  allnodes <- unique(V(all_steps)$name)
  if (sum(startsWith(allnodes, "K"))==length(allnodes)) {
    stop("all nodes are KO.")
  }
  ggraph(all_steps, layout=layout) +
    geom_edge_link(aes(filter=type %in% c("block_transition","rel")),
                       end_cap=circle(5, 'mm'),start_cap=circle(5,"mm"),
                       color="red")+
    geom_edge_link(aes(filter=!type %in% 
      c("block_transition","rel","in_block"))) + 
    geom_edge_link(aes(label=type,
                           filter=!startsWith(type,"in") & 
                             !type %in% c("block_transition","rel")),
                       angle_calc="along",
                       label_dodge = unit(2, 'mm')) + 
    geom_node_point(size=4, aes(filter=!startsWith(name,"BLOCK") &
                                  !startsWith(name,"G") &
                                !startsWith(name,"CS"))) + 
    geom_node_point(size=2, shape=21, aes(filter=startsWith(name,"BLOCK"))) + 
    geom_node_point(size=2, shape=21, aes(filter=startsWith(name,"CS") |
                                            startsWith(name,"G"))) + 
    geom_node_text(aes(label=name,
                       filter=startsWith(name,"K")),
                   repel=TRUE, size=4, bg.colour="white")+
    theme_void()
}

#' geom_node_shadowtext
#' 
#' Plot shadowtext at node position
#' 
#' @export
#' @param mapping aes mapping
#' @param data data to plot
#' @param position positional argument
#' @param show.legend whether to show legend
#' @param ... passed to `params` in `layer()` function
#' @return geom
#' @importFrom shadowtext GeomShadowText
geom_node_shadowtext <- function(mapping = NULL, data = NULL,
                           position = 'identity',
                           show.legend = NA, ...) {
  params <- list(na.rm = FALSE, ...)

  mapping <- c(mapping, aes(x=x, y=y))
  class(mapping) <- "uneval"

  layer(
    data = data, mapping = mapping, stat = StatFilter, geom = GeomShadowText,
    position = position, show.legend = show.legend, inherit.aes = FALSE,
    params = params
  )
}

#' geom_node_rect
#' 
#' add rectangular shapes to ggplot2 using GeomRect
#' @param mapping aes mapping
#' @param data data to plot
#' @param position positional argument
#' @param show.legend whether to show legend
#' @param ... passed to `params` in `layer()` function
#' @return geom
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
#' @return ggplot2 object
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
#' @return ggplot2 object
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
    plot <- plot + geom_node_rect(
      aes(filter=.data$undefined & .data$type %in% type),
      fill=object$rect_fill, color="black")

  } else {
    plot <- plot + geom_node_rect(aes(filter=.data$type %in% type),
      fill=object$rect_fill, color="black")
  }
  plot
}


#' plot_kegg_network
#' 
#' plot the output of network_graph
#' 
#' @param g graph object returned by `network()`
#' @return ggplot2 object
#' @export
plot_kegg_network <- function(g) {
  ## [TODO] Presuming G***** and CS***** is not in the symbol
  gg <- g |> as_tbl_graph() |> activate(nodes) |>
  mutate(splitn=strsplit(name,"_") |> sapply("[",1)) |>
  mutate(group=startsWith(splitn,"G") & nchar(splitn)==6,
    and_group=startsWith(splitn,"CS") & nchar(splitn)==6)
  ggraph(gg, layout="kk") +
    geom_edge_link(aes(label=type,
                       filter=!startsWith(type,"in")),
                   angle_calc="along", force_flip=FALSE,
                   label_dodge = unit(2, 'mm')) +
    geom_edge_link(aes(filter=startsWith(type,"in_and")))+ 
    geom_node_point(size=4, aes(filter=!startsWith(name,"BLOCK") &
                                  !(group)&
                                  !(and_group))) + 
    geom_node_point(size=2, shape=21, aes(filter=startsWith(name,"BLOCK"))) + 
    geom_node_point(size=2, shape=21, aes(filter=(group) |
                                            (and_group))) + 
    geom_node_text(aes(label=name,
                       filter=!(and_group) &
                        !(group) &
                         !startsWith(name,"BLOCK")),
                   repel=TRUE, size=4, bg.colour="white")+
    theme_void()
}


#' geom_kegg
#' 
#' convenient function for plotting KEGG pathway graph
#' add geom_node_rect, geom_node_text and geom_edge_link
#' @param edge_color color attribute to edge
#' @param group_color border color for group node rectangles
#' @export
#' @return ggplot2 object
geom_kegg <- function(edge_color=NULL,
                      node_label=name,
                      group_color="red") {
  structure(list(edge_color=edge_color,
                 node_label=enquo(node_label),
                 group_color=group_color),
            class = "geom_kegg")
}

#' ggplot_add.geom_kegg
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_kegg
#' @return ggplot2 object
#' @export
ggplot_add.geom_kegg <- function(object, plot, object_name) {
  plot <- plot + 
    geom_edge_link(width=0.5,
                   arrow = arrow(length = unit(1, 'mm')), 
                   start_cap = square(1, 'cm'),
                   end_cap = square(1.5, 'cm'))
  plot <- plot+ geom_node_rect(aes(filter=.data$type=="group"),
                       fill="transparent", color=object$group_color)
  plot <- plot + geom_node_rect(aes(fill=I(bgcolor),
                                filter=bgcolor!="none" & .data$type!="group"))
  plot <- plot+
    geom_node_text(aes(label=!!object$node_label,
                       filter=type!="group"), family="serif", size=2)+
    theme_void()
  
}