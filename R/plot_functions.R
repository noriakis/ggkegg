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
#' @examples
#' ## Pass multiple pathway IDs
#' multi_pathway_native(list("hsa04110","hsa03460"))
#' 
multi_pathway_native <- function(pathways, row_num=2, return_list=FALSE) {
    plen <- length(pathways)
    
    if (plen %% 2) {
        col_num <- as.integer(plen / row_num)+1; addit <- plen %% row_num
    } else {
        col_num <- plen / row_num; addit <- 0
    }
    
    tot_row <- 1
    tot_col <- 1
    miny <- 0

    ## Preallocate
    gls <- vector(mode="list", length=plen)
    for (pp in seq_len(pathways |> length())) {
        g <- pathway(pathways[pp])
        g <- g |> mutate(x=(.data$x/max(.data$x)+tot_col-1),
                           y=.data$y/min(.data$y)+miny)
        gls[[pp]] <- g
        # edf <- g |> activate("nodes") |> data.frame()
        # miny <- miny - min(edf$y)
        tot_col <- tot_col + 1
        
        if (tot_col > col_num) {
            tot_col <- 1
            tot_row <- tot_row + 1
            miny <- miny - 1
        }    
        # if (tot_row > row_num) {
        #     tot_row <- 1
        #     tot_col <- tot_col + 1
        # }
    }
    if (return_list) {return(gls)}
    Reduce(graph_join, gls)
}

#' plot_module_text
#' 
#' plot the text representation of KEGG modules
#' 
#' @param plot_list the result of `module_text()`
#' @param show_name name column to be plotted
#' @importFrom tidygraph tbl_graph
#' @import patchwork
#' @return ggplot2 object
#' @export
#' @examples
#' 
#' mo <- create_test_module()
#' 
#' ## The output of `module_text` is used for `plot_module_text()`
#' tex <- module_text(mo)
#' plt <- plot_module_text(tex)
#' 
plot_module_text <- function(plot_list, show_name="name") {
    panel_list <- lapply(seq_along(plot_list), function(concat) {
        plot_list[[concat]]$name <- plot_list[[concat]][[show_name]]
        g <- tbl_graph(nodes=plot_list[[concat]])
        ggraph(g, x=.data$x, y=1) +
         geom_node_rect(aes(filter=.data$koflag),
          fill=plot_list[[concat]][plot_list[[concat]]$koflag,]$color,
          alpha=0.5, color="black")+
         geom_node_rect(aes(filter=!.data$koflag & !.data$conflag),
          fill="transparent", color="black")+
         geom_node_text(
            aes(label=.data$name,filter=.data$koflag | .data$conflag))+
         theme_void()       
    })
    wrap_plots(panel_list, ncol=1)
}


#' plot_module_blocks
#' 
#' wrapper function for plotting network representation of 
#' module definition blocks
#' 
#' @param all_steps the result of `obtain_sequential_module_definition()`
#' @param layout ggraph layout parameter
#' @export
#' @return ggplot2 object
#' @examples
#' mo <- create_test_module()
#' ## The output of `obtain_sequential_module_definition` 
#' ## is used for `plot_module_blocks()`
#' sequential_mod <- obtain_sequential_module_definition(mo)
#' plt <- plot_module_blocks(sequential_mod)
plot_module_blocks <- function(all_steps, layout="kk") {
    allnodes <- unique(V(all_steps)$name)
    if (sum(startsWith(allnodes, "K"))==length(allnodes)) {
        stop("all nodes are KO.")
    }
    ggraph(all_steps, layout=layout) +
        geom_edge_link(aes(filter=.data$type %in% c("block_transition","rel")),
                        end_cap=circle(5, 'mm'),start_cap=circle(5,"mm"),
                        color="red")+
        geom_edge_link(aes(filter=!.data$type %in% 
                            c("block_transition","rel","in_block"))) + 
        geom_edge_link(aes(label=.data$type,
                            filter=!startsWith(.data$type,"in") & 
                            !.data$type %in% c("block_transition","rel")),
                        angle_calc="along",
                        label_dodge=unit(2, 'mm')) + 
        geom_node_point(size=4,
            aes(filter=!startsWith(.data$name,"manual_BLOCK") &
                !startsWith(.data$name,"manual_G") &
                !startsWith(.data$name,"manual_CS"))) + 
        geom_node_point(size=2, shape=21,
            aes(filter=startsWith(.data$name,"manual_BLOCK"))) + 
        geom_node_point(size=2, shape=21,
            aes(filter=startsWith(.data$name,"manual_CS") |
                    startsWith(.data$name,"manual_G"))) + 
        geom_node_text(aes(label=.data$name,
                        filter=startsWith(.data$name,"K")),
                        repel=TRUE, size=4, bg.colour="white")+
        theme_void()
}

#' geom_node_shadowtext
#' 
#' Plot shadowtext at node position, use StatFilter in ggraph
#' 
#' @export
#' @param mapping aes mapping
#' @param data data to plot
#' @param position positional argument
#' @param show.legend whether to show legend
#' @param ... passed to `params` in `layer()` function
#' @return geom
#' @importFrom shadowtext GeomShadowText
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  geom_node_shadowtext(aes(label=name))
geom_node_shadowtext <- function(mapping=NULL, data=NULL,
                           position='identity',
                           show.legend=NA, ...) {
    params <- list(na.rm=FALSE, ...)

    mapping <- c(mapping, aes(x=.data$x, y=.data$y))
    class(mapping) <- "uneval"

    layer(
        data=data, mapping=mapping, stat=StatFilter, geom=GeomShadowText,
        position=position, show.legend=show.legend, inherit.aes=FALSE,
        params=params
    )
}

#' geom_node_rect
#' 
#' Plot rectangular shapes to ggplot2 using GeomRect,
#' using StatFilter in ggraph
#' 
#' @param mapping aes mapping
#' @param data data to plot
#' @param position positional argument
#' @param show.legend whether to show legend
#' @param ... passed to `params` in `layer()` function
#' @return geom
#' @export
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  geom_node_rect()
geom_node_rect <- function(mapping=NULL, data=NULL, position='identity',
                            show.legend=NA, ...) {
    mapping1 <- mapping
    raw_mapping <- aes(xmin=.data$xmin, ymin=.data$ymin, xmax=.data$xmax, ymax=.data$ymax)
    mapping <- c(as.list(mapping1), raw_mapping[!names(raw_mapping) %in% names(mapping1)])
    class(mapping) <- "uneval"
    layer(
        data=data, mapping=mapping, stat=StatFilter, geom=GeomRect,
        position=position, show.legend=show.legend, inherit.aes=FALSE,
        params=list(na.rm=FALSE, ...)
  )
}


#' geom_node_rect_kegg
#' 
#' Wrapper function for plotting a certain type of nodes
#' with background color with geom_node_rect()
#' 
#' @param type type to be plotted (gene, map, compound ...)
#' @param rect_fill rectangular fill
#' @export
#' @return ggplot2 object
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  geom_node_rect_kegg(type="gene")
geom_node_rect_kegg <- function(type=NULL, rect_fill="grey") {
    structure(list(type=type, rect_fill=rect_fill),
        class="geom_node_rect_kegg")
}

#' ggplot_add.geom_node_rect_kegg
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_node_rect_kegg
#' @export
#' @return ggplot2 object
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  geom_node_rect_kegg(type="gene")
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
#' @param layout layout to be used, default to nicely
#' @return ggplot2 object
#' @export
#' @examples
#' ne <- create_test_network()
#' ## Output of `network_graph` must be used with plot_kegg_network
#' neg <- network_graph(ne)
#' plt <- plot_kegg_network(neg)
plot_kegg_network <- function(g, layout="nicely") {
    gg <- g |> as_tbl_graph() |> activate("nodes") |>
        mutate(splitn=strsplit(.data$name,"_") |> 
                        vapply("[",1,FUN.VALUE="character")) |>
        mutate(group=startsWith(.data$splitn,"manual_G"),
            and_group=startsWith(.data$splitn,"manual_CS"))

    ggraph(gg, layout=layout) +
        geom_edge_link(aes(label=.data$type,
                       filter=!startsWith(.data$type,"in")),
                   angle_calc="along", force_flip=FALSE,
                   label_dodge=unit(2, 'mm')) +
        geom_edge_link(aes(filter=startsWith(.data$type,"in_and")))+ 
        geom_edge_link(aes(filter=startsWith(.data$type,"in_block")),
            linetype=2)+ 
        geom_node_point(size=4,
            aes(filter=!startsWith(.data$name,"manual_BLOCK") &
                !(.data$group)&
                !(.data$and_group))) + 
        geom_node_point(size=2, shape=21,
            aes(filter=startsWith(.data$name,"manual_BLOCK"))) + 
        geom_node_point(size=2, shape=21,
            aes(filter=(.data$group) | (.data$and_group))) + 
        geom_node_text(aes(label=.data$name,
            filter=!startsWith(.data$name,"manual_")),
            repel=TRUE, size=4, bg.colour="white") +
        theme_void()
}


#' geom_node_rect_multi
#' 
#' Wrapper function for plotting multiple rects
#' with background color with geom_node_rect().
#' All columns should belong to the same scale when
#' using `asIs=FALSE`. If you need multiple scales for 
#' each element, please use `ggh4x::scale_fill_multi`
#' for each.
#' 
#' @param ... color columns
#' @param asIs treat the color as is or not
#' @export
#' @return ggplot2 object
#' @examples
#' plt <- create_test_pathway() %>% ggraph() + geom_node_rect_multi(bgcolor)
geom_node_rect_multi <- function(..., asIs=TRUE) {
  color_cols <- as.character(ensyms(...))
  structure(list(cols=color_cols, asIs=asIs),
            class="geom_node_rect_multi")}
  
#' ggplot_add.geom_node_rect_multi
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_node_rect_multi
#' @export
#' @return ggplot2 object
#' @examples
#' plt <- create_test_pathway() %>% ggraph() + geom_node_rect_multi(bgcolor)
ggplot_add.geom_node_rect_multi <- function(object, plot, object_name) {
  colnum <- length(object$cols)
  if (length(colnum)==0) {stop("Please specify at least one color column")}
  plot$data$space <- plot$data$width / colnum
  for (i in seq_len(colnum)) {
    if (object$asIs) {
      plot <- plot + 
        geom_node_rect(
          aes(xmin= .data$xmin + .data$space*!!(i-1),
              xmax= .data$xmin + .data$space*!!(i),
              fill= I(.data[[object$cols[i]]])))
    } else {
      plot <- plot + 
        geom_node_rect(
          aes(xmin= .data$xmin + .data$space*!!(i-1),
              xmax= .data$xmin + .data$space*!!(i),
              fill= .data[[object$cols[i]]]))
    }
  }
  plot
}


#' geom_kegg
#' 
#' Wrapper function for plotting KEGG pathway graph
#' add geom_node_rect, geom_node_text and geom_edge_link simultaneously
#' 
#' @param edge_color color attribute to edge
#' @param group_color border color for group node rectangles
#' @param node_label column name for node label
#' @param parallel use geom_edge_parallel() instead of geom_edge_link()
#' @export
#' @examples 
#' test_pathway <- create_test_pathway()
#' p <- ggraph(test_pathway, layout="manual", x=x, y=y)+
#' geom_kegg()
#' @return ggplot2 object
geom_kegg <- function(edge_color=NULL,
                    node_label=.data$name,
                    group_color="red",
                    parallel=FALSE) {
    structure(list(edge_color=edge_color,
                node_label=enquo(node_label),
                group_color=group_color,
                parallel=parallel),
            class="geom_kegg")
}

#' ggplot_add.geom_kegg
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_kegg
#' @return ggplot2 object
#' @export
#' @examples 
#' test_pathway <- create_test_pathway()
#' p <- ggraph(test_pathway, layout="manual", x=x, y=y)+
#' geom_kegg()
ggplot_add.geom_kegg <- function(object, plot, object_name) {
  if (object$parallel) {
        plot <- plot + 
            geom_edge_parallel(width=0.5,
                        arrow=arrow(length=unit(1, 'mm')), 
                        start_cap=square(1, 'cm'),
                        end_cap=square(1.5, 'cm'))    
    } else {
        plot <- plot +  
            geom_edge_link(width=0.5,
                        arrow=arrow(length=unit(1, 'mm')), 
                        start_cap=square(1, 'cm'),
                        end_cap=square(1.5, 'cm'))
    }

    plot <- plot + geom_node_rect(aes(filter=.data$type=="group"),
                        fill="transparent", color=object$group_color)
    plot <- plot + geom_node_rect(aes(fill=I(.data$bgcolor),
                        filter=.data$bgcolor!="none" & .data$type!="group"))
    plot <- plot +
        geom_node_text(aes(label=!!object$node_label,
                        filter=.data$type!="group"), family="serif", size=2) +
        theme_void()
  
}