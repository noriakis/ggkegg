#' ggkegg
#' 
#' main function parsing KEGG pathway data,
#' making igraph object and passing it to ggraph.
#' 
#' @param pid KEGG Pathway id e.g. hsa04110
#' @param pathway_number pathway number if passing enrichResult
#' @param layout default to "native", using KGML positions
#' @param return_igraph return the resulting igraph object
#' @param return_tbl_graph return the resulting tbl_graph object
#' (override `return_igraph` argument)
#' @param delete_undefined delete the undefined nodes from graph
#' default to FALSE, which preserves nodes but 
#' add `undefined` attribute to graph
#' @param convert_org these organism names are fetched from REST API
#' and cached, and used to convert the KEGG identifiers.
#' e.g. c("hsa", "compound")
#' @param convert_first after converting, take the first element as 
#' node name when multiple genes are listed in the node
#' @param convert_collapse if not NULL, collapse 
#' the gene names by this character
#' when multiple genes are listed in the node.
#' @param convert_reaction reaction name (graph attribute `reaction`) 
#' will be converted to reaction formula
#' @param delete_undefined delete `undefined` node specifying group,
#' should be set to `TRUE` when the layout is not from native KGML.
#' @param delete_zero_degree delete nodes with zero degree,
#' default to FALSE
#' @param module_type specify which module attributes to obtain
#' (definition or reaction)
#' @param module_definition_type `text` or `network` 
#' when parsing module definition.
#' If `text`, return ggplot object. If `network`, return `tbl_graph`.
#' @param numeric_attribute named vector for appending numeric attribute
#' @param node_rect_nudge parameter for nudging the node rect
#' @param group_rect_nudge parameter for nudging the group node rect
#' @examples
#' ## Use pathway ID to obtain `ggraph` object directly.
#' g <- ggkegg("hsa04110")
#' g + geom_node_rect()
#' @import ggraph
#' @import ggplot2
#' @importFrom tidygraph as_tbl_graph
#' @importFrom igraph induced.subgraph delete_vertex_attr
#' @importFrom methods new
#' @export
#' @return ggplot2 object
ggkegg <- function(pid,
    layout="native",
    return_igraph=FALSE,
    return_tbl_graph=FALSE,
    pathway_number=1,
    convert_org=NULL,
    convert_first=TRUE,
    convert_collapse=NULL,
    convert_reaction=FALSE,
    delete_undefined=FALSE,
    delete_zero_degree=FALSE,
    numeric_attribute=NULL,
    node_rect_nudge=0,
    group_rect_nudge=2,
    module_type="definition",
    module_definition_type="text") {
    
    if (!is.character(pid)) {
        if (attributes(pid)$class == "enrichResult") {
            org <- attributes(pid)$organism
            res <- attributes(pid)$result
            if (org != "UNKNOWN") {
                enrich_attribute <- paste0(org,
                    ":",
                    unlist(strsplit(res[pathway_number,]$geneID, "/"))
                )
            } else {
                enrich_attribute <- unlist(
                    strsplit(res[pathway_number,]$geneID, "/")
                )     
            }
            pid <- res[pathway_number,]$ID
        }
    } else {
        enrich_attribute <- NULL
    }
    

    if (is.character(pid)) {
        if (startsWith(pid, "M")) {
            mod <- module(pid)
            if (module_type == "definition") {
                if (module_definition_type == "text") {
                    plot_list <- module_text(mod, candidate_ko=enrich_attribute)
                    return(plot_module_text(plot_list))
                } else if (module_definition_type == "network") {
                    return(obtain_sequential_module_definition(mod))
                } else {
                  stop("Please specify `network` or `text`",
                  " to module_definition_type")
                }
            } else if (module_type == "reaction") {
                return(mod@reaction_graph)
            } else {
                stop("Please specify `reaction` or `definition`",
                " to module_type")
            }
        }
        if (startsWith(pid, "N")) {
            network <- network(pid)
            return(network |> network_graph() |> plot_kegg_network())
        }
    }
    ## If not module or enrichResult, return pathway
    g <- pathway(pid=pid,
                node_rect_nudge=node_rect_nudge,
                group_rect_nudge=group_rect_nudge,
                return_tbl_graph=FALSE)
  
    ## This part may be redundant, use `convert_id`
    if (!is.null(convert_org)) {
        convert_vec <- lapply(convert_org, function(co) {
               obtain_map_and_cache(co, pid)                
        }) |> unlist()
    
        V(g)$converted_name <- unlist(lapply(V(g)$name,
            function(x) {
                inc_genes <- unlist(strsplit(x, " "))
                conv_genes <- vapply(inc_genes, function(inc) {
                    convs <- convert_vec[inc]
                    if (is.na(convs)) {
                        return(x)
                    } else {
                        return(convs)
                    }
                }, FUN.VALUE="a")
                if (convert_first) {
                    conv_genes[1]
                } else {
                    paste(conv_genes, collapse=convert_collapse)
                }
            }
        ))
    }

    if (!is.null(numeric_attribute)){
        V(g)$numeric_attribute <- numeric_attribute[V(g)$name]
    }

    if (!is.null(enrich_attribute)) {
        bools <- vapply(V(g)$name, function(xx) {
            in_node <- strsplit(xx, " ") |> unlist() |> unique()
            if (length(intersect(in_node, enrich_attribute)) >= 1) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        }, FUN.VALUE=TRUE)
        V(g)$enrich_attribute <- bools
    }

    if (delete_undefined) {
        g <- induced.subgraph(g, !V(g)$name %in% "undefined")
    } else {
        V(g)$undefined <- V(g)$name %in% "undefined"
    }
    if (delete_zero_degree) {
        g <- induced.subgraph(g, degree(g)!=0)
    }

    if (convert_reaction) {
        convert_vec <- obtain_map_and_cache("reaction",NULL)
        V(g)$converted_reaction <- unlist(lapply(V(g)$reaction,
            function(x) {
                inc_genes <- unlist(strsplit(x, " "))
                conv_genes <- vapply(inc_genes, function(inc) {
                    convs <- convert_vec[inc]
                    if (is.na(convs)) {
                        return(x)
                    } else {
                        return(convs)
                    }
                }, FUN.VALUE="a")
                if (convert_first) {
                    conv_genes[1]
                } else {
                    paste(conv_genes, collapse=convert_collapse)
                }
            }
        ))
    }
    
    if (return_tbl_graph) {
        return(as_tbl_graph(g))
    }
    if (return_igraph) {
        return(g)
    }
    if (layout == "native") {
        ggraph(g, layout="manual", x=.data$x, y=.data$y)
    } else {
        g <- delete_vertex_attr(g, "x")
        g <- delete_vertex_attr(g, "y")
        ggraph(g, layout=layout)
    }
}


#' rawMap
#' 
#' given enrichResult class object,
#' return the ggplot object with raw KEGG map overlaid on
#' enriched pathway. Can be used with the function such as 
#' `clusterProfiler::enrichKEGG` and `MicrobiomeProfiler::enrichKO()`
#' 
#' @param enrich enrichResult or gseaResult class object, or list of them
#' @param pathway_number pathway number sorted by p-values
#' @param pid pathway id, override pathway_number if specified
#' @param fill_color color for genes
#' @param white_background fill background color white
#' @param how how to match the node IDs with the queries 'any' or 'all'
#' @export
#' @examples
#' if (require("clusterProfiler")) {
#'     cp <- enrichKEGG(c("1029","4171"))
#'     ## Multiple class object can be passed by list
#'     rawMap(list(cp,cp), pid="hsa04110")
#' }
#' @return ggraph with overlaid KEGG map
#' 
rawMap <- function(enrich, pathway_number=1, pid=NULL,
    fill_color="red", how="any", white_background=TRUE) {
    
    number <- length(enrich)
    if (length(fill_color) != number) {
        qqcat("Length of fill_color and enrich mismatches,",
            " taking first color\n")
        fill_color <- rep(fill_color[1], number)
    }
    if (is.list(enrich)) {
        if (is.null(pid)) {stop("Please specify pathway id.")}
    } else {
        if (attributes(enrich)$class == "enrichResult") {
            res <- attributes(enrich)$result
            if (is.null(pid)) {
                pid <- res[pathway_number, ]$ID
            }
        } else if (attributes(enrich)$class == "gseaResult") {
            res <- attributes(enrich)$result
            if (is.null(pid)) {
                pid <- res[pathway_number, ]$ID
            }         
        } else {
            stop("Please provide enrichResult")      
        }
    }
    ## For MicrobiomeProfiler
    if (startsWith(pid, "map")) {
        qqcat("Changing prefix of pathway ID from map to ko\n")
        pid <- gsub("map","ko",pid)
    }

    if (number == 1) {
        g <- pathway(pid) |> mutate(cp=append_cp(enrich, how=how, pid=pid))
        gg <- ggraph(g, layout="manual", x=.data$x, y=.data$y)+
            geom_node_rect(fill=fill_color, aes(filter=.data$cp))+
            overlay_raw_map()+theme_void()
    } else {
        g <- pathway(pid)
        for (i in seq_len(number)) {
            g <- g  |> mutate(!!paste0("cp",i) :=append_cp(enrich[[i]],
                how=how, pid=pid))
        }
        V(g)$space <- V(g)$width/number
        gg <- ggraph(g, layout="manual", x=.data$x, y=.data$y)
        nds <- g |> activate("nodes") |> data.frame()
        for (i in seq_len(number)) {
            gg <- gg +
                geom_node_rect(fill=fill_color[i],
                    data=nds[nds[[paste0("cp",i)]], ],
                    xmin=nds[nds[[paste0("cp",i)]], ]$xmin+
                        nds[nds[[paste0("cp",i)]], ]$space*(i-1)
                )
        }
        gg <- gg + overlay_raw_map()+theme_void()
    }
    if (white_background) {
        gg + theme(panel.background=element_rect(fill='white', colour='white'))
    } else {
        gg
    }
}



#' rawValue
#' 
#' given named vector of quantitative values,
#' return the ggplot object with raw KEGG map overlaid.
#' Colors can be changed afterwards.
#' 
#' @param values named vector, or list of them
#' @param pid pathway id
#' @param column column name on node table of the graph
#' @param white_background fill background color white
#' @param how how to match the node IDs with the queries 'any' or 'all'
#' @param auto_add automatically add prefix based on pathway prefix
#' @param man_graph provide manual tbl_graph
#' @param show_type type to be shown
#' typically, "gene", "ortholog", or "compound"
#' @export
#' @examples
#' ## Colorize by passing the named vector of numeric values
#' rv <- rawValue(c(1.1) |> setNames("hsa:6737"), 
#'         man_graph=create_test_pathway())
#' @return ggraph with overlaid KEGG map
#' 
rawValue <- function(values, pid=NULL, column="name", show_type="gene",
    how="any", white_background=TRUE, auto_add=FALSE, man_graph=NULL) {
    if (is.list(values)) {
        number <- length(values)
        if (auto_add) {
            pref <- gsub("[^a-zA-Z]", "", pid)
            for (i in seq_along(values)) {
                names(values[[i]]) <- paste0(pref, ":", names(values[[i]]))
            }
        }
    } else {
        number <- 1
        if (auto_add) {
            pref <- gsub("[^a-zA-Z]", "", pid)
            names(values) <- paste0(pref, ":", names(values))
        }
    }
    if (!is.null(man_graph)) {
        pgraph <- man_graph
    } else {
        pgraph <- pathway(pid)
    }
    if (number == 1) {
        g <- pgraph |> mutate(value=node_numeric(values,
            name=column, how=how))
        gg <- ggraph(g, layout="manual", x=.data$x, y=.data$y)+
            geom_node_rect(aes(fill=.data$value,
                                filter=.data$type %in% show_type))+
            overlay_raw_map()+theme_void()
    } else {
        ## Add new scales like ggh4x
        g <- pgraph
        for (i in seq_len(number)) {
            g <- g  |> mutate(!!paste0("value",i):=node_numeric(values[[i]],
                name=column,how=how))
        }
        V(g)$space <- V(g)$width/number
        gg <- ggraph(g, layout="manual", x=.data$x, y=.data$y)
        nds <- g |> activate("nodes") |> data.frame()
        nds <- nds[nds$type %in% show_type,]

        for (i in seq_len(number)) {
            nudge <- i-1
      
            gg <- gg + geom_node_rect(
                aes(fill=!!sym(paste0("value",i)),
                    filter=.data$type %in% show_type),
                    xmin=nds$xmin+nds$space*nudge,
                    xmax=nds$xmin+i*nds$space
                )
        }
        gg <- gg + overlay_raw_map()+theme_void()
    }
    if (white_background) {
        gg + theme(panel.background=element_rect(fill='white', colour='white'))
    } else {
        gg
    }
}