#' highlight_entities
#' 
#' highlight the entities in the pathway,
#' overlay raw map and return the results.
#' Note that highlighted nodes are considered to be rectangular,
#' so it is not compatible with the type like `compound`.
#' 
#' @param pathway pathway ID to be passed to `pathway()`
#' @param set vector of identifiers, or named vector of numeric values
#' @param num_combine combining function if multiple hits are obtained per node
#' @param how if `all`, if node contains multiple
#' IDs separated by `sep`, highlight if all the IDs
#' are in query. if `any`, highlight if one of the IDs
#' is in query.
#' @param name which column to search for
#' @param sep separater for node names
#' @param no_sep not separate node name
#' @param show_type entitie type, default to 'gene'
#' @param fill_color highlight color, default to 'tomato'
#' @param legend_name legend name, NULL to suppress
#' @param use_cache use cache or not
#' @param return_graph return tbl_graph instead of plot
#' @param remove_dot remove the "..." in the graphics name column
#' @param directory directroy with XML files. ignore caching when specified.
#' @return overlaid map
#' @examples
#' highlight_entities("hsa04110", c("CDKN2A"), legend_name="interesting")
#' @export
#'
highlight_entities <- function(pathway, set, how="any",
	num_combine=mean, name="graphics_name", sep=", ", no_sep=FALSE,
	show_type="gene", fill_color="tomato", remove_dot=TRUE,
	legend_name=NULL, use_cache=FALSE, return_graph=FALSE, directory=NULL) {
	graph <- pathway(pathway, use_cache=use_cache, directory=directory)
	x <- get.vertex.attribute(graph, name)
	
	if (is.null(names(set))) {## Discrete
	    vec <- vapply(seq_along(x), function(xn) {
	        if (no_sep) {
	            nn <- x[xn]
	        } else {
	            nn <- unlist(strsplit(x[xn], sep)) |> unique()
	        }
            if (remove_dot) {
                nn <- strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
            }
	        if (how == "all") {
	            if (length(intersect(nn, set)) == length(nn)) {
	                return(TRUE)
	            } else {
	                return(FALSE)
	            }
	        } else {
	            if (length(intersect(nn, set)) >= 1) {
	                return(TRUE)
	            } else {
	                return(FALSE)
	            }      
	        }
	    }, FUN.VALUE=TRUE)
	    graph <- graph |> mutate(highlight=vec)
	    if (return_graph) {return(graph)}
	    res <- ggraph(graph, layout="manual", x=.data$x, y=.data$y) + 
	        geom_node_rect(aes(filter=.data$type %in% show_type,
	            fill=.data$highlight))+
	        scale_fill_manual(values=c("grey", fill_color), name=legend_name)+
	        overlay_raw_map()+
	        theme_void()
	    if (is.null(legend_name)) {
	    	res <- res + theme(legend.position="none")    	
	    }
	} else {## Numeric
		vec <- lapply(seq_along(x), function(xn) {
			if (no_sep) {
				nn <- x[xn]
			} else {
				nn <- unlist(strsplit(x[xn], sep)) |> unique()
			}
            if (remove_dot) {
                nn <- strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
            }
            thresh <- ifelse(how=="any", 1, length(nn))
            if (length(intersect(names(set), nn)) >= thresh) {
                summed <- do.call(num_combine,
                    list(x=set[intersect(names(set), nn)]))
            } else {
                summed <- NA
            }
		}) |> unlist()
	    graph <- graph |> mutate(highlight=vec)
	    if (return_graph) {return(graph)}
	    res <- ggraph(graph, layout="manual", x=.data$x, y=.data$y) + 
	        geom_node_rect(aes(filter=.data$type %in% show_type,
	            fill=.data$highlight))+
	        scale_fill_continuous(name=legend_name)+
	        overlay_raw_map()+
	        theme_void()
	    if (is.null(legend_name)) {
	    	res <- res + theme(legend.position="none")    	
	    }
	}
	return(res)
}



#' highlight_set_nodes
#' 
#' identify if nodes are involved in specific queriy.
#' if multiple IDs are listed after separation by `sep`,
#' only return TRUE if all the IDs are in the query.
#' 
#' @param set set of identifiers
#' @param how if `all`, if node contains multiple
#' IDs separated by `sep`, highlight if all the IDs
#' are in query. if `any`, highlight if one of the IDs
#' is in query.
#' @param name which column to search for
#' @param sep separater for node names
#' @param no_sep not separate node name
#' @param remove_dot remove "..." after graphics name column
#' @export
#' @return boolean vector
#' @examples
#' graph <- create_test_pathway()
#' ## Highlight set of nodes by specifying ID
#' graph <- graph |> mutate(hl=highlight_set_nodes(c("hsa:51428")))
#' 
#' ## node column can be specified by `name` argument
#' graph <- graph |> 
#'     mutate(hl=highlight_set_nodes(c("DDX41"), name="graphics_name"))
highlight_set_nodes <- function(set, how="all",
    name="name", sep=" ", no_sep=FALSE, remove_dot=TRUE) {
    graph <- .G()
    x <- get.vertex.attribute(graph, name)
    vec <- vapply(seq_along(x), function(xn) {
        if (no_sep) {
            nn <- x[xn]
        } else {
            nn <- unlist(strsplit(x[xn], sep))
        }
        if (remove_dot) {
            nn <- strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
        }
        if (how == "all") {
            if (length(intersect(nn, set)) == length(nn)) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        } else {
            if (length(intersect(nn, set)) >= 1) {
                return(TRUE)
            } else {
                return(FALSE)
            }      
        }
    }, FUN.VALUE=TRUE)
    if (length(unique(vec))==1) {
    	cat("None of the nodes (or all the nodes) was highlighted.\n")	
    }
    vec
}


#' highlight_set_edges
#' 
#' identify if edges are involved in specific query.
#' if multiple IDs are listed after separation by `sep`,
#' only return TRUE if all the IDs are in the query.
#' 
#' @param set set of identifiers
#' @param how if `all`, if node contains multiple
#' IDs separated by `sep`, highlight if all the IDs
#' are in query. if `any`, highlight if one of the IDs
#' is in query.
#' @param name which column to search for
#' @param sep separater for node names
#' @param no_sep not separate node name
#' @export
#' @return boolean vector
#' @examples
#' graph <- create_test_pathway()
#' 
#' ## Specify edge column by `name`
#' ## In this example, edges having `degradation` value in
#' ## `subtype_name` column will be highlighted
#' graph <- graph |> activate("edges") |>
#'     mutate(hl=highlight_set_edges(c("degradation"), name="subtype_name"))
#' 
highlight_set_edges <- function(set, how="all",
    name="name", sep=" ", no_sep=FALSE) {
    graph <- .G()
    x <- get.edge.attribute(graph, name)
    vec <- vapply(seq_along(x), function(xn) {
        if (no_sep) {
            nn <- x[xn]
        } else {
            nn <- unlist(strsplit(x[xn], sep))
        }
        if (how == "all") {
            if (length(intersect(nn, set)) == length(nn)) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        } else {
            if (length(intersect(nn, set)) >= 1) {
                return(TRUE)
            } else {
                return(FALSE)
            }      
        }
    }, FUN.VALUE=TRUE)
    vec
}


#' highlight_module
#' 
#' identify if edges are involved in module reaction, and whether 
#' linked compounds are involved in the reaction. It would not be exactly 
#' the same as KEGG mapper. For instance, `R04293` involved in `M00912` 
#' is not included in KGML of `ko01100`.
#' 
#' @param graph tbl_graph
#' @param kmo kegg_module class object which stores reaction
#' @param name which column to search for
#' @param sep separator for node names
#' @param verbose show messages or not
#' @importFrom data.table :=
#' @export
#' @return boolean vector
#' @examples
#' ## Highlight module within the pathway
#' graph <- create_test_pathway()
#' mo <- create_test_module()
#' graph <- graph |> highlight_module(mo)
#' 
highlight_module <- function(graph, kmo,
                            name="name",
                            sep=" ",
                            verbose=FALSE) {
    if (attributes(kmo)$class[1] != "kegg_module") {
        stop("Please provide kegg_module class object")
    }

    edge_df <- graph |> activate("edges") |> data.frame()
    node_df <- graph |> activate("nodes") |> data.frame()

    ## First identify edges of reaction
    einds <- rep(FALSE, E(graph) |> length())
    ninds <- rep(FALSE, V(graph) |> length())

    ## Obtain each raw reaction
    rea <- kmo@reaction_each_raw
    results <- lapply(seq_len(nrow(kmo@reaction_each_raw)), function(i) {
        left <- kmo@reaction_each[i,][1] |> 
          unlist() |> as.character() |> paste0("cpd:", ...=_)
        raw_reac_string <- rea[i,][2] |> 
          unlist() |> as.character()       
        reac_list <- kmo@reaction_each[i,][2] |> unlist() |> as.character()
        right <- kmo@reaction_each[i,][3] |> 
          unlist() |> as.character() |> paste0("cpd:", ...=_)
        if (verbose) {qqcat("Checking reaction: @{raw_reac_string}\n")}

        x <- get.edge.attribute(graph, "reaction")
        ## Store edge index that meet reaction
        ind <- lapply(seq_along(x), function(xn) {
            reac <- raw_reac_string
            rls <- rep(FALSE, length(reac_list))
            names(rls) <- reac_list
            ## reactions associated with the edge
            edge_reac <- x[xn] |> strsplit(" ") |> unlist()
            if (sum(is.na(edge_reac)) != length(edge_reac)) {
                ## strip rn::
                edge_reac <- edge_reac |> gsub("rn:", "", x=_)
                for (ed in edge_reac) {
                    if (ed %in% names(rls)) {
                        rls[ed] <- TRUE
                    }
                }
                for (r in names(rls)) {
                    reac <- gsub(r, rls[r], reac)
                }
                reac <- gsub(",", "|", gsub("\\+", "&", reac))
                ## Eval boolean or length interpretation
                if (eval(parse(text=reac))) {
                    cand_node_ids <- edge_df[xn,]$orig.id
                    cand_node_ids <- cand_node_ids[!is.na(cand_node_ids)] |> 
                                        unique()
                    if (length(cand_node_ids) >= 1) {
                        for (ni in cand_node_ids) {
                            edges_ind <- node_df[node_df$orig.id %in% ni,] |> 
                                        row.names()
                            tmp_edge_df <- edge_df[edge_df$from %in% edges_ind,]
                            tmp_edge_df_2 <- edge_df[edge_df$to %in% edges_ind,]

                            node1 <- tmp_edge_df_2$from
                            node2 <- tmp_edge_df$to

                            subst <- node_df[tmp_edge_df_2$from,]$name |> 
                                    strsplit(" ") |> unlist() |> unique()
                            prod <- node_df[tmp_edge_df$to,]$name |> 
                                    strsplit(" ") |> unlist() |> unique()
                        
                            ## reversible
                            if ((length(intersect(subst,
                                    left)) == length(left) &
                                length(intersect(prod,
                                    right)) == length(right)) | 
                                (length(intersect(subst,
                                    right)) == length(right) &
                                length(intersect(prod, 
                                    left)) == length(left))) {
                                    return(list("ind"=xn,
                                        "nind"=c(node1, node2)))                             
                            }
                        }
                    }
                } else {}
            } else {} ## if edge is reaction
        }) ## each edge
        list(lapply(ind, function(x) x[["ind"]]) |> unlist(),
            lapply(ind, function(x) x[["nind"]]) |> unlist())
    })
    
    all_inds <- lapply(results, function(x) x[[1]]) |> unlist()
    nind <- lapply(results, function(x) x[[2]]) |> unlist()

    einds[all_inds] <- TRUE
    ninds[unique(as.numeric(nind))] <- TRUE

    graph |>
      activate("edges") |>
      mutate(!!kmo@ID:=einds) |>
      activate("nodes") |>
      mutate(!!kmo@ID:=ninds)
}