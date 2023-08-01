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
#' @export
#' @return boolean vector
#' @examples
#' graph <- create_test_pathway()
#' graph <- graph |> mutate(hl=highlight_set_nodes(c("hsa:51428")))
highlight_set_nodes <- function(set, how="all",
    name="name", sep=" ", no_sep=FALSE) {
    graph <- .G()
    x <- get.vertex.attribute(graph, name)
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