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
#' nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
#'                    x=c(1,1),
#'                    xmin=c(-1,-1),
#'                    xmax=c(2,2),
#'                    y=c(1,1),
#'                    ymin=c(-1,-1),
#'                    ymax=c(2,2))
#' edges <- data.frame(from=1, to=2)
#' graph <- tbl_graph(nodes, edges)
#' graph <- graph |> mutate(hl=highlight_set_nodes(c("hsa:1029")))
highlight_set_nodes <- function(set, how="all",
  name="name", sep=" ", no_sep=FALSE) {
  graph <- .G()
  x <- get.vertex.attribute(graph, name)
  vec <- NULL
  for (xn in seq_along(x)) {
    if (no_sep) {
      nn <- x[xn]
    } else {
      nn <- unlist(strsplit(x[xn], sep))
    }
    if (how=="all") {
      if (length(intersect(nn, set))==length(nn)) {
        vec <- c(vec, TRUE)
      } else {
        vec <- c(vec, FALSE)
      }
    } else {
      if (length(intersect(nn, set))>=1) {
        vec <- c(vec, TRUE)
      } else {
        vec <- c(vec, FALSE)
      }      
    }
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
#' nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
#'                    x=c(1,1),
#'                    xmin=c(-1,-1),
#'                    xmax=c(2,2),
#'                    y=c(1,1),
#'                    ymin=c(-1,-1),
#'                    ymax=c(2,2))
#' edges <- data.frame(from=1, to=2, name="K00112")
#' graph <- tbl_graph(nodes, edges)
#' graph <- graph |> activate("edges") |>
#'     mutate(hl=highlight_set_edges(c("K00112")))
#' 
highlight_set_edges <- function(set, how="all",
  name="name", sep=" ", no_sep=FALSE) {
  graph <- .G()
  x <- get.edge.attribute(graph, name)
  vec <- NULL
  for (xn in seq_along(x)) {
    if (no_sep) {
      nn <- x[xn]
    } else {
      nn <- unlist(strsplit(x[xn], sep))
    }
    if (how=="all") {
      if (length(intersect(nn, set))==length(nn)) {
        vec <- c(vec, TRUE)
      } else {
        vec <- c(vec, FALSE)
      }
    } else {
      if (length(intersect(nn, set))>=1) {
        vec <- c(vec, TRUE)
      } else {
        vec <- c(vec, FALSE)
      }      
    }
  }
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
#' nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
#'                    x=c(1,1),
#'                    xmin=c(-1,-1),
#'                    xmax=c(2,2),
#'                    y=c(1,1),
#'                    ymin=c(-1,-1),
#'                    ymax=c(2,2))
#' edges <- data.frame(from=1, to=2)
#' graph <- tbl_graph(nodes, edges)
#' mo <- create_test_module()
#' graph <- graph |> highlight_module(mo)
#' 
highlight_module <- function(graph, kmo,
                            name="name",
                            sep=" ",
                            verbose=FALSE) {
    if (attributes(kmo)$class[1]!="kegg_module") {
      stop("Please provide kegg_module class object")}

    edge_df <- graph |> activate(edges) |> data.frame()
    node_df <- graph |> activate(nodes) |> data.frame()

    ## First identify edges of reaction
    einds <- rep(FALSE, E(graph) |> length())
    ninds <- rep(FALSE, V(graph) |> length())

    all_inds <- NULL
    nind <- NULL

    ## Obtain each raw reaction
    rea <- kmo@reaction_each_raw
    for (i in seq_len(nrow(kmo@reaction_each_raw))) {
        left <- kmo@reaction_each[i,][1] |> 
          unlist() |> as.character() |> paste0("cpd:", ...=_)
        raw_reac_string <- rea[i,][2] |> 
          unlist() |> as.character() #|> paste0("rn:", ...=_)       
        reac_list <- kmo@reaction_each[i,][2] |> unlist() |> as.character()
        right <- kmo@reaction_each[i,][3] |> 
          unlist() |> as.character() |> paste0("cpd:", ...=_)
        if (verbose) {qqcat("Checking reaction: @{raw_reac_string}\n")}
        
        x <- get.edge.attribute(graph, "reaction")
        ind <- NULL ## Store edge index that meet reaction
        for (xn in seq_along(x)) {

            reac <- raw_reac_string
            
            rls <- rep(FALSE, length(reac_list))
            names(rls) <- reac_list
            
            ## reactions associated with the edge
            edge_reac <- x[xn] |> strsplit(" ") |> unlist()
            
            if (sum(is.na(edge_reac))!=length(edge_reac)) {
                ## strip rn::
                edge_reac <- edge_reac |> gsub("rn:","",x=_)
                for (ed in edge_reac) {
                    if (ed %in% names(rls)) {
                        rls[ed] <- TRUE
                    }
                }
                for (r in names(rls)) {
                    reac <- gsub(r, rls[r], reac)
                }
                reac <- gsub(",","|",gsub("\\+","&",reac))
                
                ##[TODO] eval boolean or length interpretation
                if (eval(parse(text=reac))) {
                # if (length(intersect(edge_reac,reac_list))==length(edge_reac)) {
                    cand_node_ids <- edge_df[xn,]$orig.id
                    cand_node_ids <- cand_node_ids[!is.na(cand_node_ids)] |> 
                                      unique()
                    if (length(cand_node_ids)>=1) {
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
                        # print(edge_reac)
                        # print(paste("subst",subst, "modleft", left))
                        # print(paste("prod",prod, "modright",right))
                        
                        ## reversible
                        if ((length(intersect(subst, left))==length(left) &
                          length(intersect(prod, right))==length(right)) | 
                        (length(intersect(subst, right))==length(right) &
                          length(intersect(prod, left))==length(left))) {
                            ind <- c(ind, xn)
                            nind <- c(nind, node1, node2)                     
                        }
                      }
                    }
                    
                } else {}
            } else {} ## if edge is reaction
        } ## each edge
        all_inds <- c(all_inds, ind)
    }
    einds[all_inds] <- TRUE
    ninds[unique(as.numeric(nind))] <- TRUE

    graph |>
      activate(edges) |>
      mutate(!!kmo@ID := einds) |>
      activate(nodes) |>
      mutate(!!kmo@ID := ninds)
}



## highlight_module_reaction and highlight_module_compound is deprecated.
## As we don't use them individually outside global map
#
# #' highlight_module
# #' 
# #' 
# #' Using reaction and compound IDs in the module,
# #' highlight the interactions involved.
# #' Please see `highlight_module_reaction` and `highlight_module_compound`.
# #' 
# #' @param graph graph
# #' @param kmo kegg_module class object
# #' @param name which column to search for
# #' @param sep separater for node names
# #' @export
# highlight_module <- function(graph, kmo, name="name", sep=" ") {

#   if (attributes(kmo)$class[1]=="kegg_module") {
#     graph |>
#       activate(edges) |>
#       mutate(!!kmo@ID := highlight_module_reaction(kmo, name=name, sep=sep)) |>
#       activate(nodes) |>
#       mutate(!!kmo@ID := highlight_module_compound(kmo, kmo@ID, name=name, sep=sep))
#   } else {
#     stop("please provide kegg_module class")
#   }
# }
# #' highlight_module_reaction
# #' 
# #' identify if edges are involved in module reaction.
# #' please note that all the reactions satisfying queried module's reaction
# #' will be highlighted, thus it will not be exactly the same as KEGG mapper.
# #' 
# #' @param kmo kegg_module class object which stores reaction
# #' @param name which column to search for
# #' @param sep separator for node names
# #' @param verbose show messages or not
# #' @export
# highlight_module_reaction <- function(kmo,
#                                       name="name",
#                                       sep=" ",
#                                       verbose=FALSE) {
#     graph <- .G()
#     einds <- rep(FALSE, E(graph) |> length())
#     all_inds <- NULL
    
#     if (attributes(kmo)$class[1]=="kegg_module") {
#         ## Obtain each raw reaction
#         rea <- kmo@reaction_each_raw
#         for (i in seq_len(nrow(kmo@reaction_each_raw))) {
#             left <- rea[i,][1] |> unlist() |> as.character() #|> paste0("cpd:", ...=_)
#             raw_reac_string <- rea[i,][2] |> unlist() |> as.character() #|> paste0("rn:", ...=_)       
#             reac_list <- kmo@reaction_each[i,][2] |> unlist() |> as.character()
#             right <- rea[i,][3] |> unlist() |> as.character()#|> paste0("cpd:", ...=_)
#             if (verbose) {qqcat("Checking reaction: @{raw_reac_string}\n")}
            
#             x <- get.edge.attribute(graph, "reaction")
#             ind <- NULL ## Store edge index that meet reaction
#             for (xn in seq_along(x)) {

#                 reac <- raw_reac_string
                
#                 rls <- rep(FALSE, length(reac_list))
#                 names(rls) <- reac_list
                
#                 ## reactions associated with the edge
#                 edge_reac <- x[xn] |> strsplit(" ") |> unlist()
                
#                 if (sum(is.na(edge_reac))!=length(edge_reac)) {
#                     ## strip rn::
#                     edge_reac <- edge_reac |> gsub("rn:","",x=_)
                    
#                     for (ed in edge_reac) {
#                         if (ed %in% names(rls)) {
#                             rls[ed] <- TRUE
#                         }
#                     }
#                     for (r in names(rls)) {
#                         reac <- gsub(r, rls[r], reac)
#                     }
#                     reac <- gsub(",","|",gsub("\\+","&",reac))
                    
#                     ##[TODO] eval boolean or length interpretation
#                     if (eval(parse(text=reac))) {
#                     # if (length(intersect(edge_reac,reac_list))==length(edge_reac)) {
#                         ind <- c(ind, xn)
#                     } else {}
#                 } else {} ## if edge is reaction
#             } ## each edge
#             all_inds <- c(all_inds, ind)
#         }
#         einds[all_inds] <- TRUE
#         einds
#     } else {
#         stop("please provide kegg_module class")
#     }
# }



# #' highlight_module_compound
# #' 
# #' identify if nodes are involved in compounds in module reaction.
# #' please note that all the compounds satisfying queried module's reaction
# #' will be highlighted, thus it will not be exactly the same as KEGG mapper.
# #' 
# #' @param kmo kegg_module class object which stores reaction
# #' @param label edge label
# #' @param name which column to search for
# #' @param sep separator for node names
# #' @param verbose show messages or not
# #' @importFrom tidygraph activate
# #' @importFrom dplyr mutate
# #' @export
# highlight_module_compound <- function(kmo,
#                                       label,
#                                       name="name",
#                                       sep=" ",
#                                       verbose=FALSE) {
#   graph <- .G()
#   ninds <- rep(FALSE, V(graph) |> length())
  
#   if (attributes(kmo)$class[1]=="kegg_module") {
#     return_nodes <- NULL
#     nodes <- graph |> activate(nodes) |> data.frame()
#     candidate <- graph |> activate(edges) |>
#       filter(.data[[label]] & is.na(name)) |>
#       mutate(node_to=nodes[to,"name"]) |>
#       data.frame()
#     print(candidate)
#     reacs <- NULL
#     for (i in seq(1,length(candidate$node_to),2)) {
#       reacs <- rbind(reacs, 
#                      c(candidate[i,"node_to"], ## substrate
#                        candidate[i, "reaction"],
#                        candidate[i+1, "node_to"], ## product
#                        candidate[i,"type"],
#                        candidate[i,"to"],
#                        candidate[i+1,"to"]))
#     }
#     for (i in seq_len(nrow(kmo@reaction_each))){
#       left <- kmo@reaction_each[i,][1] |> unlist() |> as.character() |> paste0("cpd:", ...=_)
#       right <- kmo@reaction_each[i,][3] |> unlist() |> as.character() |> paste0("cpd:", ...=_)
#       for (j in seq_len(nrow(reacs))) {
#         cand_left <- reacs[j,][1]
#         cand_right <- reacs[j,][3]
#         cand_id1 <- reacs[j,][5]
#         cand_id2 <- reacs[j,][6]
#         if (length(intersect(cand_left,left))==length(cand_left) & 
#             length(intersect(cand_right,right))==length(cand_right)) {
#           return_nodes <- c(return_nodes, cand_id1, cand_id2)
#         }
#         # if (length(intersect(cand_left,right))>0 & 
#         #     length(intersect(cand_right,left))>0) {
#         #   return_nodes <- c(return_nodes, cand_id1, cand_id2)
#         # }
#       }
#     }
#     ninds[unique(as.numeric(return_nodes))] <- TRUE
#     ninds
#   } else {
#     stop("please provide kegg_module class")
#   }
# }


## purrr versions
# #' highlight_module_compound
# #' 
# #' identify if nodes are involved in compounds in module reaction.
# #' please note that all the compounds satisfying queried module's reaction
# #' will be highlighted, thus it will not be exactly the same as KEGG mapper.
# #' 
# #' @param kmo kegg_module class object which stores reaction
# #' @param label edge label
# #' @param name which column to search for
# #' @param sep separator for node names
# #' @param verbose show messages or not
# #' @export
# highlight_module_compound <- function(kmo,
#                                       label,
#                                       name="name",
#                                       sep=" ",
#                                       verbose=FALSE) {
#   graph <- .G()
#   ninds <- rep(FALSE, V(graph) |> length())
  
#   if (attributes(kmo)$class[1]=="kegg_module") {
#     return_nodes <- NULL
#     nodes <- graph |> activate(nodes) |> data.frame()
#     candidate <- graph |> activate(edges) |>
#       filter(.data[[label]] & is.na(name)) |>
#       mutate(node_to=nodes[to,"name"], ) |>
#       data.frame()
#     reacs <- NULL
    
#     reac_modify <- function(i) {
#       data.frame(candidate[i,"node_to"],
#                        candidate[i, "reaction"],
#                        candidate[i+1, "node_to"],
#                        candidate[i,"type"],
#                        candidate[i,"to"],
#                        candidate[i+1,"to"])
#     }
    
    
#     reacs <- seq(1,length(candidate$node_to),2) |>
#       purrr::map(reac_modify) |>
#       purrr::list_rbind()
    
      
#     ret_node_ind <- function(i) {
#       left <- kmo@reaction_each[i,][1] |> 
#         unlist() |> as.character() |> paste0("cpd:", ...=_)
#       right <- kmo@reaction_each[i,][3] |> 
#         unlist() |> as.character() |> paste0("cpd:", ...=_)
      
#       ret_node_ind_nest <- function(j) {
#         return_nodes <- NULL
#         cand_left <- reacs[j,][1]
#         cand_right <- reacs[j,][3]
#         cand_id1 <- reacs[j,][5]
#         cand_id2 <- reacs[j,][6]
#         if (length(intersect(cand_left,left))>0 & 
#             length(intersect(cand_right,right))>0) {
#           return_nodes <- c(return_nodes, cand_id1, cand_id2)
#         }
#         if (length(intersect(cand_left,right))>0 & 
#             length(intersect(cand_right,left))>0) {
#           return_nodes <- c(return_nodes, cand_id1, cand_id2)
#         }
#         return_nodes
#       }
#       nrow(reacs) |> seq_len() |>
#         purrr::map(ret_node_ind_nest) |>
#         unlist()
#     }
      
      
#     ret_node_inds <- seq_len(nrow(kmo@reaction_each)) |>
#       purrr::map(ret_node_ind) |> unlist()
#     ninds[unique(as.numeric(ret_node_inds))] <- TRUE
#     ninds
#   } else {
#     stop("please provide kegg_module class")
#   }
# }



# #' highlight_module_reaction
# #' 
# #' identify if edges are involved in module reaction.
# #' please note that all the reactions satisfying queried module's reaction
# #' will be highlighted, thus it will not be exactly the same as KEGG mapper.
# #' 
# #' @param kmo kegg_module class object which stores reaction
# #' @param name which column to search for
# #' @param sep separator for node names
# #' @param verbose show messages or not
# #' @export
# highlight_module_reaction <- function(kmo,
#                                       name="name",
#                                       sep=" ",
#                                       verbose=FALSE) {
#   graph <- .G()
#   einds <- rep(FALSE, E(graph) |> length())
#   all_inds <- NULL
  
#   if (attributes(kmo)$class[1]=="kegg_module") {
#     each <- kmo@reaction_each
#     each_raw <- kmo@reaction_each_raw  
    
#     divide <- function(i) {
      
#       left <- each_raw[i,][1] |> unlist() |> as.character() 
#       raw_reac_string <- each_raw[i,][2] |> unlist() |> as.character()   
#       right <- each_raw[i,][3] |> unlist() |> as.character()
#       components <- each[i,][2] |> unlist() |> as.character()
      
#       rls <- rep(FALSE, length(components))
#       names(rls) <- components
      
#       check_bool <- function(x, rls) {
#         if (!is.na(x)) {
#           edge_reac <- edge_reacs[x] |> 
#             strsplit(" ") |>
#             unlist() |> 
#             gsub("rn:","",x=_)
#           rls[intersect(names(rls),edge_reac)]<-TRUE
#           raw_reac_string <- purrr::map(names(rls), function(x) {
#             raw_reac_string <- gsub(x, rls[x], raw_reac_string)
#           }) |> unlist()
#           bool <- raw_reac_string %>%
#             gsub("\\+","&",x=.) %>%
#             gsub(",","|",x=.)
#           eval(parse(text=bool))
#           }
#       }
      
#       edge_reacs <- get.edge.attribute(graph, "reaction")
      
#       ret_bool <- 
#         edge_reacs |>
#           length() |>
#           seq_len() |>
#           purrr::map(check_bool, rls) |>
#           unlist()
#       which(ret_bool)
#     }

#     inds <- nrow(kmo@reaction_each_raw) |>
#       seq_len() |>
#       purrr::map(divide) |> unlist()
    
#     einds[inds] <- TRUE
#     einds
#   } else {
#     stop("please provide kegg_module class")
#   }
# }
