#' highlight_set_nodes
#' 
#' identify if nodes are involved in specific queriy.
#' if multiple IDs are listed after separation by `sep`,
#' only return TRUE if all the IDs are in the query.
#' 
#' @param set set of identifiers
#' @param name which column to search for
#' @param sep separater for node names
#' @export
highlight_set_nodes <- function(set, name="name", sep=" ") {
  graph <- .G()
  x <- get.vertex.attribute(graph, name)
  vec <- NULL
  for (xn in seq_along(x)) {
    nn <- unlist(strsplit(x[xn], sep))
    if (length(intersect(nn, set))==length(nn)) {
      vec <- c(vec, TRUE)
    } else {
      vec <- c(vec, FALSE)
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
#' @param name which column to search for
#' @param sep separater for node names
#' @export
highlight_set_edges <- function(set, name="name", sep=" ") {
  graph <- .G()
  x <- get.edge.attribute(graph, name)
  vec <- NULL
  for (xn in seq_along(x)) {
    nn <- unlist(strsplit(x[xn], sep))
    if (length(intersect(nn, set))==length(nn)) {
      vec <- c(vec, TRUE)
    } else {
      vec <- c(vec, FALSE)
    }
  }
  vec
}



#' highlight_module
#' 
#' 
#' Using reaction and compound IDs in the module,
#' highlight the interactions involved.
#' Please see `highlight_module_reaction` and `highlight_module_compound`.
#' 
#' @param graph graph
#' @param kmo kegg_module class object
#' @param name which column to search for
#' @param sep separater for node names
#' @export
highlight_module <- function(graph, kmo, name="name", sep=" ") {

  if (attributes(kmo)$class[1]=="kegg_module") {
    graph |>
      activate(edges) |>
      mutate(!!kmo@ID := highlight_module_reaction(kmo, name=name, sep=sep)) |>
      activate(nodes) |>
      mutate(!!kmo@ID := highlight_module_compound(kmo, kmo@ID, name=name, sep=sep))
  } else {
    stop("please provide kegg_module class")
  }
}




#' highlight_module_reaction
#' 
#' identify if edges are involved in module reaction.
#' please note that all the reactions satisfying queried module's reaction
#' will be highlighted, thus it will not be exactly the same as KEGG mapper.
#' 
#' @param kmo kegg_module class object which stores reaction
#' @param name which column to search for
#' @param sep separator for node names
#' @param verbose show messages or not
#' @export
highlight_module_reaction <- function(kmo,
                                      name="name",
                                      sep=" ",
                                      verbose=FALSE) {
    graph <- .G()
    einds <- rep(FALSE, E(graph) |> length())
    all_inds <- NULL
    
    if (attributes(kmo)$class[1]=="kegg_module") {
        ## Obtain each raw reaction
        rea <- kmo@reaction_each_raw
        for (i in seq_len(nrow(kmo@reaction_each_raw))) {
            left <- rea[i,][1] |> unlist() |> as.character() #|> paste0("cpd:", ...=_)
            raw_reac_string <- rea[i,][2] |> unlist() |> as.character() #|> paste0("rn:", ...=_)       
            reac_list <- kmo@reaction_each[i,][2] |> unlist() |> as.character()
            right <- rea[i,][3] |> unlist() |> as.character()#|> paste0("cpd:", ...=_)
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
                    
                    if (eval(parse(text=reac))) {#(length(intersect(edge_reac,reac))==length(edge_reac)) {
                        ind <- c(ind, xn)
                    } else {}
                } else {} ## if edge is reaction
            } ## each edge
            all_inds <- c(all_inds, ind)
        }
        einds[all_inds] <- TRUE
        einds
    } else {
        stop("please provide kegg_module class")
    }
}



#' highlight_module_compound
#' 
#' identify if nodes are involved in compounds in module reaction.
#' please note that all the compounds satisfying queried module's reaction
#' will be highlighted, thus it will not be exactly the same as KEGG mapper.
#' 
#' @param kmo kegg_module class object which stores reaction
#' @param label edge label
#' @param name which column to search for
#' @param sep separator for node names
#' @param verbose show messages or not
#' @export
highlight_module_compound <- function(kmo,
                                      label,
                                      name="name",
                                      sep=" ",
                                      verbose=FALSE) {
  graph <- .G()
  ninds <- rep(FALSE, V(graph) |> length())
  
  if (attributes(kmo)$class[1]=="kegg_module") {
    return_nodes <- NULL
    nodes <- graph |> activate(nodes) |> data.frame()
    candidate <- graph |> activate(edges) |>
      filter(.data[[label]] & is.na(name)) |>
      mutate(node_to=nodes[to,"name"], ) |>
      data.frame()
    reacs <- NULL
    for (i in seq(1,length(candidate$node_to),2)) {
      reacs <- rbind(reacs, 
                     c(candidate[i,"node_to"],
                       candidate[i, "reaction"],
                       candidate[i+1, "node_to"],
                       candidate[i,"type"],
                       candidate[i,"to"],
                       candidate[i+1,"to"]))
    }
    for (i in seq_len(nrow(kmo@reaction_each))){
      left <- kmo@reaction_each[i,][1] |> unlist() |> as.character() |> paste0("cpd:", ...=_)
      right <- kmo@reaction_each[i,][3] |> unlist() |> as.character() |> paste0("cpd:", ...=_)
      for (j in seq_len(nrow(reacs))) {
        cand_left <- reacs[j,][1]
        cand_right <- reacs[j,][3]
        cand_id1 <- reacs[j,][5]
        cand_id2 <- reacs[j,][6]
        if (length(intersect(cand_left,left))>0 & 
            length(intersect(cand_right,right))>0) {
          return_nodes <- c(return_nodes, cand_id1, cand_id2)
        }
        if (length(intersect(cand_left,right))>0 & 
            length(intersect(cand_right,left))>0) {
          return_nodes <- c(return_nodes, cand_id1, cand_id2)
        }
      }
    }
    ninds[unique(as.numeric(return_nodes))] <- TRUE
    ninds
  } else {
    stop("please provide kegg_module class")
  }
}


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
