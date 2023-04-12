#' find_parenthesis_pairs
#' find pairs of parenthesis
#' @noRd
find_parenthesis_pairs <- function(s) {
  stack <- list()
  pairs <- list()
  for (i in 1:nchar(s)) {
    c <- substr(s, i, i)
    if (c == "(") {
      stack <- c(stack, i)
    } else if (c == ")") {
      if (length(stack) == 0) {
        stop("Mismatched parenthesis")
      }
      open <- tail(stack, 1)
      stack <- head(stack, -1)
      pairs <- c(pairs, list(c(open, i)))
    }
  }
  if (length(stack) > 0) {
    stop("Mismatched parenthesis")
  }
  pairs
}


#' highlight_set_nodes
#' 
#' identify if nodes are involved in specific queries
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

#' highlight_module
#' 
#' identify if edges are involved in module.
#' Using reaction of module.
#' 
#' @param set set of identifiers
#' @param name which column to search for
#' @param sep separater for node names
#' @export
highlight_module <- function(graph, kmo, name="name", sep=" ") {
  einds <- rep(FALSE, E(graph) |> length())
  ninds <- rep(FALSE, V(graph) |> length())
  if (attributes(kmo)$class[1]=="kegg_module") {
    rea <- kmo@reaction_each
    for (i in seq_len(nrow(kmo@reaction_each))) {
      left <- rea[i,][1] |> unlist() |> as.character() |> paste0("cpd:", ...=_)
      reac <- rea[i,][2] |> unlist() |> as.character() |> paste0("rn:", ...=_)
      right <- rea[i,][3] |> unlist() |> as.character()|> paste0("cpd:", ...=_)

      x <- get.edge.attribute(graph, "reaction")
      ind <- NULL
      for (xn in seq_along(x)) {
        edge_reac <- x[xn] |> strsplit(" ") |> unlist()
        if (sum(is.na(edge_reac))!=length(edge_reac)) {
          if (length(intersect(edge_reac,reac))==length(reac)) {
            ind <- c(ind, TRUE)
          } else {
            ind <- c(ind, FALSE)
          }
        } else {
          ind <- c(ind, FALSE)
        } ## if edge is reaction
      } ## each edge
      einds <- einds | ind


      ## Identify compounds using edge index
      cand_edges <- graph |> 
                      activate(edges) |> 
                      filter(ind) |>
                      data.frame()
      nodes <- graph |> 
                 activate(nodes) |> 
                 data.frame()
      node_ind <- cand_edges[is.na(cand_edges$name),]$to |> unique()
      cand_nodes <- nodes[ node_ind , ]
      true_nodes <- NULL
      for (node in node_ind) {
        node_comp <- nodes[node,name] |> strsplit(" ") |> unlist()
        if (sum(is.na(node_comp))!=length(node_comp)) {
          if (length(intersect(node_comp,left))==length(node_comp) |
            length(intersect(node_comp,right))==length(node_comp)) {
            ## If fullfil right or left
            true_nodes <- c(true_nodes, node)
          }
        }
      }
      ninds[true_nodes] <- TRUE
    } # each reaction

    #   x <- get.vertex.attribute(graph, name)
    #   ind <- NULL
    #   for (xn in seq_along(x)) {
    #     node_comp <- x[xn] |> strsplit(" ") |> unlist()
    #     if (sum(is.na(node_comp))!=length(node_comp)) {
    #       if (length(intersect(node_comp,left))==length(left) |
    #         length(intersect(node_comp,right))==length(right)) {
    #         ind <- c(ind, TRUE)
    #       } else {
    #         ind <- c(ind, FALSE)
    #       }
    #     } else {
    #       ind <- c(ind, FALSE)
    #     } ## if node is compound
    #   }## each node
    #   ninds <- ninds | ind
    # }# each reaction

    graph |>
      activate(nodes) |>
      mutate(!!kmo@ID := ninds) |>
      activate(edges) |>
      mutate(!!kmo@ID := einds)
  } else {
    stop("please provide kegg_module class")
  }
}


#' highlight_set_edges
#' 
#' identify if edges are involved in specific queries.
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


#' convert_id
#' 
#' convert the identifier using retrieved information
#' 
#' @param org which identifier to convert
#' @param name which column to convert
#' @param convert_column which column is parsed in 
#' obtained data frame from KEGG REST API
#' @param colon whether the original ids include colon
#' @param first_arg_comma take first argument of comma-separated
#' string, otherwise fetch all strings
#' @param first_arg_sep take first argument if multiple identifiers
#' are in the node name, otherwise parse all identifiers
#' @param sep separater to separate node names, defaul to space
#' @param divide_semicolon whether to divide string by semicolon,
#' and take the first value
#' @importFrom data.table fread
#' @export
#' 
convert_id <- function(org, name="name",
  convert_column=NULL, colon=TRUE, first_arg_comma=TRUE,
  sep=" ", first_arg_sep=TRUE, divide_semicolon=TRUE) {
  graph <- .G()
  pid <- unique(V(graph)$pathway_id)
  x <- get.vertex.attribute(graph, name)
  url <- paste0("https://rest.kegg.jp/list/",org)
  bfc <- BiocFileCache()
  path <- bfcrpath(bfc, url)
  convert <- fread(path,
                   header = FALSE,
                   sep="\t") |>
            data.frame()

  if (is.null(convert_column)) {
    if (org=="ko") {pref <- "ko:";convert_column <- 2}
    else if (org=="compound") {pref <- "cpd:"; convert_column <- 2}  
    else if (org=="glycan") {pref <- "gl:";convert_column <- 2}
    else if (org=="enzyme") {pref <- "ec:"; convert_column <- 2}
    else if (org=="reaction") {pref <- "rn:"; convert_column <- 2}
    else if (org=="pathway") {
      pref <- paste0("path:",gsub("[[:digit:]]","",pid));
      convert_column <- 2
      if (is.null(pid)) {stop("please specify pathway id")}

    }
    else {
      pref <- ""
      convert_column <- 4}
  }
  convert_vec <- convert[,convert_column]

  if (org=="pathway") {
    names(convert_vec) <- 
      paste0(pref,str_extract(convert$V1, "[[:digit:]]+"))
  } else {
    names(convert_vec) <- 
      paste0(pref,convert$V1)
  }
  if (!colon) {
    names(convert_vec) <- unlist(lapply(strsplit(names(convert_vec), ":"), "[", 2))
  }
  convs <- NULL
  for (xn in seq_along(x)) {
    if (grepl(sep,x[xn])) {
      spaced <- NULL
      for (qu in unlist(strsplit(x[xn], sep))) {
        spaced <- c(spaced, ifelse(first_arg_comma,
          strsplit(convert_vec[qu], ",")[[1]][1],
          paste0(convert_vec[qu])))
      }
      spaced <- ifelse(first_arg_sep, spaced[1],
        paste(spaced, collapse=" "))
      convs <- c(convs, spaced)
    } else {
      convs <- c(convs, ifelse(first_arg_comma,
        strsplit(convert_vec[x[xn]], ",")[[1]][1],
        convert_vec[x[xn]]))
    }
  }
  if (divide_semicolon) {
    convs <- unlist(lapply(strsplit(convs, ";"),"[", 1))
  }
  convs
}



#' obtain_map_and_cache
#' 
#' obtain list of genes, cache, and return the named vector for converting
#' 
#' @import BiocFileCache
#' @importFrom stringr str_extract str_extract_all str_pad str_locate_all
#' @noRd
obtain_map_and_cache <- function(org, pid=NULL, colon=TRUE) {
  url <- paste0("https://rest.kegg.jp/list/",org)
  bfc <- BiocFileCache()
  path <- bfcrpath(bfc, url)
  convert <- data.table::fread(path,
                               header = FALSE,
                               sep="\t")
  if (org %in% c("ko","compound")) {## KO and compound
    if (org=="compound") {pref <- "cpd"} 
    else {pref <- "ko"}
    convert_vec <- vapply(convert$V2, function(x) {
      vapply(unlist(strsplit(x, ";"))[1],
             function(x) unlist(strsplit(x,","))[1],
             FUN.VALUE="character")
      
    }, FUN.VALUE="character")
    names(convert_vec) <- paste0(pref,":",convert$V1)
  } else if (org=="reaction") {## Reaction
    pref <- "rn:"
    convert_vec <- convert$V2
    names(convert_vec) <- 
      paste0(pref,convert$V1)
  } else if (org=="pathway") {## Pathway
    pref <- paste0("path:",gsub("[[:digit:]]","",pid))
    convert_vec <- convert$V2
    names(convert_vec) <- 
      paste0(pref,str_extract(convert$V1, "[[:digit:]]+"))
  } else {## Ordinary organisms
    convert_vec <- vapply(convert$V4, function(x) {
      vapply(unlist(strsplit(x, ";"))[1],
             function(x) unlist(strsplit(x,","))[1],
             FUN.VALUE="character")
      
    }, FUN.VALUE="character")
    names(convert_vec) <- convert$V1
  }
  if (!colon) {
    names(convert_vec) <- unlist(lapply(strsplit(names(convert_vec), ":"), "[", 2))
  }
  convert_vec
}
