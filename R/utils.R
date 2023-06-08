#' find_parenthesis_pairs
#' find pairs of parenthesis
#' @noRd
find_parenthesis_pairs <- function(s) {
  stack <- list()
  pairs <- list()
  for (i in seq_len(nchar(s))) {
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


#' append_label_position
#'
#' append the label position at center of edges
#' in global map like ko01100 where line type nodes
#' are present in KGML. Add `center` column to graph edge
#' 
#' @param g graph
#' @importFrom dplyr mutate summarise group_by filter
#' @importFrom dplyr row_number n distinct ungroup
#' @importFrom stats setNames
#' @return tbl_graph
#' @examples 
#' ## For those containing nodes with the graphic type of `line`
#' gm_test <- data.frame(name="ko:K00112",type="ortholog",reaction="rn:R00112",
#'            graphics_name="K00112",fgcolor="#ff0000",bgcolor="#ffffff",
#'            graphics_type="line",coords="1,2,3,4",orig.id=1,pathway_id="test")
#' gm_test <- tbl_graph(gm_test)
#' test <- process_line(gm_test) |> append_label_position()
#' @export
append_label_position <- function(g) {
  pos <- g |>
    activate(edges) |> 
    data.frame() |>
    filter(.data$type=="line") |>
    group_by(.data$orig.id) |> 
    summarise(n=n()) |> 
    mutate(n2=n/2) |> 
    mutate(n3=as.integer(.data$n2+1))
  posvec <- pos$n3 |> setNames(pos$orig.id)
  g |> activate(edges) |> group_by(.data$orig.id) |> 
    mutate(rn=row_number()) |> ungroup() |>
    mutate(showpos=edge_numeric(name="orig.id", posvec)) |>
    mutate(center=.data$rn==.data$showpos) |>
    mutate(rn=NULL, showpos=NULL)
}

#' return_line_compounds
#' 
#' In the map, where lines are converted to edges,
#' identify compounds that are linked by the reaction.
#' Give the original edge ID of KGML (orig.id in edge table), and 
#' return the original compound node ID
#' 
#' @param g tbl_graph object
#' @param orig original edge ID
#' @return vector of original compound node IDs
#' @export
#' @examples
#' ## For those containing nodes with the graphic type of `line`
#' ## This returns no IDs as no edges are present
#' gm_test <- data.frame(name="ko:K00112",type="ortholog",reaction="rn:R00112",
#'            graphics_name="K00112",fgcolor="#ff0000",bgcolor="#ffffff",
#'            graphics_type="line",coords="1,2,3,4",orig.id=1,pathway_id="test")
#' gm_test <- tbl_graph(gm_test)
#' test <- process_line(gm_test) |> return_line_compounds(1)
return_line_compounds <- function(g, orig) {
  ndf <- g |> activate(nodes) |> data.frame()
  edf <- g |> activate(edges) |> data.frame()
  highl <- ndf[edf[edf$to %in% as.integer(ndf[ndf$orig.id %in% orig,] |> 
    row.names()),]$from,]$orig.id
  highl2 <- ndf[edf[edf$from %in% as.integer(ndf[ndf$orig.id %in% orig,] |> 
    row.names()),]$to,]$orig.id
  c(highl, highl2)
}

#' edge_numeric
#' 
#' add numeric attribute to edge of tbl_graph
#' for matrix input, use `append_node_value`
#' 
#' @param num named vector or tibble with id and value column
#' @param num_combine how to combine number when multiple hit in the same node
#' @param name name of column to match for
#' @param how `any` or `all`
#' @export
#' @return numeric vector
#' @importFrom tibble is_tibble
#' @importFrom tidygraph activate
#' @examples
#' nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
#'                     x=c(1,1),
#'                     xmin=c(-1,-1),
#'                     xmax=c(2,2),
#'                     y=c(1,1),
#'                     ymin=c(-1,-1),
#'                     ymax=c(2,2))
#' edges <- data.frame(from=1, to=2, name="K00112")
#' graph <- tbl_graph(nodes, edges)
#' graph <- graph |> activate("edges") |>
#'    mutate(num=edge_numeric(c(1.1) |> setNames("K00112")))
#' 
#' 
edge_numeric <- function(num, num_combine=mean, how="any", name="name") {
  graph <- .G()
  if (!is_tibble(num) & !is.vector(num)) {
    stop("Please provide tibble or named vector")}
  if (is_tibble(num)) {
    if (duplicated(num$id) |> unique() |> length() > 1) {
      stop("Duplicate ID found")}
    changer <- num$value
    names(changer) <- num$id
  } else {
    if (duplicated(names(num)) |> unique() |> length() > 1) {
      stop("Duplicate ID found")}
    changer <- num
  }
  x <- get.edge.attribute(graph, name)

  final_attribute <- NULL
  for (xx in x) {
    in_node <- strsplit(xx, " ") |> unlist() |> unique()
    thresh <- ifelse(how=="any", 1, length(in_node))
    if (length(intersect(names(changer), in_node)) >= thresh) {
      summed <- do.call(num_combine,
        list(x=changer[intersect(names(changer), in_node)]))
    } else {
      summed <- NA
    }
    final_attribute <- c(final_attribute,
                         summed)
  }
  final_attribute
}

#' node_numeric
#' 
#' simply add numeric attribute to node of tbl_graph
#' for matrix input, use `append_node_value`
#' 
#' @param num named vector or tibble with id and value column
#' @param num_combine how to combine number when multiple hit in the same node
#' @param name name of column to match for
#' @param how how to match the node IDs with the queries 'any' or 'all'
#' @export
#' @return numeric vector
#' @importFrom tibble is_tibble
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
#' graph <- graph |> 
#'            mutate(num=node_numeric(c(1.1) |> setNames("hsa:1029")))
#' 
#' 
node_numeric <- function(num, num_combine=mean, name="name", how="any") {
  graph <- .G()
  if (!is_tibble(num) & !is.vector(num)) {
    stop("Please provide tibble or named vector")}
  if (is_tibble(num)) {
    if (duplicated(num$id) |> unique() |> length() > 1) {
      stop("Duplicate ID found")}
    changer <- num$value
    names(changer) <- num$id
  } else {
    if (duplicated(names(num)) |> unique() |> length() > 1) {
      stop("Duplicate ID found")}
    changer <- num
  }
  x <- get.vertex.attribute(graph, name)

  final_attribute <- NULL
  for (xx in x) {
    in_node <- strsplit(xx, " ") |> unlist() |> unique()
    thresh <- ifelse(how=="any", 1, length(in_node))
    if (length(intersect(names(changer), in_node)) >= thresh) {
      summed <- do.call(num_combine,
        list(x=changer[intersect(names(changer), in_node)]))
    } else {
      summed <- NA
    }
    final_attribute <- c(final_attribute,
                         summed)
  }
  final_attribute
}


#' node_matrix
#' 
#' given the matrix representing gene as row and sample as column,
#' append the node value to node matrix and
#' return tbl_graph object
#' 
#' @param graph tbl_graph to append values to
#' @param mat matrix representing gene as row and sample as column
#' @param gene_type gene ID of matrix row
#' @param org organism ID to convert ID
#' @param org_db organism database to convert ID
#' @param num_combine function to combine multiple numeric values
#' @export
#' @importFrom AnnotationDbi select
#' @return tbl_graph
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
#' num_df <- data.frame(row.names=c("1029","4171"),
#'                      "sample1"=c(1.1,1.2),
#'                      "sample2"=c(1.1,1.2),
#'                      check.names=FALSE)
#' graph <- graph |> node_matrix(num_df, gene_type = "ENTREZID")
node_matrix <- function(graph, mat, gene_type="SYMBOL", org="hsa",
                        org_db=org.Hs.eg.db, num_combine=mean) {
  
  get_value <- function(x) {
    val <- list()
    for (xx in seq_along(x)) {
      if (x[xx]=="undefined") {val[[xx]] <- NA; next}
      vals <- strsplit(x[xx], " ") |> unlist() |> unique()
      subset_conv <- convert_df |> filter(converted %in% vals) |> data.frame()
      if (dim(subset_conv)[1]==0) {val[[xx]]<- NA; next}
      if (dim(subset_conv)[1]==1) {
        val[[xx]]<-mat[subset_conv[[gene_type]],]; next}
      val[[xx]] <- apply(mat[ subset_conv[[gene_type]],], 2, num_combine)
    }
    binded <- do.call(rbind, val)
    binded
  }

  node_df <- graph |> activate(nodes) |> data.frame()
  node_name <- node_df$name
  if (gene_type!="ENTREZID") {
    convert_df <- mat |> row.names() |> select(x=org_db,
                                          keys=_,
                                          columns="ENTREZID",
                                          keytype=gene_type)
  } else {
    convert_df <- data.frame(row.names(mat)) |> `colnames<-`(c("ENTREZID"))
  }
  
  convert_df$converted <- paste0(org, ":", convert_df[["ENTREZID"]])
  new_edges <- graph |> activate(edges) |> data.frame()
  summed <- data.frame(get_value(node_df$name))
  new_nodes <- cbind(node_df, summed)
  appended <- tbl_graph(nodes=new_nodes, edges=new_edges)
  appended
}

#' edge_matrix
#' 
#' given the matrix representing gene as row and sample as column,
#' append the edge value (sum of values of connecting nodes) to edge matrix and
#' return tbl_graph object
#' 
#' @param graph tbl_graph to append values to
#' @param mat matrix representing gene as row and sample as column
#' @param gene_type gene ID of matrix row
#' @param org organism ID to convert ID
#' @param org_db organism database to convert ID
#' @param num_combine function to combine multiple numeric values
#' @export
#' @importFrom AnnotationDbi select
#' @return tbl_graph
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
#' num_df <- data.frame(row.names=c("1029","4171"),
#'                      "sample1"=c(1.1,1.2),
#'                      "sample2"=c(1.1,1.2),
#'                      check.names=FALSE)
#' graph <- graph |> edge_matrix(num_df, gene_type = "ENTREZID")
edge_matrix <- function(graph, mat, gene_type="SYMBOL", org="hsa",
                              org_db=org.Hs.eg.db,
                              num_combine=mean) {
  
  get_value <- function(x) {
    val <- list()
    for (xx in seq_along(x)) {
      if (x[xx]=="undefined") {val[[xx]] <- NA; next}
      vals <- strsplit(x[xx], " ") |> unlist() |> unique()
      subset_conv <- convert_df |> filter(.data$converted %in% vals) |> 
                      data.frame()
      if (dim(subset_conv)[1]==0) {val[[xx]]<- NA; next}
      if (dim(subset_conv)[1]==1) {
        val[[xx]]<-mat[subset_conv[[gene_type]],]; next}
      val[[xx]] <- apply(mat[ subset_conv[[gene_type]],], 2, num_combine)
    }
    binded <- do.call(rbind, val)
    binded
  }
  
  node_df <- graph |> activate("nodes") |> data.frame()
  node_name <- node_df$name
  if (gene_type!="ENTREZID") {
    convert_df <- mat |> row.names() |> select(x=org_db,
                                                          keys=_,
                                                          columns="ENTREZID",
                                                          keytype=gene_type)
} else {
    convert_df <- data.frame(row.names(mat)) |> `colnames<-`(c("ENTREZID"))
  }
  
  convert_df$converted <- paste0(org, ":", convert_df[["ENTREZID"]])
  new_graph <- graph |> activate(edges) |>
    mutate(from_nd=node_name[.data$from], to_nd=node_name[.data$to]) |> 
    data.frame()
  summed <- data.frame(
    get_value(new_graph$from_nd) + get_value(new_graph$to_nd))
  new_edges <- cbind(new_graph, summed)
  appended <- tbl_graph(nodes=node_df, edges=new_edges)
  appended
}

#' append_cp
#' 
#' append clusterProfiler results to graph
#' 
#' @param res enrichResult class
#' @param how how to determine whether the nodes is in enrichment results
#' @param name name column to search for query
#' @param pid pathway ID, if NULL, try to infer from graph attribute
#' @return enrich_attribute column in node
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
#' if (require("clusterProfiler")) {
#'   cp <- enrichKEGG(nodes$name |>
#'                   strsplit(":") |> 
#'                   sapply("[",2))
#' }
#' graph <- graph |> mutate(cp=append_cp(cp,pid="hsa04110"))
#' @export
append_cp <- function(res, how="any", name="name", pid=NULL) {
  if (!attributes(res)$class %in% c("enrichResult","gseaResult")) {
    stop("Please provide enrichResult or gseaResult class object") }
  if (attributes(res)$class=="gseaResult") {
    gene_col <- "core_enrichment"
  } else {
    gene_col <- "geneID"
  }
  graph <- .G()
  if (is.null(pid)) {
    pid <- unique(V(graph)$pathway_id)
  }
  x <- get.vertex.attribute(graph, name)

  org <- attributes(res)$organism
  res <- attributes(res)$result

  if (org!="UNKNOWN") {
    if (org=="microbiome") {org <- "ko"; pid <- gsub("ko","map",pid)}
    enrich_attribute <- paste0(org, ":", unlist(strsplit(
      res[pid,][[gene_col]], "/")))
  } else {
    enrich_attribute <- unlist(strsplit(
      res[pid,][[gene_col]], "/"))     
  }
  bools <- NULL
  for (xx in x) {
    in_node <- strsplit(xx, " ") |> unlist() |> unique()
    if (how=="any") {
      if (length(intersect(in_node, enrich_attribute))>=1) {
        bools <- c(bools, TRUE)
      } else {
        bools <- c(bools, FALSE)
      }
    } else {
      if (length(intersect(in_node, enrich_attribute))==length(in_node)) {
        bools <- c(bools, TRUE)
      } else {
        bools <- c(bools, FALSE)
      }      
    }
  }
  bools
}



#' assign_deseq2
#' 
#' assign DESeq2 numerical values to nodes
#' 
#' @param res The result() of DESeq()
#' @param column column of the numeric attribute, default to log2FoldChange
#' @param gene_type default to SYMBOL
#' @param org_db organism database to convert ID to ENTREZID
#' @param org organism ID in KEGG
#' @param numeric_combine how to combine multiple numeric values
#' @param name column name for ID in tbl_graph nodes
#' @return numeric vector
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' @export
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
#' res <- data.frame(row.names="1029",log2FoldChange=1.2)
#' graph <- graph |> mutate(num=assign_deseq2(res, gene_type = "ENTREZID"))
assign_deseq2 <- function(res, column="log2FoldChange",
                          gene_type="SYMBOL",
                          org_db=org.Hs.eg.db, org="hsa",
                          numeric_combine=mean,
                          name="name") {
  graph <- .G()
  if (gene_type!="ENTREZID") {
    convert_df <- res |> row.names() |> select(x=org_db,
                                                keys=_,
                                                columns="ENTREZID",
                                                keytype=gene_type)
    
    nums <- data.frame(row.names(res), res[[column]]) |> 
    `colnames<-`(c(gene_type, column))
    merged <- merge(nums, convert_df, by=gene_type)
  } else {
    merged <- data.frame(row.names(res), res[[column]]) |> 
    `colnames<-`(c("ENTREZID", column))
  }
  merged$converted <- paste0(org, ":", merged[["ENTREZID"]])
  changer <- merged[[column]] |> `names<-`(merged[["converted"]])
  x <- get.vertex.attribute(graph, name)
  final_attribute <- NULL
  for (xx in x) {
    in_node <- strsplit(xx, " ") |> unlist() |> unique()
    final_attribute <- c(final_attribute,
                 do.call(numeric_combine,
                         list(x=changer[intersect(in_node, names(changer))])))
  }
  final_attribute
}



#' convert_id
#' 
#' convert the identifier using retrieved information
#' 
#' @param org which identifier to convert
#' @param name which column to convert in edge or node table
#' @param convert_column which column is parsed in 
#' obtained data frame from KEGG REST API
#' @param colon whether the original ids include colon (e.g. `ko:`)
#' If `NULL`, automatically set according to `org`
#' @param first_arg_comma take first argument of comma-separated
#' string, otherwise fetch all strings
#' @param first_arg_sep take first argument if multiple identifiers
#' are in the node name, otherwise parse all identifiers
#' @param sep separater to separate node names, defaul to space
#' @param divide_semicolon whether to divide string by semicolon,
#' and take the first value
#' @param edge if converting edges
#' @importFrom data.table fread
#' @return vector containing converted IDs
#' @export
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
#' \donttest{graph <- graph |> mutate(conv=convert_id("hsa"))}
#' 
convert_id <- function(org, name="name",
  convert_column=NULL, colon=TRUE, first_arg_comma=TRUE,
  sep=" ", first_arg_sep=TRUE, divide_semicolon=TRUE, edge=FALSE) {
  graph <- .G()
  pid <- unique(V(graph)$pathway_id)
  if (edge) {
    x <- get.edge.attribute(graph, name)
  } else {
    x <- get.vertex.attribute(graph, name)
  }
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
    names(convert_vec) <- unlist(lapply(strsplit(names(convert_vec), ":"),
      "[", 2))
  }
  convs <- NULL
  for (xn in seq_along(x)) {
    if (grepl(sep,x[xn])) {
      spaced <- NULL
      for (qu in unlist(strsplit(x[xn], sep))) {

        comma_test <- ifelse(first_arg_comma,
            strsplit(convert_vec[qu], ",")[[1]][1],
            paste0(convert_vec[qu]))

        sc_test <- ifelse(divide_semicolon,
          strsplit(comma_test, ";") |> sapply("[",1),
          comma_test)

        spaced <- c(spaced, sc_test)
      }
      spaced <- ifelse(first_arg_sep, spaced[1],
        paste(spaced, collapse=sep))
      convs <- c(convs, spaced)
    } else {
      comma_test <- ifelse(first_arg_comma,
        strsplit(convert_vec[x[xn]], ",")[[1]][1],
        convert_vec[x[xn]])
      sc_test <- ifelse(divide_semicolon,
          strsplit(comma_test, ";") |> sapply("[",1),
          comma_test)
      convs <- c(convs, sc_test)
    }
  }
  convs
  # if (divide_semicolon) {
  #   new_convs <- NULL
  #   for (i in convs) {
  #     app <- NULL
  #     for (j in strsplit(i, " ") |> unlist()) {
  #       app <- c(app, sapply(strsplit(j, ";"),"[", 1))
  #     }
  #     new_convs <- c(new_convs,
  #       paste0(app, collapse=sep))
  #   }
  #   return(new_convs)
  # } else {
  #   return(convs)
  # }
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
    names(convert_vec) <- unlist(lapply(strsplit(names(convert_vec), ":"),
      "[", 2))
  }
  convert_vec
}



#' combine_with_bnlearn
#' 
#' combine the reference KEGG pathway graph 
#' with bnlearn boot.strength output
#' 
#' @param pg reference graph (output of `pathway`)
#' @param str strength data.frame
#' @param av averaged network to plot
#' @param prefix add prefix to node name of original averaged network
#' like, `hsa:` or `ko:`.
#' @param how `any` or `all`
#' 
#' @return tbl_graph
#' @importFrom tidygraph graph_join
#' @export
#' @examples
#' if (requireNamespace("bnlearn", quietly = TRUE)) {
#'   av <- bnlearn::model2network("[1029|4171][4171]")
#'   str <- data.frame(from="4171",to="1029",strength=0.8,direction=0.7)
#'   nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
#'                       x=c(1,1),
#'                       xmin=c(-1,-1),
#'                       xmax=c(2,2),
#'                       y=c(1,1),
#'                       ymin=c(-1,-1),
#'                       ymax=c(2,2))
#'   edges <- data.frame(from=1, to=2, name="K00112")
#'   graph <- tbl_graph(nodes, edges)
#'   combined <- combine_with_bnlearn(graph, str, av, prefix="hsa:")
#' }
#' 
combine_with_bnlearn <- function(pg, str, av, prefix="ko:", how="any") {
  if (requireNamespace("bnlearn", quietly = TRUE)) {
    ## Make igraph with strength from bnlearn
    el <- av |> bnlearn::as.igraph() |> as_edgelist() |> data.frame() |>
      `colnames<-`(c("from","to"))
    g <- str |> merge(el) |> mutate(from=paste0(prefix,.data$from),
                                                 to=paste0(prefix,.data$to)) |>
                                data.frame() |> graph_from_data_frame()

    ## Merge node names with reference
    js <- NULL
    for (i in V(pg)$name) {
      if (grepl(" ",i)) {
        ref_node <- strsplit(i, " ") |> unlist()
        for (j in V(g)$name) {
          if (how=="any") {
            if (length(intersect(ref_node, j))>0) {
              js <- rbind(js, c(j, i))
            }
          } else {
            if (length(intersect(ref_node, j))==length(ref_node)) {
              js <- rbind(js, c(j, i))
            }          
          }
        }
      } else {
        js <- rbind(js, c(i, i))
      }
    }

    js <- js |> data.frame() |> `colnames<-`(c("raw","reference"))
    gdf <- as_data_frame(g)

    new_df <- NULL
    for (i in seq_len(nrow(gdf))) {
      if (gdf[i,"from"] %in% js$raw){
        new_from <- js[js[,1]==gdf[i,"from"],]$reference
        new_df <- rbind(new_df,
                        c(new_from, gdf[i,"to"], 
                          gdf[i,"strength"], gdf[i,"direction"])) |> data.frame()
      } else {
        stop("no `from` included in raw node name")
      }
    }

    gdf <- new_df |> data.frame() |> `colnames<-`(colnames(gdf))

    new_df <- NULL
    for (i in seq_len(nrow(gdf))) {
      if (gdf[i,"to"] %in% js$raw){
        new_to <- js[js[,1]==gdf[i,"to"],]$reference
        new_df <- rbind(new_df,
                        c(gdf[i,"from"], new_to, 
                          gdf[i,"strength"], gdf[i,"direction"])) |> data.frame()
      } else {
        stop("no `to` included in raw node name")    
      }
    }
    gdf <- new_df |> `colnames<-`(colnames(gdf))
    
    gdf$strength <- as.numeric(gdf$strength)
    gdf$direction <- as.numeric(gdf$direction)
    
    ## Drop duplicates
    gdf <- gdf |> distinct(.data$from, .data$to, .data$strength, .data$direction)

    joined <- graph_join(pg, gdf, by="name")
    joined    
  }
}