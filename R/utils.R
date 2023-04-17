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


#' append_cp
#' 
#' append clusterProfiler results to graph
#' 
#' @param res enrichResult class
#' @param how how to determine whether the nodes is in enrichment results
#' @param name name column to search for query
#' @return enrich_attribute column in node
#' @export
append_cp <- function(res, how="any", name="name") {
  if (attributes(res)$class!="enrichResult") { stop("Please provide enrichResult class object") }

  graph <- .G()
  pid <- unique(V(graph)$pathway_id)
  x <- get.vertex.attribute(graph, name)

  org <- attributes(res)$organism
  res <- attributes(res)$result
  if (org!="UNKNOWN") {
    enrich_attribute <- paste0(org, ":", unlist(strsplit(
      res[pid,]$geneID, "/")))
  } else {
    enrich_attribute <- unlist(strsplit(
      res[pid,]$geneID, "/"))     
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
#' @export
assign_deseq2 <- function(res, column="log2FoldChange",
                          gene_type="SYMBOL",
                          org_db=org.Hs.eg.db, org="hsa",
                          numeric_combine=mean,
                          name="name") {
  graph <- .G()
  if (gene_type!="ENTREZID") {
    convert_df <- res |> row.names() |> AnnotationDbi::select(x=org.Hs.eg.db,
                                                keys=_,
                                                columns="ENTREZID",
                                                keytype=gene_type)
    
    nums <- data.frame(row.names(res), res[[column]]) |> `colnames<-`(c(gene_type, column))
    merged <- merge(nums, convert_df, by=gene_type)
  } else {
    merged <- data.frame(row.names(res), res[[column]]) |> `colnames<-`(c("ENTREZID", column))
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
#' 
#' @return tbl_graph
#' @export
#' 
combine_with_bnlearn <- function(pg, str, av, prefix="ko:") {
  
  ## Make igraph with strength from bnlearn
  el <- av |> bnlearn::as.igraph() |> as_edgelist() |> data.frame() |>
    `colnames<-`(c("from","to"))
  g <- str |> merge(el) |> mutate(from=paste0(prefix,from),
                                               to=paste0(prefix,to)) |>
                              data.frame() |> graph_from_data_frame()
  
  ## Merge node names with reference
  js <- NULL
  for (i in V(pg)$name) {
    if (grepl(" ",i)) {
      ref_node <- strsplit(i, " ") |> unlist()
      for (j in V(g)$name) {
        if (length(intersect(ref_node, j))>0) {
          js <- rbind(js, c(j, i))
        }
      }
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
    }
  }
  gdf <- rbind(gdf, new_df |> `colnames<-`(colnames(gdf)))
  
  for (i in seq_len(nrow(gdf))) {
    if (gdf[i,"to"] %in% js$raw){
      new_to <- js[js[,1]==gdf[i,"to"],]$reference
      new_df <- rbind(new_df,
                      c(gdf[i,"from"], new_to, 
                        gdf[i,"strength"], gdf[i,"direction"])) |> data.frame()
    }
  }
  gdf <- rbind(gdf, new_df |> `colnames<-`(colnames(gdf)))
  
  gdf$strength <- as.numeric(gdf$strength)
  gdf$direction <- as.numeric(gdf$direction)
  
  joined <- graph_join(pg, gdf)
  joined
}