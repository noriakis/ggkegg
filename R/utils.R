if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' find_parenthesis_pairs
#' find pairs of parenthesis
#' @noRd
find_parenthesis_pairs <- function(s) {
    ## Preallocate
    stack <- integer(nchar(s))
    pairs <- vector(mode="list", length=nchar(s)/2)
    j <- 1
    for (i in seq_len(nchar(s))) {
        c <- substr(s, i, i)
        if (c == "(") {
            stack[i] <- 1
        } else if (c == ")") {
            if (length(which(stack==1)) == 0) {
                stop("Mismatched parenthesis")
        }
        open <- tail(which(stack==1), 1)
        stack[open] <- 0
        pairs[[j]] <- c(open, i)
        j <- j+1
        }
    }
    if (length(which(stack==1)) > 0) {
        stop("Mismatched parenthesis")
    }
    pairs[vapply(pairs, is.null, TRUE)] <- NULL
    pairs
}


#' append_label_position
#'
#' Append the label position at center of edges
#' in global map like ko01100 where line type nodes
#' are present in KGML. 
#' Add `center` column to graph edge.
#' 
#' @param g graph
#' @importFrom dplyr mutate summarise group_by filter
#' @importFrom dplyr row_number n distinct ungroup
#' @importFrom stats setNames
#' @return tbl_graph
#' @examples 
#' ## Simulate nodes containing `graphics_type` of line and `coords`
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
#' gm_test <- create_test_pathway(line=TRUE)
#' test <- process_line(gm_test) |> return_line_compounds(1)
return_line_compounds <- function(g, orig) {
    ndf <- g |> activate("nodes") |> data.frame()
    edf <- g |> activate("edges") |> data.frame()
    highl <- ndf[edf[edf$to %in% as.integer(ndf[ndf$orig.id %in% orig,] |> 
        row.names()),]$from,]$orig.id
    highl2 <- ndf[edf[edf$from %in% as.integer(ndf[ndf$orig.id %in% orig,] |> 
        row.names()),]$to,]$orig.id
    c(highl, highl2)
}

#' edge_numeric
#' 
#' add numeric attribute to edge of tbl_graph
#' 
#' @param num named vector or tibble with id and value column
#' @param num_combine how to combine number when multiple hit in the same node
#' @param name name of column to match for
#' @param sep separater for name, default to " "
#' @param remove_dot remove "..." in the name
#' @param how `any` or `all`
#' @export
#' @return numeric vector
#' @importFrom tibble is_tibble
#' @importFrom tidygraph activate
#' @examples
#' graph <- create_test_pathway()
#' graph <- graph |> activate("edges") |>
#'             mutate(num=edge_numeric(c(1.1) |>
#'             setNames("degradation"), name="subtype_name"))
edge_numeric <- function(num, num_combine=mean, how="any", name="name",
    sep=" ", remove_dot=TRUE) {
    graph <- .G()
    if (!is_tibble(num) & !is.vector(num)) {
        stop("Please provide tibble or named vector")
    }
    if (is_tibble(num)) {
        if (duplicated(num$id) |> unique() |> length() > 1) {
            stop("Duplicate ID found")
        }
        changer <- num$value
        names(changer) <- num$id
    } else {
        if (duplicated(names(num)) |> unique() |> length() > 1) {
            stop("Duplicate ID found")
        }
        changer <- num
    }

    x <- get.edge.attribute(graph, name)

    lapply(x, function(xx) {
        in_node <- strsplit(xx, sep) |> unlist() |> unique()
        if (remove_dot) {
            in_node <- lapply(in_node, function(nn) {
                strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
            }) %>% unlist()
        }
        thresh <- ifelse(how == "any", 1, length(in_node))
        if (length(intersect(names(changer), in_node)) >= thresh) {
            summed <- do.call(num_combine,
                list(x=changer[intersect(names(changer), in_node)]))
        } else {
            summed <- NA
        }
        return(summed)
    }) |> unlist()
}


#' edge_numeric_sum
#' 
#' add numeric attribute to edge of tbl_graph based on node values
#' The implementation is based on the paper by
#' Adnan et al. 2020 (https://doi.org/10.1186/s12859-020-03692-2).
#' 
#' @param num named vector or tibble with id and value column
#' @param num_combine how to combine number when multiple hit in the same node
#' @param name name of column to match for
#' @param sep separater for name, default to " "
#' @param remove_dot remove "..." in the name
#' @param how `any` or `all`
#' @export
#' @return numeric vector
#' @importFrom tibble is_tibble
#' @importFrom tidygraph activate
#' @examples
#' graph <- create_test_pathway()
#' graph <- graph |>
#'             activate("edges") |>
#'             mutate(num=edge_numeric_sum(c(1.2,-1.2) |>
#'             setNames(c("TRIM21","DDX41")), name="graphics_name"))
edge_numeric_sum <- function(num, num_combine=mean, how="any", name="name",
    sep=" ", remove_dot=TRUE) {
    graph <- .G()
  
    if (!is_tibble(num) & !is.vector(num)) {
        stop("Please provide tibble or named vector")
    }
    if (is_tibble(num)) {
        if (duplicated(num$id) |> unique() |> length() > 1) {
            stop("Duplicate ID found")
        }
        changer <- num$value
        names(changer) <- num$id
    } else {
        if (duplicated(names(num)) |> unique() |> length() > 1) {
            stop("Duplicate ID found")
        }
        changer <- num
    }
  
    node_df <- graph |> activate("nodes") |> data.frame()
    node_name <- node_df[[name]]
    new_graph <- graph |> activate(edges) |>
        mutate(from_nd=node_name[.data$from], to_nd=node_name[.data$to]) |>
        mutate(summed=edge_numeric(num, num_combine, how, name="from_nd", sep=sep, remove_dot=remove_dot)+
             edge_numeric(num, num_combine, how, name="to_nd", sep=sep, remove_dot=remove_dot)) |>
        data.frame()
    new_graph$summed
}


#' node_numeric
#' 
#' simply add numeric attribute to node of tbl_graph
#' 
#' @param num named vector or tibble with id and value column
#' @param num_combine how to combine number when multiple hit in the same node
#' @param how how to match the node IDs with the queries 'any' or 'all'
#' @param name name of column to match for
#' @param sep separater for name, default to " "
#' @param remove_dot remove "..." in the name
#' @export
#' @return numeric vector
#' @importFrom tibble is_tibble
#' @examples
#' graph <- create_test_pathway()
#' graph <- graph |>
#'             mutate(num=node_numeric(c(1.1) |> setNames("hsa:6737"))) 
#' 
node_numeric <- function(num, num_combine=mean,
    name="name", how="any", sep=" ", remove_dot=TRUE) {
    graph <- .G()
    if (!is_tibble(num) & !is.vector(num)) {
        stop("Please provide tibble or named vector")
    }
    if (is_tibble(num)) {
        if (duplicated(num$id) |> unique() |> length() > 1) {
            stop("Duplicate ID found")
        }
        changer <- num$value
        names(changer) <- num$id
    } else {
        if (duplicated(names(num)) |> unique() |> length() > 1) {
            stop("Duplicate ID found")
        }
        changer <- num
    }
    x <- get.vertex.attribute(graph, name)

    lapply(x, function(xx) {
        in_node <- strsplit(xx, sep) |> unlist() |> unique()
        if (remove_dot) {
            in_node <- lapply(in_node, function(nn) {
                strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
            }) %>% unlist()
        }
        thresh <- ifelse(how=="any", 1, length(in_node))
        if (length(intersect(names(changer), in_node)) >= thresh) {
            summed <- do.call(num_combine,
                        list(x=changer[intersect(names(changer), in_node)]))
        } else {
            summed <- NA
        }        
    }) |> unlist()
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
#' @param name name column in node data, default to node
#' @param sep separater of name, default to " "
#' @param remove_dot remove "..." in the name
#' @export
#' @return tbl_graph
#' @examples
#' 
#' ## Append data.frame to tbl_graph
#' graph <- create_test_pathway()
#' num_df <- data.frame(row.names=c("6737","51428"),
#'                     "sample1"=c(1.1,1.2),
#'                     "sample2"=c(1.5,2.2),
#'                     check.names=FALSE)
#' graph <- graph |> node_matrix(num_df, gene_type="ENTREZID")
#' 
node_matrix <- function(graph, mat, gene_type="SYMBOL", org="hsa",
                        org_db=NULL, num_combine=mean, name="name",
                        sep=" ", remove_dot=TRUE) {
    get_value <- function(x) {
        val <- lapply(seq_along(x), function(xx) {
            if (x[xx]=="undefined") {return(NA)}
            vals <- strsplit(x[xx], sep) |> unlist() |> unique()
            if (remove_dot) {
                vals <- lapply(vals, function(nn) {
                    strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
                }) %>% unlist()
            }
            subset_conv <- convert_df |>
                filter(.data$converted %in% vals) |>
                data.frame()
            if (dim(subset_conv)[1]==0) {return(NA)}
            if (dim(subset_conv)[1]==1) {
                return(mat[subset_conv[[gene_type]],])
            }
            return(apply(mat[ subset_conv[[gene_type]],], 2, num_combine))
        })
        binded <- do.call(rbind, val)
        binded
    }

    node_df <- graph |> activate("nodes") |> data.frame()
    node_name <- node_df[[name]]
    if (gene_type!="ENTREZID") {
        if (!requireNamespace("AnnotationDbi")) {
            stop("This conversion requires AnnotationDbi.")
        }
        if (is.null(org_db)) {
            stop("Please specify Annotation DB to org_db.")
        }
        convert_df <- mat %>% 
                row.names() %>%
                AnnotationDbi::mapIds(x=org_db, keys=.,
                    column="ENTREZID", keytype=gene_type) %>% 
                tibble::enframe() %>% 
                `colnames<-`(c(gene_type, "ENTREZID"))
    } else {
        convert_df <- data.frame(row.names(mat)) %>% `colnames<-`(c("ENTREZID"))
    }
  
    convert_df$converted <- paste0(org, ":", convert_df[["ENTREZID"]])
    new_edges <- graph |> activate("edges") |> data.frame()
    summed <- data.frame(get_value(node_df[[name]]))
    new_nodes <- cbind(node_df, summed)
    appended <- tbl_graph(nodes=new_nodes, edges=new_edges)
    appended
}

#' edge_matrix
#' 
#' given the matrix representing gene as row and sample as column,
#' append the edge value (sum of values of connecting nodes) to edge matrix and
#' return tbl_graph object. The implementation is based on the paper by
#' Adnan et al. 2020 (https://doi.org/10.1186/s12859-020-03692-2).
#' 
#' @param graph tbl_graph to append values to
#' @param mat matrix representing gene as row and sample as column
#' @param gene_type gene ID of matrix row
#' @param org organism ID to convert ID
#' @param org_db organism database to convert ID
#' @param num_combine function to combine multiple numeric values
#' @param name name column in node data, default to node
#' @param sep separater of name, default to " "
#' @param remove_dot remove "..." in node name
#' @export
#' @return tbl_graph
#' @examples
#' graph <- create_test_pathway()
#' num_df <- data.frame(row.names=c("6737","51428"),
#'                     "sample1"=c(1.1,1.2),
#'                     "sample2"=c(1.1,1.2),
#'                     check.names=FALSE)
#' graph <- graph %>% edge_matrix(num_df, gene_type="ENTREZID")
edge_matrix <- function(graph, mat, gene_type="SYMBOL", org="hsa",
    org_db=NULL, num_combine=mean, name="name", sep=" ", remove_dot=TRUE) {
    get_value <- function(x) {
        val <- lapply(seq_along(x), function(xx) {
            if (x[xx]=="undefined") {return(NA)}
            vals <- strsplit(x[xx], " ") %>% unlist() %>% unique()
            if (remove_dot) {
                vals <- lapply(vals, function(nn) {
                    strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
                }) %>% unlist()
            }
            subset_conv <- convert_df %>% filter(.data$converted %in% vals) %>% 
                            data.frame()
            if (dim(subset_conv)[1]==0) {
                return(NA)
            }
            if (dim(subset_conv)[1]==1) {
                return(mat[subset_conv[[gene_type]],])
            }
            return(apply(mat[ subset_conv[[gene_type]],], 2, num_combine))
        })
        binded <- do.call(rbind, val)
        binded
    }
  
    node_df <- graph %>% activate("nodes") %>% data.frame()
    node_name <- node_df$name
    if (gene_type!="ENTREZID") {
        if (!requireNamespace("AnnotationDbi")) {
            stop("This conversion requires AnnotationDbi.")
        }
        if (is.null(org_db)) {
            stop("Please specify Annotation DB to org_db.")
        }
        convert_df <- mat %>% 
                row.names() %>%
                AnnotationDbi::mapIds(x=org_db, keys=.,
                    column="ENTREZID", keytype=gene_type) %>% 
                tibble::enframe() %>% 
                `colnames<-`(c(gene_type, "ENTREZID"))
    } else {
        convert_df <- data.frame(row.names(mat)) %>% `colnames<-`(c("ENTREZID"))
    }
  
    convert_df$converted <- paste0(org, ":", convert_df[["ENTREZID"]])
    new_graph <- graph %>% activate(edges) %>%
        mutate(from_nd=node_name[.data$from], to_nd=node_name[.data$to]) %>% 
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
#' @param sep separater for name
#' @param remove_dot remove dots in the name
#' @param pid pathway ID, if NULL, try to infer from graph attribute
#' @param infer if TRUE, append the prefix to queried IDs based on pathway ID
#' @return enrich_attribute column in node
#' @examples
#' graph <- create_test_pathway()
#' nodes <- graph |> data.frame()
#' if (require("clusterProfiler")) {
#'     cp <- enrichKEGG(nodes$name |>
#'             strsplit(":") |>
#'             vapply("[", 2, FUN.VALUE="character"))
#'     ## This append graph node logical value whether the 
#'     ## enriched genes are in pathway
#'     graph <- graph |> mutate(cp=append_cp(cp, pid="hsa05322"))
#' }
#' @export
#' 
append_cp <- function(res, how="any", name="name", pid=NULL, infer=FALSE, sep=" ", remove_dot=TRUE) {
    if (!attributes(res)$class %in% c("enrichResult","gseaResult")) {
        stop("Please provide enrichResult or gseaResult class object")
    }
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

    if (name=="graphics_name") {
        ## If graphics name, use as is.
        enrich_attribute <- unlist(strsplit(res[pid,][[gene_col]], "/"))   
    } else {
        if (org!="UNKNOWN") {
            if (org=="microbiome") {org <- "ko"; pid <- gsub("ko","map",pid)}
            enrich_attribute <- paste0(org, ":", unlist(strsplit(
                                                    res[pid,][[gene_col]], "/")))
        } else {## If UNKNOWN
            ## Try to infer
            if (infer) {
                org <- gsub("[^a-zA-Z]", "", pid)
                enrich_attribute <- paste0(org, ":", unlist(strsplit(res[pid,][[gene_col]], "/")))
            } else {
                enrich_attribute <- unlist(strsplit(res[pid,][[gene_col]], "/"))
            }
        }        
    }
    bools <- vapply(x, function(xx) {
        in_node <- strsplit(xx, sep) |> unlist() |> unique()
        if (remove_dot) {
            in_node <- lapply(in_node, function(nn) {
                strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
            }) %>% unlist()
        }
        if (how=="any") {
            if (length(intersect(in_node, enrich_attribute))>=1) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        } else {
            if (length(intersect(in_node, enrich_attribute))==length(in_node)) {
                return(TRUE)
            } else {
                return(FALSE)
            }      
        } 
    }, FUN.VALUE=TRUE)
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
#' @param sep for node name
#' @param remove_dot remove dot in the name
#' @return numeric vector
#' @export
#' @examples
#' graph <- create_test_pathway()
#' res <- data.frame(row.names="6737",log2FoldChange=1.2)
#' graph <- graph |> mutate(num=assign_deseq2(res, gene_type="ENTREZID"))
assign_deseq2 <- function(res, column="log2FoldChange",
                          gene_type="SYMBOL",
                          org_db=NULL, org="hsa",
                          numeric_combine=mean,
                          name="name", sep=" ", remove_dot=TRUE) {
    graph <- .G()
    if (gene_type!="ENTREZID") {
        if (!requireNamespace("AnnotationDbi")) {
            stop("This conversion requires AnnotationDbi.")
        }
        if (is.null(org_db)) {
            stop("Please specify Annotation DB to org_db.")
        }
        convert_df <- res %>%
            row.names() %>%
            AnnotationDbi::mapIds(x=org_db, keys=.,
                column="ENTREZID", keytype=gene_type) %>%
            tibble::enframe() %>% 
            `colnames<-`(c(gene_type, "ENTREZID"))
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
    lapply(x, function(xx) {
        in_node <- strsplit(xx, sep) |> unlist() |> unique()
        if (remove_dot) {
            in_node <- lapply(in_node, function(nn) {
                strsplit(nn, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
            }) %>% unlist()
        }
        do.call(numeric_combine,
            list(x=changer[intersect(in_node, names(changer))]))       
    }) |> unlist()
}



#' convert_id
#' 
#' convert the identifier using retrieved information
#' 
#' @param org which identifier to convert
#' @param name which column to convert in edge or node table
#' @param convert_column which column is parsed in 
#' obtained data frame from KEGG REST API or local file
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
#' @param remove_dot remove dots in the name
#' @param file specify the file for conversion.
#' The column in `query_column` will be used for querying the ID in the graph.
#' @param query_column default to 1.
#' @param pref prefix for the query identifiers
#' @importFrom data.table fread
#' @return vector containing converted IDs
#' @export
#' @examples
#' library(tidygraph)
#' graph <- create_test_pathway()
#' graph <- graph %>% mutate(conv=convert_id("hsa"))
#' 
convert_id <- function(org=NULL, name="name", file=NULL, query_column=1,
    convert_column=NULL, colon=TRUE, first_arg_comma=TRUE, remove_dot=TRUE,
    pref=NULL, sep=" ", first_arg_sep=TRUE, divide_semicolon=TRUE, edge=FALSE) {
    if (is.null(org) & is.null(file)) {
        stop("Please specify org or file")
    }
    graph <- .G()
    pid <- unique(V(graph)$pathway_id)
    if (edge) {
        x <- get.edge.attribute(graph, name)
    } else {
        x <- get.vertex.attribute(graph, name)
    }
    if (is.null(file)) {
        url <- paste0("https://rest.kegg.jp/list/",org)
        bfc <- BiocFileCache()
        path <- bfcrpath(bfc, url)
        convert <- fread(path,
                    header=FALSE,
                    sep="\t") %>% data.frame()        
    } else {
        convert <- fread(file,
                    header=FALSE,
                    sep="\t") %>% data.frame()
        if (is.null(convert_column)) {
            stop("Please specify the column number for the file")
        }
        if (is.null(pref)) {
            pref <- ""
        }
    }


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
            convert_column <- 4
        }
    }
    convert_vec <- convert[,convert_column]
    if (is.null(org)) {
        names(convert_vec) <- 
            paste0(pref,convert[, query_column])        
    } else {
        if (org=="pathway") {
            names(convert_vec) <- 
                paste0(pref,str_extract(convert[, query_column], "[[:digit:]]+"))
        } else {
            names(convert_vec) <- 
                paste0(pref,convert[, query_column])
        }            
    }
    if (!colon) {
        names(convert_vec) <- unlist(
            lapply(strsplit(names(convert_vec), ":"), "[", 2)
        )
    }
    convs <- lapply(seq_along(x), function(xn) {
        if (grepl(sep,x[xn])) {
            spaced <- lapply(unlist(strsplit(x[xn], sep)), function (qu) {
                if (remove_dot) {
                    qu <- strsplit(qu, "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
                }
                comma_test <- ifelse(first_arg_comma,
                    strsplit(convert_vec[qu], ",")[[1]][1],
                    paste0(convert_vec[qu]))
                sc_test <- ifelse(divide_semicolon,
                    strsplit(comma_test, ";") |>
                    vapply("[",1,FUN.VALUE="character"),
                    comma_test)
                return(sc_test)
            }) |> unlist()
            spaced <- ifelse(first_arg_sep, spaced[1],
                paste(spaced, collapse=sep))
            return(spaced)
        } else {
            if (remove_dot) {
                qu <- strsplit(x[xn], "\\.\\.\\.") %>% vapply("[", 1, FUN.VALUE="a")
            } else {
                qu <- x[xn]
            }
            comma_test <- ifelse(first_arg_comma,
                strsplit(convert_vec[qu], ",")[[1]][1],
                convert_vec[qu])
            sc_test <- ifelse(divide_semicolon,
                strsplit(comma_test, ";") |>
                vapply("[",1,FUN.VALUE="character"),
                comma_test)
            return(sc_test)
        }
    })  
    convs |> unlist()
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
                            header=FALSE,
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
        names(convert_vec) <- paste0(pref,convert$V1)
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
        names(convert_vec) <- unlist(
            lapply(strsplit(names(convert_vec), ":"), "[", 2)
        )
  }
  convert_vec
}

#' carrow
#' 
#' make closed type arrow
#' @param length arrow length in unit()
#' @export
#' @examples
#' carrow()
#' @return arrow()
#' 
carrow <- function(length=unit(2,"mm")) {
	arrow(length=length, type="closed")
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
#' if (requireNamespace("bnlearn", quietly=TRUE)) {
#'     ## Simulating boot.strength() results
#'     av <- bnlearn::model2network("[6737|51428][51428]")
#'     str <- data.frame(from="51428",to="6737",strength=0.8,direction=0.7)
#'     graph <- create_test_pathway()
#'     combined <- combine_with_bnlearn(graph, str, av, prefix="hsa:")
#' }
#' 
combine_with_bnlearn <- function(pg, str, av, prefix="ko:", how="any") {
    if (requireNamespace("bnlearn", quietly=TRUE)) {
        ## Make igraph with strength from bnlearn
        el <- av |> bnlearn::as.igraph() |> as_edgelist() |> data.frame() |>
            `colnames<-`(c("from","to"))
        g <- str |> merge(el) |> mutate(from=paste0(prefix,.data$from),
                                        to=paste0(prefix,.data$to)) |>
                                    data.frame() |> graph_from_data_frame()

        ## Merge node names with reference
        js <- lapply(V(pg)$name, function(i) {
            if (grepl(" ",i)) {
                ref_node <- strsplit(i, " ") |> unlist()
                ret <- lapply(V(g)$name, function(j) {
                    if (how=="any") {
                        if (length(intersect(ref_node, j))>0) {
                            return(c(j, i))
                        }
                    } else {
                        if (length(intersect(ref_node, j))==length(ref_node)) {
                            return(c(j, i))
                        }          
                    }
                })
                return(do.call(rbind, ret))
            } else {
                return(c(i, i))
            }
        })

        js <- do.call(rbind, js) |>
            data.frame() |>
            `colnames<-`(c("raw","reference"))
        gdf <- as_data_frame(g)

        new_df <- lapply(seq_len(nrow(gdf)), function(i) {
            if (gdf[i,"from"] %in% js$raw){
                new_from <- js[js[,1]==gdf[i,"from"],]$reference
                return(c(new_from, gdf[i,"to"], 
                    gdf[i,"strength"], gdf[i,"direction"]))
            } else {
                stop("no `from` included in raw node name")
            }
        })

        gdf <- do.call(rbind, new_df) |>
            data.frame() |>
            `colnames<-`(colnames(gdf))

        new_df <- lapply(seq_len(nrow(gdf)), function(i) {
            if (gdf[i,"to"] %in% js$raw){
                new_to <- js[js[,1]==gdf[i,"to"],]$reference
                new_df <- return(c(gdf[i,"from"], new_to, 
                            gdf[i,"strength"], gdf[i,"direction"]))
            } else {
                stop("no `to` included in raw node name")    
            }            
        })
        gdf <- do.call(rbind, new_df) |> 
            data.frame() |>
            `colnames<-`(colnames(gdf))
    
        gdf$strength <- as.numeric(gdf$strength)
        gdf$direction <- as.numeric(gdf$direction)
    
        ## Drop duplicates
        gdf <- gdf |>
            distinct(.data$from, .data$to, .data$strength, .data$direction)

        joined <- graph_join(pg, gdf, by="name")
        joined    
    }
}

#' add_readable_edge_label
#' 
#' Add readable edge label by shadowtext.
#' Should be just under the geom_edge_* with the filtered labels.
#' 
#' @param color color mapping name in the previous layer
#' @param angle_calc only `along` or `none`
#' @param size size of the label
#' @export
#' @return ggplot2 object
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  geom_edge_link(aes(label=subtype_name, label_colour=subtype_name))+
#'  add_readable_edge_label()
add_readable_edge_label <- function(color="label_colour", angle_calc="along", size=3) {
    structure(list(color=color, angle_calc=angle_calc, size=size),
              class="add_readable_edge_label")
}

#' ggplot_add.geom_node_rect_kegg
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param ... The other arguments
#' @export ggplot_add.add_readable_edge_label
#' @export
#' @return ggplot2 object
#' @examples
#' test_pathway <- create_test_pathway()
#' plt <- ggraph(test_pathway, layout="manual", x=x, y=y) +
#'  geom_edge_link(aes(label=subtype_name, label_colour=subtype_name))+
#'  add_readable_edge_label()
ggplot_add.add_readable_edge_label <- function(object, plot, ...) {
    built <- ggplot_build(plot)
    datlen <- length(built$data)
    candidate <- built$data[[datlen]]
    if (!("label" %in% colnames(candidate))) {stop("There is no label specified in the previous layer.")}
    
    if ("xend" %in% colnames(candidate)) {
        if (object$angle_calc=="along") {
            readable <- candidate %>%
                mutate(angle2=edge_angle(x, y, xend, yend)) %>%
                group_by(group, label, .data[[object$color]]) %>%
                summarise(mx=mean(x), my=mean(y), ma=mean(angle2)) %>%
                ungroup()
        } else {
            readable <- candidate %>%
                group_by(group, label, .data[[object$color]]) %>%
                summarise(mx=mean(x), my=mean(y), ma=mean(angle)) %>%
                ungroup()  
        }
        
        plot$layers[[datlen]]$mapping$label <- NULL

    } else {## If no xend
        if (object$angle_calc=="along") {
            readable <- candidate %>%
                group_by(group, label, .data[[object$color]]) %>%
                summarise(mx=mean(x), my=mean(y),
                          x1=dplyr::first(x),
                          x2=dplyr::last(x),
                          y1=dplyr::first(y),
                          y2=dplyr::last(y)) %>%
                mutate(ma=edge_angle(x1, y1, x2, y2)) %>%
                ungroup()        
        } else {
            readable <- candidate %>%
                group_by(group, label, .data[[object$color]]) %>%
                summarise(mx=mean(x), my=mean(y), ma=mean(angle)) %>%
                ungroup()  
        }
        plot$layers[[datlen]]$mapping$label <- NULL
    }
    return(plot + shadowtext::geom_shadowtext(data=readable,
                  aes(x=mx, y=my, label=label, angle=ma, color=.data[[object$color]]),
                  size=object$size,
                  bg.colour="white"))
}