setOldClass("tbl_graph")
setClass("kegg_network",
    slots=list(
        ID="character",
        name="character",
        definition="character",
        expanded="character",
        expanded_graph="tbl_graph",
        definition_graph="tbl_graph",
        network_class="character",
        gene="character",
        metabolite="character"
    )
)

#' @importFrom GetoptLong qqcat
setMethod("show",
    signature(object="kegg_network"),
    function(object) {
        qqcat("@{object@ID}\n")
        qqcat("@{object@name}\n")
    }
)

setGeneric("get_network_attribute",
    function(x, ...) standardGeneric("get_network_attribute"))

setMethod("get_network_attribute", "kegg_network",
    function(x, attribute) attr(x, attribute))

#' KEGG network parsing function
#' 
#' parsing the network elements starting with N
#' 
#' @param nid KEGG NETWORK ID
#' @param use_cache use cache
#' @param directory directory to save raw files
#' @return list of network definition
#' @examples network("N00002")
#' @export
network <- function(nid, use_cache=FALSE, directory=NULL) {
    if (!startsWith(nid, "N")) {
        stop("Please provide a string that starts with N.")
    }
    kne <- new("kegg_network")
    kne@ID <- nid
    if (!is.null(directory)) {
        dest <- paste0(directory,"/",nid)
    } else {
        dest <- nid
    }
    if (!file.exists(dest)) {
        if (use_cache) {
            bfc <- BiocFileCache()
            dest <- bfcrpath(bfc,
            paste0("https://rest.kegg.jp/get/",nid))
        } else {
            download.file(paste0("https://rest.kegg.jp/get/",nid),
                    destfile=dest)
        }
    }
    con <- file(dest, "r")

    while ( TRUE ) {
        line <- readLines(con, n=1)
        if ( length(line) == 0 ) {
            break
        }
        if (grepl("NAME", line)) {
            name <- unlist(strsplit(line, "        "))[2]
            kne@name <- name
        }
        if (grepl("DEFINITION", line)) {
            definition <- unlist(strsplit(line, "  "))[2]
            kne@definition <- definition
        }
        if (grepl("EXPANDED", line)) {
            expanded <- unlist(strsplit(line, "  "))[3]
            kne@expanded <- expanded
        }
        if (grepl("CLASS", line)) {
            network_class <- unlist(strsplit(line, "  "))[2]
            kne@network_class <- network_class
        }
    }
    close(con)
    kne@expanded_graph <- convert_expanded_to_graph(kne)
    kne@definition_graph <- convert_definition_to_graph(kne)
    kne
}

#' @noRd
convert_expanded_to_graph <- function(kne) {
    sp <- kne@expanded |> strsplit(" ") |> unlist()
    edges <- lapply(seq(1,length(sp), 2), function(i) {
        if (i!=length(sp)) {
            left <- sp[i]
            edge <- sp[i+1]
            right <- sp[i+2]
            return(c(left, right, edge))
        } else {}       
    })
    edges <- do.call(rbind, edges) |> data.frame() |>
        `colnames<-`(c("from","to","type"))
    return(as_tbl_graph(edges))
}

#' @noRd
convert_definition_to_graph <- function(kne) {
    sp <- kne@definition |> strsplit(" ") |> unlist()
    edges <- lapply(seq(1,length(sp), 2), function(i) {
        if (i!=length(sp)) {
            left <- sp[i]
            edge <- sp[i+1]
            right <- sp[i+2]
            return(c(left, right, edge))
        } else {}       
    })
    edges <- do.call(rbind, edges) |> data.frame() |>
        `colnames<-`(c("from","to","type"))
    return(as_tbl_graph(edges))
}



#' network_graph
#' 
#' obtain tbl_graph of KEGG network
#' 
#' @param kne network object
#' @param type definition or expanded
#' @return tbl_graph
#' @examples
#' ne <- create_test_network()
#' neg <- network_graph(ne)
#' @export
#' 
network_graph <- function (kne, type="definition") {
    if (type=="definition") {
        raw_nodes <- kne@definition_graph |> activate("nodes") |> data.frame()
        raw_edges <- kne@definition_graph |> activate("edges") |> data.frame()
    } else {
        raw_nodes <- kne@expanded_graph |> activate("nodes") |> data.frame()
        raw_edges <- kne@expanded_graph |> activate("edges") |> data.frame()
    }

    res <- lapply(seq_along(raw_nodes$name), function(nn) {
        bln <- paste0("manual_BLOCK",nn,"_",kne@ID)
        ## In NETWORK definition, "-" is included in gene symbol
        ## Also, names like `Ca2+` is present, manually curate them
        input_string <- gsub("Ca2\\+","Ca",raw_nodes$name[nn])
        gra <- module_graph(input_string, skip_minus=TRUE)
        if (is.character(gra)) {
            # blocks <- rbind(blocks, c(gra, bln))
        } else {
            es <- as_data_frame(gra)
            es[,1] <- ifelse(startsWith(es[,1],"manual_CS"),
                paste0(es[,1],"_",nn,"_",kne@ID) ,es[,1])
            es[,2] <- ifelse(startsWith(es[,2],"manual_CS"),
                paste0(es[,2],"_",nn,"_",kne@ID) ,es[,2])
            es[,1] <- ifelse(startsWith(es[,1],"manual_G"),
                paste0(es[,1],"_",nn,"_",kne@ID) ,es[,1])
            es[,2] <- ifelse(startsWith(es[,2],"manual_G"),
                paste0(es[,2],"_",nn,"_",kne@ID) ,es[,2])

            vs <- data.frame(V(gra)$name, bln)
            vs[,1] <- ifelse(startsWith(vs[,1],"manual_CS"),
                paste0(vs[,1],"_",nn,"_",kne@ID) ,vs[,1])
            vs[,1] <- ifelse(startsWith(vs[,1],"manual_G"),
                paste0(vs[,1],"_",nn,"_",kne@ID) ,vs[,1])
            vs <- do.call(rbind, lapply(vs[,1], function(j) {
                c(j,bln,"in_block")
            })) |> data.frame() |> `colnames<-`(c("from","to","type"))
            list(rbind(es, vs), nn)
        }
    })
    edges <- do.call(rbind, lapply(res, function(x) x[[1]]))
    name_change <- lapply(res, function(x) x[[2]]) |> unlist()
    nns <- lapply(res, function(x) x[[2]]) |> unlist()

    name_change <- paste0("manual_BLOCK",name_change,"_",kne@ID)
    names(name_change) <- as.character(nns)

    new_edges_from <- NULL
    new_edges_to <- NULL

    new_edges_from <- lapply(raw_edges$from, function(i) {
        if (i %in% names(name_change)) {
            as.character(name_change[as.character(i)])
        } else {
            raw_nodes$name[i]
        }
    }) |> unlist()

    new_edges_to <- lapply(raw_edges$to, function(i) {
        if (i %in% names(name_change)) {
            as.character(name_change[as.character(i)])
        } else {
            raw_nodes$name[i]
        }
    }) |> unlist()

    raw_edges$from <- new_edges_from
    raw_edges$to <- new_edges_to
    raw_edges$subtype <- "reference"
    if (!is.null(edges)) {
        edges$subtype <- "manual"
    }
    all_edges <- rbind(raw_edges |>
        `colnames<-`(c("from","to","type","subtype")), edges)
    g <- as_tbl_graph(all_edges, directed=TRUE)
    g <- g |> activate("nodes") |>
        mutate(network_name=kne@name, network_ID=kne@ID)
    g
}

#' create_test_network
#' @return test network
#' @export
#' @examples create_test_network()
create_test_network <- function() {
    ne <- new("kegg_network")
    ne@ID <- "test"
    ne@name <- "test network"
    ne@definition <- "DDX41 -> IRF3"
    ne@definition_graph <- convert_definition_to_graph(ne)
    ne
}
