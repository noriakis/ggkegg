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
        ))

#' @importFrom GetoptLong qqcat
setMethod("show",
  signature(object="kegg_network"),
  function(object) {
    qqcat("@{object@ID}\n")
    qqcat("@{object@name}\n")
  })


#' KEGG network parsing function
#' network
#' parsing the network elements starting with N
#' @return list of network definition
#' @export
network <- function(nid) {
  if (!startsWith(nid, "N")) {stop("Please provide a string that starts with N.")}
  kne <- new("kegg_network")
  kne@ID <- nid
  if (!file.exists(nid)) {
    download.file(paste0("https://rest.kegg.jp/get/",nid),
                  destfile=nid)
  }
  con = file(nid, "r")

  while ( TRUE ) {
    line = readLines(con, n = 1)
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
  edges <- NULL
  for (i in seq(1,length(sp), 2)) {
    if (i!=length(sp)) {
      left <- sp[i]
      edge <- sp[i+1]
      right <- sp[i+2]
      edges <- rbind(edges, c(left,right,edge))
    } else {}
  }
  edges <- edges |> data.frame() |>
    `colnames<-`(c("from","to","type"))
  return(as_tbl_graph(edges))
}

#' @noRd
convert_definition_to_graph <- function(kne) {
  sp <- kne@definition |> strsplit(" ") |> unlist()
  edges <- NULL
  for (i in seq(1,length(sp), 2)) {
    if (i!=length(sp)) {
      left <- sp[i]
      edge <- sp[i+1]
      right <- sp[i+2]
      edges <- rbind(edges, c(left,right,edge))
    } else {}
  }
  edges <- edges |> data.frame() |>
    `colnames<-`(c("from","to","type"))
  return(as_tbl_graph(edges))
}