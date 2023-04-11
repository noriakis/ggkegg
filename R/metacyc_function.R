#'
#' parseMetaCycPathwayReactions
#' 
#' parse MetaCyc "pathways.dat"
#' and return the nested list of:
#' 
#' UNIQUE-ID
#' REACTION-LAYOUT
#' COMMON-NAME
#' SPECIES
#' TAXONOMIC-RANGE
#' 
#' @param file path to pathways.dat
#' @examples
#' file <- "pathways.dat"
#' \dontrun{parseMetaCycPathwayReactions(file, candSp="all")}
#' @export
#' 
parseMetaCycPathwayReactions <- function(file) {
  flg <- FALSE
  con = file(file, "r")
  allList <- list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (startsWith(line,"UNIQUE-ID - ")) {
      pwy <- gsub("UNIQUE-ID - ","",line)
      flg <- TRUE
    }
    if (flg) {
      if (startsWith(line, "REACTION-LAYOUT")) {
        reaclayout <- gsub("REACTION-LAYOUT - ","",line)
        allList[[pwy]][["REACTION-LAYOUT"]] <- c(allList[[pwy]][["REACTION-LAYOUT"]],
                                                 reaclayout)
      }
      if (startsWith(line, "COMMON-NAME")) {
        commn <- gsub("COMMON-NAME - ","",line)
        allList[[pwy]][["COMMON-NAME"]] <- commn
      }
      if (startsWith(line, "SPECIES - ")) {
        spec <- gsub("SPECIES - ","",line)
        allList[[pwy]][["SPECIES"]] <- c(allList[[pwy]][["SPECIES"]],spec)
        
      }
      if (startsWith(line, "TAXONOMIC-RANGE - ")) {
        taxr <- gsub("TAXONOMIC-RANGE - ","",line)
        allList[[pwy]][["TAXONOMIC-RANGE"]] <- taxr
      }
      if (startsWith(line,"//")) {
        flg <- FALSE
      }
    }
  }
  close(con)
  return(allList)
}

#' parse_reaction_layout
#' 
#' parse reaction layout of metacyc pathway
#' and return tbl_graph object
#' 
#' @export
parse_reaction_layout <- function(reacs) {
  eds <- NULL
  for (r in reacs) {
    f <- unlist(strsplit(r, " \\("))
    edge_label <- gsub("\\(","",f[1])
    left <- gsub("\\)","",strsplit(f[2]," ")[[1]][2])
    direction <- f[3]
    right <- gsub("\\)","",strsplit(f[4]," ")[[1]][2])
    
    if (grepl("L2R",direction)) {
      eds <- rbind(eds, c(left, right, edge_label))
    } else if (grepl("R2L",direction)) {
      eds <- rbind(eds, c(right, left, edge_label))
    } else {
      stop("cannot understand direction")
    }
  }
  eds <- eds |> data.frame() |> `colnames<-`(c("from","to","reaction"))
  as_tbl_graph(eds)
}