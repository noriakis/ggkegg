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
