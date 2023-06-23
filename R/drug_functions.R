#' drug
#' KEGG drug parsing function
#' @param did KEGG drug ID
#' @return list of components in drug
#' @noRd
drug <- function(did) {
  if (!startsWith(did, "D")) {
    stop("Please provide a string that starts with N.")
  }
  if (!file.exists(did)) {
    download.file(paste0("https://rest.kegg.jp/get/",did),
                  destfile=did)
  }
  current_sub_id <- NULL
  con <- file(did, "r")
  content_list <- list()
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (!startsWith(line, " ")) {
      current_id <- strsplit(line, " ") |> vapply("[", 1, FUN.VALUE="character")
      current_sub_id <- NULL
    } else {
      if (substr(line, 1, 12)!="            ") {
        current_sub_id <- gsub(" ","", substr(line, 1, 12))
      }
    }
    ## Treating subclass as class here
    if (!current_id %in% c("REFERENCE","///")) {
      content <- substr(line, 13, nchar(line))
      if (!is.null(current_sub_id)) {
        content_list[[current_sub_id]] <- 
          c(content_list[[current_sub_id]], content) 
      } else {
        content_list[[current_id]] <- c(content_list[[current_id]], content) 
      }
    }
  }
  close(con)
  content_list$ENTRY <- substr(content_list$ENTRY, 1, 6)
  content_list
}
