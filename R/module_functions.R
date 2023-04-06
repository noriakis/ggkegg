#' obtain_module
#' @noRd
obtain_module <- function(mid) {
  module <- list()
  if (!file.exists(mid)) {
    download.file(paste0("https://rest.kegg.jp/get/",mid),
                  destfile=mid)
  }
  con = file(mid, "r")
  reac <- FALSE
  reacs <- NULL
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (grepl("NAME", line)) {
      name <- unlist(strsplit(line, "        "))[2]
      module[["name"]] <- name
    }
    if (grepl("DEFINITION", line)) {
      definition <- unlist(strsplit(line, "  "))[2]
      module[["definition"]] <- definition
    }
  
    if (grepl("COMPOUND", line)) {reac <- FALSE}
    if (grepl("REACTION", line)) {reac <- TRUE}
    if (reac) {
        if (grepl("REACTION", line)) {
          reacs <- c(reacs, unlist(strsplit(line, "    "))[2])
        } else {
          reacs <- c(reacs, unlist(strsplit(line, "            "))[2])
        }
    }
  }
  module[["reaction"]] <- reacs
  close(con)
  module
}

#' parse_module
#' @noRd
parse_module <- function(mod, type="reaction") {
  if (type=="reaction") {
    reac <- NULL
    for (rea in mod$reaction) {
      left <- unlist(strsplit(rea, "->"))[1]
      right <- unlist(strsplit(rea, "->"))[2]
      if (grepl("\\+",right)) {
        right <- unlist(strsplit(right, "\\+"))
      }
      right <- gsub(" ","", right)
      
      # left
      left2 <- gsub(" ", "", unlist(strsplit(left, "  "))[2])
      left1 <- gsub(" ", "", unlist(strsplit(left, "  "))[1])
      for (l2 in unlist(strsplit(left2, ","))) {
        for (ll2 in unlist(strsplit(l2, "\\+"))) {
          for (l1 in unlist(strsplit(left1, ","))) {
            reac <- rbind(reac, c(ll2, l1))
          }
        }
      }
      for (l1 in unlist(strsplit(left1, ","))) {
        for (r in right) {
          reac <- rbind(reac, c(l1, r))
        }       
      }
    }
  } else if (type=="definition") {
    divide_string <- function(input_string) {
      steps <- c()
      current_step <- ""
      open_parentheses <- 0
      
      for (i in 1:nchar(input_string)) {
        char <- substr(input_string, i, i)
        
        if (char == "(") {
          open_parentheses <- open_parentheses + 1
        } else if (char == ")") {
          open_parentheses <- open_parentheses - 1
        } else if (char == " " & open_parentheses == 0) {
          steps <- append(steps, current_step)
          current_step <- ""
          next
        }
        
        current_step <- paste0(current_step, char)
      }
      
      steps <- append(steps, current_step)
      return(steps)
    }
    
    result <- divide_string(mod$definition)
    
    message("Currently returning list of the KO and steps")
    node_list <- list(step=result)
    return(node_list)
  }
}