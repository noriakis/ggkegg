#' module_text
#' Obtain textual representation of modules for all the steps
#' @export
module_text <- function(def, candidate_ko=NULL, paint_colour="tomato", convert=NULL) {
  plot_list <- list()
  for (step in seq_along(def$step)) {
    input_string <- def$step[step]
    ppos <- NULL
    for (i in find_parenthesis_pairs(input_string)) {
      ppos <- rbind(ppos, c(i[[1]], i[[2]], i[[2]]-i[[1]]))
    }
    if (!is.null(ppos)) {
      posmat <- ppos[order(ppos[,3]),]
      if (is.vector(posmat)) {dfs <- data.frame(t(posmat))} else {
        dfs <- data.frame(posmat)
      }
      posmat <- dfs |> `colnames<-`(c("xmin","xmax","length"))
      posmat$name <- paste0("G",str_pad(1:nrow(posmat),5,pad="0"))
      ul <- sort(unique(posmat$length))
      he <- (1:length(ul))+1
      names(he) <- ul
      posmat$height <- he[as.character(posmat$length)]/2
      posmat$text <- apply(posmat, 1, function(row) substr(input_string, row["xmin"], row["xmax"]))
      posmat$rawtext <- apply(posmat, 1, function(row) substr(input_string, row["xmin"], row["xmax"]))
      posmat$size <- posmat$length
      posmat$x <- (as.numeric(posmat$xmin)+as.numeric(posmat$xmax))/2
      # posmat$length <- NULL
    }
    
    ## All-KO
    kopos <- NULL
    for (i in unlist(def$ko_in_step[step])) {
      findko <- str_locate_all(input_string, i)
      for (ff in findko) {
        for (rn in seq_len(nrow(ff))) {
          kopos <- rbind(kopos,
                         c(i, ff[rn, 1], ff[rn, 2]))
        }
      }
    }
    kopos <- data.frame(kopos) |> `colnames<-`(c("name","xmin","xmax"))
    kopos$x <- (as.numeric(kopos$xmin)+as.numeric(kopos$xmax))/2
    kopos$height <- 0.5
    kopos$size <- 6
    
    locate_and_append <- function(concat, loc) {
      if (grepl(loc, input_string)) {
        put <- gsub("\\\\", "", loc)
        if (put==" ") {put <- "->"}
        if (put==",") {put <- "alt"}
        findx <- str_locate_all(input_string, loc)
        for (pos in findx) {
          for (rn in seq_len(nrow(pos))) {
            concat <- rbind(concat, c(put, pos[rn, 1], pos[rn, 1], pos[rn, 1], 0.5, 1))
          }
        }
      }
      concat
    }
  
  
  
    if (!is.null(ppos)) {
      concat <- rbind(posmat[,c("name","xmin","xmax","x","height","size")], kopos)
    } else {
      concat <- kopos
    }
    
    concat <- locate_and_append(concat, "\\+")
    concat <- locate_and_append(concat, "\\-")
    concat <- locate_and_append(concat, " ")
    concat <- locate_and_append(concat, ",")
    
    concat$xmin <- as.numeric(concat$xmin)
    concat$xmax <- as.numeric(concat$xmax)
    concat$x <- as.numeric(concat$x)
    concat$height <- as.numeric(concat$height)
    
    concat$ymin <- 1 - concat$height
    concat$ymax <- 1 + concat$height
    
    bgcol <- NULL
    for (nm in concat$name) {
      if (nm %in% candidate_ko) {
        bgcol <- c(bgcol, paint_colour)
      } else {
        bgcol <- c(bgcol, "transparent")
      }
    }
    concat$color <- bgcol

    if (!is.null(convert)) {
      conv <- convert[concat$name]
      conv[is.na(conv)] <- concat$name[is.na(conv)]
      concat$converted_name <- conv
    }
    concat$koflag <- startsWith(concat$name, "K")
    concat$conflag <- concat$name %in% c("->","+","-","alt")
    plot_list[[step]] <- concat
  }
  plot_list
}

#' obtain_sequential_module_definition
#' 
#' Given module definition and step number,
#' Recursively obtain graphical represencation of step and 
#' connect them by pseudo-nodes representing steps.
#' @export
obtain_sequential_module_definition <- function(def, step=NULL) {
  
  if (is.null(step)) {cand_step <- def$step}
  all_steps <- NULL
  orders <- NULL

  for (i in seq_along(cand_step)) {
    if (def$num_in_step[i]!=1) {
      for (ko in def$ko_in_step[[i]]) {
        all_steps <- rbind(all_steps, c(ko, paste0("STEP",i),"instep"))
        all_steps <- all_steps |> data.frame() |> `colnames<-`(c("from","to","type"))
      }
      plotg <- as_data_frame(get_module_graph(def$step[i]))
      ## Need to change naming
      frm <- plotg$from
      frm[!startsWith(frm,"K")] <- paste0(frm[!startsWith(frm,"K")],"_",i)
      to <- plotg$to
      to[!startsWith(to,"K")] <- paste0(to[!startsWith(to,"K")],"_",i)
      plotg$from <- frm
      plotg$to <- to

      all_steps <- rbind(all_steps, plotg)
      orders <- c(orders, paste0("STEP",i))
    } else {
      ## Some steps have only "-K*****"
      orders <- c(orders, def$step[i])
    }
  } ## For each step, obtain the plot
  
  ords <- NULL
  for (i in seq_along(orders)) {
    if (i!=length(orders)) {
      ords <- rbind(ords, c(orders[i], orders[i+1],"step_transition"))
    }
  }
  if (!is.null(ords)) {
    all_steps <- rbind(all_steps, 
                     data.frame(ords)|>
                       `colnames<-`(c("from","to","type")))
  } else {
    ## Only one step
    all_steps <- subset(all_steps, all_steps$type!="instep")

  }
  return(list(all_steps=all_steps, definition=def))
}


#' get_module_graph
#' obtain graphical representation of module definition
#' @return igraph object
#' @noRd
get_module_graph <- function(input_string) {
  ppos <- NULL
  for (i in find_parenthesis_pairs(input_string)) {
    ppos <- rbind(ppos, c(i[[1]], i[[2]], i[[2]]-i[[1]]))
  }
  if (!is.null(ppos)) {
    posmat <- ppos[order(ppos[,3]),]
    if (is.vector(posmat)) {dfs <- data.frame(t(posmat))} else {
      dfs <- data.frame(posmat)
    }
    posmat <- dfs |> `colnames<-`(c("xmin","xmax","length"))
    posmat$name <- paste0("G",str_pad(1:nrow(posmat),5,pad="0"))
    ul <- sort(unique(posmat$length))
    he <- (1:length(ul))+1
    names(he) <- ul
    posmat$height <- he[as.character(posmat$length)]/2
    posmat$text <- apply(posmat, 1, function(row) substr(input_string, row["xmin"], row["xmax"]))
    posmat$rawtext <- apply(posmat, 1, function(row) substr(input_string, row["xmin"], row["xmax"]))

    converted_string <- input_string
    for (i in seq_along(posmat$text)) {
      if (i < nrow(posmat)) {
        converted_string <- gsub(posmat$text[i], posmat$name[i], converted_string, fixed = TRUE)
        posmat$text[(i+1):nrow(posmat)] <- gsub(posmat$text[i], posmat$name[i], posmat$text[(i+1):nrow(posmat)], fixed=TRUE)
      }
    }
  } else {
    converted_string <- input_string
  }
  
  # Process "+" or "-" or " " between comma
  cssnum <- 1
  retcss <- function(converted_string, cssnum) {
    css <- NULL
    for (i in gsub("\\)","",gsub("\\(","",unlist(strsplit(converted_string, ","))))){
      if (nchar(i)!=6) {
        css <- rbind(css, c(i,paste0("CS",str_pad(cssnum,4,pad="0"))))
        cssnum <- cssnum + 1
      }
    }
    if (!is.null(css)) {
      css <- data.frame(css) |> `colnames<-`(c("text","name"))
    }
    list(css=css, num=cssnum)
  }
  
  retcss1 <- retcss(converted_string, cssnum)
  
  css <- retcss1$css
  cssnum <- retcss1$num
  
  for (i in seq_along(css$text)) {
    converted_string <- gsub(css$text[i], css$name[i], converted_string, fixed = TRUE)
  }
  
  converted_string
  
  altnodes <- gsub("\\)","",gsub("\\(","",unlist(strsplit(converted_string,","))))
  rels <- NULL
  for (i in altnodes) {
    for (j in altnodes) {
      if (i!=j) {
        rels <- rbind(rels, c(i,j, "alt"))
      }
    }
  }
  
  if (!is.null(rels)) {
    rels <- rels |> data.frame() |> `colnames<-`(c("from","to","type"))
    
    relg <- graph_from_data_frame(rels,directed = FALSE)
    g <- simplify(relg,edge.attr.comb = "first")
    
    noparen <- function(x) gsub("\\)","",gsub("\\(","",x))
    
    gs <- list()
    for (k in seq_along(posmat$name)) {
      if (k==nrow(posmat)) {next}
      i <- posmat$name[k]
      tmpg <- NULL
      tmp <- subset(posmat, posmat$name==i)$text
      checkcss <- retcss(tmp, cssnum) # check plus spacee
      if (is.null(checkcss$css)) {## No space or plus in text
        commas <- unlist(strsplit(noparen(tmp), ","))
        for (co in commas) {
          for (co2 in commas) {
            if (co!=co2) {
              tmpg <- rbind(tmpg, c(co,co2,"alt"))
              tmpg <- rbind(tmpg, c(i, co,"ingroup"))
            }
          }
        }
      } else {
        noparentmp <- noparen(tmp)
        tmpcss <- retcss(noparentmp, cssnum)
        css2 <- tmpcss$css
        css <- rbind(css, css2)
        cssnum <- tmpcss$num
        for (j in seq_along(css2$text)) {
          noparentmp <- gsub(css2$text[j], css2$name[j], noparentmp, fixed=TRUE)
        }
        commas <- unlist(strsplit(noparentmp, ","))
        if (length(commas)!=1) {
          for (co in commas) {
            for (co2 in commas) {
              if (co!=co2) {
                tmpg <- rbind(tmpg, c(co,co2,"alt"))
                tmpg <- rbind(tmpg, c(i, co,"ingroup"))
              }
            }
          }
        } else {## No comma
          tmpg <- rbind(tmpg, c(commas, i, "ingroup"))
        }##CSS
      }
      tmp_g <- data.frame(tmpg) |> `colnames<-`(c("from","to","type"))
      
      el <- simplify(graph_from_data_frame(tmpg, directed=FALSE), edge.attr.comb = "first")
      tmp_g <- data.frame(as_data_frame(el)) |> `colnames<-`(c("from","to","type"))
      gs[[i]] <- tmp_g
    }
    
    all_g <- as_data_frame(g) |> `colnames<-`(c("from","to","type"))
    
    for (i in names(gs)) {
      all_g <- rbind(all_g, gs[[i]])
    }
  }
  
  ## Plus, space, minus
  ## There are no comma left
  return_certain_sep <- function(text, sep) {
    plls <- NULL
    for (tex in text) {
      pl <- unlist(strsplit(tex, sep))
      for (p in pl) {
        for (q in pl) {
          if (p!=q)
            plls <- rbind(plls, c(p,q,gsub("\\\\","",sep)))
        }
      }
    }
    plls
  }
  
  
  parse_css <- function (row) {
    text <- row$text
    csname <- row$name
    gra <- NULL
    space_sep <- unlist(strsplit(text, " ")) ## NEED TO BE DIRECTED!
    if (length(space_sep)>1) {
      for (i in seq_along(space_sep)) {
        if (i==length(space_sep)) {next}
        gra <- rbind(gra, c(space_sep[i], space_sep[i+1], "rel"))
      }
    } else {
      gra <- return_certain_sep(space_sep, "\\+")
    }
    for (i in unique(c(gra[,1], gra[,2]))) {
      gra<-rbind(gra, c(as.character(csname), i, "incss"))
    }
    return(gra)
  }
  
  ## Plus signs are used to represent a complex and 
  ## a minus sign denotes a non-essential component in the complex.
  if (!is.null(css)) {
    cssparsed <- NULL
    for (i in seq_len(nrow(css))) {
      cssparsed <- rbind(cssparsed, parse_css(css[i,]))
    }
    
    cssparsed <- data.frame(cssparsed) |>
      `colnames<-`(c("from","to","type"))
  } else {
    cssparsed<-NULL
  }
  
  ## Minus
  if (!is.null(cssparsed)) {
    for (i in seq_len(nrow(cssparsed))) {
      row <- cssparsed[i,]
      if (grepl("-", row$from)) {
        mi <- unlist(strsplit(row$from,"-"))
        for (j in seq_along(mi)) {
          if (j==length(mi)) {next}
          ## Append the relationship recursively
          cssparsed <- rbind(cssparsed, c(mi[j],mi[j+1],"-"))
        }
        ## Leave the original row the first occurrence
        cssparsed[i, "from"] <- mi[1]
      }
      if (grepl("-", row$to)) {
        mi <- unlist(strsplit(row$to,"-"))
        for (j in seq_along(mi)) {
          if (j==length(mi)) {next}
          ## Append the relationship recursively
          cssparsed <- rbind(cssparsed, c(mi[j],mi[j+1],"-"))
        }
        ## Leave the original row the first occurrence
        cssparsed[i, "to"] <- mi[1]
      }
    }
  }
  
  norel <- subset(cssparsed, cssparsed$type!="rel")
  reledges <- subset(cssparsed, cssparsed$type=="rel") ## Need to preserve order
  if (!is.null(rels)) {
    plotg <- simplify(graph_from_data_frame(rbind(all_g, norel), directed = FALSE),edge.attr.comb = "first")
    plotg <- graph_from_data_frame(rbind(as_data_frame(plotg), reledges), directed=TRUE)
  } else {
    plotg <- simplify(graph_from_data_frame(norel, directed = FALSE),edge.attr.comb = "first")
    plotg <- graph_from_data_frame(rbind(as_data_frame(plotg), reledges), directed=TRUE)
  }
  return(plotg)
}



#' obtain_module
#' @export
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
#' @export
#' @noRd
parse_module <- function(mod, type="reaction") {
  if (type=="reaction") {

    ## Try to represent reaction as nodes
    ## This fails because some steps use the same reaction
    # reac <- NULL
    # for (rea in mod$reaction) {
    #   left <- unlist(strsplit(rea, "->"))[1]
    #   right <- unlist(strsplit(rea, "->"))[2]
    #   if (grepl("\\+",right)) {
    #     right <- unlist(strsplit(right, "\\+"))
    #   }
    #   right <- gsub(" ","", right)
      
    #   # left
    #   left2 <- gsub(" ", "", unlist(strsplit(left, "  "))[2])
    #   left1 <- gsub(" ", "", unlist(strsplit(left, "  "))[1])
    #   for (l2 in unlist(strsplit(left2, ","))) {
    #     for (ll2 in unlist(strsplit(l2, "\\+"))) {
    #       for (l1 in unlist(strsplit(left1, ","))) {
    #         reac <- rbind(reac, c(ll2, l1))
    #       }
    #     }
    #   }
    #   for (l1 in unlist(strsplit(left1, ","))) {
    #     for (r in right) {
    #       reac <- rbind(reac, c(l1, r))
    #     }       
    #   }
    # }

    ## Try to represent reaction as edge label
    reac <- NULL
    for (rea in mod$reaction) {
      left <- unlist(strsplit(rea, "->"))[1]
      right <- unlist(strsplit(rea, "->"))[2]
      if (grepl("\\+",right)) {
        right <- unlist(strsplit(right, "\\+"))
      }
      right <- gsub(" ","", right)
      
      # left
      ## Reaction
      left2 <- gsub(" ", "", unlist(strsplit(left, "  "))[2])
      left1 <- gsub(" ", "", unlist(strsplit(left, "  "))[1])
      for (l2 in unlist(strsplit(left2, ","))) {
        l2sp <- unlist(strsplit(l2, "\\+"))
        for (ll2 in seq_along(l2sp)) {
          for (r1 in seq_along(right)) {
            reac <- rbind(reac, c(l2sp[r1], right[r1], left1))
          }
        }
      }
      reac
    }

    return(reac)
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
    pattern <- "K\\d{5}"
    matches <- str_extract_all(mod$definition, pattern)
    num_step <- NULL
    ko_in_step <- list()
    for (i in seq_along(result)) {
      ko_in_step[[i]] <- unlist(str_extract_all(result[i], pattern))
      num_step <- c(num_step, length(unlist(str_extract_all(result[i], pattern))))
    }
    message("Currently returning list of the KO and steps")
    node_list <- list(step=result, KO=unlist(matches), num_in_step=num_step, ko_in_step=ko_in_step)
    return(node_list)
  }
}