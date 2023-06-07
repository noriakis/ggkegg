setOldClass("tbl_graph")
setClass("kegg_module",
         slots=list(
           ID="character",
           name="character",
           definition_raw="list",
           definitions="list",
           definition_components="character",
           orthology="character",
           kegg_class="character",
           pathway="character",
           reaction="character",
           compound="character",
           rmodule="character",
           reaction_each="tbl_df",
           reaction_each_raw="tbl_df",
           reaction_graph="tbl_graph",
           reaction_components="character"
         ))

#' @importFrom GetoptLong qqcat
setMethod("show",
          signature(object="kegg_module"),
          function(object) {
            qqcat("@{object@ID}\n")
            qqcat("@{object@name}\n")
          })

# setOldClass("tbl_graph")
# setClass("kegg_module",
#   slots=list(
#         ID="character",
#         name="character",
#         definition_raw="list",
#         definitions="list",
#         # definition_block="vector",
#         # definition_kos="character",
#         # definition_num_in_block="vector",
#         # definition_ko_in_block="list",
#         definition_components="character",

#         reaction="character",
#         reaction_each="tbl_df",
#         reaction_each_raw="tbl_df",
#         reaction_graph="tbl_graph",
#         reaction_components="character"
#         ))

# #' @importFrom GetoptLong qqcat
# setMethod("show",
#   signature(object="kegg_module"),
#   function(object) {
#     qqcat("@{object@ID}\n")
#     qqcat("@{object@name}\n")
#   })
# #' KEGG module parsing function
# #' module
# #' @return list of module definition and reaction
# #' @export
# module <- function(mid) {
#   kmo <- new("kegg_module")
#   kmo@ID <- mid
#   if (!file.exists(mid)) {
#     download.file(paste0("https://rest.kegg.jp/get/",mid),
#                   destfile=mid)
#   }
#   con = file(mid, "r")
#   reac <- FALSE
#   defflg <- FALSE
#   defnum <- 1
#   reacs <- NULL
#   definitions <- list()
#   while ( TRUE ) {
#     line = readLines(con, n = 1)
#     if ( length(line) == 0 ) {
#       break
#     }
#     if (gsub(" ","",line)=="") {next} ## If blank line, skip
#     ## happens in some modules describing parallel reactions
#     if (grepl("NAME", line)) {
#       name <- unlist(strsplit(line, "        "))[2]
#       kmo@name <- name
#     }
#     if (grepl("ORTHOLOGY", line)) {defflg <- FALSE}
#     if (grepl("DEFINITION", line)) {defflg <- TRUE}
#     if (defflg) {
#       if (grepl("DEFINITION", line)) {
#         definitions[[defnum]] <- unlist(strsplit(line, "  "))[2]
#       } else {
#         definitions[[defnum]] <- unlist(strsplit(line, "            "))[2]
#       }
#       defnum <- defnum + 1
#     }
#     if (grepl("COMPOUND", line)) {reac <- FALSE}
#     if (grepl("REACTION", line)) {reac <- TRUE}
#     if (reac) {
#         if (grepl("REACTION", line)) {
#           reacs <- c(reacs, unlist(strsplit(line, "    "))[2])
#         } else {
#           reacs <- c(reacs, unlist(strsplit(line, "            "))[2])
#         }
#     }
#   }
#   if (!is.null(definitions)) {
#       kmo@definition_raw <- definitions
#   }
#   if (!is.null(reacs)) {
#     kmo@reaction <- reacs
#   }
#   close(con)
#   pattern <- "K\\d{5}"
#   kos <- paste0("ko:",unlist(str_extract_all(kmo@definition_raw |> unlist(), pattern)))
#   pattern <- "C\\d{5}"
#   cos <- paste0("cpd:",unlist(str_extract_all(kmo@reaction, pattern)))
#   pattern <- "R\\d{5}"
#   rns <- paste0("rn:",unlist(str_extract_all(kmo@reaction, pattern)))
#   kmo@definition_components <- c(kos)
#   kmo@reaction_components <- c(cos, rns)
#   kmo <- parse_module(kmo)
#   kmo
# }

#' module
#' KEGG module parsing function
#' @param mid KEGG module ID
#' @return list of module definition and reaction
#' @examples \donttest{module("M00003")}
#' @export
module <- function(mid) {
  kmo <- new("kegg_module")
  kmo@ID <- mid
  if (!file.exists(mid)) {
    download.file(paste0("https://rest.kegg.jp/get/",mid),
                  destfile=mid)
  }
  con = file(mid, "r")
  content_list <- list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (!startsWith(line, " ")) {
      current_id <- strsplit(line, " ") |> sapply("[", 1)
    }
    if (!current_id %in% c("REFERENCE","///")) {
      content <- substr(line, 13, nchar(line))
      content_list[[current_id]] <- c(content_list[[current_id]], content) 
    }
  }
  close(con)
  kmo@ID <- substr(content_list$ENTRY, 1, 6)
  kmo@name <- content_list$NAME
  if (!is.null(content_list$DEFINITION)) {
      kmo@definition_raw <- as.list(content_list$DEFINITION)
  }
  if (!is.null(content_list$ORTHOLOGY)) {
      kmo@orthology <- content_list$ORTHOLOGY
  }  
  if (!is.null(content_list$CLASS)) {
    kmo@kegg_class <- content_list$CLASS
  }  
  if (!is.null(content_list$PATHWAY)) {
    kmo@pathway <- content_list$PATHWAY
  }
  if (!is.null(content_list$REACTION)) {
    reacs <- content_list$REACTION
    kmo@reaction <- reacs[reacs!=""]
  }
  if (!is.null(content_list$COMPOUND)) {
    kmo@compound <- content_list$COMPOUND
  }
  if (!is.null(content_list$RMODULE)) {
    kmo@rmodule <- content_list$RMODULE
  }
  pattern <- "K\\d{5}"
  kos <- paste0("ko:",unlist(str_extract_all(kmo@definition_raw |> unlist(), pattern)))
  pattern <- "C\\d{5}"
  cos <- paste0("cpd:",unlist(str_extract_all(kmo@reaction, pattern)))
  pattern <- "R\\d{5}"
  rns <- paste0("rn:",unlist(str_extract_all(kmo@reaction, pattern)))
  kmo@definition_components <- c(kos)
  kmo@reaction_components <- c(cos, rns)
  kmo <- parse_module(kmo)
  kmo
}


#' module_text
#' Obtain textual representation of module definition for all the blocks
#' @param kmo module object
#' @param name name of definition
#' @param candidate_ko KO to highlight
#' @param paint_colour color to highlight
#' @param convert named vector converting the KO to gene name
#' @export
#' @examples
#' mo <- create_test_module()
#' tex <- module_text(mo)
#' @return textual description of module definitions
module_text <- function(kmo, name="1", candidate_ko=NULL, paint_colour="tomato", convert=NULL) {
  kmo <- kmo@definitions[[name]]
  plot_list <- list()
  for (block in seq_along(kmo$definition_block)) {

    input_string <- kmo$definition_block[block]
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
      posmat$name <- paste0("G",str_pad(seq_len(nrow(posmat)),5,pad="0"))
      ul <- sort(unique(posmat$length))
      he <- (seq_len(length(ul)))+1
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
    for (i in unlist(kmo$definition_ko_in_block[block])) {
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
        if (put==" ") {put <- "and"}
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
    concat$conflag <- concat$name %in% c("and","+","-","alt")
    plot_list[[block]] <- concat
  }
  plot_list
}

#' module_completeness
#' 
#' @param kmo module object
#' @param query vector of KO
#' @param name name of definitions when multiple definitions are present
#' @examples
#' test_complete <- module_completeness(create_test_module(), c("K00112"))
#' @export
#' @return tibble
module_completeness <- function(kmo, query, name="1") {
  if (length(kmo@definitions)>1) {message("Multiple definitions found, taking specified definition");
    message(paste0("  number of definitions: ",length(kmo@definitions)));
  }
  kmo <- kmo@definitions[[name]] ## Take first definition
  complete <- NULL
  pres <- NULL
  pres_ratio <- NULL
  alls <- NULL
  for (i in seq_along(kmo$definition_ko_in_block)) {
    if (sum(identical(kmo$definition_ko_in_block[[i]],
      character(0)))!=0) {message("No KO found in definition block. Returning NULL");
      message(paste0("  ", paste(kmo$definition_block, collapse=" ")));
      return(NULL)}
    present <- kmo$definition_ko_in_block[[i]] %in% query
    names(present) <- kmo$definition_ko_in_block[[i]]
    bool <- gsub("\\+","&",gsub(" ", "&", gsub(",", "|", kmo$definition_block[i])))
    for (j in names(present)) {
      bool <- gsub(j, present[j], bool)
    }
    alls <- c(alls, length(kmo$definition_ko_in_block[[i]]))
    complete <- c(complete, eval(parse(text=bool)))
    pres <- c(pres, sum(present))
    pres_ratio <- c(pres_ratio, sum(present) / kmo$definition_num_in_block[i])
  }
  tibble(block=kmo$definition_block,
         all_num=alls,
         present_num=pres,
         ratio=pres_ratio,
         complete=complete)
}

#' obtain_sequential_module_definition
#' 
#' Given module definition and block number,
#' Recursively obtain graphical represencation of block and 
#' connect them by pseudo-nodes representing blocks.
#' @param kmo module object
#' @param name name of definition when multiple definitions are present
#' @param block specify if need to parse specific block
#' @export
#' @examples
#' mo <- create_test_module()
#' sequential_mod <- obtain_sequential_module_definition(mo)
#' @return list of module definitions
obtain_sequential_module_definition <- function(kmo, name="1", block=NULL) {
  kmo <- kmo@definitions[[name]]
  if (is.null(block)) {cand_step <- kmo$definition_block} else {
    cand_step <- block
  }
  all_steps <- NULL
  orders <- NULL

  for (i in seq_along(cand_step)) {
    if (kmo$definition_num_in_block[i]!=1) {
      for (ko in kmo$definition_ko_in_block[[i]]) {
        all_steps <- rbind(all_steps, c(ko, paste0("BLOCK",i),"in_block"))
        all_steps <- all_steps |> data.frame() |> `colnames<-`(c("from","to","type"))
      }
      plotg <- as_data_frame(module_graph(kmo$definition_block[i]))
      ## Need to change naming
      frm <- plotg$from
      frm[!startsWith(frm,"K")] <- paste0(frm[!startsWith(frm,"K")],"_",i)
      to <- plotg$to
      to[!startsWith(to,"K")] <- paste0(to[!startsWith(to,"K")],"_",i)
      plotg$from <- frm
      plotg$to <- to

      all_steps <- rbind(all_steps, plotg)
      orders <- c(orders, paste0("BLOCK",i))
    } else {
      ## Some steps have only "-K*****"
      orders <- c(orders, kmo$definition_block[i])
    }
  } ## For each step, obtain the plot
  
  ords <- NULL
  for (i in seq_along(orders)) {
    if (i!=length(orders)) {
      ords <- rbind(ords, c(orders[i], orders[i+1],"block_transition"))
    }
  }
  if (!is.null(ords)) {
    all_steps <- rbind(all_steps, 
                     data.frame(ords)|>
                       `colnames<-`(c("from","to","type")))
  } else {
    ## Only one step
    all_steps <- subset(all_steps, all_steps$type!="in_block")
  }
  return(as_tbl_graph(all_steps))
}


#' module_graph
#' obtain graphical representation of module definition
#' @return igraph object
#' @noRd
module_graph <- function(input_string, skip_minus=FALSE) {
  ## [TODO] need verbose to identify what this function does
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
    posmat$name <- paste0("G",str_pad(seq_len(nrow(posmat)),5,pad="0"))
    ul <- sort(unique(posmat$length))
    he <- (seq_len(length(ul)))+1
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
  css <- NULL
  
  if (!skip_minus) {
    if (grepl("\\+",converted_string) | 
        grepl("-",converted_string) | 
        grepl(" ",converted_string)) {
      retcss1 <- retcss(converted_string, cssnum)
      css <- retcss1$css
      cssnum <- retcss1$num
      
      for (i in seq_along(css$text)) {
        converted_string <- gsub(css$text[i], css$name[i], converted_string, fixed = TRUE)
      }
   }
  } else {
    if (grepl("\\+",converted_string) |
        grepl(" ",converted_string)) {
      retcss1 <- retcss(converted_string, cssnum)
      css <- retcss1$css
      cssnum <- retcss1$num
      
      for (i in seq_along(css$text)) {
        converted_string <- gsub(css$text[i], css$name[i], converted_string, fixed = TRUE)
      }
    }
  }
  
  

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
              tmpg <- rbind(tmpg, c(i, co,"in_group"))
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
                tmpg <- rbind(tmpg, c(i, co,"in_group"))
              }
            }
          }
        } else {## No comma
          tmpg <- rbind(tmpg, c(commas, i, "in_group"))
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
      gra<-rbind(gra, c(as.character(csname), i, "in_and"))
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
  if (!skip_minus) {
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
  }
  
  # norel <- subset(cssparsed, cssparsed$type!="rel")
  # reledges <- subset(cssparsed, cssparsed$type=="rel") ## Need to preserve order
  if (!is.null(rels)) {
    if (!is.null(cssparsed)) {
      plotg <- simplify(graph_from_data_frame(rbind(all_g, cssparsed), directed = FALSE),
                        edge.attr.comb = "first")
    } else {
      plotg <- simplify(graph_from_data_frame(rbind(all_g), directed = FALSE),edge.attr.comb = "first")
    }
  } else {
    if (!is.null(cssparsed)) {
      plotg <- simplify(graph_from_data_frame(rbind(cssparsed), directed = FALSE),
                        edge.attr.comb = "first")
    } else {
      plotg <- input_string
    }
    # plotg <- simplify(graph_from_data_frame(norel, directed = FALSE),edge.attr.comb = "first")
    # plotg <- graph_from_data_frame(rbind(as_data_frame(plotg), reledges), directed=TRUE)
  }
  return(plotg)
}



#' parse_module
#' @importFrom dplyr tibble
#' @noRd
parse_module <- function(kmo) {

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
    each_reacs <- NULL
    each_reacs_raw <- NULL
    for (rea in kmo@reaction) {
      left <- unlist(strsplit(rea, "->"))[1]
      right <- unlist(strsplit(rea, "->"))[2]
      right_raw <- right
      if (grepl("\\+",right)) {
        right <- unlist(strsplit(right, "\\+"))
      }
      right <- gsub(" ","", right)
      Cpattern <- "C\\d{5}"
      Rpattern <- "R\\d{5}"
      left_Rmatches <- str_extract_all(left, Rpattern) |> unlist() |> tibble()
      left_Cmatches <- str_extract_all(left, Cpattern) |> unlist() |> tibble()
      right_Cmatches <- str_extract_all(right, Cpattern) |> unlist() |> tibble()
      each_reacs <- rbind(each_reacs, c(left_Cmatches, left_Rmatches, right_Cmatches))
      
      # left
      ## Reaction
      left2 <- gsub(" ", "", unlist(strsplit(left, "  "))[2])
      left1 <- gsub(" ", "", unlist(strsplit(left, "  "))[1])
      if (is.na(left2)) {
        message("Some modules cannot be parsed properly by the delimiter '  ', changing the split parameter")
        message(paste0("  ",left))
        left2 <- gsub(" ", "", unlist(strsplit(left, " "))[2])
        left1 <- gsub(" ", "", unlist(strsplit(left, " "))[1])        
      }
      each_reacs_raw <- rbind(each_reacs_raw,
        c(left2, left1, right_raw |> gsub(" ","",x=_)))
      for (l2 in unlist(strsplit(left2, ","))) {
        l2sp <- unlist(strsplit(l2, "\\+"))
        for (ll2 in seq_along(l2sp)) {
          for (r1 in seq_along(right)) {
            ## Orders important or not
            reac <- rbind(reac, c(l2sp[ll2], right[r1], left1))
          }
        }
      }
      reac
    }
    if (!is.null(reac)) {
      each <- as_tibble(each_reacs)
      names(each) <- c("left","reaction","right")
      eachraw <- as_tibble(each_reacs_raw)
      names(eachraw) <- c("left","reaction","right")

      reac <- reac |> data.frame() |> `colnames<-`(c("from","to","reaction"))
      kmo@reaction_graph <- as_tbl_graph(reac)
      kmo@reaction_each <- each
      kmo@reaction_each_raw <- eachraw
    }

    divide_string <- function(input_string) {
      steps <- c()
      current_step <- ""
      open_parentheses <- 0
      
      for (i in seq_len(nchar(input_string))) {
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
    
    definitions <- list()
    for (defnum in seq_along(kmo@definition_raw)) {
      result <- divide_string(kmo@definition_raw[[defnum]])
      result <- result[result!=""]
      if ("--" %in% result) {
        message("Found '--' sep, ignoring");
        message(paste0("  ",kmo@definition_raw[[defnum]]));
        result <- result[result!="--"]
      }
      pattern <- "K\\d{5}"
      matches <- str_extract_all(kmo@definition_raw[[defnum]], pattern)
      num_step <- NULL
      ko_in_step <- list()
      for (i in seq_along(result)) {
        ko_in_step[[i]] <- unlist(str_extract_all(result[i], pattern))
        num_step <- c(num_step, length(unlist(str_extract_all(result[i], pattern))))
      }
      ind <- defnum |> as.character()
      definitions[[ind]][["definition_block"]] <- result
      definitions[[ind]][["definition_kos"]] <- unlist(matches)
      definitions[[ind]][["definition_num_in_block"]] <- num_step
      definitions[[ind]][["definition_ko_in_block"]] <- ko_in_step
    }

    
    kmo@definitions <- definitions
    # kmo@definition_kos <- unlist(matches)
    # kmo@definition_num_in_block <- num_step
    # kmo@definition_ko_in_block <- ko_in_step
    
    return(kmo)
}

#' module_abundance
#' weighted mean abundance of fraction of present KO in the block
#' @param mod_id module ID
#' @param vec KO-named vector of abundance without prefix `ko:`
#' @param num definition number when multiple definitions are present
#' @param calc calculation of final results, mean or weighted_mean
#' @export
#' @importFrom stats weighted.mean
#' @examples \donttest{module_abundance("M00003",c(1.2) |> setNames("K00927"))}
#' @return numeric value
module_abundance <- function(mod_id, vec, num=1, calc="weighted_mean") {
  mod <- module(mod_id)
  ko_abun <- NULL
  for (kos in mod@definitions[[num]]$definition_ko_in_block) {
    if (length(intersect(kos, names(vec))) >= 1) {
      mean_kos <- vec[intersect(kos, names(vec))] |> mean()
    } else {
      mean_kos <- 0
    }
    ko_abun <- c(ko_abun, mean_kos)
  }
  comp <- module_completeness(mod, names(vec))
  comp$abundance <- ko_abun
  comp$strict <- comp$complete * comp$abundance
  if (calc=="mean") {
    mean(comp$strict)
  } else if (calc=="weighted_mean") {
    weighted.mean(comp$strict, w=comp$ratio)
  } else {
    stop("calc must be mean or weighted_mean")
  }
}

#' pathway_abundance
#' @param id pathway id
#' @param vec named vector of abundance
#' @param num number of module definition
#' @return numeric value
#' @examples \donttest{pathway_abundance("ko00270", c(1.2) |> `setNames`("K00927"))}
#' @export
pathway_abundance <- function(id, vec, num=1) {
  pway <- pathway_info(id)
  mods <- pway$MODULE |> strsplit(" ") |> sapply("[", 1) |> unique()
  abuns <- NULL
  for (mod in mods) {
    abuns <- c(abuns, module_abundance(mod_id=mod, num=num, vec=vec))
  }
  tibble(
    module=mods,
    abundance=abuns
  )
}


#' create_test_module
#' @export
#' @examples create_test_module()
#' @return return a test module to use in examples
create_test_module <- function() {
  ## Completely random module
  mo <- new("kegg_module")
  mo@ID <- "test"
  mo@name <- "test module"
  mo@reaction_each <- tibble(left=list("C00112"),reaction=list("R00112"),
    right=list("C00224"))
  mo@reaction_each_raw <- tibble(left="C00112",reaction="R00112",right="C00224")
  mo@definition_raw <- list(c("K00112+K00224"))
  mo@definitions <- list("1"=list("definition_block"="K00112+K00224",
                                  "definition_kos"=c("K00112","K00224"),
                                  "definition_num_in_block"=2,
                                  "definition_ko_in_block"=list(c("K00112","K00224"))))
  mo
}