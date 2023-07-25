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
    )
)

#' @importFrom GetoptLong qqcat
setMethod("show",
	signature(object="kegg_module"),
	function(object) {
	    qqcat("@{object@ID}\n")
        qqcat("@{object@name}\n")
    }
)

#' module
#' KEGG module parsing function
#' @param mid KEGG module ID
#' @param use_cache use cache
#' @param directory directory to save raw files
#' @return list of module definition and reaction
#' @examples \dontrun{module("M00003")}
#' @export
module <- function(mid, use_cache=FALSE, directory=NULL) {
	kmo <- new("kegg_module")
	kmo@ID <- mid

	if (!is.null(directory)) {
    	dest <- paste0(directory,"/",mid)
  	} else {
    	dest <- mid
  	}
	if (!file.exists(dest)) {
    	if (use_cache) {
      		bfc <- BiocFileCache()
      		dest <- bfcrpath(bfc,
        		paste0("https://rest.kegg.jp/get/",mid))
    	} else {
      		download.file(paste0("https://rest.kegg.jp/get/",mid),
                destfile=dest)
    	}
  	}
	con <- file(dest, "r")
	content_list <- list()
	while ( TRUE ) {
    	line <- readLines(con, n = 1)
    	if ( length(line) == 0 ) {
      		break
    	}
	    if (!startsWith(line, " ")) {
    		current_id <- strsplit(line, " ") |>
    		    vapply("[", 1, FUN.VALUE="character")
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
	kos <- paste0("ko:",
		unlist(str_extract_all(kmo@definition_raw |> unlist(), pattern)))
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
module_text <- function(kmo, name="1", candidate_ko=NULL,
	paint_colour="tomato", convert=NULL) {
    kmo <- kmo@definitions[[name]]
    plot_list <- lapply(seq_along(kmo$definition_block),
        function(x) {
            input_string <- kmo$definition_block[x]
            ppos_list <- lapply(find_parenthesis_pairs(input_string),
                function(y) {
                    c(y[[1]], y[[2]], y[[2]]-y[[1]]) 
                })
            ppos <- do.call(rbind, ppos_list)
            if (!is.null(ppos)) {
                posmat <- ppos[order(ppos[,3]),]
                if (is.vector(posmat)) {dfs <- data.frame(t(posmat))} else {
                  dfs <- data.frame(posmat)
                }
                posmat <- dfs |> `colnames<-`(c("xmin","xmax","length"))
                posmat$name <- paste0("manual_G",
                	str_pad(seq_len(nrow(posmat)),5,pad="0"))
                ul <- sort(unique(posmat$length))
                he <- (seq_len(length(ul)))+1
                names(he) <- ul
                posmat$height <- he[as.character(posmat$length)]/2
                posmat$text <- apply(posmat, 1,
                	function(row)
                		substr(input_string, row["xmin"], row["xmax"]))
                posmat$rawtext <- apply(posmat, 1,
                	function(row)
                		substr(input_string, row["xmin"], row["xmax"]))
                posmat$size <- posmat$length
                posmat$x <- (as.numeric(posmat$xmin)+as.numeric(posmat$xmax))/2
            }

            kopos <- lapply(unlist(kmo$definition_ko_in_block[x]),
                function(z) {
                    findko <- str_locate_all(input_string, z)
                    unlist(lapply(findko, function(i) {
                        unlist(lapply(seq_len(nrow(i)), function(j) {
                            c(z, i[j, 1], i[j, 2])
                        }))
                    }))
                })
            kopos <- do.call(rbind, kopos)
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
                    concat <- rbind(concat,
                    	c(put, pos[rn, 1], pos[rn, 1], pos[rn, 1], 0.5, 1))
                  }
                }
              }
              concat
            }

            if (!is.null(ppos)) {
              concat <- rbind(
              		posmat[,c("name","xmin","xmax","x","height","size")], kopos
              	)
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

            concat$color <- lapply(concat$name, function(k) {
                if (k %in% candidate_ko) {
                  paint_colour
                } else {
                  "transparent"
                }
            }) |> unlist()

            if (!is.null(convert)) {
                conv <- convert[concat$name]
                conv[is.na(conv)] <- concat$name[is.na(conv)]
                concat$converted_name <- conv                
            }
            concat$koflag <- startsWith(concat$name, "K")
            concat$conflag <- concat$name %in% c("and","+","-","alt")
            concat
        })
    plot_list
}

#' module_completeness
#' 
#' This converts module definitions consisting of KO identifiers to the expression
#' by converting `+` and ` ` to `AND`, and `,` to `OR`. After that, KO IDs specified
#' by `query` is inserted to expression by `TRUE` or `FALSE`, and is evaluated.
#' Please feel free to contact the bug, or modules that cannot be calculated.
#' (Module definitions consisting of module IDs [M*] cannot be calculated)
#' 
#' Below is quoted from https://www.genome.jp/kegg/module.html
#' 
#' `A space or a plus sign, representing a connection in the pathway or the molecular complex,
#' is treated as an AND operator and a comma, used for alternatives, is treated as an OR operator.
#' A minus sign designates an optional item in the complex.`
#' 
#' @param kmo module object
#' @param query vector of KO
#' @param name name of definitions when multiple definitions are present
#' @examples
#' test_complete <- module_completeness(create_test_module(), c("K00112"))
#' @export
#' @return tibble
module_completeness <- function(kmo, query, name="1") {
	if (length(kmo@definitions)>1) {
		message("Multiple definitions found, taking specified definition");
    	message(paste0("  number of definitions: ",length(kmo@definitions)));
  	}
	kmo <- kmo@definitions[[name]] ## Take first definition by default
  	results <- lapply(seq_along(kmo$definition_ko_in_block), function (i) {
  		if (sum(identical(kmo$definition_ko_in_block[[i]], character(0)))!=0) {
  			message("No KO found in definition block. Returning NULL")
		    message(paste0("  ", paste(kmo$definition_block, collapse=" ")))
      		return(NULL)
  		}
	    present <- kmo$definition_ko_in_block[[i]] %in% query
	    names(present) <- kmo$definition_ko_in_block[[i]]
	    bool <- gsub("\\+","&",
	    	gsub(" ", "&", gsub(",", "|", kmo$definition_block[i])))
	    for (j in names(present)) {
	      bool <- gsub(j, present[j], bool)
	    }
	    list(length(kmo$definition_ko_in_block[[i]]),
	    	eval(parse(text=bool)),
	    	sum(present),
	    	sum(present) / kmo$definition_num_in_block[i])
  	})
  	
	tibble(block=kmo$definition_block,
        all_num=lapply(results, function(x) x[[1]]) |> unlist(),
        present_num=lapply(results, function(x) x[[3]]) |> unlist(),
        ratio=lapply(results, function(x) x[[4]]) |> unlist(),
        complete=lapply(results, function(x) x[[2]]) |> unlist())
}

#' obtain_sequential_module_definition
#' 
#' Given module definition and block number,
#' Recursively obtain graphical represencation of block and 
#' connect them by pseudo-nodes representing blocks.
#' 
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
	if (is.null(block)) {
		cand_step <- kmo$definition_block
	} else {
		cand_step <- block
	}
	all_steps <- NULL
	orders <- NULL
	all_steps <- lapply(seq_along(cand_step), function(i) {
		if (kmo$definition_num_in_block[i]!=1) {
			tmp <- lapply(kmo$definition_ko_in_block[[i]], function(ko) {
				return(c(ko, paste0("manual_BLOCK",i),"in_block"))			
			})
			do.call(rbind, tmp)
		}
	})

	all_steps <- do.call(rbind, all_steps) |>
		data.frame() |>
		`colnames<-`(c("from","to","type"))
	
  	plotg <- lapply(seq_along(cand_step), function(i) {
		plotg <- as_data_frame(module_graph(kmo$definition_block[i]))  		
		## Need to change naming
		frm <- plotg$from
		frm[!startsWith(frm,"K")] <- paste0(frm[!startsWith(frm,"K")],"_",i)
		to <- plotg$to
		to[!startsWith(to,"K")] <- paste0(to[!startsWith(to,"K")],"_",i)
		plotg$from <- frm
		plotg$to <- to
		plotg
  	})
  	plotg <- do.call(rbind, plotg)

	all_steps <- rbind(all_steps, plotg)
	
	orders <- lapply(seq_along(cand_step), function(i) {
		if (kmo$definition_num_in_block[i]!=1) {
			return(paste0("manual_BLOCK",i))
		} else {
			return(kmo$definition_block[i])	
		}
	}) |> unlist()

	ords <- lapply(seq_along(orders), function (i) {
	    if (i!=length(orders)) {
	      return(c(orders[i], orders[i+1],"block_transition"))
	    }		
	})
	ords <- do.call(rbind, ords)

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
	ppos <- lapply(find_parenthesis_pairs(input_string), function (y) {
		return(c(y[[1]], y[[2]], y[[2]]-y[[1]]))
  	})
	ppos <- do.call(rbind, ppos)
	if (!is.null(ppos)) {
	    posmat <- ppos[order(ppos[,3]),]
	    if (is.vector(posmat)) {dfs <- data.frame(t(posmat))} else {
	    	dfs <- data.frame(posmat)
	    }
	    posmat <- dfs |> `colnames<-`(c("xmin","xmax","length"))
	    posmat$name <- paste0("manual_G",
	    	str_pad(seq_len(nrow(posmat)),5,pad="0"))
	    ul <- sort(unique(posmat$length))
	    he <- (seq_len(length(ul)))+1
	    names(he) <- ul
	    posmat$height <- he[as.character(posmat$length)]/2
	    posmat$text <- apply(posmat, 1,
	    	function(row) substr(input_string, row["xmin"], row["xmax"]))
	    posmat$rawtext <- apply(posmat, 1,
	    	function(row) substr(input_string, row["xmin"], row["xmax"]))
	    
	    converted_string <- input_string
	    for (i in seq_along(posmat$text)) {
	    	if (i < nrow(posmat)) {
	        	converted_string <- gsub(posmat$text[i],
	        		posmat$name[i],
	        		converted_string,
	        		fixed = TRUE)
	        	posmat$text[(i+1):nrow(posmat)] <- gsub(posmat$text[i],
	        		posmat$name[i],
	        		posmat$text[(i+1):nrow(posmat)],
	        		fixed=TRUE)
	      	}
	    }
	} else {
	    posmat <- NULL
	    converted_string <- input_string
	}

	# Process "+" or "-" or " " between comma
	# Need css number not to be reset, so preallocate
    cssnum <- 1
    retcss <- function(converted_string, cssnum) {
    	alloc <- gsub("\\)","",
    		gsub("\\(","",unlist(strsplit(converted_string, ","))))
    	css <- vector(mode="list", length=length(alloc))
    	j <- 1
    	for (i in alloc) {
			if (nchar(i)!=6) {
        		cssnum <- cssnum + 1
        		css[[j]] <- c(i,paste0("manual_CS",str_pad(cssnum,4,pad="0")))
				j <- j + 1
			}    		
    	}
    	css <- do.call(rbind, css)
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
				converted_string <- gsub(css$text[i],
					css$name[i],
					converted_string,
					fixed = TRUE)
			}
		}
	} else {
		if (grepl("\\+",converted_string) |
		    grepl(" ",converted_string)) {
			
			retcss1 <- retcss(converted_string, cssnum)
			css <- retcss1$css
			cssnum <- retcss1$num
	  
			for (i in seq_along(css$text)) {
				converted_string <- gsub(css$text[i],
					css$name[i],
					converted_string,
					fixed = TRUE)
			}
		}
	}
	

	altnodes <- gsub("\\)","",
		gsub("\\(","",unlist(strsplit(converted_string,","))))
	rels <- lapply(altnodes, function(x) {
		ret1 <- lapply(altnodes, function(y) {
			if (x!=y) {
				return(c(x,y,"alt"))
			}
		})
		do.call(rbind, ret1)
	})
	rels <- do.call(rbind, rels)

	if (!is.null(rels)) {
    	rels <- rels |> data.frame() |> `colnames<-`(c("from","to","type"))    
	    relg <- graph_from_data_frame(rels,directed = FALSE)
    	g <- simplify(relg,edge.attr.comb = "first")
    
    	noparen <- function(x) gsub("\\)","",gsub("\\(","",x))
    	
    	gs <- lapply(seq_along(posmat$name), function(k) {
			if (k==nrow(posmat)) {return(NULL)}
	      	i <- posmat$name[k]
			tmp <- subset(posmat, posmat$name==i)$text
			checkcss <- retcss(tmp, cssnum) # check plus space

			if (is.null(checkcss$css)) {## No space or plus in text
		        commas <- unlist(strsplit(noparen(tmp), ","))
		        tmpg <- lapply(commas, function(co) {
		        	cos <- lapply(commas, function(co2) {
		        		if (co!=co2) {
		        			return(rbind(c(co, co2, "alt"),
		        				c(i, co, "in_group")))
		        		}
		        	})
		        	do.call(rbind, cos)
		        })
		        tmpg <- do.call(rbind, tmpg)
		    } else {
	        	noparentmp <- noparen(tmp)
		        tmpcss <- retcss(noparentmp, cssnum)
		        css2 <- tmpcss$css
		        css <- rbind(css, css2)
		        cssnum <- tmpcss$num

		        for (j in seq_along(css2$text)) {
					noparentmp <- gsub(css2$text[j],
						css2$name[j],
						noparentmp,
						fixed=TRUE)
		        }

		        commas <- unlist(strsplit(noparentmp, ","))
	    	    if (length(commas)!=1) {
			        tmpg <- lapply(commas, function(co) {
			        	cos <- lapply(commas, function(co2) {
			        		if (co!=co2) {
			        			return(rbind(c(co, co2, "alt"),
			        				c(i, co, "in_group")))
			        		}
			        	})
			        	do.call(rbind, cos)
			        })
			        tmpg <- do.call(rbind, tmpg)
        		} else {## No comma
          			tmpg <- c(commas, i, "in_group")
        		}##CSS		    	
		    }
	        tmp_g <- data.frame(tmpg) |> `colnames<-`(c("from","to","type"))      
	        el <- simplify(graph_from_data_frame(tmpg, directed=FALSE),
	        			    edge.attr.comb = "first")
	        tmp_g <- data.frame(as_data_frame(el)) |>
	            `colnames<-`(c("from","to","type"))
	        list(tmp_g, css)
    	})
	    all_g <- as_data_frame(g) |> `colnames<-`(c("from","to","type"))
        all_g <- rbind(all_g, do.call(rbind, lapply(gs, function(x) x[[1]])))
        css <- rbind(css, do.call(rbind, lapply(gs, function(x) x[[2]])))
  	}

	## Plus, space, minus
	## There are no comma left
    return_certain_sep <- function(text, sep) {
		pl <- unlist(strsplit(text, sep))
		plls <- lapply(pl, function(x) {
				tmp_plls <- lapply(pl, function(y) {
					if (x!=y) {
						return(c(x,y,gsub("\\\\","",sep)))
					}
				})
				do.call(rbind, tmp_plls)
    		})
    	do.call(rbind, plls)
  	}

	parse_css <- function (row) {
    	text <- row$text
    	csname <- row$name
    	space_sep <- unlist(strsplit(text, " ")) ## NEED TO BE DIRECTED!
    	if (length(space_sep)>1) {
     		gra <- lapply(seq_along(space_sep), function(i) {
		        if (i==length(space_sep)) {} else {
	                return(c(space_sep[i], space_sep[i+1], "rel"))        	
		        }
     		})
     		gra <- do.call(rbind, gra)
		} else {
	        gra <- return_certain_sep(space_sep, "\\+")
    	}
	    for (i in unique(c(gra[,1], gra[,2]))) {
	      gra <- rbind(gra, c(as.character(csname), i, "in_and"))
	    }
	    return(gra)
  	}

	## Plus signs are used to represent a complex and 
	## a minus sign denotes a non-essential component in the complex.
    if (!is.null(css)) {
    	cssparsed <- lapply(seq_len(nrow(css)), function(x) {
    		parse_css(css[x,])
    	})
    	cssparsed <- do.call(rbind, cssparsed)
    	cssparsed <- data.frame(cssparsed) |>
    		`colnames<-`(c("from","to","type"))
  	} else {
    	cssparsed <- NULL
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
  
    if (!is.null(rels)) {
	    if (!is.null(cssparsed)) {
    	    plotg <- simplify(graph_from_data_frame(rbind(all_g, cssparsed),
    	    			directed = FALSE),
                        edge.attr.comb = "first")
    	} else {
      		plotg <- simplify(graph_from_data_frame(rbind(all_g),
      					directed = FALSE),
      					edge.attr.comb = "first")
    	}
    } else {
  	    if (!is.null(cssparsed)) {
      		plotg <- simplify(graph_from_data_frame(rbind(cssparsed),
      					directed = FALSE),
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

    reac_list <- lapply(kmo@reaction, function(rea) {
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
		
		reac <- lapply(unlist(strsplit(left2, ",")), function(l2) {
	        l2sp <- unlist(strsplit(l2, "\\+"))
	        tmp2 <- lapply(seq_along(l2sp), function(ll2) {
	        	tmp <- lapply(seq_along(right), function(r1) {
	        		return(c(l2sp[ll2], right[r1], left1))
	        	})
	        	do.call(rbind, tmp)
	        })
	        do.call(rbind, tmp2)
		})
		reac <- do.call(rbind, reac)
		list("each_reacs"=c(left_Cmatches, left_Rmatches, right_Cmatches),
			"each_reacs_raw"=c(left2, left1, right_raw |> gsub(" ","",x=_)),
			"reac"=reac)
    })
    reac <- as_tibble(do.call(rbind, lapply(reac_list, function(x) x[["reac"]])))
    if (!is.null(reac)) {
	    each <- as_tibble(do.call(rbind, lapply(reac_list, function(x) x[["each_reacs"]])))
	    names(each) <- c("left","reaction","right")
	    eachraw <- as_tibble(do.call(rbind, lapply(reac_list, function(x) x[["each_reacs_raw"]])))
	    names(eachraw) <- c("left","reaction","right")
		reac <- reac |> data.frame() |> `colnames<-`(c("from","to","reaction"))
		kmo@reaction_graph <- as_tbl_graph(reac)
		kmo@reaction_each <- each
		kmo@reaction_each_raw <- eachraw
    }

    ## Preallocate here
    divide_string <- function(input_string) {
		steps <- vector(mode="list", length=nchar(input_string))
		j <- 1
		current_step <- ""
		open_parentheses <- 0
      
		for (i in seq_len(nchar(input_string))) {
			char <- substr(input_string, i, i)
        
	        if (char == "(") {
				open_parentheses <- open_parentheses + 1
			} else if (char == ")") {
				open_parentheses <- open_parentheses - 1
			} else if (char == " " & open_parentheses == 0) {
				steps[[j]] <- current_step
				j <- j + 1
				current_step <- ""
				next
        	}        
    		current_step <- paste0(current_step, char)
      	}
		steps[[j]] <- current_step
		j <- j + 1
        steps[sapply(steps, is.null)] <- NULL
		return(steps)
    }
	## Preallocate   
    
    definitions <- lapply(seq_along(kmo@definition_raw), function(defnum) {
		result <- divide_string(kmo@definition_raw[[defnum]])
		result <- result[result!=""]
		if ("--" %in% result) {
			message("Found '--' sep, ignoring");
			message(paste0("  ",kmo@definition_raw[[defnum]]));
			result <- result[result!="--"]
      	}
		pattern <- "K\\d{5}"
		matches <- str_extract_all(kmo@definition_raw[[defnum]], pattern)
		ko_in_step <- lapply(seq_along(result), function(i) {
			unlist(str_extract_all(result[i], pattern))
		})
		num_step <- lapply(seq_along(result), function(i) {
			length(unlist(str_extract_all(result[i], pattern)))
		}) |> unlist()
		list("definition_block"=unlist(result),
			"definition_kos"=unlist(matches),
			"definition_num_in_block"=num_step,
			"definition_ko_in_block"=ko_in_step)
    })
    names(definitions) <- as.character(seq_len(length(definitions)))
    kmo@definitions <- definitions
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
#' @examples \dontrun{module_abundance("M00003",c(1.2) |> setNames("K00927"))}
#' @return numeric value
module_abundance <- function(mod_id, vec, num=1, calc="weighted_mean") {
	mod <- module(mod_id)
	ko_abun <- lapply(mod@definitions[[num]]$definition_ko_in_block, function(kos) {
	    if (length(intersect(kos, names(vec))) >= 1) {
	      mean_kos <- vec[intersect(kos, names(vec))] |> mean()
	    } else {
	      mean_kos <- 0
	    }
	    return(mean_kos)		
	}) |> unlist()
	
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
#' @examples \dontrun{pathway_abundance("ko00270", c(1.2) |> `setNames`("K00927"))}
#' @export
pathway_abundance <- function(id, vec, num=1) {
	pway <- pathway_info(id)
	mods <- pway$MODULE |> strsplit(" ") |> vapply("[", 1, FUN.VALUE="character") |> unique()
	abuns <- lapply(mods, function(mod) {
		module_abundance(mod_id=mod, num=num, vec=vec)
	}) |> unlist()
	tibble(
    	module=mods,
    	abundance=abuns
  	)
}


#' create_test_module
#' 
#' Test kegg_module for examples and vignettes.
#' The module has no biological meanings.
#' 
#' @export
#' @examples create_test_module()
#' @return return a test module to use in examples
create_test_module <- function() {
	## Completely random module
	mo <- new("kegg_module")
	mo@ID <- "test"
	mo@name <- "test module"
	mo@reaction_each <- tibble(left=list("C00065"),reaction=list("R00586"),
		right=list("C00979"))
	mo@reaction_each_raw <- tibble(left="C00065",reaction="R00586",right="C00979")
	mo@definition_raw <- list(c("(K00174+K00175,K00382) (K01902+K01903,K01899"))
	mo@definitions <- list("1"=list("definition_block"=c("K00174+K00175,K00382","K01902+K01903,K01899"),
    	"definition_kos"=c("K00174","K00175","K00382","K01902","K01903","K01899"),
    	"definition_num_in_block"=c(3,3),
    	"definition_ko_in_block"=list(c("K00174","K00175","K00382"),c("K01902","K01903","K01899"))))
	mo
}
