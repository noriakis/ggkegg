#' pathway
#' 
#' KEGG pathway parsing function
#'
#' @param pid pathway id
#' @param directory directory to download KGML
#' @param use_cache whether to use BiocFileCache
#' @param add_pathway_id add pathway id to graph, default to TRUE
#' needed for the downstream analysis
#' @param group_rect_nudge nudge the position of group node
#' default to add slight increase to show the group node
#' @param node_rect_nudge nudge the position of all node
#' @param invert_y invert the y position to match with R graphics
#' @param return_image return the image URL
#' @param return_tbl_graph return tbl_graph object, if FALSE, return igraph
#' @return tbl_graph by default
#' @importFrom igraph graph_from_data_frame
#' @import igraph
#' @importFrom tidygraph .G
#' @importFrom XML xmlParse xmlApply
#' @importFrom tibble as_tibble
#' @importFrom utils download.file head tail
#' @examples pathway("hsa04110")
#' @export
pathway <- function(pid,
    directory=NULL,
    use_cache=FALSE,
    group_rect_nudge=2,
    node_rect_nudge=0,
    invert_y=TRUE,
    add_pathway_id=TRUE,
    return_tbl_graph=TRUE,
    return_image=FALSE) {
    
    ## Specification of KGML format is available at:
    ## https://www.genome.jp/kegg/xml/docs/

    file_name <- paste0(pid,".xml")
    if (!is.null(directory)) {
        file_name <- paste0(directory,"/",file_name)
    }
    if (!file.exists(file_name)) {
        if (use_cache) {
            bfc <- BiocFileCache()
            file_name <- bfcrpath(bfc,
                paste0("https://rest.kegg.jp/get/",pid,"/kgml"))
        } else {
            download.file(url=paste0("https://rest.kegg.jp/get/",pid,"/kgml"),
                            destfile=file_name)
        }
    }

    xml <- xmlParse(file_name)
    node_sets <- getNodeSet(xml, "//entry")
    
    ## Preallocate
    all_nodes <- vector(mode="list", length=length(node_sets))
    grs <- vector(mode="list", length=length(node_sets))
    rev_grs <- vector(mode="list", length=length(node_sets))

    node_names <- c("id","name","type","reaction",
                   "graphics_name",
                   "x","y","width","height","fgcolor","bgcolor",
                   "graphics_type","coords")

    pwy <- getNodeSet(xml, "//pathway")[[1]]

    pwy_name <- xmlAttrs(pwy)["name"]
    pwy_org <- xmlAttrs(pwy)["org"]
    pwy_number <- xmlAttrs(pwy)["number"]
    pwy_title <- xmlAttrs(pwy)["title"]
    pwy_image <- xmlAttrs(pwy)["image"]
    pwy_link <- xmlAttrs(pwy)["link"]

    if (return_image) return(pwy_image)
    
    ni <- 1
    for (node in node_sets) {
        id <- xmlAttrs(node)["id"]
        name <- xmlAttrs(node)["name"]
        type <- xmlAttrs(node)["type"]
        reac <- xmlAttrs(node)["reaction"]

        gls <- getNodeSet(node, "graphics")

        ## Preallocate
        mult_coords <- vector(mode="list",
            length=length(xmlApply(gls, function(x) xmlAttrs(x)["coords"])))
        for (gl in gls) {
            glname <- xmlAttrs(gl)["name"]
            gltype <- xmlAttrs(gl)["type"]

            ## If multiple graphics, take the last
            ## parameters and append only the multiple coordinates
            ## Otherwise graph will have duplicate 'original' ID
            
            glcoords <- xmlAttrs(gl)["coords"]
            mult_coords <- c(mult_coords, glcoords)

            x <- as.numeric(xmlAttrs(gl)["x"])
            if (invert_y) {
                y <- -1*as.numeric(xmlAttrs(gl)["y"])
            } else {
                y <- as.numeric(xmlAttrs(gl)["y"])
            }
            
            w <- as.numeric(xmlAttrs(gl)["width"])
            h <- as.numeric(xmlAttrs(gl)["height"])
            fg <- xmlAttrs(gl)["fgcolor"]
            bg <- xmlAttrs(gl)["bgcolor"]
            
            if (type=="group") {
                for (comp in xmlElementsByTagName(node,"component")) {
                    grs[[as.character(id)]] <- 
                        c(grs[[as.character(id)]],
                            as.character(xmlAttrs(comp)["id"]))
                    rev_grs[[as.character(xmlAttrs(comp)["id"])]] <- 
                        c(rev_grs[[as.character(xmlAttrs(comp)["id"])]],
                            as.character(id))
                }
            }   
        }
        all_nodes[[ni]] <- c(id, name, type, reac,
                        glname, x, y, w, h, fg, bg, gltype,
                        paste0(mult_coords %>% unlist(), collapse="|")) %>%
                        setNames(node_names)
        ni <- ni + 1
    }
    
    all_nodes[vapply(all_nodes, is.null, TRUE)] <- NULL
    grs[vapply(grs, is.null, TRUE)] <- NULL
    rev_grs[vapply(rev_grs, is.null, TRUE)] <- NULL

    kegg_nodes <- dplyr::bind_rows(all_nodes) %>% data.frame() %>%
        `colnames<-`(node_names)

    kegg_nodes$x <- as.numeric(kegg_nodes$x)
    kegg_nodes$y <- as.numeric(kegg_nodes$y)
    kegg_nodes$width <- as.numeric(kegg_nodes$width)
    kegg_nodes$height <- as.numeric(kegg_nodes$height)

    kegg_nodes$xmin <- kegg_nodes$x-kegg_nodes$width/2-node_rect_nudge
    kegg_nodes$xmax <- kegg_nodes$x+kegg_nodes$width/2+node_rect_nudge
    kegg_nodes$ymin <- kegg_nodes$y-kegg_nodes$height/2-node_rect_nudge
    kegg_nodes$ymax <- kegg_nodes$y+kegg_nodes$height/2+node_rect_nudge
  
    kegg_nodes[kegg_nodes$type=="group",]$xmin <- 
        kegg_nodes[kegg_nodes$type=="group",]$xmin-group_rect_nudge
    kegg_nodes[kegg_nodes$type=="group",]$ymin <- 
        kegg_nodes[kegg_nodes$type=="group",]$ymin-group_rect_nudge
    kegg_nodes[kegg_nodes$type=="group",]$xmax <- 
        kegg_nodes[kegg_nodes$type=="group",]$xmax+group_rect_nudge
    kegg_nodes[kegg_nodes$type=="group",]$ymax <- 
        kegg_nodes[kegg_nodes$type=="group",]$ymax+group_rect_nudge
  
    kegg_nodes$orig.id <- kegg_nodes$id ## Store ID as orig.id
  
    rel_sets <- getNodeSet(xml, "//relation")
    ## Preallocate
    all_rels <- vector(mode="list", length=length(rel_sets))
    ei <- 1
    rel_names <- c("entry1","entry2","type",
            "subtype_name","subtype_value")
    for (rel in rel_sets) {
        entry1 <- xmlAttrs(rel)["entry1"]
        entry2 <- xmlAttrs(rel)["entry2"]
        rel_type <- xmlAttrs(rel)["type"]
        # rel_subtype <- xmlAttrs(rel[["subtype"]])["name"]
        rel_subtypes <- xmlElementsByTagName(rel,"subtype")
        if (length(rel_subtypes)!=0) {
            for (rs in rel_subtypes) {
                all_rels[[ei]] <- c(entry1, entry2, rel_type,
                    xmlAttrs(rs)["name"], xmlAttrs(rs)["value"]) %>%
                setNames(rel_names)
                ei <- ei + 1
            }
        } else {
            all_rels[[ei]] <- c(entry1, entry2, rel_type, NA, NA) %>%
            setNames(rel_names)
            ei <- ei + 1
        }
    }

    if (length(all_rels) != 0) {
        kegg_edges <- dplyr::bind_rows(all_rels) %>% data.frame() %>%
            `colnames<-`(c("entry1","entry2","type",
                "subtype_name","subtype_value"))
    } else {
        kegg_edges <- NULL
    }

    gr_rels <- lapply(names(grs), function(gr_name) {
        tmp_rel <- lapply(grs[[gr_name]], function(comp_name) {
            ## Pad other values by `in_group`
            return(c(gr_name, comp_name, "in_group", "in_group", "in_group"))         
        })
        do.call(rbind, tmp_rel)
    })
    gr_rels <- do.call(rbind, gr_rels)
    

    if (length(getNodeSet(xml, "//reaction"))!=0) {
        kegg_reac <- get_reaction(xml)
        if (!is.null(kegg_edges)) {kegg_edges$reaction <- NA
            kegg_edges$reaction_id <- NA}
        kegg_edges <- rbind(kegg_edges, kegg_reac)
    }

    ## Append grouping
    if (!is.null(kegg_edges)) {
        if (!is.null(gr_rels)) {
            gr_rels <- gr_rels %>% 
                data.frame() %>% 
                `colnames<-`(c("entry1","entry2","type",
                    "subtype_name","subtype_value"))
            if ("reaction" %in% colnames(kegg_edges)) {
                gr_rels$reaction <- "in_group"
            }
            kegg_edges <- rbind(kegg_edges, gr_rels)
        }
    }

    if (!is.null(kegg_edges)) {
        g <- graph_from_data_frame(kegg_edges, vertices=kegg_nodes)
    } else {
        g <- tbl_graph(nodes=kegg_nodes)
    }


    if (add_pathway_id) {
        V(g)$pathway_id <- pid
        E(g)$pathway_id <- pid
    }
    if (return_tbl_graph) {
        return(as_tbl_graph(g))
    } else {
        return(g)
    }
}
parse_kgml <- pathway

#' process_line
#' 
#' process the KGML containing graphics type of `line`, like
#' global maps e.g. ko01100. Recursively add nodes and edges 
#' connecting them based on `coords` properties in KGML.
#' 
#' We cannot show directed arrows, as coords are not ordered to show direction.
#' 
#' @param g graph
#' @param invert_y whether to invert the position, default to TRUE
#' should match with `pathway` function
#' @param verbose show progress
#' @importFrom tidygraph bind_nodes bind_edges
#' @export
#' @return tbl_graph
#' @examples 
#' ## For those containing nodes with the graphic type of `line`,
#' ## parse the coords attributes to edges.
#' gm_test <- create_test_pathway(line=TRUE)
#' test <- process_line(gm_test)
process_line <- function(g, invert_y=TRUE, verbose=FALSE) {

    df <- as_tbl_graph(g)
    name_col_node <- c("name","x","y","type","original_name","node.orig.id")
    name_col_edge <- c("from","to","type","name",
        "bgcolor","fgcolor","reaction","orig.id")
    results <- lapply(seq_along(V(g)$name), function(i) {
        if (V(g)$graphics_type[i]=="line") {
            raw_name <- V(g)$name[i]
            bgcol <- V(g)$bgcolor[i]
            fgcol <- V(g)$fgcolor[i]
            reac <- V(g)$reaction[i]
            origid <- V(g)$orig.id[i]
            rawco <- V(g)$coords[i]
            
            if (grepl("\\|",rawco)) {
                rawcos <- strsplit(rawco, "\\|") %>% unlist()
            } else {
                rawcos <- rawco
            }
            
            lapply(seq_along(rawcos), function(rc) {
                co <- unlist(strsplit(rawcos[rc], ","))
                lapply(seq_len(length(co)), function(h) {
                    if (is.na(co[h+2])) {return(NULL)}
                    if (h %% 2 == 0) {return(NULL)}
                    ## Assign unique identifiers each node
                    list(
                            c(paste0(raw_name,"_",i,"_",rc,"_",h),
                                co[h], co[h+1], "line", raw_name, origid) %>%
                                setNames(name_col_node),
                            c(paste0(raw_name,"_",i,"_",rc,"_",h+1),
                                co[h+2], co[h+3], "line", raw_name, origid)%>%
                                setNames(name_col_node),
                            c(paste0(raw_name,"_",i,"_",rc,"_",h),
                                paste0(raw_name,"_",i,"_",rc,"_",h+1),
                                "line", raw_name, bgcol, fgcol, reac, origid) %>%
                                setNames(name_col_edge)        
                        )                   
                })
            })
        }
    })
    
    results <- results %>% unlist(recursive=FALSE)
    results <- results %>% unlist(recursive=FALSE)
    results[vapply(results, is.null, TRUE)] <- NULL

    cos <- do.call(rbind, lapply(results, function(x) {
        rbind(x[[1]],x[[2]])
    })) %>% data.frame() %>% `colnames<-`(name_col_node)
    eds <- do.call(rbind, lapply(results, function(x) {
        x[[3]]
    })) %>% data.frame() %>% `colnames<-`(name_col_edge)

    
    cos$x <- as.numeric(cos$x);
    if (invert_y) {
        cos$y <- -1 * as.numeric(cos$y)
    } else {
        cos$y <- as.numeric(cos$y)
    }
    
    df_add <- df %>% bind_nodes(cos) %>% bind_edges(eds)
    df_add %>% activate("nodes") %>%
        mutate(original_name=vapply(seq_len(length(.data$original_name)),
            function(x){ 
                if(is.na(.data$original_name[x])) {
                    .data$name[x]  
                }  else {
                    .data$original_name[x]
                }
            },
            FUN.VALUE="character"))
}

#' process_reaction
#' 
#' process the kgml of global maps
#' e.g. in ko01100
#' 
#' Typically, `process_line` function is used to draw relationships 
#' as in the original KGML positions, however, the `coords` properties
#' is not considering the direction of reactions (substrate -> product),
#' thus if it is preferred, `process_reaction` is used to populate
#' new edges corresponding to `substrate -> product` and `product -> substrate`
#' if the reaction is reversible.
#' 
#' @param g graph
#' @param single_edge discard one edge when edge type is `reversible`
#' @param keep_no_reaction keep edges not related to reaction
#' @importFrom tidygraph bind_nodes bind_edges
#' @export
#' @return tbl_graph
#' @examples
#' gm_test <- create_test_pathway(line=TRUE)
#' test <- process_reaction(gm_test)
#' 
process_reaction <- function(g, single_edge=FALSE, keep_no_reaction=TRUE) {
    ## This is perhaps dirty ways to obtain edges. Perhaps directly
    ## parsing substrate -> product would be reasonable with
    ## assigning "reversible" and "irreversible"

    ## Obtain raw nodes
    nds <- g %>% activate("nodes") %>% data.frame()

    ## Obtain raw edges
    eds <- g %>% activate("edges") %>% data.frame()
    no_reacs <- eds[is.na(eds$reaction_id),]
    reacs <- eds$reaction_id %>% unique()
    reacs <- reacs[!is.na(reacs)]
    ## Prepare new edges
    
    new_eds <- lapply(reacs, function(reac_id) {
        konm <- nds[nds$orig.id %in% reac_id,]$name %>% unique()
        konm <- ifelse(is.null(konm), NA, konm)
        in_reacs <- eds[eds$reaction_id %in% reac_id, ]
        reac_name <- in_reacs$reaction %>% unique()
        row.names(in_reacs) <- seq_len(nrow(in_reacs))
        reac_type <- in_reacs$type %>% unique()
        
        subst_ind <- which(in_reacs$subtype_name == "substrate")
        prod_ind <- which(in_reacs$subtype_name == "product")

        eds <- lapply(subst_ind, function(subst) {
        	lapply(prod_ind, function(prod) {
        		fr <- in_reacs[subst, ]$from
        		to <- in_reacs[prod, ]$to
        		reac_info <- nds[in_reacs[subst, ]$to, ]
        		if (reac_type=="irreversible") {
        			return(c(fr, to, reac_type, reac_name,
        				konm, reac_info$bgcolor %>% unique(),
        				reac_info$fgcolor %>% unique()))
        		} else if (reac_type=="reversible") {
                    if (single_edge) {
                        return(rbind(
                            c(fr, to, reac_type,
                              reac_name, konm,
                              reac_info$bgcolor %>% unique(),
                              reac_info$fgcolor %>% unique())
                            ))                          
                    } else {                        
                        return(rbind(
                            c(fr, to, reac_type,
                              reac_name, konm,
                              reac_info$bgcolor %>% unique(),
                              reac_info$fgcolor %>% unique()),
                            c(to, fr, reac_type,
                              reac_name, konm,
                              reac_info$bgcolor %>% unique(),
                              reac_info$fgcolor %>% unique())
                            ))
                    }       			
        		} else {
                    stop("Unknown reaction type detected")
        		}
        	})
        })
        return(eds)
    })

    new_eds <- unlist(new_eds, recursive=FALSE)
    new_eds <- do.call(rbind, unlist(new_eds, recursive=FALSE)) %>%
        data.frame() %>%
        `colnames<-`(c("from","to","type","reaction",
            "name","bgcolor","fgcolor"))
    
    new_eds <- new_eds[!duplicated(new_eds),]
    new_eds$from <- as.integer(new_eds$from)
    new_eds$to <- as.integer(new_eds$to)
    if (keep_no_reaction) {
        if (dim(no_reacs)[1]!=0) {## If the no-reaction row is present
            all_columns <- union(colnames(no_reacs), colnames(new_eds))
            for (coln in all_columns) {
                if (!coln %in% colnames(new_eds)) {
                    new_eds[[coln]] <- NA
                }
                if (!coln %in% colnames(no_reacs)) {
                    no_reacs[[coln]] <- NA
                }

            }
            new_eds <- rbind(no_reacs, new_eds)            
        }
    }
    new_g <- tbl_graph(nodes=nds, edges=new_eds)
    new_g
}


#' get_reaction
#' 
#' Parse the reaction in KGML.
#' Used internally in pathway().
#' 
#' @noRd
#' @importFrom XML xmlAttrs getNodeSet xmlElementsByTagName
get_reaction <- function(xml) {
    rea_sets <- getNodeSet(xml, "//reaction")
    all_reas <- lapply(rea_sets, function(rea) {
        id <- xmlAttrs(rea)["id"]
        name <- xmlAttrs(rea)["name"]
        type <- xmlAttrs(rea)["type"]
        subs <- xmlElementsByTagName(rea,"substrate")
        prod <- xmlElementsByTagName(rea,"product")
        ## Looking for `alt` tag
        ## Multiple products or substrates are to be expected
        lapply(subs, function(ss) {
            lapply(prod, function(pp) {
                return(c(id, name, type,
                        xmlAttrs(ss)["id"], xmlAttrs(ss)["name"],
                        xmlAttrs(pp)["id"], xmlAttrs(pp)["name"]))          
            })
        })
    })
    all_reas <- unlist(all_reas, recursive=FALSE)
    all_reas <- do.call(rbind, unlist(all_reas, recursive=FALSE)) %>%
        data.frame() %>% 
        `colnames<-`(c("id","reac_name",
                    "type","substrate_id","substrate_name",
                    "product_id","product_name"))

    ## Perhaps this parsing would lead to wrong interpretation
    ## But for preserving Compound -> KO edges, this function
    ## adds edges of: 
    ##     substrate -> ID (KO) (type: type, reaction: reaction)
    ##     ID (KO) -> product (type: type, reaction: reaction)
    ## Later used in `process_reaction()`.
    ## Changed this layout to drop duplicates by distinct()
    
    rsp_rels <- lapply(seq_len(nrow(all_reas)), function(i) {
        lapply(unlist(strsplit(all_reas[i,"id"], " ")), function(j) {
            return(
                rbind(
                    c(all_reas[i,"substrate_id"], j, all_reas[i,"type"],
                        "substrate", NA, all_reas[i, "reac_name"],
                        all_reas[i, "id"]),
                    c(j, all_reas[i,"product_id"], all_reas[i,"type"],
                        "product", NA, all_reas[i, "reac_name"],
                        all_reas[i, "id"])
                    )
                )
        })
    })


    rsp_rels <- do.call(rbind, unlist(rsp_rels, recursive=FALSE)) %>%
        data.frame() %>%
        dplyr::distinct() %>%
        `colnames<-`(c("entry1","entry2","type",
            "subtype_name","subtype_value","reaction","reaction_id"))
    rsp_rels
}


#' pathway_info
#' 
#' obtain the list of pathway information
#' @param pid KEGG Pathway id
#' @param use_cache whether to use cache
#' @param directory directory of file
#' @return list of orthology and module contained in the pathway
#' @examples pathway_info("hsa04110")
#' @export
pathway_info <- function(pid, use_cache=FALSE, directory=NULL) {
    if (!is.null(directory)){
        dest <- paste0(directory, "/", pid)
    } else {
        dest <- pid
    }
    if (!file.exists(pid)) {
        if (use_cache) {
            bfc <- BiocFileCache()
            dest <- bfcrpath(bfc, paste0("https://rest.kegg.jp/get/",pid))  
        } else {
            download.file(paste0("https://rest.kegg.jp/get/",pid),
                destfile=dest)      
        }
    }

    con <- file(dest, "r")
    content_list <- list()
    while ( TRUE ) {
        line <- readLines(con, n=1)
        if ( length(line) == 0 ) {
            break
        }
        if (!startsWith(line, " ")) {
            current_id <- strsplit(line, " ") %>%
                vapply("[", 1, FUN.VALUE="character")
        }
        if (!current_id %in% c("REFERENCE","///")) {
            content <- substr(line, 13, nchar(line))
            content_list[[current_id]] <- c(content_list[[current_id]], content) 
        }
    }
    close(con)
    content_list$ENTRY <- strsplit(content_list$ENTRY, " ") %>%
        vapply("[", 1, FUN.VALUE="character")
    content_list
}


#' create_test_pathway
#' 
#' As downloading from KEGG API is not desirable
#' in vignettes or examples, return the `tbl_graph`
#' with two nodes and two edges.
#' @param line return example containing graphics type line
#' @examples create_test_pathway()
#' @export
#' @return tbl_graph
create_test_pathway <- function(line=FALSE) {

    if (line) {
        gm_test <- data.frame(name=c("cpd:C99998","cpd:C99999","ko:K00224"),
            type=c("compound","compound","ortholog"),
            graphics_type=c("circle","circle","line"),
            graphics_name=c("C99998","C99999","K00224"),
            coords=c(NA, NA, "1,2,3,4,5"),
            reaction=c(NA,NA,"rn:R99999"),
            orig.id=c(1,2,3),
            fgcolor=c("#ff0000","#ff0000","#ff0000"),
            bgcolor=c("#ffffff","#ffffff","#ffffff"))

        gm_test_edges <- rbind(
            data.frame(from=1,to=3,reaction="rn:R99999",
                subtype_name="substrate",
                type="irreversible",reaction_id="1"),
            data.frame(from=3,to=2,reaction="rn:R99999",
                subtype_name="product",
                type="irreversible", reaction_id="1"))
        gm_test <- tbl_graph(gm_test, gm_test_edges)
        return(gm_test)
    } else {
        ddx <- data.frame(
            name="hsa:51428",
            type="gene",
            reaction=NA,
            graphics_name="DDX41",
            x=500, y=-400,
            width=20,height=9,
            bgcolor="#BFFFBF",
            pathway_id="test"
        )
      
        trim <- data.frame(
            name="hsa:6737",
            type="gene",
            reaction=NA,
            graphics_name="TRIM21",
            x=560, y=-400,
            width=20,height=9,
            bgcolor="#BFFFBF",
            pathway_id="test"
        )
      
        nodes <- rbind(trim, ddx)
        nodes$xmin <- nodes$x-nodes$width/2
        nodes$ymin <- nodes$y-nodes$height/2
        nodes$xmax <- nodes$x+nodes$width/2
        nodes$ymax <- nodes$y+nodes$height/2
      
        edges <- rbind(c(from=1, to=2,
            subtype_name="degradation",pathway_id="test"),
                       c(from=1, to=2,
            subtype_name="ubiquitination",pathway_id="test")) %>%
            data.frame()
        edges$from <- as.integer(edges$from)
        edges$to <- as.integer(edges$to)
        tbl_graph(nodes, edges)        
    }
}

