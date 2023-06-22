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
#' @importFrom XML xmlParse
#' @importFrom tibble as_tibble
#' @importFrom utils download.file head tail
#' @examples \donttest{pathway("hsa04110")}
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
  ## Specification of KGML format
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
  all_nodes <- list()
  grs <- list()
  rev_grs <- list()

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
    mult_coords <- NULL
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
            c(grs[[as.character(id)]], as.character(xmlAttrs(comp)["id"]))
          rev_grs[[as.character(xmlAttrs(comp)["id"])]] <- 
            c(rev_grs[[as.character(xmlAttrs(comp)["id"])]], as.character(id))
        }
      }
    }
    all_nodes[[ni]] <- c(id, name, type, reac,
                        glname, x, y, w, h, fg, bg, gltype,
                        paste0(mult_coords, collapse="|")) |>
                        setNames(node_names)
    ni <- ni + 1
  }
  kegg_nodes <- dplyr::bind_rows(all_nodes) |> data.frame() |>
    `colnames<-`(node_names)

  kegg_nodes$x <- as.numeric(kegg_nodes$x)
  kegg_nodes$y <- as.numeric(kegg_nodes$y)
  kegg_nodes$width <- as.numeric(kegg_nodes$width)
  kegg_nodes$height <- as.numeric(kegg_nodes$height)

  kegg_nodes$xmin <- kegg_nodes$x-kegg_nodes$width/2-node_rect_nudge
  kegg_nodes$xmax <- kegg_nodes$x+kegg_nodes$width/2+node_rect_nudge
  kegg_nodes$ymin <- kegg_nodes$y-kegg_nodes$height/2-node_rect_nudge
  kegg_nodes$ymax <- kegg_nodes$y+kegg_nodes$height/2+node_rect_nudge
  
  kegg_nodes[kegg_nodes$type=="group",]$xmin <- kegg_nodes[kegg_nodes$type=="group",]$xmin-group_rect_nudge
  kegg_nodes[kegg_nodes$type=="group",]$ymin <- kegg_nodes[kegg_nodes$type=="group",]$ymin-group_rect_nudge
  kegg_nodes[kegg_nodes$type=="group",]$xmax <- kegg_nodes[kegg_nodes$type=="group",]$xmax+group_rect_nudge
  kegg_nodes[kegg_nodes$type=="group",]$ymax <- kegg_nodes[kegg_nodes$type=="group",]$ymax+group_rect_nudge
  
  # kegg_nodes$xmin <- kegg_nodes$x-rect_width#kegg_nodes$width
  # kegg_nodes$xmax <- kegg_nodes$x+rect_width#kegg_nodes$width
  # kegg_nodes$ymin <- kegg_nodes$y-rect_height#kegg_nodes$height
  # kegg_nodes$ymax <- kegg_nodes$y+rect_height#kegg_nodes$height
  kegg_nodes$orig.id <- kegg_nodes$id ## Store ID as orig.id
  
  rel_sets <- getNodeSet(xml, "//relation")
  all_rels <- NULL
  ei <- 1
  for (rel in rel_sets) {
    entry1 <- xmlAttrs(rel)["entry1"]
    entry2 <- xmlAttrs(rel)["entry2"]
    rel_type <- xmlAttrs(rel)["type"]
    # rel_subtype <- xmlAttrs(rel[["subtype"]])["name"]
    rel_subtypes <- xmlElementsByTagName(rel,"subtype")
    for (rs in rel_subtypes) {
      all_rels[[ei]] <- c(entry1, entry2, rel_type,
        xmlAttrs(rs)["name"], xmlAttrs(rs)["value"]) |>
      setNames(c("entry1","entry2","type","subtype_name","subtype_value"))
      ei <- ei + 1
    }
  }
  if (!is.null(all_rels)) {
    kegg_edges <- dplyr::bind_rows(all_rels) |> data.frame() |>
      `colnames<-`(c("entry1","entry2","type","subtype_name","subtype_value"))
  } else {
    kegg_edges <- NULL
  }

  gr_rels <- NULL
  for (gr_name in names(grs)) {
    for (comp_name in grs[[gr_name]]) {
      ## Pad other values by `in_group`
      gr_rels <- rbind(gr_rels, 
        c(gr_name, comp_name, "in_group", "in_group", "in_group"))
    }
  }
  ## Include grouping edge
  # for (i in seq_len(nrow(kegg_edges))) {
  #   if (kegg_edges[i,"entry1"] %in% names(grs)) {
  #     for (tmp_node in grs[[kegg_edges[i,"entry1"]]]) {
  #       kegg_edges <- rbind(kegg_edges, c(tmp_node,
  #                                         kegg_edges[i,"entry2"],
  #                                         kegg_edges[i,"type"],
  #                                         kegg_edges[i,"subtype"]))
  #     }
  #   }
  #   if (kegg_edges[i,"entry2"] %in% names(grs)) {
  #     for (tmp_node in grs[[kegg_edges[i,"entry2"]]]) {
  #       kegg_edges <- rbind(kegg_edges, c(kegg_edges[i,"entry1"],
  #                                         tmp_node,
  #                                         kegg_edges[i,"type"],
  #                                         kegg_edges[i,"subtype"]))
  #     }
  #   }
  # }
  if (length(getNodeSet(xml, "//reaction"))!=0) {
    kegg_reac <- get_reaction(xml)
    if (!is.null(kegg_edges)) {kegg_edges$reaction <- NA}
    kegg_edges <- rbind(kegg_edges, kegg_reac)
  }

  ## Append grouping
  if (!is.null(kegg_edges)) {
    if (!is.null(gr_rels)) {
      gr_rels <- gr_rels |> data.frame() |> `colnames<-`(c("entry1","entry2","type",
        "subtype_name","subtype_value"))
      if ("reaction" %in% colnames(kegg_edges)) {
        gr_rels$reaction <- "in_group"
      }
      kegg_edges <- rbind(kegg_edges, gr_rels)
    }
  }

  if (!is.null(kegg_edges)) {
    g <- graph_from_data_frame(kegg_edges, vertices = kegg_nodes)
  } else {
    # message("No edges in the map, returning node data")
    g <- tbl_graph(nodes=kegg_nodes)
    # return(kegg_nodes)
  }

  ## Assign grouping
  # group <- NULL
  # for (i in V(g)$orig.id) {
  #   if (i %in% names(rev_grs)) {
  #     group <- c(group, paste(as.character(rev_grs[[i]]), collapse=","))
  #   } else {
  #     group <- c(group, NA)
  #   }
  # }
  # V(g)$group <- unlist(group)


  if (add_pathway_id) {
    V(g)$pathway_id <- pid
    E(g)$pathway_id <- pid
  }
  if (return_tbl_graph) {
    return(as_tbl_graph(g))
  } else {
    return(g)
    # return(as.igraph(g))
  }
}
parse_kgml <- pathway

#' process_line
#' 
#' process the KGML containing graphics type of `line`, like
#' global maps e.g. ko01100. Recursively add nodes and edges 
#' connecting them based on `coords` properties in KGML.
#' 
#' @param g graph
#' @param invert_y whether to invert the position, default to TRUE
#' should match with `pathway` function
#' @param verbose show progress
#' @importFrom tidygraph bind_nodes bind_edges
#' @export
#' @return tbl_graph
#' @examples 
#' ## For those containing nodes with the graphic type of `line`
#' gm_test <- data.frame(name="ko:K00112",type="ortholog",reaction="rn:R00112",
#'            graphics_name="K00112",fgcolor="#ff0000",bgcolor="#ffffff",
#'            graphics_type="line",coords="1,2,3,4",orig.id=1,pathway_id="test")
#' gm_test <- tbl_graph(gm_test)
#' test <- process_line(gm_test)
process_line <- function(g, invert_y=TRUE, verbose=FALSE) {
  ## [TODO] speed up
  ## [TODO] add verbose
  ## [TODO] add positional argument to coords to show arrow
  ## We cannot do this, as coords are not ordered to show direction
  df <- as_tbl_graph(g)

  cos <- list()
  eds <- list()
  j <- 1
  e <- 1
  name_col_node <- c("name","x","y","type","original_name")
  name_col_edge <- c("from","to","type","name",
    "bgcolor","fgcolor","reaction","orig.id")

  for (i in seq_along(V(g)$name)) {

    if (V(g)$graphics_type[i]=="line") {

      raw_name <- V(g)$name[i]
      bgcol <- V(g)$bgcolor[i]
      fgcol <- V(g)$fgcolor[i]
      reac <- V(g)$reaction[i]
      origid <- V(g)$orig.id[i]
      rawco <- V(g)$coords[i]

      ## reversible or irreversible
      # rev <- E(g)[E(g)$reaction==reac]$type |> unique()
      # rev <- ifelse(identical(rev,character(0)),NA,rev)

      if (grepl("\\|",rawco)) {
        rawcos <- strsplit(rawco, "\\|") |> unlist()
      } else {
        rawcos <- rawco
      }

      for (rc in rawcos) {
        co <- unlist(strsplit(rc, ","))
        q <- 1
        if (verbose) {
          GetoptLong::qqcat("@{raw_name}, cooridnates components: @{length(co)}\n")
        }
        for (h in seq_len(length(co))) {
          if (is.na(co[q+2])) {
            # eds[[e-1]]["position"] <- "End"
            break
          }
          if (verbose) {
            GetoptLong::qqcat("  @{paste(co[q], co[q+1])} -> @{paste(co[q+2], co[q+3])}\n")
          }
          cos[[j]] <- c(paste0(raw_name,"_",j), co[q], co[q+1], "line",raw_name) |>
          setNames(name_col_node)
          cos[[j+1]] <- c(paste0(raw_name,"_",j+1), co[q+2], co[q+3], "line",raw_name)|>
          setNames(name_col_node)
          eds[[e]] <- c(paste0(raw_name,"_",j), paste0(raw_name,"_",j+1),
                              "line",raw_name,bgcol,fgcol,reac,origid
                              )|>
          setNames(name_col_edge)
          # if (h==1) {
          #   eds[[e]]["position"] <- "Start"
          # }
          e <- e+1
          j <- j+2
          q <- q+2
        }
      }
    }
  }

  cos <- dplyr::bind_rows(cos) |> data.frame() |> 
    `colnames<-`(name_col_node)
  cos$x <- as.numeric(cos$x);
  if (invert_y) {
    cos$y <- -1 * as.numeric(cos$y)
  } else {
    cos$y <- as.numeric(cos$y)
  }
  eds <- dplyr::bind_rows(eds) |> data.frame() |> 
    `colnames<-`(name_col_edge)

  df_add <- df |> bind_nodes(cos) |> bind_edges(eds)
  df_add |> activate(nodes) |>
    mutate(original_name=vapply(seq_len(length(original_name)),
     function(x){ if(is.na(original_name[x])) name[x] else original_name[x]},
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
#' @importFrom tidygraph bind_nodes bind_edges
#' @export
#' @return tbl_graph
#' @examples
#' gm_test <- rbind(data.frame(name="cpd:C99998",type="compound",
#'            graphics_name="C99998",fgcolor="#ff0000",bgcolor="#ffffff"),
#'            data.frame(name="cpd:C99999",type="compound",
#'                       graphics_name="C99999",fgcolor="#ff0000",bgcolor="#ffffff"),
#'            data.frame(name="ko:K00224",type="ortholog",
#'                       graphics_name="K00224",fgcolor="#ff0000",bgcolor="#ffffff")
#'            )
#' gm_test_edges <- rbind(data.frame(from=1,to=3,reaction="rn:R99999",subtype_name="substrate",
#'                             type="irreversible"),
#'                        data.frame(from=3,to=2,reaction="rn:R99999",subtype_name="product",
#'                                   type="irreversible"))
#' gm_test <- tbl_graph(gm_test, gm_test_edges)
#' test <- process_reaction(gm_test)
#' 
process_reaction <- function(g) {

  ## Obtain raw nodes
  nds <- g |> activate(nodes) |> data.frame()

  ## Obtain raw edges
  eds <- g |> activate(edges) |> data.frame()

  ## Prepare new edges
  new_eds <- NULL
  k <- 1
  for (i in eds$reaction |> unique()) {
      konm <- nds[nds$reaction %in% i,]$name
      konm <- ifelse(is.null(konm),NA,konm)
      tmp <- eds[eds$reaction %in% i,]
      if (tmp$type |> unique()=="irreversible") {
          fs <- tmp[tmp$subtype_name=="substrate",]$from
          tos <- tmp[tmp$subtype_name=="product",]$to
          for (cfs in fs) {
              for (ctos in tos) {
                  new_eds[[k]] <- c(cfs, ctos, "irreversible",
                    tmp$reaction |> unique(), konm)
                  k <- k + 1
              }
          }
      } else {
          fs <- tmp[tmp$subtype_name=="substrate",]$from
          tos <- tmp[tmp$subtype_name=="product",]$to
          for (cfs in fs) {
              for (ctos in tos) {
                  new_eds[[k]] <- c(cfs, ctos, "reversible",
                    tmp$reaction |> unique(), konm)
                                    # tmp$pathway_id |> unique(), tmp$name |> unique(),
                                    # tmp$bgcolor |> unique(),
                                    # tmp$fgcolor |> unique(),
                                    # tmp$orig.id |> unique())
                  k <- k + 1
              }
          } 
          for (ctos in tos) {
              for (cfs in fs) {
                  new_eds[[k]] <- c(ctos, cfs, "reversible",
                    tmp$reaction |> unique(), konm)
                  k <- k + 1
              }
          }    
      }
  }

  new_eds <- do.call(rbind, new_eds) |> data.frame() |>
  `colnames<-`(c("from","to","type","reaction","name"))

  new_eds <- new_eds[!duplicated(new_eds),]
  new_eds$from <- as.integer(new_eds$from)
  new_eds$to <- as.integer(new_eds$to)

  new_g <- tbl_graph(nodes=nds, edges = new_eds)
}


#' get_reaction
#' parse the reaction in KGML
#' @noRd
#' @importFrom XML xmlAttrs getNodeSet xmlElementsByTagName
get_reaction <- function(xml) {
  rea_sets <- getNodeSet(xml, "//reaction")
  all_reas <- NULL
  for (rea in rea_sets) {
    id <- xmlAttrs(rea)["id"]
    name <- xmlAttrs(rea)["name"]
    type <- xmlAttrs(rea)["type"]
    subs <- xmlElementsByTagName(rea,"substrate")
    prod <- xmlElementsByTagName(rea,"product")
    ## Looking for `alt` tag
    ## Multiple products or substrates are to be expected
    for (ss in subs) {
      for (pp in prod) {
        all_reas <- rbind(all_reas, c(id, name, type,
                                      xmlAttrs(ss)["id"], xmlAttrs(ss)["name"],
                                      xmlAttrs(pp)["id"], xmlAttrs(pp)["name"]))
      }
    }
  }
  all_reas <- all_reas |> data.frame() |> `colnames<-`(c("id","reac_name",
                                             "type","substrate_id","substrate_name",
                                             "product_id","product_name"))

  ## Maybe this parsing will lead to wrong interpretation
  ## substrate -> ID (KO) (type: type, reaction: reaction)
  ## ID (KO) -> product (type: type, reaction: reaction)
  rsp_rels <- NULL
  for (i in seq_len(nrow(all_reas))) {
    for (j in unlist(strsplit(all_reas[i,"id"], " "))) {
      rsp_rels <- rbind(rsp_rels,
      c(all_reas[i,"substrate_id"], j, all_reas[i,"type"], "substrate", NA, all_reas[i, "reac_name"]),
      c(j, all_reas[i,"product_id"], all_reas[i,"type"], "product", NA, all_reas[i, "reac_name"]))
    }
  }

  ## substrate -> product, type: type, subtype: ID (KO), reaction: reaction
  # rsp_rels <- NULL
  # for (i in seq_len(nrow(all_reas))) {
  #   for (j in unlist(strsplit(all_reas[i,"id"], " "))) {
  #     rsp_rels <- rbind(rsp_rels,
  #     c(all_reas[i,"substrate_id"], all_reas[i,"product_id"], all_reas[i,"type"], j, all_reas[i, "reac_name"]))
  #   }
  # }


  rsp_rels <- data.frame(rsp_rels) |> 
    `colnames<-`(c("entry1","entry2","type","subtype_name","subtype_value","reaction"))
  rsp_rels
}


#' pathway_info
#' 
#' obtain the list of pathway information
#' @param pid KEGG Pathway id
#' @param use_cache whether to use cache
#' @param directory directory of file
#' @return list of orthology and module contained in the pathway
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
  pway <- list()
  con = file(dest, "r")
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
  content_list$ENTRY <- strsplit(content_list$ENTRY, " ") |> sapply("[", 1)
  content_list
}



#' return_pathway_example
#' 
#' As downloading from KEGG API is not desirable
#' in vignettes or examples, return the `tbl_graph`
#' with two nodes and two edges.
#' @examples return_pathway_example
#' @export
#' @return tbl_graph
return_pathway_example <- function() {
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
  nodes$xmin = nodes$x-nodes$width/2
  nodes$ymin = nodes$y-nodes$height/2
  nodes$xmax = nodes$x+nodes$width/2
  nodes$ymax = nodes$y+nodes$height/2
  
  edges <- rbind(c(from=1, to=2, subtype_name="degradation",pathway_id="test"),
                 c(from=1, to=2, subtype_name="ubiquitination",pathway_id="test")) |>
    data.frame()
  edges$from <- as.integer(edges$from)
  edges$to <- as.integer(edges$to)
  tbl_graph(nodes, edges)
}

