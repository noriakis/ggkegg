#' pathway
#' 
#' KEGG pathway parsing function
#'
#' @param pid pathway id
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
           group_rect_nudge=2,
           node_rect_nudge=0,
           invert_y=TRUE,
           add_pathway_id=TRUE,
           return_tbl_graph=TRUE,
           return_image=FALSE) {
  ## Specification of KGML format
  ## https://www.genome.jp/kegg/xml/docs/

  file_name <- paste0(pid,".xml")
  if (!file.exists(file_name)) {
    download.file(url=paste0("https://rest.kegg.jp/get/",pid,"/kgml"),
                  destfile=file_name)
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
#' process the kgml containing graphics type of `line`
#' e.g. in ko01100
#' @param g graph
#' @param invert_y whether to invert the position, default to TRUE
#' should match with `pathway` function
#' @param verbose show progress
#' @importFrom tidygraph bind_nodes bind_edges
#' @export
#' @return tbl_graph
#' @examples 
#' ## For those containing nodes with the graphic type of `line`
#' \donttest{g <- pathway("ko01100") |> process_line()}
process_line <- function(g, invert_y=TRUE, verbose=FALSE) {
  ## [TODO] speed up
  ## [TODO] add verbose
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
          if (is.na(co[q+2])) {break}
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
#' @return list of orthology and module contained in the pathway
#' @export
pathway_info <- function(pid) {
  if (!file.exists(pid)) {
    download.file(paste0("https://rest.kegg.jp/get/",pid),
                  destfile=pid)
  }
  pway <- list()
  con = file(pid, "r")
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