#' @param rect_width manual specification of rectangle
#' @param rect_height manual specification of rectangle
#' @param add_pathway_id add pathway id to graph
#' @param return_tbl_graph return tbl_graph object
#' @export
parse_kgml <- function(pid,
                       rect_width=25,
                       rect_height=10,
                       convert_org=NULL,
                       convert_collapse=NULL,
                       convert_first=TRUE,
                       group_rect_nudge=2,
                       node_rect_nudge=0,
                       invert_y=TRUE,
                       add_pathway_id=TRUE,
                       return_tbl_graph=TRUE) {
  ## Specification of KGML format
  ## https://www.genome.jp/kegg/xml/docs/

  file_name <- paste0(pid,".xml")
  if (!file.exists(file_name)) {
    download.file(url=paste0("https://rest.kegg.jp/get/",pid,"/kgml"),
                  destfile=file_name)
  }

  xml <- xmlParse(file_name)
  node_sets <- getNodeSet(xml, "//entry")
  all_nodes <- NULL
  grs <- list()
  rev_grs <- list()
  for (node in node_sets) {
    id <- xmlAttrs(node)["id"]
    name <- xmlAttrs(node)["name"]
    type <- xmlAttrs(node)["type"]
    reac <- xmlAttrs(node)["reaction"]
    gl <- node[["graphics"]]
    glname <- xmlAttrs(gl)["name"]
    gltype <- xmlAttrs(gl)["type"]
    glcoords <- xmlAttrs(gl)["coords"]
    x <- as.numeric(xmlAttrs(gl)["x"])
    if (invert_y) {
      y <- -1*as.numeric(xmlAttrs(gl)["y"])
    } else {
      y <- as.numeric(xmlAttrs(gl)["y"])
    }
    w <- as.numeric(xmlAttrs(gl)["width"])
    h <- as.numeric(xmlAttrs(gl)["height"])
    bg <- xmlAttrs(gl)["bgcolor"]
    if (type=="group") {
      for (comp in xmlElementsByTagName(node,"component")) {
        grs[[as.character(id)]] <- 
          c(grs[[as.character(id)]], as.character(xmlAttrs(comp)["id"]))
        rev_grs[[as.character(xmlAttrs(comp)["id"])]] <- 
          c(rev_grs[[as.character(xmlAttrs(comp)["id"])]], as.character(id))
      }
    }
    all_nodes <- rbind(all_nodes, c(id, name, type, reac,
                                    glname, x, y, w, h, bg, gltype, glcoords))
  }
  kegg_nodes <- all_nodes |> data.frame() |>
    `colnames<-`(c("id","name","type","reaction",
                   "graphics_name",
                   "x","y","width","height","bgcolor","graphics_type","coords"))

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
  for (rel in rel_sets) {
    entry1 <- xmlAttrs(rel)["entry1"]
    entry2 <- xmlAttrs(rel)["entry2"]
    rel_type <- xmlAttrs(rel)["type"]
    rel_subtype <- xmlAttrs(rel[["subtype"]])["name"]
    all_rels <- rbind(all_rels, c(entry1, entry2, rel_type, rel_subtype))
  }
  if (!is.null(all_rels)) {
    kegg_edges <- all_rels |> data.frame() |>
      `colnames<-`(c("entry1","entry2","type","subtype"))
  } else {
    kegg_edges <- NULL
  }
  ## Include grouping
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
    kegg_edges <- rbind(kegg_edges, kegg_reac)
  }
  if (!is.null(kegg_edges)) {
    g <- graph_from_data_frame(kegg_edges, vertices = kegg_nodes)
  } else {
    message("No edges in the map, returning node data")
    return(kegg_nodes)
  }

  ## Assign grouping
  group <- NULL
  for (i in V(g)$orig.id) {
    if (i %in% names(rev_grs)) {
      group <- c(group, paste(as.character(rev_grs[[i]]), collapse=","))
    } else {
      group <- c(group, NA)
    }
  }
  V(g)$group <- unlist(group)

  convert_vec <- NULL
  if (!is.null(convert_org)) {
    for (co in convert_org) {
      convert_vec <- c(convert_vec,
                       obtain_map_and_cache(co, pid))
    }
    V(g)$converted_name <- unlist(lapply(V(g)$name,
                                         function(x) {
                                           inc_genes <- unlist(strsplit(x, " "))
                                           conv_genes <- NULL
                                           for (inc in inc_genes) {
                                             convs <- convert_vec[inc]
                                             if (is.na(convs)) {
                                               conv_genes <- c(conv_genes, x)
                                             } else {
                                               conv_genes <- c(conv_genes, convs)
                                             }
                                           }
                                           if (convert_first) {
                                             conv_genes[1]
                                           } else {
                                             paste(conv_genes, collapse=convert_collapse)
                                           }
                                         }))
  }
  if (add_pathway_id) {
    V(g)$pathway_id <- pid
  }
  if (return_tbl_graph) {
    as_tbl_graph(g)
  } else {
    return(g)
  }
}


#' process_line
#' process the kgml containing graphics type of `line`
#' e.g. in ko01100
#' @export
process_line <- function(g, invert_y=TRUE) {
  df <- as_tbl_graph(g)

  cos <- NULL
  eds <- NULL

  for (i in seq_along(V(g)$name)) {
    j <- 0
    if (V(g)$graphics_type[i]=="line") {
      co <- unlist(strsplit(V(g)$coords[i], ","))
      for (q in seq(1,length(co),4)) {
        cos <- rbind(cos, c(paste0(V(g)$name[i],"_",j), co[q], co[q+1], "line",V(g)$name[i]))
        cos <- rbind(cos, c(paste0(V(g)$name[i],"_",j+1), co[q+2], co[q+3], "line",V(g)$name[i]))
        eds <- rbind(eds, c(paste0(V(g)$name[i],"_",j), paste0(V(g)$name[i],"_",j+1),
                            "line",V(g)$name[i]))
        j <- j +2
      }
    }
  }
  cos <- cos |> data.frame() |> `colnames<-`(c("name","x","y","type","original_name"))
  cos$x <- as.numeric(cos$x);
  if (invert_y) {
    cos$y <- -1 * as.numeric(cos$y)
  } else {
    cos$y <- as.numeric(cos$y)
  }
  eds <- eds |> data.frame() |> `colnames<-`(c("from","to","type","name"))
  df_add <- df |> bind_nodes(cos) |> bind_edges(eds)
  df_add |> activate(nodes) |>
  mutate(original_name=vapply(1:length(original_name),
   function(x){ if(is.na(original_name[x])) name[x] else original_name[x]},
   FUN.VALUE="character"))
}


#' @noRd
get_reaction <- function(xml) {
  rea_sets <- getNodeSet(xml, "//reaction")
  all_reas <- NULL
  for (rea in rea_sets) {
    id <- xmlAttrs(rea)["id"]
    name <- xmlAttrs(rea)["name"]
    type <- xmlAttrs(rea)["type"]
    subs <- xmlAttrs(xmlElementsByTagName(rea,"substrate")[[1]],"substrate")
    prod <- xmlAttrs(xmlElementsByTagName(rea,"product")[[1]],"product")
    ## Looking for `alt` tag
    all_reas <- rbind(all_reas, c(id, name, type,
                                  subs["id"], subs["name"],
                                  prod["id"], prod["name"]))
  }
  all_reas <- all_reas |> data.frame() |> `colnames<-`(c("id","reac_name",
                                             "type","substrate_id","substrate_name",
                                             "product_id","product_name"))
  rsp_rels <- NULL
  for (i in seq_len(nrow(all_reas))) {
    for (j in unlist(strsplit(all_reas[i,"id"], " "))) {
      rsp_rels <- rbind(rsp_rels,
      c(j, all_reas[i,"substrate_id"], all_reas[i,"type"], "substrate"),
      c(j, all_reas[i,"product_id"], all_reas[i,"type"], "product"))
    }
  }
  rsp_rels <- data.frame(rsp_rels) |> 
    `colnames<-`(c("entry1","entry2","type","subtype"))
  rsp_rels
}