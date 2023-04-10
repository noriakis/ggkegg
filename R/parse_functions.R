#' @param rect_width manual specification of rectangle
#' @param rect_height manual specification of rectangle
#' @export
parse_kgml <- function(file_name,
                       pid=NULL,
                       rect_width=25,
                       rect_height=10,
                       convert_org=NULL,
                       convert_collapse=NULL,
                       convert_first=TRUE,
                       group_rect_nudge=2,
                       node_rect_nudge=0,
                       invert_y=TRUE) {
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
    x <- as.numeric(xmlAttrs(gl)["x"])
    if (invert_y) {
      y <- -1*as.numeric(xmlAttrs(gl)["y"])
    } else {
      y <- as.numeric(xmlAttrs(gl)["y"])
    }
    w <- as.numeric(xmlAttrs(gl)["width"])
    h <- as.numeric(xmlAttrs(gl)["height"])
    if (type=="group") {
      for (comp in xmlElementsByTagName(node,"component")) {
        grs[[as.character(id)]] <- 
          c(grs[[as.character(id)]], as.character(xmlAttrs(comp)["id"]))
        rev_grs[[as.character(xmlAttrs(comp)["id"])]] <- 
          c(rev_grs[[as.character(xmlAttrs(comp)["id"])]], as.character(id))
      }
    }
    all_nodes <- rbind(all_nodes, c(id, name, type, reac,
                                    glname, x, y, w, h))
  }
  kegg_nodes <- all_nodes |> data.frame() |>
    `colnames<-`(c("id","name","type","reaction",
                   "graphics_name",
                   "x","y","width","height"))

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
  kegg_edges <- all_rels |> data.frame() |>
    `colnames<-`(c("entry1","entry2","type","subtype"))
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
  g <- graph_from_data_frame(kegg_edges, vertices = kegg_nodes)

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
  return(g)
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

#' @import BiocFileCache
#' @importFrom stringr str_extract
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
    convert_vec <- vapply(convert$V2, function(x) {
      vapply(unlist(strsplit(x, ";"))[1],
             function(x) unlist(strsplit(x,","))[1],
             FUN.VALUE="character")
      
    }, FUN.VALUE="character")
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
