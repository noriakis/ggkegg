#' ggkegg
#' 
#' main function parsing KEGG pathway data,
#' making igraph object and passing it to ggraph.
#' 
#' @param pid KEGG Pathway id e.g. hsa04110
#' @param layout default to "native", using KGML positions
#' @param return_igraph return the resulting igraph object
#' @param return_tidy_graph return the resulting tidygraph object
#' (override `return_igraph` argument)
#' @param delete_undefined delete the undefined nodes from graph
#' default to FALSE, which preserves nodes but 
#' add `undefined` attribute to graph
#' @param convert_org these organism names are fetched from REST API
#' and cached, and used to convert the KEGG identifiers.
#' e.g. c("hsa", "compound")
#' @param convert_first after converting, take the first element as 
#' node name when multiple genes are listed in the node
#' @param convert_collapse if not NULL, collapse the gene names by this character
#' when multiple genes are listed in the node.
#' @param convert_reaction reaction name (graph attribute `reaction`) 
#' will be converted to reaction formula
#' @param delete_undefined delete `undefined` node specifying group,
#' should be set to `TRUE` when the layout is not from native KGML.
#' @param delete_zero_degree delete nodes with zero degree, default to FALSE
#' @param module_type specify which module attributes to obtain
#' (definition or reaction)
#' @param module_definition_type `text` or `network` when parsing module definition.
#' If `text`, return ggplot object. If `network`, return `tbl_graph`.
#' @import igraph ggraph ggplot2 tidygraph
#' @export
ggkegg <- function(pid,
                   layout="native",
                   return_igraph=FALSE,
                   return_tbl_graph=FALSE,
                   pathway_number=1,
                   convert_org=NULL,
                   convert_first=TRUE,
                   convert_collapse=NULL,
                   convert_reaction=FALSE,
                   delete_undefined=FALSE,
                   delete_zero_degree=FALSE,
                   numeric_attribute=NULL,
                   node_rect_nudge=0,
                   group_rect_nudge=2,
                   module_type="definition",
                   module_definition_type="text") {
  enrich_attribute <- NULL
  if (!is.character(pid)) {
    if (attributes(pid)$class=="enrichResult") {
      org <- attributes(pid)$organism
      res <- attributes(pid)$result
      if (org!="UNKNOWN") {
        enrich_attribute <- paste0(org, ":", unlist(strsplit(
          res[pathway_number,]$geneID, "/")))
      } else {
        enrich_attribute <- unlist(strsplit(
          res[pathway_number,]$geneID, "/"))     
      }
      pid <- res[pathway_number,]$ID
    }
  }

  if (is.character(pid)) {
    if (startsWith(pid, "M")) {
      mod <- obtain_module(pid)
      def <- parse_module(mod, module_type)
      if (module_type=="definition") {
        if (module_definition_type=="text") {
          plot_list <- module_text(def, candidate_ko = enrich_attribute)
          return(plot_module_text(plot_list))
        } else if (module_definition_type=="network") {
          return(obtain_sequential_module_definition(def))
        } else {
          stop("Please specify `network` or `text` to module_definition_type")
        }
      } else if (module_type=="reaction") {
        return(def)
      } else {
        stop("Please specify `reaction` or `definition` to module_type")
      }
    }
  }

  g <- parse_kgml(pid=pid, convert_org=convert_org,
                  convert_first=convert_first, convert_collapse=convert_collapse,
                  node_rect_nudge=node_rect_nudge, group_rect_nudge=group_rect_nudge)

  if (!is.null(numeric_attribute)){
    V(g)$numeric_attribute <- numeric_attribute[V(g)$name]
  }

  if (!is.null(enrich_attribute)) {
    V(g)$enrich_attribute <- V(g)$name %in% enrich_attribute
  }

  if (delete_undefined) {
    g <- induced.subgraph(g, !V(g)$name %in% "undefined")
  } else {
    V(g)$undefined <- V(g)$name %in% "undefined"
  }
  if (delete_zero_degree) {
    g <- induced.subgraph(g, degree(g)!=0)
  }

  if (convert_reaction) {
    convert_vec <- obtain_map_and_cache("reaction",NULL)
    V(g)$converted_reaction <- unlist(lapply(V(g)$reaction,
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
  if (return_tbl_graph) {
    return(as_tbl_graph(g))
  }
  if (return_igraph) {
    return(g)
  }
  if (layout=="native") {
    ggraph(g, x=x, y=y)
  } else {
    g <- delete_vertex_attr(g, "x")
    g <- delete_vertex_attr(g, "y")
    ggraph(g, layout=layout)
  }
}