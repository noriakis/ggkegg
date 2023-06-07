#' overlay_raw_map
#' 
#' Overlay the raw KEGG pathway image on ggraph
#' 
#' @param pid pathway ID
#' @param transparent_colors make these colors transparent to overlay
#' Typical choice of colors would be: "#CCCCCC", "#FFFFFF","#BFBFFF","#BFFFBF", "#7F7F7F", "#808080",
#' "#ADADAD","#838383"
#' @param clip clip the both end of x- and y-axis by one dot
#' @param adjust adjust the x-axis location by 0.5 in data coordinates
#' @import magick
#' @return ggplot2 object
#' @export
#' @examples
#' nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
#'                    x=c(1,1),
#'                    xmin=c(-1,-1),
#'                    xmax=c(2,2),
#'                    y=c(1,1),
#'                    ymin=c(-1,-1),
#'                    ymax=c(2,2),
#'                    pathway_id="hsa04110")
#' edges <- data.frame(from=1, to=2, name="K00112")
#' graph <- tbl_graph(nodes, edges)
#' \donttest{ggraph(graph) + overlay_raw_map()}
#'
overlay_raw_map <- function(pid=NULL,
                            transparent_colors=c("#FFFFFF",
                              "#BFBFFF","#BFFFBF","#7F7F7F",
                              "#808080"),
                            adjust=TRUE, clip=FALSE) {
  structure(list(pid=pid,
                 transparent_colors=transparent_colors,
                 adjust=adjust,
                 clip=clip),
            class = "overlay_raw_map")
}

#' ggplot_add.overlay_raw_map
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.overlay_raw_map
#' @return ggplot2 object
#' @importFrom grDevices as.raster
#' @export
#' @examples
#' nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
#'                    x=c(1,1),
#'                    xmin=c(-1,-1),
#'                    xmax=c(2,2),
#'                    y=c(1,1),
#'                    ymin=c(-1,-1),
#'                    ymax=c(2,2),
#'                    pathway_id="hsa04110")
#' edges <- data.frame(from=1, to=2, name="K00112")
#' graph <- tbl_graph(nodes, edges)
#' \donttest{ggraph(graph) + overlay_raw_map()}
#'
ggplot_add.overlay_raw_map <- function(object, plot, object_name) {
  if (is.null(object$pid)) {
    infer <- plot$data$pathway_id |> unique()
    object$pid <- infer[!is.na(infer)]
  }
  ## Return the image URL, download and cache
  url <- paste0(as.character(pathway(object$pid,
                                     return_image=TRUE)))
  bfc <- BiocFileCache()
  path <- bfcrpath(bfc, url)
  
  ## Load, transparent and rasterize
  magick_image <- image_read(path)
  img_info <- image_info(magick_image)
  w <- img_info$width
  h <- img_info$height
  
  for (col in object$transparent_colors) {
    magick_image <- magick_image |> 
      image_transparent(col)
  }
  
  ras <- as.raster(magick_image)


  xmin=0
  xmax=w
  ymin=-1*h
  ymax=0

  if (object$clip) {
    ras <- ras[seq_len(nrow(ras)-1),
    seq_len(ncol(ras)-1)]
  }
  if (object$adjust) {
    xmin <- xmin - 0.5
    xmax <- xmax - 0.5
  }
  plot + 
    annotation_raster(ras, xmin=xmin, ymin=ymin,
      xmax=xmax, ymax=ymax, interpolate = TRUE)+
    coord_fixed(xlim = c(xmin,xmax), ylim=c(ymin,ymax))
}
