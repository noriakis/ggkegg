#' overlay_raw_map
#' 
#' Overlay the raw KEGG pathway image on ggraph
#' 
#' @param pid pathway ID
#' @param transparent_colors make these colors transparent to overlay
#' Typical colors would be: "#CCCCCC", "#FFFFFF","#BFBFFF","#BFFFBF"
#' @import magick
#' @export
overlay_raw_map <- function(pid,
                            transparent_colors=c("#FFFFFF","#BFBFFF","#BFFFBF")) {
  structure(list(pid=pid,
                 transparent_colors=transparent_colors),
            class = "overlay_raw_map")
}

#' ggplot_add.overlay_raw_map
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.overlay_raw_map
#' @export
ggplot_add.overlay_raw_map <- function(object, plot, object_name) {
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
  plot + 
    annotation_raster(ras, xmin=0, ymin=0, xmax=w, ymax=-1*h,
                    interpolate = TRUE)+
    coord_fixed(xlim = c(0,w), ylim=c(-1*h,0))
}
