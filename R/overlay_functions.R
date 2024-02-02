#' overlay_raw_map
#' 
#' Overlay the raw KEGG pathway image on ggraph
#' 
#' @param pid pathway ID
#' @param directory directory to store images if not use cache
#' @param transparent_colors make these colors transparent to overlay
#' Typical choice of colors would be:
#' "#CCCCCC", "#FFFFFF","#BFBFFF","#BFFFBF", "#7F7F7F", "#808080",
#' "#ADADAD","#838383","#B3B3B3"
#' @param clip clip the both end of x- and y-axis by one dot
#' @param adjust adjust the x-axis location by 0.5 in data coordinates
#' @param adjust_manual_x adjust the position manually for x-axis
#' Override `adjust`
#' @param adjust_manual_y adjust the position manually for y-axis
#' Override `adjust`
#' @param use_cache whether to use BiocFileCache()
#' @param interpolate parameter in annotation_raster()
#' @param high_res Use high resolution (2x) image for the overlay
#' @import magick
#' @return ggplot2 object
#' @export
#' @examples
#' ## Need `pathway_id` column in graph 
#' ## if the function is to automatically infer
#' graph <- create_test_pathway() |> mutate(pathway_id="hsa04110")
#' ggraph(graph) + overlay_raw_map()
#'
overlay_raw_map <- function(pid=NULL, directory=NULL,
                            transparent_colors=c("#FFFFFF",
                                "#BFBFFF","#BFFFBF","#7F7F7F",
                                "#808080"),
                            adjust=TRUE,
                            adjust_manual_x=NULL,
                            adjust_manual_y=NULL,
                            clip=FALSE,
                            use_cache=TRUE,
                            interpolate=TRUE,
                            high_res=FALSE) {
    structure(list(pid=pid,
                    transparent_colors=transparent_colors,
                    adjust=adjust,
                    clip=clip,
                    adjust_manual_x=adjust_manual_x,
                    adjust_manual_y=adjust_manual_y,
                    directory=directory,
                    use_cache=use_cache,
                    interpolate=interpolate,
                    high_res=high_res),
            class="overlay_raw_map")
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
#' ## Need `pathway_id` column in graph 
#' ## if the function is to automatically infer
#' graph <- create_test_pathway() |> mutate(pathway_id="hsa04110")
#' ggraph(graph) + overlay_raw_map()
#'
ggplot_add.overlay_raw_map <- function(object, plot, object_name) {
    if (is.null(object$pid)) {
        infer <- plot$data$pathway_id |> unique()
        object$pid <- infer[!is.na(infer)]
        if (object$high_res) {
        	## Convert to reference ID
        	cur_id <- object$pid
        	object$pid <- paste0("map",
        		regmatches(cur_id, gregexpr("[[:digit:]]+", cur_id)) %>% unlist())
        }
    }
    if (!grepl("[[:digit:]]", object$pid)) {
        warning("Looks like not KEGG ID for pathway")
        return(1)
    }
    ## Return the image URL, download and cache
    ## From 1.1.10
    url <- paste0("https://rest.kegg.jp/get/",object$pid,"/image")
    if (object$high_res) {
    	if (!startsWith(object$pid, "map")) {
    		stop("High resolution image can be obtained for the reference pathway.")
    	}
    	url <- paste0(url, "2x")
    }
    if (object$use_cache) {
        bfc <- BiocFileCache()
        path <- bfcrpath(bfc, url)    
    } else {
        path <- paste0(object$pid, ".png")
        if (!is.null(object$directory)) {
            path <- paste0(object$directory,"/",path)
        }
        download.file(url=url, destfile=path, mode='wb')
    }
  
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


    xmin <- 0
    xmax <- w
    ymin <- -1*h
    ymax <- 0

    if (object$clip) {
        ras <- ras[seq_len(nrow(ras)-1),
                    seq_len(ncol(ras)-1)]
    }
    if (!is.null(object$adjust_manual_x)) {
        object$adjust <- FALSE
        xmin <- xmin + object$adjust_manual_x
        xmax <- xmax + object$adjust_manual_x
    }
    if (!is.null(object$adjust_manual_y)) {
        object$adjust <- FALSE
        ymin <- ymin + object$adjust_manual_y
        ymax <- ymax + object$adjust_manual_y
    }
    if (object$adjust) {
        xmin <- xmin - 0.5
        xmax <- xmax - 0.5
        # ymin <- ymin - 0.5
        # ymax <- ymax - 0.5
    }
    p <- plot + 
        annotation_raster(ras, xmin=xmin, ymin=ymin,
            xmax=xmax, ymax=ymax, interpolate=object$interpolate)+
        coord_fixed(xlim=c(xmin,xmax), ylim=c(ymin,ymax))
    attr(p, "original_width") <- w
    attr(p, "original_height") <- h
    return(p)
}


#' ggkeggsave
#' @param filename file name of the image
#' @param plot plot to be saved
#' @param dpi dpi, passed to ggsave
#' @param wscale width scaling factor for pixel to inches
#' @param hscale height scaling factor fo pixel to inches
#' @return save the image
#' @export
ggkeggsave <- function(filename, plot, dpi=300, wscale=90, hscale=90) {
	ggsave(filename, plot, dpi=dpi, width=attr(plot, "original_width")/wscale,
		height=attr(plot, "original_height")/hscale, units="in")
}