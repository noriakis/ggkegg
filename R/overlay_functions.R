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
#' @param adjust adjust the x- and y-axis location by 0.5 in data coordinates
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
                            adjust=FALSE,
                            adjust_manual_x=NULL,
                            adjust_manual_y=NULL,
                            clip=FALSE,
                            use_cache=TRUE,
                            interpolate=FALSE,
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
    xmax <- w-1
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
        ymin <- ymin - 0.5
        ymax <- ymax - 0.5
    }
    p <- plot + 
        annotation_raster(ras, xmin=xmin, ymin=ymin,
            xmax=xmax, ymax=ymax, interpolate=object$interpolate)+ 
        scale_x_continuous(expand=c(0,0), limits=c(0,w-1)) +
        scale_y_continuous(expand=c(0,0), limits=c(-1*h+1, 0))
    attr(p, "original_width") <- w
    attr(p, "original_height") <- h
    return(p)
}


#' ggkeggsave
#' 
#' save the image respecting the original width and height of the image.
#' Only applicable for the ggplot object including `overlay_raw_map` layers.
#' 
#' @param filename file name of the image
#' @param plot plot to be saved
#' @param dpi dpi, passed to ggsave
#' @param wscale width scaling factor for pixel to inches
#' @param hscale height scaling factor fo pixel to inches
#' @return save the image
#' @export
#' 
ggkeggsave <- function(filename, plot, dpi=300, wscale=90, hscale=90) {
	ggsave(filename, plot, dpi=dpi, width=attr(plot, "original_width")/wscale,
		height=attr(plot, "original_height")/hscale, units="in")
}


#' output_overlay_image
#' 
#' The function first exports the image, combine it with the original image.
#' Note that if the legend is outside the pathway image, the result will not 
#' show it correctly. Place the legend inside the panel by adding the theme 
#' such as theme(legend.position=c(0.5, 0.5)).
#' 
#' If the legend must be placed outside the image, the users can set 
#' with_legend_image to TRUE. This will create another legend only image
#' and concatenate it with the pathway image. legend_space option can be 
#' specified to control the spacing for the legend. If need to append horizontal
#' legend, enable legend_horiz option.
#' 
#' By default, unlink option is enabled which means the function will delete
#' the intermediate files.
#' 
#' 
#' @param gg ggraph object
#' @param with_legend if legend (group-box) is in gtable, output them
#' @param use_cache use BiocFileCache for caching the image
#' @param high_res use 2x resolution image
#' @param res resolution parameter passed to saving the ggplot2 image
#' @param out output file name
#' @param directory specify if you have already downloaded the image
#' @param transparent_colors transparent colors
#' @param unlink unlink the intermediate image
#' @param with_legend_image append legend image instead of using gtable
#' @param legend_horiz append legend to the bottom of the image
#' @param legend_space legend spacing specification (in pixel)
#' @export
#' @importFrom grDevices dev.off png
#' @import gtable
#' @return output the image
#' @examples
#' \dontrun{
#'     ouput_overlay_image(ggraph(pathway("hsa04110")))
#' } 
#' 
#' 
output_overlay_image <- function(gg, with_legend=TRUE,
    use_cache=TRUE, high_res=FALSE, res=72, out=NULL, directory=NULL,
    transparent_colors=c("#FFFFFF", "#BFBFFF","#BFFFBF","#7F7F7F", "#808080"),
    unlink=TRUE, with_legend_image=FALSE, legend_horiz=FALSE, legend_space=100
) {
    pid <- gg$data$pathway_id %>% unique()
    if (length(pid)>1) {stop("Only one pathway is supported.")}
    url <- paste0("https://rest.kegg.jp/get/",pid,"/image")
    if (high_res) {
        ## Convert to reference ID
        cur_id <- pid
        pid <- paste0("map", regmatches(cur_id, gregexpr("[[:digit:]]+", cur_id)) %>% unlist())

        ## sanity check
        if (!startsWith(pid, "map")) {
            stop("High resolution image can be obtained for the reference pathway.")
        }
        url <- paste0("https://rest.kegg.jp/get/",pid,"/image")
        url <- paste0(url, "2x")
    }
    if (use_cache) {
        bfc <- BiocFileCache()
        path <- bfcrpath(bfc, url)    
    } else {
        path <- paste0(pid, ".png")
        if (!is.null(directory)) {
            path <- paste0(directory,"/",path)
        }
        download.file(url=url, destfile=path, mode='wb')
    }
    magick_image <- image_read(path)
    info <- image_info(magick_image)
    for (col in transparent_colors) {
        magick_image <- magick_image %>%
            image_transparent(col)
    }

    ## Modify original gg to align with the image
    gg <- gg + scale_x_continuous(expand=c(0,0), limits=c(0,info$width-1)) +
        scale_y_continuous(expand=c(0,0), limits=c(-1*info$height+1, 0))

    ## Obtain grob and get panel
    ggGrob <- ggplotGrob(gg)
    legendGrob <- NULL
    panelGrob <- gtable::gtable_filter(ggGrob, "panel")
    if (length(gtable::gtable_filter(ggGrob, "guide-box"))!=0) {
        legendGrob <- gtable::gtable_filter(ggGrob, "guide-box")
    }

    ## Export grob
    timestamp <- as.numeric(Sys.time())
    ggname <- paste0(pid, "_", timestamp, ".png")
    png(ggname, width=info$width, height=info$height, res=res, units="px")
        grid::grid.draw(panelGrob)
        if (!with_legend_image & with_legend & !is.null(legendGrob)) {
            grid::grid.draw(legendGrob)
        }
    dev.off()
    if (with_legend_image & !is.null(legendGrob)) {
        ggLegendName <- paste0(pid, "_legend_", timestamp, ".png")
        if (legend_horiz) {
            lw <- info$width
            lh <- legend_space
        } else {
            lw <- legend_space
            lh <- info$height
        }
        png(ggLegendName, width=lw, height=lh, res=res, units="px")
        grid::grid.draw(legendGrob)
        dev.off()
        from_gg_legend <- image_read(ggLegendName)
    }
    
    from_gg <- image_read(ggname)
    
    if (unlink) {
        unlink(ggname)
        if (with_legend_image & !is.null(legendGrob)) {
            unlink(ggLegendName)
        }
    }
    
    flat <- image_flatten(c(from_gg, magick_image))
    if (with_legend_image & !is.null(legendGrob)) {
        if (legend_horiz) {
            flat <- image_append(c(flat, from_gg_legend), stack=TRUE)
        } else {
            flat <- image_append(c(flat, from_gg_legend))
        }
    }
    if (is.null(out)) {
        out <- paste0(pid, "_ggkegg.png")
    }
    image_write(flat, out)
}