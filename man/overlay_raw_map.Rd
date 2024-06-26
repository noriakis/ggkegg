% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlay_functions.R
\name{overlay_raw_map}
\alias{overlay_raw_map}
\title{overlay_raw_map}
\usage{
overlay_raw_map(
  pid = NULL,
  directory = NULL,
  transparent_colors = c("#FFFFFF", "#BFBFFF", "#BFFFBF"),
  adjust = FALSE,
  adjust_manual_x = NULL,
  adjust_manual_y = NULL,
  clip = FALSE,
  use_cache = TRUE,
  interpolate = TRUE,
  high_res = FALSE,
  fix_coordinates = TRUE
)
}
\arguments{
\item{pid}{pathway ID}

\item{directory}{directory to store images if not use cache}

\item{transparent_colors}{make these colors transparent to overlay
Typical choice of colors would be:
"#CCCCCC", "#FFFFFF","#BFBFFF","#BFFFBF", "#7F7F7F", "#808080",
"#ADADAD","#838383","#B3B3B3"}

\item{adjust}{adjust the x- and y-axis location by 0.5 in data coordinates}

\item{adjust_manual_x}{adjust the position manually for x-axis
Override `adjust`}

\item{adjust_manual_y}{adjust the position manually for y-axis
Override `adjust`}

\item{clip}{clip the both end of x- and y-axis by one dot}

\item{use_cache}{whether to use BiocFileCache()}

\item{interpolate}{parameter in annotation_raster()}

\item{high_res}{Use high resolution (2x) image for the overlay}

\item{fix_coordinates}{fix the coordinate (coord_fixed)}
}
\value{
ggplot2 object
}
\description{
Overlay the raw KEGG pathway image on ggraph
}
\examples{
## Need `pathway_id` column in graph 
## if the function is to automatically infer
graph <- create_test_pathway() |> mutate(pathway_id="hsa04110")
ggraph(graph) + overlay_raw_map()

}
