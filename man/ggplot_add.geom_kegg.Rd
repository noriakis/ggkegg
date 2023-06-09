% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{ggplot_add.geom_kegg}
\alias{ggplot_add.geom_kegg}
\title{ggplot_add.geom_kegg}
\usage{
\method{ggplot_add}{geom_kegg}(object, plot, object_name)
}
\arguments{
\item{object}{An object to add to the plot}

\item{plot}{The ggplot object to add object to}

\item{object_name}{The name of the object to add}
}
\value{
ggplot2 object
}
\description{
ggplot_add.geom_kegg
}
\examples{
nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
                   x=c(1,1),
                   xmin=c(-1,-1),
                   xmax=c(2,2),
                   y=c(1,1),
                   ymin=c(-1,-1),
                   ymax=c(2,2),
                   type=c("gene","gene"),
                   bgcolor=c("red","blue"))
edges <- data.frame(from=1, to=2)
graph <- tbl_graph(nodes, edges)
p <- ggraph(graph, layout="manual", x=x, y=y)+
geom_kegg()
}
