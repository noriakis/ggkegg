% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{append_cp}
\alias{append_cp}
\title{append_cp}
\usage{
append_cp(res, how = "any", name = "name", pid = NULL)
}
\arguments{
\item{res}{enrichResult class}

\item{how}{how to determine whether the nodes is in enrichment results}

\item{name}{name column to search for query}

\item{pid}{pathway ID, if NULL, try to infer from graph attribute}
}
\value{
enrich_attribute column in node
}
\description{
append clusterProfiler results to graph
}
\examples{
nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
                   x=c(1,1),
                   xmin=c(-1,-1),
                   xmax=c(2,2),
                   y=c(1,1),
                   ymin=c(-1,-1),
                   ymax=c(2,2))
edges <- data.frame(from=1, to=2)
graph <- tbl_graph(nodes, edges)
if (require("clusterProfiler")) {
  cp <- enrichKEGG(nodes$name |>
                  strsplit(":") |> 
                  vapply("[",2,FUN.VALUE="character"))
}
graph <- graph |> mutate(cp=append_cp(cp,pid="hsa04110"))
}
