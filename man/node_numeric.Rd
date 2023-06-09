% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{node_numeric}
\alias{node_numeric}
\title{node_numeric}
\usage{
node_numeric(num, num_combine = mean, name = "name", how = "any")
}
\arguments{
\item{num}{named vector or tibble with id and value column}

\item{num_combine}{how to combine number when multiple hit in the same node}

\item{name}{name of column to match for}

\item{how}{how to match the node IDs with the queries 'any' or 'all'}
}
\value{
numeric vector
}
\description{
simply add numeric attribute to node of tbl_graph
for matrix input, use `append_node_value`
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
graph <- graph |> 
           mutate(num=node_numeric(c(1.1) |> setNames("hsa:1029")))


}
