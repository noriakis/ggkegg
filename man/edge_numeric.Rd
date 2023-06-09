% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{edge_numeric}
\alias{edge_numeric}
\title{edge_numeric}
\usage{
edge_numeric(num, num_combine = mean, how = "any", name = "name")
}
\arguments{
\item{num}{named vector or tibble with id and value column}

\item{num_combine}{how to combine number when multiple hit in the same node}

\item{how}{`any` or `all`}

\item{name}{name of column to match for}
}
\value{
numeric vector
}
\description{
add numeric attribute to edge of tbl_graph
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
edges <- data.frame(from=1, to=2, name="K00112")
graph <- tbl_graph(nodes, edges)
graph <- graph |> activate("edges") |>
   mutate(num=edge_numeric(c(1.1) |> setNames("K00112")))


}
