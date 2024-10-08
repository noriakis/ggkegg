% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{edge_numeric}
\alias{edge_numeric}
\title{edge_numeric}
\usage{
edge_numeric(
  num,
  num_combine = mean,
  how = "any",
  name = "name",
  sep = " ",
  remove_dot = TRUE
)
}
\arguments{
\item{num}{named vector or tibble with id and value column}

\item{num_combine}{how to combine number when multiple hit in the same node}

\item{how}{`any` or `all`}

\item{name}{name of column to match for}

\item{sep}{separater for name, default to " "}

\item{remove_dot}{remove "..." in the name}
}
\value{
numeric vector
}
\description{
add numeric attribute to edge of tbl_graph
}
\examples{
graph <- create_test_pathway()
graph <- graph |> activate("edges") |>
            mutate(num=edge_numeric(c(1.1) |>
            setNames("degradation"), name="subtype_name"))
}
