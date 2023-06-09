% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metacyc_function.R
\name{parseMetaCycPathwayReactions}
\alias{parseMetaCycPathwayReactions}
\title{parseMetaCycPathwayReactions}
\usage{
parseMetaCycPathwayReactions(file)
}
\arguments{
\item{file}{path to pathways.dat}
}
\value{
nested list of MetaCyc pathway descriptions
}
\description{
parse MetaCyc "pathways.dat"
and return the nested list of:
}
\details{
UNIQUE-ID
REACTION-LAYOUT
COMMON-NAME
SPECIES
TAXONOMIC-RANGE
}
\examples{
file <- "pathways.dat"
\dontrun{parseMetaCycPathwayReactions(file, candSp="all")}
}
