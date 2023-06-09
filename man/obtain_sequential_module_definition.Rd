% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module_functions.R
\name{obtain_sequential_module_definition}
\alias{obtain_sequential_module_definition}
\title{obtain_sequential_module_definition}
\usage{
obtain_sequential_module_definition(kmo, name = "1", block = NULL)
}
\arguments{
\item{kmo}{module object}

\item{name}{name of definition when multiple definitions are present}

\item{block}{specify if need to parse specific block}
}
\value{
list of module definitions
}
\description{
Given module definition and block number,
Recursively obtain graphical represencation of block and 
connect them by pseudo-nodes representing blocks.
}
\examples{
mo <- create_test_module()
sequential_mod <- obtain_sequential_module_definition(mo)
}
