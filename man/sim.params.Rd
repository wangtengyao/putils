% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{sim.params}
\alias{sim.params}
\title{Simulation parameter data frame generation}
\usage{
sim.params(...)
}
\arguments{
\item{...}{each argument should be of the form of tag = vector, meaning the variable named 'tag' takes values in 'vector'.}
}
\description{
create a dataframe of all possible parameter combinations in lexicographic order (if tags are supplied, use tag for column names)
}
\details{
A sample usage is sim.params(tag1 = vec1, tag2 = vec2, tag3 = vec3).
}
