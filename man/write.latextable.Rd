% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{write.latextable}
\alias{write.latextable}
\title{Write Latex Table}
\usage{
write.latextable(x, d = NA, s = NA, no.rounding = numeric(), file = "")
}
\arguments{
\item{d}{if specified, all floats will be rounded to d decimal places}

\item{s}{if specified (and d is unspecified), then all floats will be rounded
to s signif}

\item{no.rounding}{column indices not to apply rounding}

\item{file}{filename for output, default to screen output}
}
\description{
Write Latex Table
}
\examples{
df <- data.frame(name=c('Alpha', 'Beta', 'Gamma', 'Delta'),
                 size=c(100L,200L,300L,400L), score=c(23.091,19.978,1119.229, 0.03089))
write.latextable(df, s=3)
}
