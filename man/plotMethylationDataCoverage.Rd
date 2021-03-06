% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/coverage.R
\name{plotMethylationDataCoverage}
\alias{plotMethylationDataCoverage}
\title{Plot methylation data coverage}
\usage{
plotMethylationDataCoverage(methylationData1, methylationData2 = NULL, breaks,
  regions = NULL, conditionsNames = NULL, context = "CG",
  proportion = TRUE, labels = NULL, col = NULL, pch = c(1, 0, 16, 2, 15,
  17), lty = c(4, 1, 3, 2, 6, 5), contextPerRow = FALSE)
}
\arguments{
\item{methylationData1}{the methylation data in condition 1
(see \code{\link{methylationDataList}}).}

\item{methylationData2}{the methylation data in condition 2
(see \code{\link{methylationDataList}}). This is optional.}

\item{breaks}{a \code{numeric} vector specifing the different values for the
thresholds when computing the coverage.}

\item{regions}{a \code{\link{GRanges}} object with the regions where to
compute the coverage. If \code{NULL}, the coverage is computed genome-wide.}

\item{conditionsNames}{a vector of character with the names of the conditions
for \code{methylationData1} and \code{methylationData2}.}

\item{context}{the context in which the DMRs are computed (\code{"CG"},
\code{"CHG"} or \code{"CHH"}).}

\item{proportion}{a \code{logical} value indicating whether proportion or
counts will be plotted.}

\item{labels}{a \code{vector} of \code{character} used to add a subfigure
character to the plot. If \code{NULL} nothing is added.}

\item{col}{a \code{character} vector with the colors. It needs to contain a
minimum of 2 colors per condition. If not or if \code{NULL}, the defalut
colors will be used.}

\item{pch}{the R symbols used to plot the data. It needs to contain a minimum
of 2 symbols per condition. If not or if \code{NULL}, the defalut symbols
will be used.}

\item{lty}{the line types used to plot the data. It needs to contain a
minimum of 2 line types per condition. If not or if \code{NULL}, the defalut
line types will be used.}

\item{contextPerRow}{a \code{logical} value indicating if the each row
represents an individual context. If \code{FALSE}, each column will represent
an individual context.}
}
\value{
Invisibly returns \code{NULL}
}
\description{
This function plots the coverage for the bisulfite sequencing data.
}
\details{
This function plots the proportion of cytosines in a specific
context that have at least a certain number of reads (x-axis)
}
\examples{
# load the methylation data
data(methylationDataList)

# plot the coverage in CG context
par(mar=c(4, 4, 3, 1)+0.1)
plotMethylationDataCoverage(methylationDataList[["WT"]],
                           methylationDataList[["met1-3"]],
                           breaks = c(1,5,10,15), regions = NULL,
                           conditionsNames = c("WT","met1-3"),
                           context = c("CG"), proportion = TRUE,
                           labels = LETTERS, col = NULL,
                           pch = c(1,0,16,2,15,17), lty = c(4,1,3,2,6,5),
                           contextPerRow = FALSE)

\dontrun{
# plot the coverage in all three contexts
plotMethylationDataCoverage(methylationDataList[["WT"]],
                           methylationDataList[["met1-3"]],
                           breaks = 1:15, regions = NULL,
                           conditionsNames = c("WT","met1-3"),
                           context = c("CG", "CHG", "CHH"),
                           proportion = TRUE, labels = LETTERS, col = NULL,
                           pch = c(1,0,16,2,15,17), lty = c(4,1,3,2,6,5),
                           contextPerRow = FALSE)
}
}
\author{
Nicolae Radu Zabet
}
\seealso{
\code{\link{computeMethylationDataCoverage}},
\code{\link{methylationDataList}}
}

