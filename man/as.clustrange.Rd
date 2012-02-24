\name{as.clustrange}
\alias{as.clustrange}
\alias{as.clustrange.twins}
\alias{as.clustrange.hclust}
\alias{plot.clustrange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build a clustrange object to compare different clustering solutions.
}
\description{
Convert a 
}
\usage{
as.clustrange(object, diss, weights=NULL, ...)
\method{as.clustrange}{twins}(object, diss, weights=NULL, ncluster, ...) 
\method{as.clustrange}{hclust}(object, diss, weights=NULL, ncluster, ...) 
\method{plot}{clustrange}(x, stat="noCH", legendpos="bottomright", 
    norm="none", withlegend=TRUE, lwd=1, col=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The object to convert such as a data.frame.}
  \item{diss}{A dissimilarity matrix or a dist object (see \code{\link{dist}}).}
  \item{weights}{Optional numerical vector containing weights.}
  \item{ncluster}{Integer. Maximum number of cluster. The range will include all clustering solution starting from two to \code{ncluster}.}
  \item{x}{A \code{clustrange} object to be plotted.}
  \item{stat}{Character. The list of statistics to plot or "noCH" to plot all statistics except "CH" and "CHsq" or "all" for all statistics. See \code{\link{wcClusterQuality}} for a list of possible values. It is also possible to use "RHC" to plot the quality measure 1-HC. Unlike HC, RHC should be maximized as all other quality measures.}
  \item{legendpos}{Character. legend position, see \code{\link{legend}}.}
  \item{norm}{Character. Normalization method of the statistics can be one of "none" (no normalization), "range" (given as (value -min)/(max-min), "zscore" (adjusted by mean and standard deviation) or "zscoremed" (adjusted by median and median of the difference to the median).}
  \item{withlegend}{Logical. If \code{FALSE}, the legend is not plotted.}
  \item{lwd}{Numeric. Line width, see \code{\link{par}}.}
  \item{col}{A vector of line colors, see \code{\link{par}}. If \code{NULL}, a default set of color is used.}
  \item{\dots}{Additionnal parameters passed to/from methods.}
}
\details{
\code{as.clustrange} convert objects to \code{clustrange} objects. \code{clustrange} objects contains a list of clustering solution with associated statistics and can be used to find the optimal clustering solution.

If \code{object} is a \code{data.frame} or a \code{matrix}, each column should be a clustering solution to be evaluated.

If \code{object} is an \code{hclust} or \code{twins} objects (i.e. hierarchical clustering output, see \code{\link{hclust}}, \code{\link{diana}} or \code{\link{agnes}}), the function compute all clustering solution ranging from two to \code{ncluster} and compute the associated statistics.
}
\value{
An object of class \code{clustrange} with the following elements:
	\describe{
		\item{\code{clustering}:}{A \code{data.frame} of all clustering solutions.}
		\item{\code{stats}:}{A \code{matrix} containing the clustering statistics of each cluster solution.}
	}
}

\seealso{
See also \code{\link{wcKMedRange}}, \code{\link{wcClusterQuality}}.
}


\examples{
data(mvad)
## Aggregating state sequence
aggMvad <- wcAggregateCases(mvad[, 17:86], weights=mvad$weight)

## Creating state sequence object
mvad.seq <- seqdef(mvad[aggMvad$aggIndex, 17:86], weights=aggMvad$aggWeights)

## COmpute distance using Hamming distance
diss <- seqdist(mvad.seq, method="HAM")

## Ward clustering
wardCluster <- hclust(as.dist(diss), method="ward", members=aggMvad$aggWeights)

## Computing clustrange from Ward clustering
wardRange <- as.clustrange(wardCluster, diss=diss, 
		weights=aggMvad$aggWeights, ncluster=15)

## Plot all statistics (standardized)
plot(wardRange, stat="all", norm="zscoremed", lwd=3)

## Plot HC, RHC and ASW
plot(wardRange, stat=c("HC", "RHC", "ASWw"), norm="zscore", lwd=3)

}