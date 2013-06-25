\name{broccoli}
\alias{broccoli}
\docType{data}
\title{The Effects of Drought Stress on Leaf Development in a \emph{Brassica oleracea} population}
\description{The effect of drought stress on \emph{Brassica oleracea} should be investigated, selecting drought stress resistant out of a population of different DH genotypes. The study was carried out on 48 DH lines developed from F1 plants of a cross between the rapid cycling chinese kale (\emph{Brassica oleracea} var. \emph{alboglabra} (L.H. Bailey) Musil) and broccoli (\emph{Brassica oleracea} var. \emph{italica} Plenck). 2 stress treatments (not watered and a watered control) are randomly assigned to 4 plants per genotype (2 per treatment) resulting in 192 plants in total. For the genotypes 5, 17, 31, 48, additional 12 plants (6 per treatment) are included into the completely randomized design, which results in a total of 240 plants. For each plant the length of the youngest leaf at the beginning of the experiment is measured daily for a period of 16 days. For the additional 12 plants of the 4 genotypes the leaf water potential was measured as a secondary endpoint (omitted here); due to these destructive measurements some dropouts occur.}
\usage{data(broccoli)}
\format{
  A data frame with 3689 observations on the following 5 variables.
  \describe{
    \item{\code{LeafLength}}{Length of the youngest leaf [cm]}
    \item{\code{ID}}{Plant identifier for 240 plants}
    \item{\code{Stress}}{Drought stress treatment with 2 levels (control/drought)}
    \item{\code{Genotype}}{Genotype ID with 48 levels}
    \item{\code{Day}}{Day of repeated measurement (1,2,...,16)}
  }
}
\references{
Uptmoor, R.; Osei-Kwarteng, M.; Guertler, S. & Stuetzel, H. Modeling the Effects of Drought Stress on Leaf Development in a Brassica oleracea Doubled Haploid Population Using Two-phase Linear Functions. Journal of the American Society for Horticultural Science, 2009, 134, 543-552.
}



\keyword{datasets}
