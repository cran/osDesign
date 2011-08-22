\name{infants0709}
\alias{infants0709}
\docType{data}
\title{
Infant vital statistic data from North Carolina
}
\description{
Individual-level infant mortality data on 387,705 births in the U.S. state of North Carolina, between 2007 and 2009, inclusive.
}
\usage{data(infants0709)}
\format{
  A data frame with 387,705 observations on the following 17 variables.
  \describe{
    \item{\code{county}}{North Carolina county in which birth occurred}
    \item{\code{year}}{Year of birth}
    \item{\code{sex}}{Infant's gender}
    \item{\code{race}}{Race of mother/child}
    \item{\code{mage}}{Age of mother}
    \item{\code{plural}}{Plurality}
    \item{\code{weeks}}{Completed weeks of gestation (calculated)}
    \item{\code{weight}}{Birth weight group}
    \item{\code{gained}}{Weight gained by mother during pregnancy, in pounds (up to 98 lbs; values of 98 represent 98 or more pounds gained during pregnancy)}
    \item{\code{cignum}}{Number of cigarettes smoked by mother per day, up to 97.  98 represents a smoking mother, where the number of cigarettes smoked per day is unknown.}
    \item{\code{nonwhite}}{Logical vector: race variable not equal to 1}
    \item{\code{smoker}}{Logical vector: cignum not equal to 0}
    \item{\code{olderm}}{Logical vector: mother's age >= 35}
    \item{\code{lowbw}}{Logical vector: infant's birth weight <= 2500 grams}
    \item{\code{lowgain}}{Logical vector: mother gained less than 15 pounds during pregnancy}
    \item{\code{female}}{Logical vector: infant's sex is female}
    \item{\code{preemie}}{Logical vector: infant born before completing gestational week 37}
  }
}
\details{
The \code{sex} variable is coded as follows:
  1 = Male
  2 = Female

The \code{race} variable is coded as follows:
  0 = Other non-white
  1 = White
  2 = Black
  3 = American Indian
  4 = Chinese
  5 = Japanese
  6 = Hawaiin
  7 = Filipino
  8 = Other Asian or Pacific Islander

The \code{plural} variable is coded as follows:
  1 = Singleton
  2 = Twins
  3 = Triplets
  4 = Quadruplets
  5 = Quintuplets of higher
  9 = Unknown  
  
The \code{weight} variable is coded as follows:  
  0 = 500 grams or less
  1 = 501 to 1000 grams
  2 = 1001 to 1500 grams
  3 = 1501 to 2000 grams
  4 = 2001 to 2500 grams
  5 = 2501 to 3000 grams
  6 = 3001 to 3500 grams
  7 = 3501 to 4000 grams
  8 = 4001 to 4500 grams
  9 = 4501 grams or more
}
\source{
The data were compiled by the North Carolina State Center for Health Statistics (http://www.irss.unc.edu/).
}
\keyword{datasets}
