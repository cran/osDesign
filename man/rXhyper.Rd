\name{rXhyper}
\Rdversion{1.7}
\alias{rXhyper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random generation for the multivariate hypergeometric distribution
}
\description{
Generates random observations from a multivariate hypergeometric distribution.
}
\usage{
rXhyper(theta, data, number = 1)
}
\arguments{
  \item{theta}{ Odds ratio}
  \item{data}{ Population margins}
  \item{number}{ Number of random observations to be drawn}
}
%\details{
%}
\value{
Returns a vector.
}
\author{
Sebastien Haneuse
}
\seealso{
\code{\link{rhyper}}.
}
\examples{
##
    mdata = list(M0 = 50, M1 = 50, N0 = 70, N1 = 30)
    rXhyper(1.2, mdata, 1 )
}
