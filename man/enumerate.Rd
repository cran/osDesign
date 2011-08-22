\name{enumerate}
\alias{enumerate}
\alias{enumerate.count}
\alias{enumerate.window}
\title{
  Enumerate Function
}
\description{
  enumerate() generates a matrix of vectors that meet the margin totals MM and NN.  
  enumerate.count() gives the total number of vectors that meet the margin totals MM and NN.
}
\usage{
enumerate(MM, NN)
enumerate.count(MM, NN)
}
\arguments{
  \item{MM}{
MM is a matrix of margin totals.  Rows are groups, columns are margin totals
}
  \item{NN}{
NN is a matrix of outcome margin totals.  Rows are groups, columns are margin totals.  NN is always a K x 2 matrix, where K is the number of groups
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
enumerate returns a matrix
enumerate.count returns a numver
}
\references{
Smoot, E., and S. Haneuse. "On the Analysis of Hybrid Designs that Combine Group- and Individual-Level Data." Biometrics (in press, 2014).
}
\author{
  G. Malecha, E. Smoot   
}
%\note{
%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
data(infants0709, package = "osDesign")

#  Get marginal totals for low birth weight and smoking status by county
MM = table(infants0709$county, infants0709$smoker)
NN = table(infants0709$county, infants0709$lowbw)

#  Determine the number of possible solutions to margin totals for county 48
enumerate.count(MM[48,], NN[48,])

#  Generate matrix with all possible solutions to margin totals for county 48
enumerate(MM[48,], NN[48,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
