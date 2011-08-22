\name{hybdes}
\alias{hybdes}
\alias{hyblik}
\title{
  Hybrid Design MLE and likelihood
}
\description{
  hybdes() computes the MLE for a Hybrid Design.  
  hyblik() computes the likelihood for a Hybrid Design at a specified parameter vector
}
\usage{
hybdes(MM, NN, cc, ntrue = 0, aprx='binom', start.mle=NA, group.int=FALSE, betafct =
 function(x){return(x[1] + c(0,x[-1]) )}, print.level = 0, iterlim = 100)
hyblik(beta.matrix, MM, NN, cc, aprx = 'binom', ntrue = 0, group.int=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MM}{
MM is a matrix of margin totals.  Rows are groups, columns are margin totals
}
  \item{NN}{
NN is a matrix of outcome margin totals.  Rows are groups, columns are margin totals.  
NN is always a K x 2 matrix, where K is the number of groups
}
  \item{cc}{
cc is a list of case-control data.  Each element is a table with 
exposure (rows) and outcome (columns)
}
  \item{ntrue}{
The number of groups that should be calculated using the true hybrid likelihood, rather than an approximation.
}
  \item{aprx}{
Type of approximation to use when calculating the hybrid likelihood.  Default is the binomial approximation.
}
  \item{start.mle}{
Starting value for the Newton-Raphson algorithm used to determine Hybrid Design MLE.
}
  \item{group.int}{
A logical indicator of whether or not groups should be treated as having different intercept parameters.
}
\item{betafct}{
A function used to specify the model of interest by reparamterizing the hybrid likelihood.
betafct() takes in group-specific parameters associated with each level of the exposure variable.
The default function corresponds to a model with an intercept parameter and log-odds-ratio parameters relating levels of X to the baseline level, X = 0 (i.e. column 1 of MM).
}
\item{print.level}{ Argument passed into nlm()}
\item{iterlim}{ Argument passed into nlm()}
  \item{beta.matrix}{
Parameter values for likelihood calculation; used only in hyblik().  This should be entered in the form of a matrix, with one row per group and one column per parameter.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  \item{mle }{MLE of the hybrid design}
  \item{start.mle }{Result of clogit function (stratified case-control MLE) }
%% ...
}
\references{
Smoot, E., and S. Haneuse. "On the Analysis of Hybrid Designs that Combine Group- and Individual-Level Data." Biometrics (in press, 2014).
}
\author{
E. Smoot
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#hybdes(MM, NN, cc, approx='NA')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
