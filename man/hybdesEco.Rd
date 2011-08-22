\name{hybdesEco}
\alias{hybdesEco}
\title{
  Hybrid Design in the Pure Ecological setting -- MLE
}
\description{
  Computes the MLE for a Hybrid Design in the pure ecological setting, with two binary covariates, Z and W, in the model.  hybdesEco is specific to the model 
  logit(P(Y=1)) = b0 + b1*Z + b2*W
}
\usage{
hybdesEco(MM.Z, MM.W, NN, cc, aprx = "binom", start.mle = NA, group.int = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MM.Z}{
MM.Z is a matrix of margin totals for covariate Z.  Rows are groups, columns are margin totals
}
  \item{MM.W}{
MM.W is a matrix of margin totals for covariate W.  Rows are groups, columns are margin totals
}
  \item{NN}{
NN is a matrix of outcome margin totals.  Rows are groups, columns are margin totals.  NN is always a K x 2 matrix, where K is the number of groups
}
  \item{cc}{
cc is a list of case-control data.  Each element is a table with joint Z-W exposure (rows) and outcome (columns)
}
  \item{aprx}{
Type of approximation to use when calculating the hybrid likelihood.  Default is the binomial approximation.
}
  \item{start.mle}{
Starting value for the Newton-Raphson algorithm used to determine Hybrid Design MLE.
}
  \item{group.int}{
A logical indicator of whether or not groups should be treated as having different intercept parameters.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Smoot, E., and S. Haneuse. "On the Analysis of Hybrid Designs that Combine Group- and Individual-Level Data." Biometrics (in press, 2014).
}
\author{
E. Smoot
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
