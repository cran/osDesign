\name{logit}
\alias{logit}
\alias{expit}
\title{
  Logit and inverse-logit functions
}
\description{
  logit() performs the logit function
  expit() performs the inverse logit function
}
\usage{
logit(p)
expit(x)
}
\arguments{
 \item{p}{A vector of probabilities on the interval [0,1] to be transformed into numbers on the real line.}
\item{x}{A vector of numbers on the real line to be transformed into probabilities.}
}
\details{
The \code{logit} function is log(p / (1 - p)).
The \code{expit} function is exp(x)/(1 + exp(x)).
}
\value{
\code{logit} returns a vector of real numbers.
\code{expit} returns a vector of probabilities.
}
%\references{
% ~put references to the literature/web site here ~
%}
%\author{
%  ~~who you are~~
%}
%\note{
%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
