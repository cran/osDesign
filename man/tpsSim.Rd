\name{tpsSim}
\Rdversion{1.7}
\alias{tpsSim}
\title{
Simulation function for two-phase study designs.
}
\description{
  Monte Carlo based evaluation of operating characteristics for estimators of the components of a logistic regression model, based on the two-phase and case-control study designs (Breslow and Chatterjee, 1999; Prentice and Pykle, 1979).
}
\usage{
tpsSim(B=1000, betaTruth, X, N, strata, expandX="all", etaTerms=NULL,
       nII0=NULL, nII1=NULL, nII=NULL, nCC=NULL,
       alpha=0.05, threshold=c(-Inf, Inf), digits=1, betaNames=NULL,
       referent=2, monitor=NULL, cohort=TRUE, NI=NULL, returnRaw=FALSE)
}
\arguments{
  \item{B}{
    The number of datasets generated by the simulation.
  }
  \item{betaTruth}{
    Regression coefficients from the logistic regression model.
  }
  \item{X}{
    Design matrix for the logistic regression model. The first column should correspond to intercept. For each exposure, the baseline group should be coded as 0, the first level as 1, and so on.
  }
  \item{N}{
    A numeric vector providing the sample size for each row of the design matrix, \code{X}.
  }
  \item{strata}{
    A list of numeric vectors indicating which columns of the design matrix, \code{X}, are used to form the phase I stratification schemes. \code{strata=1} specifies the intercept and is, therefore, equivalent to a case-control study. \code{strata=0} indicates all possible stratified two-phase sampling schemes. \code{strata=list(2,3)} indicates 2 two-phase designs (that stratify on the 2nd and 3rd columns, separately) are to be considered.
  }
  \item{expandX}{
    Character vector indicating which columns of \code{X} to expand as a series of dummy variables. Useful when at least one exposure is continuous (and should not be expanded). Default is `all'; other options include `none' or character vector of column names. See Details, below.
  }
  \item{etaTerms}{
    Character vector indicating which columns of \code{X} are to be included in the model. See Details, below.
  }
  \item{nII0}{
    A numeric vector of sample sizes at phase II for controls. The length must correspond to the number of unique values for the \code{strata} specification.
  }
  \item{nII1}{
    A numeric vector of sample sizes at phase II for cases. The length must correspond to the number of unique values for the \code{strata} specification.
  }
  \item{nII}{
    A pair of numbers providing the sample sizes for controls and cases at phase II. This is only used when simulating all stratified two-phase sampling schemes (i.e., \code{strata=0}).
  }
  \item{nCC}{
    A pair of sample sizes at phase II for controls and cases in a case-control design. If left \code{NULL}, the values case-control sample sizes are taken as the sums of \code{n1} and \code{n0}, respectively.
  }
  \item{alpha}{
    Type I error rate assumed for the evaluation of coverage probabilities and power.
  }
  \item{threshold}{
    An interval that specifies truncation of the Monte Carlo sampling distribution of each estimator.
  }
  \item{digits}{
    Integer indicating the precision to be used for the output.
  }
  \item{betaNames}{
    An optional character vector of names for the regression coefficients, 
    \code{betaTruth}.
  }
  \item{referent}{
    An numeric value specifying which estimator is taken as the referent (denominator) for the relative uncertainty calculation. 1=CD, 2=CC, 3=WL, 4=PL, 5=ML (see Details below).
  }
  \item{monitor}{
    Numeric value indicating how often \code{tpsSim()} reports real-time progress on the simulation, as the \code{B} datasets are generated and evaluated. The default of \code{NULL} indicates no output.
  }
  \item{cohort}{
    Logical flag. TRUE indicates phase I is drawn as a cohort; FALSE indicates phase I is drawn as a case-control sample.
  }
  \item{NI}{
    A pair of integers providing the outcome-specific phase I sample
    sizes when the phase I data are drawn as a case-control sample. The
    first element corresponds to the controls and the second to the
    cases.
  }
\item{returnRaw}{
   Logical indicator of whether or not the raw coefficient and standard error estimates for each of the design/estimator combinations should be returned.
  }
}
\details{
A simulation study is performed to evaluate the operating
characteristics of various designs/estimators for \code{betaTruth}:
\itemize{
  \item{}{(a) complete data maximum likelihood (CD)}
  \item{}{(b) case-control maximum likelihood (CC)}
  \item{}{(c) two-phase weighted likelihood (WL)}
  \item{}{(d) two-phase pseudo- or profile likelihood (PL)}
  \item{}{(e) two-phase maximum likelihood (ML)}
}
  The operating characteristics are evaluated using the Monte Carlo
  sampling distribution of each estimator. The latter is generated using
  the following steps:
  \itemize{
    \item{}{(i) Specify the (joint) marginal exposure distribution of underlying population, using \code{X} and \code{N}.} 
    \item{}{(ii) Simulate outcomes for all sum(\code{N}) individuals in the population, based on an underlying logistic regression model specified via \code{betaTruth}.}
    \item{}{(iii) Evaluate the CD estimator on the basis of the complete data.}
    \item{}{(iv) Sample either (a) \code{ccDesign} controls and cases or (b) sum(\code{n0}) controls and sum(\code{n1}) cases, (without regard to the \code{strata} variable) and evaluate the CC estimator.}
    \item{}{(v) Stratify the population according to outcome and the \code{strata} argument, to form the phase I data.}
    \item{}{(vi) Sample \code{n0} controls and \code{n1} cases from their respective phase I strata.}
    \item{}{(vii) Evaluate the WL, PL and ML estimators.}
    \item{}{(viii) Repeat steps (ii)-(vii) \code{B} times.}
  }
  Both the CD and CC estimators are evaluated using the generic \code{\link{glm}} function. The three two-phase estimators are based on the \code{\link{tps}} function.
  
  The correspondence between \code{betaTruth} and \code{X}, specifically the ordering of elements, is based on successive use of \code{\link{factor}} to each column of \code{X} which is expanded via the \code{expandX} argument. Each exposure that is expanded must conform to a {0, 1, 2, ...} integer-based coding convention.
  
  The \code{etaTerms} argument is useful when only certain columns in \code{X} are to be included in the model. In the context of the two-phase design, this might be the case if phase I stratifies on some surrogate exposure and a more detailed/accurate measure is to be included in the main model.
  
  When evaluating operating characteristics, some simulated datasets may result in unusually large or small estimates. Particularly, when the the case-control/phase II sample sizes are small. In some settings, it may be desirable to truncate the Monte Carlo sampling distribution prior to evaluating operating characteristics. The \code{threshold} argument indicates the interval beyond which point estimates are ignored. The default is such that all \code{B} datasets are kept.

NOTE: In some settings, the current implementation of the ML estimator returns point estimates that do not satisfy the phase I and/or phase II constraints. If this is the case a warning is printed and the 'fail' elements of the returned list is set to TRUE. An example of this is phenomenon is given the help file for \code{\link{tps}}. When this occurs, \code{tpsSim()} considers ML estimation for the particular dataset to have failed.

}
\value{
  \code{tpsSim()} returns an object of class 'tpsSim', a list containing all the input arguments, as well list \code{results} with the following components:
    \item{betaMean}{
    	Mean of the Monte Carlo sampling distribution for each regression coefficient estimator.
    }
    \item{betaMeanBias}{
    	Bias based on the mean, calculated as \code{betaMean} - \code{betaTruth}. 
    }
    \item{betaMeanPB}{
    	Percent bias based on mean, calculated as ((\code{betaMean} - \code{betaTruth}) / \code{betaTruth}) x 100. If a regression coefficient is
      zero, percent bias is not calculated and an NA is returned.
    }
    \item{betaMedian}{
    	Median of the Monte Carlo sampling distribution for each regression coefficient estimator.
    }  
    \item{betaMedianBias}{
    	Bias based on the median, calculated as \code{betaMedian} - \code{betaTruth}. 
    }
    \item{betaMedianPB}{
    	Percent bias based on median, calculated as ((\code{betaMedian} - \code{betaTruth}) / \code{betaTruth}) x 100. If a regression coefficient is
      zero, median percent bias is not calculated and an NA is returned.
    }
    \item{betaSD}{
    	Standard deviation of the Monte Carlo sampling distribution for each regression coefficient estimator.
    }
    \item{betaMSE}{
      Mean squared error of the Monte Carlo sampling distribution for each regression coefficient
      estimator.
    }
    \item{seMean}{
    	Mean of the Monte Carlo sampling distribution for the standard error estimates reported by glm().
    }
    \item{seRatio}{
    	Ratio of the mean reported standard error to the standard deviation of the Monte Carlo sampling
    	distribution for each regression coefficient estimator. The ratio is multiplied by 100.    
    }
    \item{betaCP}{
    	Coverage probability for Wald-based confidence intervals, evaluated on the basis of an \code{alpha} type I error rate.
    }
    \item{betaPower}{
    	Power against the null hypothesis that the regression coefficient is zero for a Wald-based test with an \code{alpha} type I error rate.
    }
    \item{betaRU}{
    	The ratio of the standard deviation of the Monte Carlo sampling
    	distribution for each estimator to the standard deviation of the Monte Carlo sampling distribution
    	for the estimator corresponding to \code{refDesign}. The ratio is multiplied by 100.
    }
Also returned is an object \code{failed} which is a vector consisting of the number of datasets excluded from the power calculations (i.e. set to \code{NA}), for each simulation performed. For the evaluation of general operating characteristics, the four reasons are: (1) lack of convergence indicated by \code{NA} point estimates returned by \code{\link{glm}} or \code{\link{tps}}; (2) lack of convergence indicated by \code{NA} standard error point estimates returned by \code{\link{glm}} or \code{\link{tps}}; (3) exclusion on the basis of the \code{threshold} argument; and (4) for the ML estimator only, the phase I and/or phase II constraints are not satisfied.
}
\note{
	A generic print method provides formatted output of the results.
}
\references{
Prentice, R. and Pyke, R. (1979) "Logistic disease incidence models and case-control studies." Biometrika 66:403-411.

Breslow, N. and Chatterjee, N. (1999) "Design and analysis of two phase studies with binary outcome applied to Wilms tumour prognosis." Applied Statistics 48:457-468.

Haneuse, S. and Saegusa, T. and Lumley, T. (2011) "osDesign: An R Package for the Analysis, Evaluation, and Design of Two-Phase and Case-Control Studies." Journal of Statistical Software, 43(11), 1-29.
}
\author{
  Sebastien Haneuse, Takumi Saegusa
}
\examples{
##
data(Ohio)

## Design matrix that forms the basis for model and 
## phase I strata specification
##
XM <- cbind(Int=1, Ohio[,1:3])      ## main effects only
XI <- cbind(XM, SbyR=XM[,3]*XM[,4]) ## interaction between sex and race

## 'True' values for the underlying logistic model
##
fitM <- glm(cbind(Death, N-Death) ~ factor(Age) + Sex + Race, data=Ohio,
            family=binomial)
fitI <- glm(cbind(Death, N-Death) ~ factor(Age) + Sex * Race, data=Ohio,
            family=binomial)
##
betaNamesM <- c("Int", "Age1", "Age2", "Sex", "Race")
betaNamesI <- c("Int", "Age1", "Age2", "Sex", "Race", "SexRace")

## Two-phase design stratified by age
##   * sample 50 from each of 6 phase I strata
##   * show primary output (% bias, 95% CP, relative uncertainty)
##
\donttest{
ocAge <- tpsSim(B=1000, betaTruth=fitM$coef, X=XM, N=Ohio$N, strata=2,
                nII0=c(50,50,50), nII1=c(50,50,50), betaNames=betaNamesM,
                monitor=100)
ocAge}

## All possible balanced two-phase designs
##  * 250 controls and 250 cases
##  * only show the relative uncertainty output
##
\donttest{
ocAll <- tpsSim(B=1000, betaTruth=fitM$coef, X=XM, N=Ohio$N, strata=0,
                nII=c(250, 250), betaNames=betaNamesM,
                monitor=100)
ocAll}

## Two-phase design stratified by race
##  * balanced solely on outcome
##  * only show the relative uncertainty output
##
\donttest{
ocRace <- tpsSim(B=1000, betaTruth=fitI$coef, X=XI, N=Ohio$N, strata=4,
                 nII0=c(200, 50), nII1=c(200, 50), betaNames=betaNamesI,
                 monitor=100)
ocRace}

## Comparison of two case-control designs
##  * 240 controls and 260 cases
##  * 240 controls and 260 cases
##  * only show the relative uncertainty output
##
\donttest{
ocCC <- tpsSim(B=1000, betaTruth=fitM$coef, X=XM, N=Ohio$N, strata=1,
               nII0=240, nII1= 260, nCC=c(200,300),
               betaNames=betaNamesM,
               monitor=100)
ocCC}

## Illustration of setting where one of the covariates is continuous
##  * restrict to black and white children born in 2003
##  * dichotomize smoking, mothers age, weight gain during pregnancy and weight weight
##  * note the use of 'etaTerms' to restrict to specific variables (the majority of which
##    are created)
##  * note the use of 'strata=list(11,12)' to simultaneously investigate stratification by
##     - 11th column in XM: derived 'smoker' variable
##     - 12th column in XM: derived 'teen' variable
##
## Warning: takes a long time!
##
\donttest{
data(infants)
##
infants <- infants[infants$year == 2003,]
##
infants$race[!is.element(infants$race, c(1,2))] <- NA  ## White/Black = 0/1
infants$race <- infants$race - 1
infants <- na.omit(infants)
##
infants$smoker  <- as.numeric(infants$cignum > 0)
infants$teen    <- as.numeric(infants$mage < 20)
infants$lowgain <- as.numeric(infants$gained < 20)
infants$lbw     <- as.numeric(infants$weight < 2500)
infants$weeks   <- (infants$weeks - 36) / 4       ## estimate a 4-week contrast
##
fitM  <- glm(death ~ smoker + teen + race + male+ lowgain + lbw + weeks,
             data=infants,
             family=binomial)
betaM <- fitM$coef
XM    <- cbind(Int=1, infants)
etaM  <- c("Int", "smoker", "teen", "race", "male", "lowgain", "lbw", "weeks")
##
tpsSim(B=1000, betaTruth=fitM$coef, X=XM, N=rep(1, nrow(XM)), strata=list(11,12),
       expand="none", etaTerms=etaM, nII=c(1000,1000),
       threshold=c(-20,20),
       monitor=100)}
}

