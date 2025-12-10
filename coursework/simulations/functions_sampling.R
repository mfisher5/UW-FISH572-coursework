#' Function to draw a single random catchability from a gamma distribution with mean (mm) and variance (vv)
#' @author Allan Hicks, Stan Kotwicki
#' @param mm mean of the distribution
#' @param vv variance of the distribution
#' @return a vector of random draws for catchability
#' @export
#' @examples
#' catchability.gamma(10,.8,.01)

catchability.gamma=function(nn,mm,vv){
  shape=mm^2/vv
  scale=vv/mm
  catchability=rgamma(nn,shape=shape, scale=scale)
  return(catchability)
}


#' samples simulated data and applies observation error and catchability
#' @description applies unbiased lognormal observation error and catchability to produce observed values
#' I^{\hat}=qIe^{\epsilon}
#' @author Allan Hicks
#' @par n number of samples to produce
#' @par median of the simulated (population) value
#' @par cv coefficient of variation of observation error
#' @par catchabilityParams a list of mean and variance for gamma catchability. Must have names 'mean' and 'var'. var less than or equal to 0 fixes catchability at the mean parameter. NULL sets catchability to 1.
#' @examples
#' observationsLN.fn(3, 21, 0.2, list(mean=1.0,var=0.0000001))
#' observationsLN.fn(3, 21, 0.2, list(mean=1.0,var=0.01))
#' observationsLN.fn(3, 21, 0.2, list(mean=1.0,var=0.0))
#' observationsLN.fn(3, 21, 0.2, list(mean=0.8,var=-1))
#' observationsLN.fn(3, 21, 0.2, NULL)

observationsLN.fn <- function(n, median, cv, catchabilityParams=list(mean=1,var=0.1)) {
  logSigma <- sqrt(log(cv^2+1))
  samp <- rlnorm(n, log(median), logSigma)
  if(is.null(catchabilityParams)) {
    qq <- 1
  } else {
    if(catchabilityParams$var<=0) {
      qq <- catchabilityParams$mean
    } else {
      qq <- catchability.gamma(n, catchabilityParams$mean, catchabilityParams$var)			
    }
  }
  out <- data.frame(observation=qq*samp, catchability=qq, samp=samp)
  return(out)
}


#' samples simulated data from a grid of X and Y locations and applies observation error and catchability to produce a sample
#' @description applies unbiased lognormal observation error and catchability to produce observed values
#' @author Allan Hicks
#' @par dat a simulated data.frame with columns for X, Y, and the simulated value for that cell. This is the sample and unobserved grid cells have been removed.
#' @par obsCV observation coefficient of variation (set to zero if obs error already applied)
#' @par catchabilityPars a list of mean and variance for gamma catchability. Must have names 'mean' and 'var'. var less than or equal to 0 fixes catchability at the mean parameter. NULL sets catchability to 1.
#' @par varName the name of the variable in the dataframe that is to have error/catchability applied
#' @examples
#' dat <- data.frame(X=c(1,1,2,3,3,3),Y=c(1,2,2,1,2,3),eta=c(1:6))
#' obsCV <- 0.1
#' catchabilityPars <- list(mean=1,var=0.01)
#' sampleGrid.fn(xy, dat, obsCV, catchabilityPars)

sampleGrid.fn <- function(dat, obsCV, catchabilityPars,varName="eta_scaled") {
  out <- observationsLN.fn(n=nrow(dat), median=dat[,varName], obsCV, catchabilityPars)
  return(cbind(dat,out))
}

