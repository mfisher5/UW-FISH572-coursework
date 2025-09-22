#' @description A filled contour plot of the simulated population
#' @author Allan Hicks
#' @par dat a simulated data.frame with columns for X, Y, eta (the simulated value for that cell), etc 
#' @par valueName name of the column to be plotted as shades of colors
#' @par ... plotting parameters
#' @examples
#' dat <- data.frame(X=c(1,1,1,2,2,2,3,3,3),Y=c(1,2,3,1,2,3,1,2,3),eta=c(1:9))
#' filledContourPlot.fn(dat, valueName="eta",nlevels=9)

filledContourPlot.fn <- function(dat, valueName="eta_scaled",...) {
  X <- as.numeric(names(table(dat$X)))
  Y <- as.numeric(names(table(dat$X)))
  
  val <- matrix(NA,nrow=length(Y),ncol=length(X))
  for(xx in 1:nrow(val)) {
    ind <- dat$X==xx
    val[xx,dat$Y[ind]] <- dat[ind,valueName]
  }
  # val <- matrix(dat.yr1$eta,100,100) #alternative way to do it, but make sure it maps correctly to X and Y
  
  filled.contour(x=X,y=Y,val,...)#nlevels = 30,color.palette=colorRampPalette(c("darkseagreen1", "darkgreen")))
  
  dots <- list(...)
  if(!("xlab" %in% names(dots))) mtext("X",side=1,outer=T,line=-2)
  if(!("ylab" %in% names(dots))) mtext("Y",side=2,outer=T,line=-2)
  
  invisible(list(x=X,y=Y,z=val))
}


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
  out <- list(observed=qq*samp, catchability=qq, samp=samp)
  return(out)
}


#' samples simulated data from a grid of X and Y locations and applies observation error and catchability to produce a sample
#' @description applies unbiased lognormal observation error and catchability to produce observed values
#' @author Allan Hicks
#' @par xy data.frame of locations along the X and Y axes with names 'X' and 'Y'
#' @par dat a simulated data.frame with columns for X, Y, and eta (the simulated value for that cell) that will be sampled from
#' @par obsCV observation coefficient of variation
#' @par catchabilityPars a list of mean and variance for gamma catchability. Must have names 'mean' and 'var'. var less than or equal to 0 fixes catchability at the mean parameter. NULL sets catchability to 1.
#' @par varName the name of the variable in the dataframe that is to error applied
#' @examples
#' dat <- data.frame(X=c(1,1,1,2,2,2,3,3,3),Y=c(1,2,3,1,2,3,1,2,3),eta=c(1:9))
#' xy <- data.frame(X=c(1,1,2,3,3,3),Y=c(2,3,2,1,2,3))
#' obsCV <- 0.1
#' catchabilityPars <- list(mean=1,var=0.01)
#' sampleGrid.fn(xy, dat, obsCV, catchabilityPars)

sampleGrid.fn <- function(xy, dat, obsCV, catchabilityPars,varName="eta_scaled") {
  gridLocs <- paste(dat$X,dat$Y)
  sampLocs <- paste(xy$X,xy$Y)
  out <- dat[gridLocs%in%sampLocs,]
  
  out$observed <- observationsLN.fn(n=length(out$eta), median=out[,varName], obsCV, catchabilityPars)$observed
  return(out)
}


###################################
# use the functions

# load/source the above functions

#load data (use the files with a suffix of 's')
#eta is the true value
#eta_scaled is the true value scaled to have the same mean as sim set 20
#setwd("C:\\Users\\stan.kotwicki\\Work\\projects\\2022\\UW_class\\code\\Fish572-main\\Simulated Data-20221223T205153Z-001\\Simulated Data")
dat <- readRDS(paste0("./Simulated Data-20221223T205153Z-001/Simulated Data/","sim_dat_1s.RDS"))
dat.yr1 <- dat[dat$year==1,]


###############################################
#Explore data
#plot locations (a grid of 100 X 100)
plot(dat.yr1$X,dat.yr1$Y,pch=20,cex=0.5,xlab="X",ylab="Y",las=1)

#plot eta_scaled
filledContourPlot.fn(dat.yr1, valueName="eta_scaled", nlevels = 30,color.palette=colorRampPalette(c("darkseagreen1", "darkgreen")))


#############################################
#example using sampleGrid

#all grid locations
xyAll <- paste(dat.yr1$X,dat.yr1$Y)

#sample from grid locations
xySample <- sample(xyAll, 10, replace=F)
tmp <-  t(as.data.frame(strsplit(xySample,"\\s+")))
rownames(tmp) <- NULL
xy <- data.frame(X=as.numeric(tmp[,1]),Y=as.numeric(tmp[,2]))

obsCV <- 0.1
catchabilityPars <- list(mean=1,var=0.0)
sampleGrid.fn(xy, dat.yr1, obsCV, catchabilityPars)


####################################################
# loop to show the difference between years for a specific simulated dataset
#setwd("C:\\IPHC\\OneDrive - International Pacific Halibut Commission\\UW-SAFS\\Courses\\FISH572-Survey\\Simulations\\Data")
file=paste0("./Simulated Data-20221223T205153Z-001/Simulated Data/","sim_dat_1s.RDS")
dat <- readRDS(file)
for(yr in unique(dat$year)) {
  dat.yr<- dat[dat$year==yr,]
  mean.eta=mean(dat.yr$eta_scaled)
  filledContourPlot.fn(dat.yr, nlevels = 30,valueName="eta_scaled",zlim=c(0,max(dat$eta_scaled)),main=paste("Year", yr, file, "Mean.Eta", round(mean.eta,4)))
  #aa=locator(n=1)
  #if the margins are not right, close all graphics windows and try again
}


####################################################
# loop to show the differences between all 20 simulated datasets
#setwd("C:\\IPHC\\OneDrive - International Pacific Halibut Commission\\UW-SAFS\\Courses\\FISH572-Survey\\Simulations\\Data")
#setwd("./Simulated Data-20221223T205153Z-001/Simulated Data")
for(ii in 1:20) {
  # file=paste0("sim_dat_",ii,".RDS")
  file=paste0("./Simulated Data-20221223T205153Z-001/Simulated Data/","sim_dat_",ii,"s.RDS")
  dat <- readRDS(file)
  dat.yr1 <- dat[dat$year==2,]
  mean.eta=mean(dat.yr1$eta_scaled)
  filledContourPlot.fn(dat.yr1, nlevels = 30,valueName="eta_scaled",zlim=c(0,120),
                       main=ii)
  aa=locator(n=1)
  #if the margins are not right, close all graphics windows and try again
}
