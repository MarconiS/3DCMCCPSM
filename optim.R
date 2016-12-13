library(truncnorm)
library(ineq)
library(hydroGOF)

setwd("/Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/src")
curr.dir <- getwd()

target.dir <- "/Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/input/DEHainich"

running.string <- "./3D_CMCC_Forest_Model -i /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/input/DEHainich -o /Users/sergiomarconi//Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/output_6.1/debug_output_6.1/debug_output -b /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/output_6.1/daily_output_6.1/daily_output -f /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/output_6.1/monthly_output_6.1/monthly_output -e /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/output_6.1/annual_output_6.1/annual_output -n /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/output_6.1/soil_output_6.1/soil_output -d /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/software/3D-CMCC-Forest-Model/input/DEHainich/Hainich_input_1ha_2000.txt -m /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/input/DEHainich/111_111_2000.txt -s /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/input/DEHainich/Hainich_site.txt -c /Users/sergiomarconi/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/input/DEHainich/Hainich_settings.txt"


eclipseoutput.dir <- "~/Documents/Classes/Statistics\ for\ Biological\ sciences/3D-CMCC-Forest-Model/output_6.1/daily_output_6.1"
observed <-read.table(file="/Users/sergiomarconi/Documents/Statistics\ for\ Biological\ sciences/Optimization/HaiNEE2000i.csv", header=TRUE)


setwd(target.dir)
target.file <- data.frame(read.table(file="FagussylvaticaOptim.txt", header=FALSE))
params.list <- c("ALPHA","K","SLA","SLA_RATIO", "MAXAGE","GROWTHSTART", "R0CTEM","F0CTEM","RES0CTEM","CN_LEAVES", "CN_FINE_ROOTS", "LEAF_FALL_FRAC_GROWING", "FINE_ROOT_LEAF",  "SAP_LEAF", "HMAX", "CRB", "DENMAX", "MINDAYLENGTH", "RHOMIN", "MAXINTCPTN", "MAXCOND","CN_LIVE_WOODS", "HMAX", "TBB", "TRHO")

nparams <- length(params.list)
row.numpars <- rep(0,nparams)
for(i in 1: nparams){	
  parm <- params.list[i]
  row.numpars[i] <- which(target.file[,1]==parm, arr.ind=TRUE)
}
for(runs in 1:1)
{
  #new.values <- target.file[row.numpars,2]*2 # just to try
  #fill with the specific param set
  new.fill <- read.delim("/Users/sergiomarconi/Shared/git/3D-CMCC-FEM/parameters/FagussylvaticaPrior.txt", header=FALSE)
  
  new.values<- rep(0,nparams)
  for(i in 1: nparams){	
    #rtrunknorm because the values of the parameters are all positive
    new.values[i] <-rtruncnorm(n=1, a=0, b=(new.fill[i,2]+2*new.fill[i,3]), mean=new.fill[i,2], sd=new.fill[i,3])
  }
  
  target.file[row.numpars,2] <- new.values
  write.table(x=target.file, file="/Users/sergiomarconi/Shared/git/3D-CMCC-FEM/software/3D-CMCC-Forest-Model/input/DEhainich/FagussylvaticaOptim.txt", quote=FALSE, row.names=FALSE, col.names = FALSE)
  
  setwd(curr.dir)
  system(running.string)
  
  predicted <- read.table(file=paste(eclipseoutput.dir, "/daily_output_Hainich_DNDC.txt", sep=""), header=TRUE)
  
  if(is.na.data.frame(mean(predicted$NEE)))
  {
    out.val <- c(cor(observed, predicted$NEE),"NA", "NA", "NA")
  }else{
    lse <- lsfit(observed,predicted$NEE, intercept = FALSE, tolerance = 1e-07, yname = NULL)
    out.val <- c(cor(observed, predicted$NEE),NSE.data.frame(predicted$NEE, observed),
                 mae.data.frame(predicted$NEE, observed))
  }
  
  FF <- as.matrix(t(c(new.values, out.val)))
  write.table(FF, file = "/Users/sergiomarconi/Documents/Statistics\ for\ Biological\ sciences/Optimization/OptimHainich_NEE2000.csv", sep = ",", 
              col.names = FALSE, append=TRUE, row.names = FALSE)
}
###### Sensitivity analisys
#used NSE as driving proxy; could do the same with the ranking score once calculated 
# ranking score calculated on excel for convenience; trivial to implement in the code if required
OptimNEE2000 <- read.csv("~/Documents/Statistics for Biological sciences/Optimization/OptimHainich_NEE2000.csv")
for(i in 1 : 27)
{
  png(paste("/Users/sergiomarconi/Documents/Statistics\ for\ Biological\ sciences/Optimization/Results/",i,"_HaiSens",".png", sep=""))
  plot(OptimHainich_NEE2000[c(i, 27)], pch = ".", ylim = c(-9, 1))
  dev.off()
}

#at this point is important to run the model with preset of parameter, whose outputs may be stored into predicted. That to give an idea of the improvements associated to the method. Done manually even though trivial to do automatically (change running.string to point to another settings file addressing to the classical species parameterization)
#predicted: after changing settings file

#target.file[row.numpars,2] <- new.values
#write.table(x=target.file, file="/Users/sergiomarconi/Shared/git/3D-CMCC-FEM/software/3D-CMCC-Forest-Model/input/DEhainich/FagussylvaticaOptim.txt", quote=FALSE, row.names=FALSE, col.names = FALSE)  

setwd(curr.dir)
system(running.string)
predicted <- read.table(file=paste(eclipseoutput.dir, "/daily_output_Hainich_DNDC.txt", sep=""), header=TRUE)

#predicted2
new.values <-c(0.067620181,	0.575746923,	69.84829306,	2.526925344,	221.5844218, 684.2302099,	0.019739911,	0.918312209,	0.186772948,	51.97827487,	92.78528923, 44.45703033	, 0.077708988, 0.027068638,	1805.537747,	26.94124062,	0.048033345, 0.33432972,	0.009427448,	0.266655572)

target.file[row.numpars,2] <- new.values
write.table(x=target.file, file="/Users/sergiomarconi/Shared/git/3D-CMCC-FEM/software/3D-CMCC-Forest-Model/input/DEhainich/FagussylvaticaOptim.txt", quote=FALSE, row.names=FALSE, col.names = FALSE)

setwd(curr.dir)
system(running.string)
predicted2 <- read.table(file=paste(eclipseoutput.dir, "/daily_output_Hainich_DNDC.txt", sep=""), header=TRUE)

png(paste("/Users/sergiomarconi/Documents/Statistics\ for\ Biological\ sciences/Optimization/Results/",i,"_HaiSeries",".png", sep=""))
plot(observed, pch =".", xlim = c(0,365),col = "black", ylim = c(-6,9))
par(new = TRUE)
lines(observed, xlim = c(0,365), ylim = c(-6,9), xlab = "", ylab = "")
par(new = TRUE)
plot(predicted$NEE, pch =".", col="green", xlim = c(0,365), ylim = c(-6,9), xlab = "", ylab = "")
par(new = TRUE)
lines(predicted$NEE, col="green", xlim = c(0,365), ylim =c(-6,9), xlab = "", ylab = "")
par(new = TRUE)
plot(predicted2$NEE, pch =".", col="red", xlim = c(0,365), ylim = c(-6,9), xlab = "", ylab = "")
par(new = TRUE)
lines(predicted2$NEE, col="red", xlim = c(0,365), ylim =c(-6,9), xlab = "", ylab = "")
par(new = TRUE)
dev.off()

y<-predicted2$NEE
png(paste("/Users/sergiomarconi/Documents/Statistics\ for\ Biological\ sciences/Optimization/Results/",i,"_HaiScatter",".png", sep=""))
plot(observed$NEE_st_fANN, predicted2$NEE , col = "red", xlim = c(-6,9), ylim = c(-6,9))
abline(lm(observed$NEE_st_fANN~predicted2$NEE), col="red",xlab = "", ylab = "") # regression line (y~x) 
par(new = TRUE)
plot(y, predicted2$NEE , pch =".", xlim =c(-6,9), ylim = c(-6,9))
abline(lm(y~predicted2$NEE),xlab = "", ylab = "") # regression line (y~x) 
par(new = TRUE)
plot(observed$NEE_st_fANN, predicted$NEE , col = "green",  xlim = c(-6,9), ylim = c(-6,9),xlab = "", ylab = "")
abline(lm(observed$NEE_st_fANN~predicted$NEE), col="green",xlab = "", ylab = "") # regression line (y~x) 
par(new = TRUE)
dev.off()

