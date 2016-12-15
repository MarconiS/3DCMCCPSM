library(truncnorm)
library(ineq)
library(hydroGOF)

pr.dir <- "/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model"
sites <- c("ITCastelporziano", "ITCollelongo", "ITRenon", "BEBrasschaat", "FIHyytiala")

setwd(paste(pr.dir, "/src", sep=""))
curr.dir <- getwd()
build = TRUE
if(build){
  system("sh buildtool.sh")
}

for (j in sites) {
  
  target.dir <- paste("/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/", j, sep="")
  par1 <- paste("./3D_CMCC_Forest_Model -i /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/", j, sep="")
  par2 <- " -o /Users/sergiomarconi//Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/debug_output_6.1/debug_output"
  par3 <- " -b /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/daily_output_6.1/daily_output"
  par4 <- " -f /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/monthly_output_6.1/monthly_output"
  par5 <- " -e /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/annual_output_6.1/annual_output"
  par6 <- " -n /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/soil_output_6.1/soil_output"
  par7 <- paste(" -d /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/",j, "/input.txt", sep="")
  par8 <- paste(" -m ",pr.dir, "/input/",j,"/111_111_2000.txt", sep="")
  par9 <- paste(" -s /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/",j,"/site.txt", sep="")
  par10 <- paste(" -c /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/",j,"/settings.txt", sep="")
  
  running.string <- paste(par1, par2, par3,par4,par5,par6,par7,par8,par9,par10, sep="")
  
  eclipseoutput.dir <- "~/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/daily_output_6.1"
  observed <-read.table(file=paste("/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/Optimization/",j,"NEE.csv", sep=""), header=TRUE)
  
  
  setwd(target.dir)
  target.file <- data.frame(read.table(file="speciesOptim.txt", header=FALSE))
  params.list <- as.data.frame.array(read.csv(paste(pr.dir, "/input/PrPar.csv", sep=""), header = FALSE, stringsAsFactors = FALSE)) 
  
  nparams <- length(params.list[,1])
  row.numpars <- rep(0,nparams)
  for(i in 1: nparams){	
    parm <- params.list[i,]
    row.numpars[i] <- which(target.file[,1]==parm, arr.ind=TRUE)
  }
  for(runs in 1:1)
  {
    #new.values <- target.file[row.numpars,2]*2 # just to try
    #fill with the specific param set
    new.fill <- read.delim("Prior.txt", header=FALSE)
    
    new.values<- rep(0,nparams)
    for(i in 1: nparams){	
      #rtrunknorm because the values of the parameters are all positive
      new.values[i] <-rtruncnorm(n=1, a=0.001, b=(new.fill[i,2]+2*new.fill[i,3]), mean=new.fill[i,2], sd=new.fill[i,3])
      if(params.list[i,]=="SWPOPEN" || params.list[i,] == "SWPCLOSE"){
        new.values[i] = -new.values[i]
      }
    }
    
    target.file[row.numpars,2] <- new.values
    write.table(x=target.file, file=paste(getwd(),"/speciesOptim.txt",sep=""), quote=FALSE, row.names=FALSE, col.names = FALSE)
    
    setwd(curr.dir)
    system(running.string)
    predicted <- read.delim(file=paste(eclipseoutput.dir, "/daily_output_", j,"_DNDC.txt", sep=""), header=TRUE)
    
    if(is.na.data.frame(mean(predicted$NEE)))
    {
      out.val <- c(cor(observed, predicted$NEE),"NA", "NA", "NA")
    }else{
      lse <- lsfit(observed,predicted$NEE, intercept = FALSE, tolerance = 1e-07, yname = NULL)
      out.val <- c(cor(observed, predicted$NEE),NSE.data.frame(predicted$NEE, observed),
                   mae.data.frame(predicted$NEE, observed))
    }
    
    FF <- as.matrix(t(c(new.values, out.val)))
    write.table(FF, file = paste("/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/Optimization",j,"_boot.csv", sep=""), sep = ",", 
                col.names = FALSE, append=TRUE, row.names = FALSE)
  }
} 
  
  

