library(truncnorm)
library(ineq)
library(hydroGOF)

# setting environment -----------------------------------------------------

pr.dir <- "/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model"
sites <- c("ITCastelporziano", "ITCollelongo", "ITRenon", "BEBrasschaat", "FIHyytiala")
sites <- c("ITCastelporziano")

setwd(paste(pr.dir, "/src", sep=""))
curr.dir <- getwd()
build = FALSE
if(build){
  system("sh buildtool.sh")
}

# defining function -------------------------------------------------------

nsecalc <- function(guess){
  setwd(target.dir)
  new.values <- exp(guess);
  #fill with the specific param set
  target.file[row.numpars,2] <- new.values
  write.table(x=target.file, file=paste(getwd(),"/speciesOptim.txt",sep=""), quote=FALSE, row.names=FALSE, col.names = FALSE)
  setwd(curr.dir)
  system(running.string)
  predicted <- read.delim(file=paste(eclipseoutput.dir, "/daily_output_", j,"_DNDC.txt", sep=""), header=TRUE)
  
  if(is.na.data.frame(mean(predicted$NEE)))
  {
    out.val <- 1000000000000
  }else{
    row.index <- rep(0,length(observed$DoY))
    for(i in 1: length(observed$DoY)){	
      row.index[i] <- which(as.numeric(predicted[,1]==observed[i,1]) * 
                              as.numeric(predicted[,2]==observed[i,2])==1 , arr.ind=TRUE)
    }
    train.obs <- observed$NEE_st_fANN
    train.pred <- predicted[row.index,]$NEE
    out.val <- -NSE.data.frame(train.obs,train.pred)
  }
  return(out.val)
}


# main --------------------------------------------------------------------


for (j in sites) {
  
  j = sites
  target.dir <- paste("/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/", j, sep="")
  par1 <- paste("./3D_CMCC_Forest_Model -i /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/", j, sep="")
  par2 <- " -o /Users/sergiomarconi//Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/debug_output_6.1/debug_output"
  par3 <- " -b /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/daily_output_6.1/daily_output"
  par4 <- " -f /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/monthly_output_6.1/monthly_output"
  par5 <- " -e /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/annual_output_6.1/annual_output"
  par6 <- " -n /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/soil_output_6.1/soil_output"
  par7 <- paste(" -d /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/",j, "/input.txt", sep="")
  par8 <- readChar(paste(pr.dir, "/input/",j,"/trainMet.txt", sep=""), file.info(paste(pr.dir, "/input/",j,"/trainMet.txt", sep=""))$size)
  par9 <- paste(" -s /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/",j,"/site.txt", sep="")
  par10 <- paste(" -c /Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/input/",j,"/settings.txt", sep="")
  
  running.string <- paste(par1, par2, par3,par4,par5,par6,par7,par8,par9,par10, sep="")
  
  eclipseoutput.dir <- "~/Documents/Classes/StatisticsForBiologicalSciences/3D-CMCC-Forest-Model/output_6.1/daily_output_6.1"
  observed <-read.table(file=paste("/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/Optimization/",j,"NEE.csv", sep=""), header=TRUE, sep=",")
  
  
  setwd(target.dir)
<<<<<<< HEAD
  target.file <<- data.frame(read.table(file="speciesOptim.txt", header=FALSE))
=======
  target.file <- data.frame(read.table(file="speciesOptim.txt", header=FALSE))
>>>>>>> b2b9ac53d5e94312128b4f16099280040be1f771
  params.list <- as.data.frame.array(read.csv(paste(target.dir, "/PrPar.csv", sep=""), header = FALSE, stringsAsFactors = FALSE)) 
  
  nparams <- length(params.list[,1])
  row.numpars <- rep(0,nparams)
  for(i in 1: nparams){	
    parm <- params.list[i,1]
    row.numpars[i] <- which(target.file[,1]==parm, arr.ind=TRUE)
  }
  
  #fill with the specific param set
  new.fill <- params.list[,2:3]
  #randomizer 
  for(runs in 1:1000)
  {
<<<<<<< HEAD
    setwd(target.dir)
=======
    #fill with the specific param set
    new.fill <- params.list[,2:3]
    
>>>>>>> b2b9ac53d5e94312128b4f16099280040be1f771
    new.values<- rep(0,nparams)
    for(i in 1: nparams){	
      #rtrunknorm because the values of the parameters are all positive
      new.values[i] <-rtruncnorm(n=1, a=0.001, b=(new.fill[i,1]+2*new.fill[i,2]), mean=new.fill[i,1], sd=new.fill[i,2])
      if(params.list[i,1]=="SWPOPEN" || params.list[i,1] == "SWPCLOSE"){
        new.values[i] = -new.values[i]
      }
    }
    
    target.file[row.numpars,2] <- new.values
    write.table(x=target.file, file=paste(getwd(),"/speciesOptim.txt",sep=""), quote=FALSE, row.names=FALSE, col.names = FALSE)
    target.file[row.numpars,]
    
    setwd(curr.dir)
    system(running.string)
    predicted <- read.delim(file=paste(eclipseoutput.dir, "/daily_output_", j,"_DNDC.txt", sep=""), header=TRUE)
    
    if(is.na.data.frame(mean(predicted$NEE)))
    {
      out.val <- c("NA", "NA", "NA")
    }else{
      row.index <- rep(0,length(observed$DoY))
      for(i in 1: length(observed$DoY)){	
        row.index[i] <- which(as.numeric(predicted[,1]==observed[i,1]) * 
<<<<<<< HEAD
                                as.numeric(predicted[,2]==observed[i,2])==1 , arr.ind=TRUE)
=======
                                 as.numeric(predicted[,2]==observed[i,2])==1 , arr.ind=TRUE)
>>>>>>> b2b9ac53d5e94312128b4f16099280040be1f771
      }
      
      train.obs <- observed$NEE_st_fANN
      train.pred <-predicted[row.index,]$NEE
<<<<<<< HEAD
      #lse <- lsfit(train.obs,train.pred, intercept = FALSE, tolerance = 1e-07, yname = NULL)
=======
      lse <- lsfit(train.obs,train.pred, intercept = FALSE, tolerance = 1e-07, yname = NULL)
>>>>>>> b2b9ac53d5e94312128b4f16099280040be1f771
      out.val <- c(cor(train.obs,train.pred),NSE.data.frame(train.obs,train.pred),
                   mae.data.frame(train.obs,train.pred))
    }
    
    FF <- as.matrix(t(c(new.values, out.val)))
    write.table(FF, file = paste("/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/Optimization",j,"_boot.csv", sep=""), sep = ",", 
                col.names = FALSE, append=TRUE, row.names = FALSE)
  }
  #GA
  #fill with the specific param set
  best.random <- read.delim(file=paste("/Users/sergiomarconi/Documents/Classes/StatisticsForBiologicalSciences/Optimization",j,"_boot.csv", sep=""), sep=",", header=TRUE)
  best.random <- na.omit(best.random)
  best.initial <- best.random[which(best.random$NSE==max(best.random$NSE)),]
  new.fill <- as.vector(t(best.initial[,1:45]))

  target.file[row.numpars,2] <- new.fill
  setwd(target.dir)
  write.table(x=target.file, file=paste(getwd(),"/speciesOptim.txt",sep=""), quote=FALSE, row.names=FALSE, col.names = FALSE)
  setwd(curr.dir)
  system(running.string)
  
  trial.values <- new.fill
  trial.run <- nsecalc(trial.values)
  trial.optim <- optim(par=log(trial.values), fn=nsecalc, method="Nelder-Mead")
  best.params <- exp(trial.optim$par)
  best.params
  best.nse <- -trial.optim$value
  print(best.nse)
  
  target.file[row.numpars,2] <- best.params
  setwd(target.dir)
  write.table(x=target.file, file=paste(getwd(),"/speciesOptim.txt",sep=""), quote=FALSE, row.names=FALSE, col.names = FALSE)

  #run parameterization on test set
  par8 <- readChar(paste(pr.dir, "/input/",j,"/testMet.txt", sep=""), file.info(paste(pr.dir, "/input/",j,"/testMet.txt", sep=""))$size)
  running.string <- paste(par1, par2, par3,par4,par5,par6,par7,par8,par9,par10, sep="")
  setwd(curr.dir)
  system(running.string)
  
} 



