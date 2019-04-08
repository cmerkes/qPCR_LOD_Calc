## Before getting started, please read the notes and personalize lines 17, 20, 23, 26,
##   and 33.  There are interactive codes, it is designed to be run using R studio 
##   where you can click on "source" after modifying lines 17, 20, 23, 26, and 33.

## This is a script to calculate Limit of Detection for eDNA assays.
## It requires you to assemble your data for many replicates of a dilution
##   curve of low-copy standards into a single csv file with the following headings:
##   Target, Cq, SQ (Note: headings are not case-sensitive, but must be spelled the same.)
##   Additional columns may be included for readability but are not used in calculations.
##   SQ (Starting Quantity) should be the expected values for your standards.

## Load packages:
library(ggplot2)
library(drc)

## Set your working directory to where your csv file is saved (MODIFY AS NEEDED):
setwd("C:/Users/cmerkes/Desktop/Data/LoD/")
#setwd("C:/Users/cmerkes/Desktop/Data/LoD/emy_data/")

## Read in your data file (MODIFY FILE NAME AS NEEDED):
DAT <- read.csv("GEDWG_LOD_DATA3.csv")
#DAT <- read.csv("Helbing-LOD_LOQ_Data-022019CCH.csv")
#DAT <- read.csv("Emy-data.csv")

## Define your CV threshold for LoQ:
LOQ.Threshold <- 0.35

## Define which logarithmic function to use for LoD model:
LOD.FCT <- "Best"
## Selecting "Best" will signal the code to automatically select the best fitting
##   model choice. Run the function getMeanFunctions() to print the list of all choices.
## Example: LOD.FCT <- W2.4()
##   This example will use the Weibull type II, 4 parameter function.

## Define which model to use for LoQ model:
LOQ.FCT <- "Best"
## Selecting "Best" will signal the code to automatically select the model with lowest
##   residual standard error. Change to "Decay" to use exponential decay model, "Linear"
##   to use linear model, "Pn" to use an nth-order polynomial model where n is numerical.
##   Example: "P2" will use a 2nd order polynomial model, "P3" will use 3rd order, and etc.
##   Selecting "Best" will test polynomial models up to 6th-order.

## Create an analysis log file:
write(paste0("Analysis started: ",date(),"\n\n"),file="Analysis Log.txt")

## Check the data:
if(sum(colnames(DAT)=="Target")!=1) { #Is there a "Target" column?
  A <- grep("target",colnames(DAT),ignore.case=T)
  if(length(A)==1) { colnames(DAT)[A] <- "Target" } #Rename target column if it is mispelled but can be identified and there is only 1.
  if(length(A)!=1) { write("There is a problem with the 'Target' column.\n\n",file="Analysis Log.txt",append=T) } #Add error message to analysis log.
  if(length(A)>1) { cat("ERROR: multiple 'Target' columns detected.",colnames(DAT)[A],sep="\n") }
  if(length(A)==0) { print("ERROR: cannot detect 'Target' column.") }
}
if(sum(colnames(DAT)=="Cq")!=1) { #Is there a "Cq" column?
  A <- grep("cq|ct|cycle",colnames(DAT),ignore.case=T)
  if(length(A)==1) { colnames(DAT)[A] <- "Cq" } #Rename cq column if it is mispelled but can be identified and there is only 1.
  if(length(A)!=1) { write("There is a problem with the 'Cq' column.\n\n",file="Analysis Log.txt",append=T) } #Add error message to analysis log.
  if(length(A)>1) { cat("ERROR: multiple 'Cq' columns detected.",colnames(DAT)[A],sep="\n") }
  if(length(A)==0) { print("ERROR: cannot detect 'Cq' column.") }
}
if(sum(colnames(DAT)=="SQ")!=1) { #Is there a "SQ" column?
  A <- grep("sq|copies|starting|quantity",colnames(DAT),ignore.case=T)
  if(length(A)==1) { colnames(DAT)[A] <- "SQ" } #Rename SQ column if it is mispelled but can be identified and there is only 1.
  if(length(A)!=1) { write("There is a problem with the 'SQ' column.\n\n",file="Analysis Log.txt",append=T) } #Add error message to analysis log.
  if(length(A)>1) { cat("ERROR: multiple 'SQ' columns detected.",colnames(DAT)[A],sep="\n") }
  if(length(A)==0) { print("ERROR: cannot detect 'SQ' column.") }
}

## Ensure data is in the proper format:
DAT$Target <- as.factor(DAT$Target)
DAT$Cq <- suppressWarnings(as.numeric(as.character(DAT$Cq))) #Non-numerical values (i.e. negative wells) will be converted to NAs
DAT$SQ <- suppressWarnings(as.numeric(as.character(DAT$SQ))) #Non-numerical values (i.e. NTC) will be converted to NAs
if(sum(is.na(DAT$SQ))>0) {
  write(paste0("WARNING: ",sum(is.na(DAT$SQ))," data points excluded without a valid starting quantity (SQ)!\nHere is a sample of the data being excluded:\n"),
        file="Analysis Log.txt",append=T)
  suppressWarnings(write.table(head(DAT[is.na(DAT$SQ),]),file="Analysis Log.txt",append=T,
                               sep="\t",eol="\n",row.names=F,col.names=T))
  write("\n",file="Analysis Log.txt",append=T)
  print(paste0("WARNING: ",sum(is.na(DAT$SQ))," data points excluded without a valid starting quantity (SQ)!"))
  print(head(DAT[is.na(DAT$SQ),]))
}

## Check for wild outliers that the user should go back and review:
Targets <- unique(DAT$Target)
## Get matchups of all standards and markers used:
for(i in 1:length(Targets)) {
  if(i==1) {
    Standards <- unique(DAT$SQ[DAT$Target==Targets[i]&!is.na(DAT$SQ)])
    Target <- rep(as.character(Targets[i]),length(Standards))
  }
  else {
    Standards <- c(Standards,unique(DAT$SQ[DAT$Target==Targets[i]&!is.na(DAT$SQ)]))
    Target <- c(Target,rep(as.character(Targets[i]),
                           length(unique(DAT$SQ[DAT$Target==Targets[i]&!is.na(DAT$SQ)]))))
  }
}
OUTS <- data.frame(Target=Target,Standard=Standards,Outliers=NA)
## Identify any wells where the Cq value is more than 10% away from the median for
##   that standard.
for(i in 1:nrow(OUTS)) {
  MED <- median(DAT$Cq[DAT$SQ==OUTS$Standard[i]&DAT$Target==OUTS$Target[i]],na.rm=T)
  A <- which(DAT$SQ==OUTS$Standard[i]&DAT$Target==OUTS$Target[i]&DAT$Cq<0.9*MED&!is.na(DAT$Cq))
  B <- which(DAT$SQ==OUTS$Standard[i]&DAT$Target==OUTS$Target[i]&DAT$Cq>1.1*MED&!is.na(DAT$Cq))
  if(length(c(A,B))>0) {
    OUTS$Outliers[i] <- paste(c(A,B),collapse=",")
  }
}
## If any outliers are detected, export the raw data as csv and make a note in
##   the analysis log.
if(sum(!is.na(OUTS$Outliers))>0) {
  OUT.ROW <- paste(OUTS$Outliers[!is.na(OUTS$Outliers)],collapse=",")
  OUT.ROW2 <- unlist(strsplit(OUT.ROW,split=","))
  write.csv(DAT[OUT.ROW2,],file="Potential-Outliers.csv",row.names=F)
  write("Potential outliers have been detected. Please review the data exported as
Potential-Outliers.csv, and determine if any data points need to be excluded
or adjusted due to false positives or poorly normalized baselines.",
        file="Analysis Log.txt",append=T)
  write("\n",file="Analysis Log.txt",append=T)
}

## Generate standard curves using all data and calculate copy estimates for each
##   replicate using the curves:
curve.list <- ""
DAT$Copy.Estimate <- rep(NA,nrow(DAT))
DAT$Mod <- rep(0,nrow(DAT))
for(i in 1:length(Targets)) {
  STDS <- data.frame(S=unique(DAT$SQ[DAT$Target==Targets[i]]),R=NA)
  ## Calculate detection rates for each standard:
  for(j in 1:nrow(STDS)) {
    STDS$R[j] <- sum(!is.na(DAT$Cq)&DAT$SQ==STDS$S[j]&DAT$Target==Targets[i],na.rm=T)/sum(DAT$SQ==STDS$S[j]&DAT$Target==Targets[i],na.rm=T)
  }
  ## Only use standards with 50% or greater detection rates for linear regression:
  if(sum(STDS$R>=0.5,na.rm=T)>2) {
    STDS2 <- STDS$S[STDS$R>=0.5&!is.na(STDS$R)&!is.na(STDS$S)]
  }
  ## If there are not at least 3 standards with 50% or greater detection, use the top 3:
  if(sum(STDS$R>=0.5,na.rm=T)<3) {
    STDS2 <- STDS$S[order(STDS$R,decreasing=T)][1:3]
  }
  ## Identify the 2nd and 3rd quartiles of each used standard for inclusion in the
  ##   standard curve calculations
  for(j in 1:length(STDS2)) {
    D <- DAT$Cq[DAT$Target==Targets[i]&DAT$SQ==STDS2[j]]
    DAT$Mod[DAT$Target==Targets[i]&DAT$SQ==STDS2[j]&DAT$Cq>=quantile(D,na.rm=T)[2]&DAT$Cq<=quantile(D,na.rm=T)[4]&!is.na(DAT$SQ)] <- 1
  }
  if(length(unique(DAT$SQ[DAT$Target==Targets[i]]))!=length(STDS2)) {
    ToWrite <- paste0("These standards not included in ",Targets[i],
                      " standard curve regression for copy estimate calculations, because they detected below 50%: ",
                      paste(setdiff(unique(DAT$SQ[DAT$Target==Targets[i]]),STDS2),collapse=", "),"\n\n")
    write(ToWrite,file="Analysis Log.txt",append=T)
  }
  assign(paste0("curve",i),lm(Cq~log10(SQ),data=DAT[DAT$Target==Targets[i]&DAT$Mod==1,]))
  curve.list <- c(curve.list,paste0("curve",i))
  Intercept <- coef(get(curve.list[i+1]))[1]
  Slope <- coef(get(curve.list[i+1]))[2]
  DAT$Copy.Estimate[DAT$Target==Targets[i]] <- 10^((DAT$Cq[DAT$Target==Targets[i]]-Intercept)/Slope)
}

## Summarize the data:
DAT2 <- data.frame(Standards=Standards,Target=Target,Reps=NA,Detects=NA,Cq.mean=NA,
                   Cq.sd=NA,Copy.CV=NA,Cq.CV=NA)
## Fill in replicate counts, positive detect counts, mean Cq values, standard
##   deviations of Cq values, and coefficient of variation of copy estimates for
##   each standard and marker combination:
for(i in 1:nrow(DAT2)) {
  DAT2$Reps[i] <- sum(DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i],na.rm=T)
  DAT2$Detects[i] <- sum(!is.na(DAT$Cq)&DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i],na.rm=T)
  DAT2$Cq.mean[i] <- mean(DAT$Cq[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]],na.rm=T)
  DAT2$Cq.sd[i] <- sd(DAT$Cq[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]],na.rm=T)
  DAT2$Copy.CV[i] <- sd(DAT$Copy.Estimate[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]],na.rm=T)/mean(DAT$Copy.Estimate[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]],na.rm=T)
  DAT2$Cq.CV[i] <- sqrt(2^(DAT2$Cq.sd[i]^2*log(2))-1)
}
## Calculate positive detection rate for each standard and marker combination:
DAT2$Rate <- DAT2$Detects/DAT2$Reps
write("Data Summary:",file="Analysis Log.txt",append=T)
suppressWarnings(write.table(DAT2,file="Analysis Log.txt",append=T,sep="\t",eol="\n",
                             row.names=F,col.names=T))

## Determine the lowest standard with 95% or greater detection:
for(i in 1:length(Targets)) {
  A <- min(DAT2$Standards[DAT2$Rate>=0.95&DAT2$Target==Targets[i]])
  ToWrite <- paste0("For ",Targets[i],", the lowest standard with 95% or greater detection is: ",A," copies/reaction.")
  ToWrite2 <- ""
  if(length(which(DAT2$Rate<0.95&DAT2$Target==Targets[i]))>0) {
    B <- max(DAT2$Standards[DAT2$Rate<0.95&DAT2$Target==Targets[i]])
    if(B>A) {
      ToWrite2 <- paste0("WARNING: For ",Targets[i],", ",B," copies/reaction standard detected at lower rate than ",A," copies/reaction standard.\nPlease retest.")
    }
  }
  if(length(which(DAT2$Rate<0.95&DAT2$Target==Targets[i]))==0) {
    ToWrite2 <- paste0("WARNING: LoD cannot be determined for ",Targets[i],", because it is lower than the lowest standard you tested.\nReport as <",A," copies/reaction, or retest with lower concentrations.")
  }
  write(paste0("\n\n",ToWrite,"\n"),file="Analysis Log.txt",append=T)
  if(ToWrite2!="") { write(paste0(ToWrite2,"\n\n"),file="Analysis Log.txt",append=T) }
  cat(ToWrite,ToWrite2,sep="\n")
}

## Determine LoD and LoQ by modeling, and summarize each assay:
## NOTE: LoD is now determined by dose-response modeling. Probit modeling code remains,
##         but has been converted to comments.
DAT$Detect <- as.numeric(!is.na(DAT$Cq))
#LOD.list <- ""
LOD.list2 <- ""
LOD.list3 <- ""
LOQ.list <- ""
DAT3 <- data.frame(Assay=Targets,R.squared=NA,Slope=NA,Intercept=NA,Low.95=NA,
                   LOD=NA,LOQ=NA,rep2.LOD=NA,rep3.LOD=NA,rep4.LOD=NA,rep5.LOD=NA,rep8.LOD=NA)
LOD.FCTS <- list(LL.2(),LL.3(),LL.3u(),LL.4(),LL.5(),W1.2(),W1.3(),W1.4(),W2.2(),W2.3(),
                 W2.4(),AR.2(),AR.3(),MM.2(),MM.3())
for(i in 1:length(Targets)) {
  ## Check input suitability for probit or dose-response modeling:
  if(sum(DAT2$Rate[DAT2$Target==Targets[i]]!=1&DAT2$Rate[DAT2$Target==Targets[i]]!=0)==0) {
    ToWrite <- paste0("WARNING: For ",Targets[i],", all standards detected fully or failed fully.  Therefore, the LoD model will not converge.")
    write(paste0(ToWrite,"\n\n"),file="Analysis Log.txt",append=T)
    print(ToWrite)
  }
  if(sum(DAT2$Rate[DAT2$Target==Targets[i]]!=1&DAT2$Rate[DAT2$Target==Targets[i]]!=0)==1) {
    ToWrite <- paste0("WARNING: For ",Targets[i],", only 1 standard detected in the informative range (not 0% and not 100%).  Therefore, the LoD model results will be less reliable.")
    write(paste0(ToWrite,"\n\n"),file="Analysis Log.txt",append=T)
    print(ToWrite)
  }
  ## Define probit model:
  #assign(paste0("LOD.mod",i),glm(Detect~SQ,data=DAT[DAT$Target==Targets[i],],
  #                               family=binomial(link="probit")))
  #LOD.list <- c(LOD.list,paste0("LOD.mod",i))
  ## Define LOQ model using lowest residual standard error selection:
  if(LOQ.FCT=="Best") {
    ## Remove previous marker LOQ models from environment if they exist:
    suppressWarnings(rm(LOQ1,LOQ2,LOQ3,LOQ4,LOQ5,LOQ6,LOQ7))
    tryCatch({ #skip if model cannot be determined.
      LOQ1 <- nls(Cq.CV~SSasymp(log10(Standards),Asym,R0,lrc),
                  data=DAT2[DAT2$Target==Targets[i],])
    }, error=function(e) {
      e
      cat("ERROR: decay LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
    tryCatch({ #skip if model cannot be determined.
      LOQ2 <- lm(Cq.CV~log10(Standards),data=DAT2[DAT2$Target==Targets[i],])
    }, error=function(e) {
      e
      cat("ERROR: linear LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
    tryCatch({ #skip if model cannot be determined.
      LOQ3 <- lm(Cq.CV~poly(log10(Standards),2),data=DAT2[DAT2$Target==Targets[i],])
    }, error=function(e) {
      e
      cat("ERROR: 2nd polynomial LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
    tryCatch({ #skip if model cannot be determined.
      LOQ4 <- lm(Cq.CV~poly(log10(Standards),3),data=DAT2[DAT2$Target==Targets[i],])
    }, error=function(e) {
      e
      cat("ERROR: 3rd polynomial LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
    tryCatch({ #skip if model cannot be determined.
      LOQ5 <- lm(Cq.CV~poly(log10(Standards),4),data=DAT2[DAT2$Target==Targets[i],])
    }, error=function(e) {
      e
      cat("ERROR: 4th polynomial LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
    tryCatch({ #skip if model cannot be determined.
      LOQ6 <- lm(Cq.CV~poly(log10(Standards),5),data=DAT2[DAT2$Target==Targets[i],])
    }, error=function(e) {
      e
      cat("ERROR: 5th polynomial LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
    tryCatch({ #skip if model cannot be determined.
      LOQ7 <- lm(Cq.CV~poly(log10(Standards),6),data=DAT2[DAT2$Target==Targets[i],])
    }, error=function(e) {
      e
      cat("ERROR: 6th polynomial LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
    ## Determine which models were able to be determined:
    A <- sapply(c("LOQ1","LOQ2","LOQ3","LOQ4","LOQ5","LOQ6","LOQ7"),exists)
    B <- names(A)[A==T]
    ## If at least 1 LOQ model was determined, select the one with the lowest
    ##   residual standard error:
    if(length(B)>0) {
      LOQ.res <- rep(NA,length(B))
      for(j in 1:length(B)) {
        LOQ.res[j] <- summary(get(B[j]))$sigma
      }
      C <- which(LOQ.res==min(LOQ.res,na.rm=T))
      assign(paste0("LOQ.mod",i),get(B[C]))
      LOQ.list <- c(LOQ.list,paste0("LOQ.mod",i))
    }
  }
  ## Define LOQ model by exponential decay modeling:
  if(LOQ.FCT=="Decay") {
    tryCatch({ #skip if model cannot be determined.
      assign(paste0("LOQ.mod",i),nls(Cq.CV~SSasymp(log10(Standards),Asym,R0,lrc),
                  data=DAT2[DAT2$Target==Targets[i],]))
      LOQ.list <- c(LOQ.list,paste0("LOQ.mod",i))
    }, error=function(e) {
      e
      cat("ERROR: decay LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
  }
  ## Define LOQ model by linear modeling:
  if(LOQ.FCT=="Linear") {
    tryCatch({ #skip if model cannot be determined.
      assign(paste0("LOQ.mod",i),lm(Cq.CV~log10(Standards),
                                    data=DAT2[DAT2$Target==Targets[i],]))
      LOQ.list <- c(LOQ.list,paste0("LOQ.mod",i))
    }, error=function(e) {
      e
      cat("ERROR: linear LOQ model cannot be defined for ",as.character(Targets[i]),sep="")
    })
  }
  ## Define LOQ model by polynomial modeling:
  if(substr(LOQ.FCT,1,1)=="P") {
    Z <- as.numeric(substr(LOQ.FCT,2,nchar(LOQ.FCT)))
    tryCatch({ #skip if model cannot be determined.
      assign(paste0("LOQ.mod",i),lm(Cq.CV~poly(log10(Standards),Z),
                                    data=DAT2[DAT2$Target==Targets[i],]))
      LOQ.list <- c(LOQ.list,paste0("LOQ.mod",i))
    }, error=function(e) {
      e
      cat("ERROR: ",Z,"-order polynomial LOQ model cannot be defined for ",
          as.character(Targets[i]),sep="")
    })
  }
  ## Signal undetermined model with NA:
  if(length(LOQ.list)<i+1) {
    LOQ.list <- c(LOQ.list,NA)
  }
  ## Define the logarithmic model for LOD using user-selected function:
  if(is.list(LOD.FCT)==T) {
    tryCatch({ #skip if model cannot be determined.
      assign(paste0("LOD.mod2",i),drm(Detect~SQ,data=DAT[DAT$Target==Targets[i],],fct=LOD.FCT))
      LOD.list2 <- c(LOD.list2,paste0("LOD.mod2",i))
      LOD.list3 <- c(LOD.list3,LOD.FCT$name)
    }, error=function(e) {
      e
      cat("ERROR: LOD model cannot be defined for ",as.character(Targets[i]),sep="")
    })
  }
  ## Define the logarithmic model with function automatically selected:
  if(is.character(LOD.FCT)) {
    if(LOD.FCT=="Best") {
      tryCatch({ #skip if model cannot be determined.
        ## Pull out data for specific assay:
        TEMP.DAT <- DAT[DAT$Target==Targets[i],]
        ## Define a model to start with:
        LOD.mod <- drm(Detect~SQ,data=TEMP.DAT,fct=W2.4())
        ## Test all available models and select the best one:
        LOD.FCT2 <- row.names(mselect(LOD.mod,LOD.FCTS))[1]
        LOD.FCT3 <- getMeanFunctions(fname=LOD.FCT2)
        assign(paste0("LOD.mod2",i),drm(Detect~SQ,data=DAT[DAT$Target==Targets[i],],fct=LOD.FCT3[[1]]))
        LOD.list2 <- c(LOD.list2,paste0("LOD.mod2",i))
        LOD.list3 <- c(LOD.list3,LOD.FCT2)
      }, error=function(e) {
        e
        cat("ERROR: LOD model cannot be defined for ",as.character(Targets[i]),sep="")
      })
    }
  }
  ## Signal undetermined model with NA:
  if(length(LOD.list2)<i+1) {
    LOD.list2 <- c(LOD.list2,NA)
    LOD.list3 <- c(LOD.list3,NA)
  }
  ## Populate summary data:
  DAT3$R.squared[i] <- summary(get(curve.list[i+1]))$r.squared
  DAT3$Slope[i] <- coef(get(curve.list[i+1]))[2]
  DAT3$Intercept[i] <- coef(get(curve.list[i+1]))[1]
  DAT3$Low.95[i] <- min(DAT2$Standards[DAT2$Rate>=0.95&DAT2$Target==Targets[i]])
  ## Only get LOD values if the LOD model is defined:
  if(!is.na(LOD.list2[i+1])) {
    DAT3$LOD[i] <- ED(get(LOD.list2[i+1]),0.95,type="absolute")[1]
    DAT3$rep2.LOD[i] <- ED(get(LOD.list2[i+1]),1-sqrt(0.05),type="absolute")[1]
    DAT3$rep3.LOD[i] <- ED(get(LOD.list2[i+1]),1-0.05^(1/3),type="absolute")[1]
    DAT3$rep4.LOD[i] <- ED(get(LOD.list2[i+1]),1-0.05^0.25,type="absolute")[1]
    DAT3$rep5.LOD[i] <- ED(get(LOD.list2[i+1]),1-0.05^0.2,type="absolute")[1]
    DAT3$rep8.LOD[i] <- ED(get(LOD.list2[i+1]),1-0.05^0.125,type="absolute")[1]
    ## Residual code using probit method:
    #DAT3$LOD[i] <- (qnorm(0.95)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
    #DAT3$rep2.LOD[i] <- (qnorm(0.50)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
    #DAT3$rep3.LOD[i] <- (qnorm(1/3)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
    #DAT3$rep4.LOD[i] <- (qnorm(0.25)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
    #DAT3$rep5.LOD[i] <- (qnorm(0.2)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
    #DAT3$rep8.LOD[i] <- (qnorm(0.125)-coef(get(LOD.list[i+1]))[1])/coef(get(LOD.list[i+1]))[2]
  }
  ## Generate prediction data for LoQ:
  ## Only get LOQ if LOQ model is determined:
  if(!is.na(LOQ.list[i+1])) {
    newData <- data.frame(Standards = seq(1, 10000))
    newData$Cq.CV <- predict(get(LOQ.list[i+1]), newData)
    ## Determine what type of LOQ model is used and calculate LOQ accordingly:
    ## For exponential decay:
    if(as.character(get(LOQ.list[i+1])$call)[1]=="nls") {
      ## Look up lowest modeled standard below the CV threshold:
      DAT3$LOQ[i] <- min(newData$Standards[newData$Cq.CV<=LOQ.Threshold])
      ## Unless... If the background variation exceeds the CV threshold, adjust threshold:
      ## Determine the highest standard used:
      A <- max(DAT2$Standards[DAT2$Target==Targets[i]])
      if(min(newData$Cq.CV[newData$Standards<=A])>LOQ.Threshold) {
        ## Set the adjusted threshold to 1.5x the lowest simulated Cq.CV 
        ##   within the range of data tested:
        B <- min(newData$Cq.CV[newData$Standards<=A])
        DAT3$LOQ[i] <- min(newData$Standards[newData$Cq.CV<=B*1.5])
        ## Make a note of the adjusted threshold in the analysis log:
        ToWrite <- paste0("Note: All standards tested for ",Targets[i],
                          " yielded higher Cq.CV values than the user-defined CV threshold of ",
                          LOQ.Threshold,". The CV threshold has been adjusted to ",
                          B*1.5," for the LOQ of this marker.")
        write(paste0(ToWrite,"\n\n"),file="Analysis Log.txt",append=T)
      }
    }
    if(as.character(get(LOQ.list[i+1])$call)[1]=="lm") {
      ## For polynomial:
      if(grepl("poly",as.character(get(LOQ.list[i+1])$call)[2])==T) {
        ## Determine the highest standard used:
        A <- max(DAT2$Standards[DAT2$Target==Targets[i]])
        ## Adjust if the tested range does not cross below the CV threshold:
        if(min(DAT2$Cq.CV[DAT2$Target==Targets[i]],na.rm=T)>LOQ.Threshold) {
          B <- min(DAT2$Cq.CV[DAT2$Target==Targets[i]],na.rm=T)*1.5
          ## Make a note of the adjusted threshold in the analysis log:
          ToWrite <- paste0("Note: All standards tested for ",Targets[i],
                            " yielded higher Cq.CV values than the user-defined CV threshold of ",
                            LOQ.Threshold,". The CV threshold has been adjusted to ",
                            B," for the LOQ of this marker.")
          write(paste0(ToWrite,"\n\n"),file="Analysis Log.txt",append=T)
        }
        else {
          B <- LOQ.Threshold
        }
        ## Look up highest modeled standard below the CV threshold:
        C <- max(newData$Standards[newData$Cq.CV<=B&newData$Standards<=A])
        ## Look up the highest modeled standard above the CV threshold...
        ##   and also below the highest standard below the CV threshold.
        ##   This captures the farthest right crossing point on a downward slope.
        D <- max(newData$Standards[newData$Cq.CV>B&newData$Standards<C])
        ## LOQ is D + 1 to get back less than or equal to the CV threshold.
        DAT3$LOQ[i] <- D+1
      }
      # For linear:
      else {
        ## Look up lowest modeled standard below the CV threshold:
        DAT3$LOQ[i] <- min(newData$Standards[newData$Cq.CV<=LOQ.Threshold])
      }
    }
    ## If modeled LOQ is calculated to be below the 95% LOD, set LOD as LOQ:
    if(is.na(DAT3$LOD[i])==F) {
      if(DAT3$LOQ[i]<DAT3$LOD[i]) {
        DAT3$LOQ[i] <- DAT3$LOD[i]
      }
    }
    ## If modeled LOQ is calculated to be below the lowest standard tested,
    ##   set the lowest standard as the LOQ:
    if(DAT3$LOQ[i]<min(DAT2$Standards[DAT2$Target==Targets[i]])) {
      DAT3$LOQ[i] <- min(DAT2$Standards[DAT2$Target==Targets[i]])
    }
  }
}
write("Assay summary:",file="Analysis Log.txt",append=T)
write("\nR.squared: The R-squared value of linear regression of all standards Cq-values vs log10 of the starting quantities.",file="Analysis Log.txt",append=T)
write("Slope: The slope of the linear regression.",file="Analysis Log.txt",append=T)
write("Intercept: The y-intercept of the linear regression.",file="Analysis Log.txt",append=T)
write("\nLow.95: The lowest standard with at least 95% positive detection.",file="Analysis Log.txt",append=T)
write("LOD: The 95% limit of detection as determined by probit modeling.",file="Analysis Log.txt",append=T)
write(paste0("LOQ: The limit of quantification as determined by decay modeling, using the user-selected CV threshold of: ",LOQ.Threshold),file="Analysis Log.txt",append=T)
write("\nrep2.LOD: The effective limit of detection if analyzing in 2 replicates.",file="Analysis Log.txt",append=T)
write("rep3.LOD: The effective limit of detection if analyzing in 3 replicates.",file="Analysis Log.txt",append=T)
write("rep4.LOD: The effective limit of detection if analyzing in 4 replicates.",file="Analysis Log.txt",append=T)
write("rep5.LOD: The effective limit of detection if analyzing in 5 replicates.",file="Analysis Log.txt",append=T)
write("rep8.LOD: The effective limit of detection if analyzing in 8 replicates.\n\n",file="Analysis Log.txt",append=T)
write.csv(DAT3,file="Assay summary.csv",row.names=F)

## Plot Cq value vs Standard Concentration standard curves:
DAT$Mod[DAT$Mod==0] <- "Excluded"
DAT$Mod[DAT$Mod==1] <- "Modeled"
for(i in 1:length(Targets)) {
  ggOut <- ggplot(data=DAT[DAT$Target==Targets[i]&is.na(DAT$SQ)==F,],
                  aes(x=SQ,y=Cq,color=factor(Mod),shape=factor(Mod),size=factor(Mod))) + 
    geom_jitter(width=0.1,alpha=0.75) + 
    scale_shape_manual("",values=c(3,20),guide=F) +
    scale_size_manual("",values=c(1,3)) +
    scale_x_log10() +
    scale_color_manual("",values=c("blue", "black")) +
    xlab("Standard Concentrations (Copies / Reaction)") +
    ylab("Cq-value") +
    geom_abline(intercept=coef(get(curve.list[i+1]))[1],
                slope=coef(get(curve.list[i+1]))[2]) +
    geom_vline(xintercept=DAT3$LOD[i],linetype=2) +
    geom_vline(xintercept=DAT3$LOQ[i],colour="red") +
    annotate("text",y=max(DAT$Cq[DAT$Target==Targets[i]],na.rm=T)*0.99,
             x=DAT3$LOD[i]*0.8,angle=90,label="LOD") +
    annotate("text",y=max(DAT$Cq[DAT$Target==Targets[i]],na.rm=T)*0.94,color="red",
             x=DAT3$LOQ[i]*0.8,angle=90,label="LOQ") +
    theme_bw() + theme(legend.justification=c(1,1),legend.position=c(1,0.99)) +
    ggtitle(paste0("Standard curve for: ",Targets[i])) +
    theme(plot.title=element_text(hjust=0.5,size=20),
          axis.title=element_text(size=16)) +
    theme(legend.title=element_blank(),
          legend.text=element_text(size=11)) +
    annotate("text",y=min(DAT$Cq[DAT$Target==Targets[i]&is.na(DAT$SQ)==F],na.rm=T)*1.05,
             x=min(DAT$SQ[DAT$Target==Targets[i]&is.na(DAT$SQ)==F],na.rm=T)*1.01,hjust=0,
             label=(paste0("R-squared: ",DAT3$R.squared[i],"\ny = ",DAT3$Slope[i],"x + ",DAT3$Intercept[i])))
  print(ggOut)
  readline(prompt="Press [Enter] for next plot.")
  print("Calculating... Please wait.")
}

## Plot the LOD models for each assay:
for(i in 1:length(Targets)) {
  if(!is.na(LOD.list2[i+1])) {
    DAT4 <- rbind(ED(get(LOD.list2[i+1]),0.95,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-sqrt(0.05),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^(1/3),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.25,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.2,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.125,interval="delta",type="absolute"))
    if(substr(LOD.list3[i+1],1,3)=="LL2") {
      DAT4 <- exp(DAT4)
    }
    DAT4 <- data.frame(DAT4,LoD=c("1rep.LOD","2rep.LOD","3rep.LOD","4rep.LOD",
                                  "5rep.LOD","8rep.LOD"),
                       Assay=rep(Targets[i],nrow(DAT4)))
    DAT4$Assay <- as.character(DAT4$Assay)
    if(i==1) {
      LOD.CI <- DAT4
    }
    if(i>1) {
      LOD.CI <- rbind(LOD.CI,DAT4)
    }
    if(sum(!is.na(DAT4[,3])&DAT4[,3]<=0)>0) { #Unable to plot negative lower limits, converting any lower limit values to 0.0001
      DAT4[!is.na(DAT4[,3])&DAT4[,3]<=0,3] <- 0.0001
    }
    plot(get(LOD.list2[i+1]),main=paste0("LoD Plot for: ",Targets[i]),
         ylab="Detection Probability",xlab="Standard concentrations (Copies / Reaction)",
         xlim=c(min(DAT4[,1:4],na.rm=T),max(DAT$SQ,na.rm=T)))
    LODS <- sum(!is.na(DAT4[,1]))
    COLS <- c(rgb(0.8,0.47,0.65),rgb(0,0.45,0.7),rgb(0.94,0.89,0.26),
              rgb(0.84,0.37,0),rgb(0,0.62,0.45),rgb(0.90,0.62,0))
    PNTS <- c(15,16,17,18,25,3)
    YS <- c(0.95,1-sqrt(0.05),1-0.05^(1/3),1-0.05^0.25,1-0.05^0.2,1-0.05^0.125)
    LODS2 <- c("Limit of Detection","2 Replicates LoD","3 Replicates LoD",
               "4 Replicates LoD","5 Replicates LoD","8 Replicates LoD")
    if(LODS<6) {
      LODS2[(LODS+1):6] <- gsub("licates LoD","s: Insufficient Data",LODS2[(LODS+1):6])
    }
    points(x=DAT4[1:LODS,1],y=YS[1:LODS],pch=PNTS,col=COLS,cex=1.2)
    for(j in 1:LODS) {
      lines(x=DAT4[j,3:4],y=rep(YS[j],2),col=COLS[j],lwd=2)
      lines(x=rep(DAT4[j,3],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=2,col=COLS[j])
      lines(x=rep(DAT4[j,4],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=2,col=COLS[j])
    }
    legend("bottomright",legend=LODS2,pch=PNTS,col=COLS,text.col=COLS)
    Pval <- modelFit(get(LOD.list2[i+1]))[[5]][2]
    mtext(paste0("FCT used: ",LOD.list3[i+1],"    Lack of fit test: p = ",Pval),side=3)
  }
  if(is.na(LOD.list2[i+1])) {
    plot(DAT2$Rate[DAT2$Target==Targets[i]]~log10(DAT2$Standards[DAT2$Target==Targets[i]]),
         ylim=c(0,1),ylab="Detection Probability",
         xlab=expression("Log of standard concentrations (Log"[10]*"Copies / Reaction)"),
         main=paste0("LoD for: ",Targets[i]," unsolvable"))
  }
  readline(prompt="Press [Enter] for next plot.")
  print("Calculating... Please wait.")
}
LOD.CI <- LOD.CI[,c(6,5,1,3,4,2)]
write.csv(LOD.CI,file="LOD_confint.csv",row.names=F)

## Plot the LoQ models of each assay:
for(i in 1:length(Targets)) {
  if(is.na(LOQ.list[i+1])==F) {
    ## Re-generate prediction data for the model:
    newData <- data.frame(Standards = seq(1, 10000))
    newData$Cq.CV <- predict(get(LOQ.list[i+1]), newData)
    ## Define LOQ polygon coordinates:
    PDAT <- data.frame(x=c(min(DAT2$Standards[DAT2$Target==Targets[i]],na.rm=T),
                           min(DAT2$Standards[DAT2$Target==Targets[i]],na.rm=T),
                           DAT3$LOQ[DAT3$Assay==Targets[i]],
                           DAT3$LOQ[DAT3$Assay==Targets[i]]),
                       y=c(min(c(DAT2$Cq.CV[DAT2$Target==Targets[i]],newData$Cq.CV[newData$Standards<=max(DAT2$Standards[DAT2$Target==Targets[i]])&newData$Standards>=min(DAT2$Standards[DAT2$Target==Targets[i]])]),na.rm=T)*0.9,
                           newData$Cq.CV[newData$Standards==DAT3$LOQ[DAT3$Assay==Targets[i]]],
                           newData$Cq.CV[newData$Standards==DAT3$LOQ[DAT3$Assay==Targets[i]]],
                           min(c(DAT2$Cq.CV[DAT2$Target==Targets[i]],newData$Cq.CV[newData$Standards<=max(DAT2$Standards[DAT2$Target==Targets[i]])&newData$Standards>=min(DAT2$Standards[DAT2$Target==Targets[i]])]),na.rm=T)*0.9))
    if(DAT3$LOQ[DAT3$Assay==Targets[i]]!=floor(DAT3$LOQ[DAT3$Assay==Targets[i]])) {
      PDAT$y[2:3] <- LOQ.Threshold
    }
  }

  Decay.Plot <- ggplot(DAT2[DAT2$Target==Targets[i],], aes(x= Standards, y = Cq.CV)) +
    geom_point(size=2) +
    scale_x_continuous(trans = 'log10') +
    ylab("Coefficient of variation for Cq-Values") +
    xlab("Standard concentrations (Copies / Reaction)") +
    geom_vline(xintercept=DAT3$LOD[DAT3$Assay==Targets[i]],color="red") +
    annotate("text",y=max(DAT2$Cq.CV[DAT2$Target==Targets[i]],na.rm=T)*0.99,
             x=DAT3$LOD[i]*0.8,angle=90,label="LOD",color="red") +
    theme(legend.position="none") +
    theme(plot.title=element_text(hjust=0.5))
  
  if(is.na(LOQ.list[i+1])==F) {
    if(DAT3$LOQ[DAT3$Assay==Targets[i]]<=min(DAT2$Standards[DAT2$Target==Targets[i]])) {
      PDAT$x[3:4] <- NA
      Decay.Plot <- Decay.Plot + 
        annotate("text",y=max(DAT2$Cq.CV[DAT2$Target==Targets[i]],na.rm=T)*0.99,
                 x=median(DAT2$Standards[DAT2$Target==Targets[i]]),
                 label="LOQ may be outside tested range.",hjust=0)
    }
    Decay.Plot <- Decay.Plot + geom_polygon(data=PDAT,aes(x=x,y=y,alpha=0.5))

    if(as.character(get(LOQ.list[i+1])$call)[1]=="nls") {
      Decay.Plot <- Decay.Plot + 
        stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE) +
        ggtitle(paste0("Exponential Decay LOQ model for: ",Targets[i]))
    }

    if(as.character(get(LOQ.list[i+1])$call)[1]=="lm") {
      if(grepl("poly",as.character(get(LOQ.list[i+1])$call)[2])==T) {
        B <- length(get(LOQ.list[i+1])$coefficients)-1
        Decay.Plot <- Decay.Plot +
          stat_smooth(method = "lm", formula = y ~ poly(x,B),se=F) +
          ggtitle(paste0(B,"-order polynomial LOQ model for: ",Targets[i]))
      }
      else {
        Decay.Plot <- Decay.Plot +
          stat_smooth(method = "lm", formula = y ~ x,se=F) +
          ggtitle(paste0("Linear LOQ model for: ",Targets[i]))
      }
    } 
  }
  
  if(is.na(LOQ.list[i+1])==T) {
    Decay.Plot <- Decay.Plot + 
      ggtitle(paste0("LOQ model for: ",Targets[i]," not solvable."))
  }
  
  print(Decay.Plot)
  readline(prompt="Press [Enter] for next plot.")
  print("Calculating... Please wait.")
}
