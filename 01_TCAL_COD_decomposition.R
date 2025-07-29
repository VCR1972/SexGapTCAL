
cO <- c(2450)
names <- c("USA")
names2 <- c("United States")
P=1

options(digits = 4)

library(data.table)
library(reshape2)

Causes<-c("Total", #1,
          "Neoplasms", #2
          "Neoplasms (without gynecologic cancer)", #3
          "Neoplasms (without breast cancer)", #4
          "Cardiovascular disease", #5
          "Respiratory disease", #6
          "External causes") #7


library(RColorBrewer)
mypalette<-brewer.pal(8,"YlGnBu")
mypalette2<-rev(brewer.pal(8,"YlOrRd"))
WildColors<-c(mypalette[5:8],"white","white",mypalette2[c(6,4,2,1)])
WildColors<-rev(c(rev(WildColors[1:4]),"white","white",WildColors[7:10]))
# the different levels of the Z values
levels<- c(-1,-0.1,-0.01,-0.001,-0.0001,0,.0001,.001,.01,.1,1)
options(scipen=10)

customAxis <- function() { 
  n <- length(levels) 
  y <- seq(min(levels), max(levels), length.out=n)
  rect(y[1:(n-1)], 0, y[2:n], 1, col=WildColors)
  axis(1, at=y, labels=levels) 
}

source("Functions/Filled.contour3.R")
source("Functions/Filled.legend_flip.R")



MakeLetter <- function(a, where="topleft", cex=1)
{legend(where, pt.cex=0, bty="n", title=a, cex=cex, legend=NA)}

#### Functions to Calculate TCAL diff #####

CALDecompFunction<-function(Mx1,Mx2,Y1,Y2,Name1,Name2){
  CALlx<-c()
  CALlx1<-c()
  CALlx2<-c()
  PxCh<-c()
  
  YM<-Y2-Y1
  
  for (x in 1:111){
    if (x <(YM+1)){
      px1<-c()
      px2<-c()
      for (z in 1:x){
        px1<-c(px1,Mx1[z,YM-x+z])
        px2<-c(px2,Mx2[z,YM-x+z])
      }
      pxCH<-c(log(px2/px1),rep(0,111-x)) 
      
      lx1<-prod(px1)
      lx2<-prod(px2)
      
      
    }
    if (x >(YM)){
      px1<-c()
      px2<-c()
      for (z in (x-YM+1):x){
        px1<-c(px1,Mx1[z,YM-x+z])
        px2<-c(px2,Mx2[z,YM-x+z])
      }
      
      px1<-c(rep(1,(x-YM)),px1)
      px2<-c(rep(1,(x-YM)),px2)
      pxCH<-c(log(px2/px1),rep(0,111-x))	 	
      
      lx1<-prod(px1)
      lx2<-prod(px2)
      
    }
    CALlx1<-c(CALlx1,lx1)
    CALlx2<-c(CALlx2,lx2)
    
    PxCh<-cbind(PxCh,pxCH)
    
  }
  CALlx<- t(matrix(rep((CALlx1+ CALlx2)/2,111),111))
  
  PxCh[is.na(PxCh)]<-0
  
  ## as Guillot calculates this plus a one for l(0)
  A1<-sum(c(1,CALlx1))+.5
  A2<-sum(c(1,CALlx2))+.5
  A3<-sum(CALlx2)-sum(CALlx1)
  A4<-sum(PxCh* CALlx)
  
  #print(rbind(c(paste("CAL-",Name1),paste("CAL-",Name2),"Diff","est-Diff"),round(c(A1,A2,A3,A4),2)))
  return(PxCh* CALlx)
}

#### Functions to Calculate TCAL diff by Causes #####

CALDecompFunctionCause<-function(Mx1,Mx2,Y1,Y2,CAU1,CAU2,Z){
  
  CALlx<-c()
  CALlx1<-c()
  CALlx2<-c()
  PxCh<-c()
  PxChc<-c()
  
  YM<-Y2-Y1
  
  for (x in 1:111){
    if (x <(YM+1)){
      px1<-c()
      px2<-c()
      
      px1c<-c()
      px2c<-c()
      for (z in 1:x){
        px1<-c(px1,Mx1[z,YM-x+z])
        px2<-c(px2,Mx2[z,YM-x+z])
        
        px1c<-c(px1c,Mx1[z,YM-x+z]^CAU1[z,YM-x+z])
        px2c<-c(px2c,Mx2[z,YM-x+z]^CAU2[z,YM-x+z])
      }
      pxCH<-c(log(px2/px1),rep(0,111-x)) 
      pxCHc<-c(log(px2c/px1c),rep(0,111-x)) 
      
      lx1<-prod(px1)
      lx2<-prod(px2)
      
      
    }
    if (x >(YM)){
      px1<-c()
      px2<-c() 
      
      px1c<-c()
      px2c<-c()
      for (z in (x-YM+1):x){
        px1<-c(px1,Mx1[z,YM-x+z])
        px2<-c(px2,Mx2[z,YM-x+z])
        
        px1c<-c(px1c,Mx1[z,YM-x+z]^CAU1[z,YM-x+z])
        px2c<-c(px2c,Mx2[z,YM-x+z]^CAU2[z,YM-x+z])
        
      }
      
      px1<-c(rep(1,(x-YM)),px1)
      px2<-c(rep(1,(x-YM)),px2)
      
      px1c<-c(rep(1,(x-YM)),px1c)
      px2c<-c(rep(1,(x-YM)),px2c)
      
      pxCH<-c(log(px2/px1),rep(0,111-x))	 	
      pxCHc<-c(log(px2c/px1c),rep(0,111-x)) 
      
      lx1<-prod(px1)
      lx2<-prod(px2)
      
    }
    CALlx1<-c(CALlx1,lx1)
    CALlx2<-c(CALlx2,lx2)
    
    PxCh<-cbind(PxCh,pxCH)
    PxChc<-cbind(PxChc,pxCHc)
    
  }
  CALlx<- t(matrix(rep((CALlx1+ CALlx2)/2,111),111))
  
  PxCh[is.na(PxCh)]<-0
  PxChc[is.na(PxChc)]<-0
  
  ## as Guillot calculates this plus a one for l(0)
  A1<-sum(c(1,CALlx1))+.5
  A2<-sum(c(1,CALlx2))+.5
  A3<-sum(CALlx2)-sum(CALlx1)
  A4<-sum(PxCh* CALlx)
  A5<-sum(PxChc* CALlx)
  
  
  # print(rbind(c(paste("CAL-",Name1),paste("CAL-",Name2),"Diff","est-Diff","Cause i"),round(c(A1,A2,A3,A4,A5),2)))
  
  if(Z==1){return(PxChc* CALlx)}
  if(Z==2){return(round(c(A1,A2,A3,A4,A5),2))}
}

#### Prepare the smoothed data as data.frame ####

read.ICD <- function(cO,names){
  
  B1<-t(read.csv(paste("Causes/","Country",cO,"Cause1","m",".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
  B2<-t(read.csv(paste("Causes/","Country",cO,"Cause1","f",".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
  
  for (j in 2:7){
    B1C<-t(read.csv(paste("Causes/","Country",cO,"Cause",j,"m",".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
    colnames(B1C)[2:112] <- as.character(c(0:110))
    colnames(B1C)[1] <- "Year"
    longB1C <- reshape2::melt(as.data.frame(B1C),id.vars="Year")
    colnames(longB1C)[2] <- "Age"
    longB1C$Cause <- j
    longB1C$Country <- "male"
    B2C<-t(read.csv(paste("Causes/","Country",cO,"Cause",j,"f",".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
    colnames(B2C)[2:112] <- as.character(c(0:110))
    colnames(B2C)[1] <- "Year"
    longB2C <- reshape2::melt(as.data.frame(B2C),id.vars="Year")
    colnames(longB2C)[2] <- "Age"
    longB2C$Cause <- j
    longB2C$Country <- "female"
    if(j==2){doble <- rbind(longB1C,longB2C)}
    if(j>2){doble <- rbind(doble,rbind(longB1C,longB2C))}
  }
  return(doble)
}

#### Functions to simulate counts by causes ####

count2rows2 <- function(cou){
  n <- sum(cou)
  p <- cou/n
  if(n>0){
    newdat <- 
      suppressWarnings(
        cbind(c(2:7), Freq = rmultinom(1, n, p)))
  }
  if(n==0){
    newdat <- cbind(c(2:7), Freq = rep(0, 12))
  }
  return(newdat[,2])
}

#### Functions to simulate the TCAL Cause specific diff ####

CALDecompFunctionCauseboot<-function(P,R){
  
  Names<-names[P]
  Names2<-names2[P]
  CO<-cO[P]
  
  data <- read.ICD(CO, Names)#Reading data
  # First caclulating TCAL by cause
  
  NAME<-paste("external/HMD/",Names,".mltper_1x1",".txt",sep="")
  A1<-read.table(NAME,header=TRUE,fill=TRUE,skip=1)
  
  NAME<-paste("external/HMD/",Names,".fltper_1x1",".txt",sep="")
  A2<-read.table(NAME,header=TRUE,fill=TRUE,skip=1)
  
  B1<-t(read.csv(paste("Causes/","Country",CO,"Cause1","m",".txt",sep=""),
                 header=TRUE,fill=TRUE,skip=0)[,-1])
  B2<-t(read.csv(paste("Causes/","Country",CO,"Cause1","f",".txt",sep=""),
                 header=TRUE,fill=TRUE,skip=0)[,-1])
  
  Y2<-2020
  Y1<-Y2-65
  
  if(Names=="PRT"){Y2<-2019}
  if(Names=="NZL_NP"){Y2<-2018}
  
  data <- data[data$Year %in% c(Y1:Y2),]
  #Ensuring we are working on the common range of years
  
  A2<-A2[(A2$Year>(Y1-1))&(A2$Year<(Y2+1)),]
  A1<-A1[(A1$Year>(Y1-1))&(A1$Year<(Y2+1)),]
  
    
  b1<-B1[(B1[,1]>(Y1-1))&(B1[,1]<(Y2+1)),-1]
  b2<-B2[(B2[,1]>(Y1-1))&(B2[,1]<(Y2+1)),-1]

  qx1<-matrix(1-A1$qx,111)
  qx2<-matrix(1-A2$qx,111)
  
  #   CALDecompFunctionCause(qx1,qx2,Y2,cty1,cty2,bb1c,bb2c)
  Repl <- matrix(NA,nrow=R,ncol=7)
  #Now the bootstrap Conf. Inf.
  
  for(i in 1:R){
    data$Strat <- paste(data$Age,data$Year,data$Country,sep=":")
    a <- tapply(data$value,data$Strat,count2rows2)#Resampling
    b <- as.data.frame(do.call(rbind, a))
    b$id <- row.names(b)
    c <- reshape2::melt(b, id=c("id"))
    c$Cause <- as.numeric(gsub("V","",(c$variable)))+1
    newdata <- data
    newdata <- merge(data,c,by.x=c("Strat","Cause"),by.y=c("id","Cause"))
    newdata <- newdata[order(newdata$Country,newdata$Cause,newdata$Year,newdata$Age),]
    
    ##I have to create from data matrix bb1c and bb2c, which are matrices having (for country 1 and 2) the share of deaths, age-year-country-specific, with respect to total deaths. TOTAL1 and TOTAL2 have the total n. of deaths for each combination of age class (rows) and year(cols)
    
    ## Now the same matrices should be created having only the deaths (for country 1 and 2 ) due to specific cause of death
    
    for(j in range(data$Cause)[1]:range(data$Cause)[2]){
      SPEC1 <- t(matrix(newdata$value.y[newdata$Country=="male"&newdata$Cause==j],111))
      SPEC2 <- t(matrix(newdata$value.y[newdata$Country=="female"&newdata$Cause==j],111))
      
      ## Now CAU1 and CAU2 ive the share of deaths due to specific cause. Note that these are still in age classes, while in the original code data had been split into one year age classes (and qx are in one year age classes)
      CAU1 <- t(SPEC1/b1)
      CAU2 <- t(SPEC2/b2)
      print(c(i,j))
      
      Repl[i,j] <- sum(CALDecompFunctionCause(qx1,qx2,Y1,Y2,CAU1,CAU2,1),na.rm = T)
    }
  }
  return(Repl)
}

#### Plotting functions ####

CALDecompPlot <- function(P,CI.CAL){
  
  Names<-names[P]
  Names2<-names2[P]
  CO<-cO[P]
  
  CONTI <- c()
  

  gap <- read.table(paste("gap_CI.txt",sep=""))
  gap <- gap[gap$Country==Names,]

  panel <- paste(Names2,".pdf",sep = "")
  
  pdf(panel, width = 12, height = 8)
  plot.new()
  
  for (CC in c(1,2,3,4,5,6,7)){
    
    femalenametb <- paste("external/HMD/",Names[1],".fltper_1x1.txt",sep="")
    malenametb <- paste("external/HMD/",Names[1],".mltper_1x1.txt",sep="")
    
    A1<-read.table(malenametb,header=TRUE,fill=TRUE,skip=1)
    A2<-read.table(femalenametb,header=TRUE,fill=TRUE,skip=1)
    
    B1<-t(read.csv(paste("Causes/","Country",CO[1],"Cause1","m",".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
    B1C<-t(read.csv(paste("Causes/","Country",CO[1],"Cause",CC,"m",".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
    
    B2<-t(read.csv(paste("Causes/","Country",CO[1],"Cause1","f",".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
    B2C<-t(read.csv(paste("Causes/","Country",CO[1],"Cause",CC,"f",".txt",sep=""),header=TRUE,fill=TRUE,skip=0)[,-1])
    
    Y2<-2020
    Y1<-Y2-65
    
    Year1<-Y1-1
    Year2<-Y2+1
    
    A2<-A2[(A2$Year>(Y1-1))&(A2$Year<(Y2+1)),]
    A1<-A1[(A1$Year>(Y1-1))&(A1$Year<(Y2+1)),]
    
    b1<-B1[(B1[,1]>(Y1-1))&(B1[,1]<(Y2+1)),-1]
    b2<-B2[(B2[,1]>(Y1-1))&(B2[,1]<(Y2+1)),-1]
      
    b1c<-B1C[(B1C[,1]>(Y1-1))&(B1C[,1]<(Y2+1)),-1]
    b2c<-B2C[(B2C[,1]>(Y1-1))&(B2C[,1]<(Y2+1)),-1]
    
    bb1c<-t(b1c/b1)
    bb2c<-t(b2c/b2)
    
    qx1<-matrix(1-A1$qx,111)
    qx2<-matrix(1-A2$qx,111)
    
    #Total
    CALlxDecomp<-CALDecompFunction(qx1,qx2,Y1,Y2,"Males","Females")
    
    # # Cause-Deleted
    # CALlxDecompcause<-CALDecompFunctionCause(qx1,qx2,Y2,1-bb1c,1-bb2c,1)
    
    CALlxDecompcause<-CALDecompFunctionCause(qx1,qx2,Y1,Y2,bb1c,bb2c,1)
    # Causes-specific 
    
    if(CC!=1){
      CI.lower <- round(quantile(CI.CAL[,CC],(1-0.95)/2),2)
      CI.upper <- round(quantile(CI.CAL[,CC],1-(1-0.95)/2),2)
      CI.mean <- round(mean(CI.CAL[,CC]),2)
    }
    if(CC==1){
      CI.lower <- round(gap$low,2)
      CI.upper <- round(gap$up,2)
      CI.mean <- round(gap$mean,2)
    }
    
    
    CONTI <- rbind(CONTI,c(CALDecompFunctionCause(qx1,qx2,Y1,Y2,bb1c, bb2c,2),
                           CI.mean,CI.lower,CI.upper))
    
    CALlxD<-matrix(0,111,111)
    CALlxDS<-CALlxD
    CALlxDS2<-CALlxD
    CALlxD2<-CALlxD
    
    Age<-c(0:110)
    
    YEARS<-c((Y2-110):Y2)
    
    for (y in 1:111){
      for (x in 1:y){
        CALlxD[x,(111-y+x)]<-CALlxDecomp[x,y]			
        CALlxDS[x,(111-y+x)]<-sum(CALlxDecomp[(1:x),y])
        CALlxD2[x,(111-y+x)]<-CALlxDecompcause[x,y]			
        CALlxDS2[x,(111-y+x)]<-sum(CALlxDecompcause[(1:x),y])}}
    
    if (CC == 1){
      par(new = "TRUE",plt = c(0.05,0.425,0.2,0.9),las = 1,cex.axis = 1)
      filled.contour3(YEARS,Age,t(CALlxDS2),
                      xlim = c(1950,2020),levels=levels, 
                      col=WildColors,key.axes=customAxis(),
                      ylab="Age-contribution",cex.lab=1.2,
                      frame.plot = F,plot.axes = F)
      segments(x0=2020,y0=0,y1=110)
      axis(1)
      axis(2)
      par(xpd = NA)
      MakeLetter(paste(Causes[CC], ", \n gap = ", 
                       round(gap$mean,2),
                       " [",round(gap$low,2),",",round(gap$up,2),"]", 
                       sep = ""))}
    
    if (CC == 2){
      par(new = "TRUE",plt = c(0.525,0.95,0.2,0.9),las = 1,cex.axis = 1)
      filled.contour3(YEARS,Age,t(CALlxDS2),
                      xlim = c(1950,2020),levels=levels,
                      col=WildColors,key.axes=customAxis(),
                      ylab="Age-contribution",cex.lab=1.2,
                      frame.plot = F,plot.axes = F)
      segments(x0=2020,y0=0,y1=110)
      axis(1)
      axis(2)
      MakeLetter(paste(Causes[CC], ", \n gap = ", 
                       CONTI[2,5], 
                       " [",CONTI[2,7],",",CONTI[2,8],"]", 
                       sep = ""))}
    
    par(new = "TRUE",plt = c(0.1,0.9,0.05,0.1),las = 1,cex.axis = 1)

    filled.legend(YEARS,Age,t(CALlxDS2),levels=levels,
                  col=WildColors,key.axes = customAxis(),
                  las=1,cex.lab=1.2)
    
    # mtext("Males",# mtext("Males",# mtext("Males",1,0.5,adj=.9,cex=1.2)
    # mtext("Females",3,0.5,adj=.9,cex=1.2)
    # mtext(names2[P],3,3,adj=0.35,cex=1.8)
  }
  dev.off()
  return(CONTI)
}


CI.CAL <- CALDecompFunctionCauseboot(1,10)
result.CAL <- CALDecompPlot(1,CI.CAL)
 