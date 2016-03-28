## This script performs the analyses reported in the main text and reproduces (in rough form) figures 1 and 3. It begins by reading in local ranges of genera from the raw data file ("Local.Ranges.csv", which can be downloaded either at the GitHub page or at Data Dryad) and processessing these data to extract global first and last appearances of each genus and extimates of the geographic distribution of the genus in each timeslice (lines 25-242). It then uses the distribution of each genus on each paleocontinental block and modeled SST estimates from Herrmann et al. 2004 to estimate the thermal tolerance range of each genus on each paleoncontinental block and then determine whether it would continue to be able to access coastal habitats within its thermal tolerance ranges in a modeled climate state shift (greenhouse to icehouse and, for comparison, icehouse to greenhouse)(lines 243-573).  Paleocontinental blocks are defined based on paleolatitudinal reconstructions from Cocks et al., 2011, as figured in Rasmussen and Harper 2011. Blocks that were joined or so close that the distance would not have presented a substantial dispersal barrier (for example Avalonia and Baltica) were considered as single blocks.  See supplemental materials for further information about blocks used. The data frame output by this script is then used as input for analyses of extinction patterns within each timeslice and comparisons among alternative latest Katian (Katian 5) models (lines 577-906).

## Script takes 10-20 minutes to run depending on processing power

## load all packages to be used
library(gdata)
library(plyr)
library(miscTools)
library(reshape2)
library(ggplot2)
library(fields)
library(PBSmapping)
library(AICcmodavg)
library(dismo)
library(MuMIn)
library(gbm)
library(rpart)
library(rpart.plot)

## modify this line to set working directory
setwd(dir="~/Documents/Ord.Brach.Extinctions.Code")
BrachData <- read.csv("Local.Ranges.csv",header=TRUE,stringsAsFactors=FALSE)

## choose whether to use Genera or Subgenera
BrachData$Genus.Use <- BrachData$Operative_Genus.Subgenus

## remove genera first appearing in Llandovery to avoid edge effects
BrachData <- BrachData[(BrachData$Appearance != "Llandovery"),]

## choose binning factor for lat measures (determines how coarsely range metrics are binned)
bin <- 1

## remove unresolved occurrences, unnamed occurrences, inarticulates, and genera with uncertain higher taxonomy
BrachData <- drop.levels(subset(BrachData,BrachData$interval_resolved == 1 & BrachData$Genus_Modified != 1 & BrachData$Old.Class != "Inarticulate" & BrachData$Order != "Indet."))
BrachData <- BrachData[!is.na(BrachData$Genus.Use), ]

## make function to expand occurrences ranges 
expand <- function(df){
  Interval <- seq(as.numeric(df$FAD_interval),as.numeric(df$LAD_interval),by=1)
  Operative_Genus <- rep(df$Genus.Use,(as.numeric(df$LAD_interval)-as.numeric(df$FAD_interval)+1))
  strat.ranges <- data.frame(Operative_Genus,Interval)
  expanded <- merge(strat.ranges,df)
}
## apply function to each occurrence

Expanded.BrachData <- ddply(BrachData,.(Occurrence_num),expand,.progress="text")

## add column of interval IDs, coarser than original
Interval <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
Interval_Top <- c(456,452,448,447.5,447,446.5,446.2,445.8,445.5,445.3,445,444.7,444.05,443.4, 442.1,440.5,439.0,438.0,437.0,436.0,433.4,430.8,428)
Interval_Bottom <- c(460.9,456,452,448,447.5,447,446.5,446.2,445.8,445.5,445.3,445,444.7,444.05,443.4, 442.1,440.5,439.0,438.0,437.0,436.0,433.4,430.8)
Combined_Intervals <- data.frame(Interval,Interval_Top,Interval_Bottom)
Expanded.BrachData <- merge(Expanded.BrachData,Combined_Intervals,by="Interval")


## determine where apparent section truncations occur (no brachiopod occurrences in an area in which there were brachiopods reported in the previous interval). Suppress to save time unless needed for resampling routines

Truncations <- function(df) {
  packages <- list(2)
  for(i in 1:23) {
    Interval <- i
    present <- ifelse(min(df$FAD_interval) <= i & max(df$LAD_interval) >= i,1,0)
    packages[[i]] <- data.frame(Interval,present)}
  packages <- do.call("rbind",packages)
  return(packages)
}

loc.packages <- ddply(Expanded.BrachData,.(Loc_ID),each(Truncations),.progress="text")
nexts <- loc.packages[2:nrow(loc.packages),]
colnames(nexts) <- c("Loc_ID2","Interval2","present2")
nexts[nrow(nexts)+1,] <- c(NA,NA,NA)
rownames(nexts) <- 1:nrow(nexts)
loc.packages2 <- data.frame(loc.packages,nexts)
loc.packages2$Truncation <- ifelse(loc.packages2$Loc_ID == loc.packages2$Loc_ID2 & loc.packages2$present== 1 & loc.packages2$present2 == 0,1,0)
Packages <- loc.packages2[c("Loc_ID","Interval","Truncation")]

## merge in package truncations to original
Expanded.BrachData <- merge(Expanded.BrachData,Packages,by=c("Loc_ID","Interval"),all.y=FALSE)


## tabulate genus richness in each time-space combination, mean depth of each time-space combination
Local.Genera <- function(df) length(df$Genus.Use)
Loc.Depth <- function(df) mean(df$mean_depth)
Deep.Only <- function(df) {
  df$shallow <- ifelse(df$min_depth < 3,1,0)
  meand <- mean(df$shallow)
  deep.only <- ifelse(meand > 0,0,1)
  return(deep.only)
}

Loc.Gen <- ddply(Expanded.BrachData,.(Interval_Top,Palaeolongitude,Palaeolatitude),each(Local.Genera,Loc.Depth,Deep.Only))
Expanded.BrachData <- merge(Expanded.BrachData,Loc.Gen,by=c("Interval_Top","Palaeolongitude","Palaeolatitude"))
Expanded.BrachData$North.Hemisphere <- ifelse(Expanded.BrachData$Palaeolatitude > 0,1,0)

#round lats
Expanded.BrachData$Palaeolatitude <- round(Expanded.BrachData$Palaeolatitude,0)
Expanded.BrachData$AbsPalaeolatitude <- abs(Expanded.BrachData$Palaeolatitude)
Expanded.BrachData$Palaeolongitude <- round(Expanded.BrachData$Palaeolongitude,0)

##TODO: check whether still needed
write.table(Expanded.BrachData,"Brachiopod_occurrences_expanded.csv",sep=",")

## create global lat-long grid
LatSeq <- seq(-90,90,by=10)
LongSeq <- seq(-180,180,by=30)

polys<- makeGrid (x=round(LongSeq,0), y=round(LatSeq,0),projection="UTM")

## range.buffer: set to determine how many degree latitude should be added to the range apparent range.  Assumes that observed paleolatitudinal ranges are always underestimates of true paleolatitudinal ranges and adds a set buffer to both the north and south range boundaries

range.buffer <- 2

## make function to extract range and richness attributes
ranges <- function(df){
  #df <- Expanded.BrachData[23,]
  if(length(na.omit(df$Palaeolatitude))>1){
    events <- data.frame(EID=1:length(df$Palaeolatitude), X = as.numeric(as.vector(df$Palaeolongitude)),Y=as.numeric(as.vector(df$Palaeolatitude)))
    events <- as.EventData(na.omit(events), projection='UTM')
    fc <- findCells(events, polys)
    richness <- length(unique(df$Genus.Species))
    eac <-length(unique(subset(fc,select=c(PID,SID)))[,1])
    maxdist <- max(rdist.earth(events[,2:3]))
    gcd <-ifelse(maxdist > 10, maxdist, 10)
    MaxLat <- max(df$Palaeolatitude) + range.buffer
    MinLat <- min(df$Palaeolatitude) - range.buffer
    MeanAbsLat <- round(mean(abs(df$Palaeolatitude))/bin,0)*bin
    LatRange <- MaxLat - MinLat
    Single_occurrence <- 0
    MaxAbsLat <- max(abs(df$Palaeolatitude)) + range.buffer
    MinAbsLat <- min(abs(df$Palaeolatitude)) - range.buffer
    MaxAbsLat <- ifelse(MaxAbsLat > 90,90,MaxAbsLat)
    MinAbsLat <- ifelse(MinAbsLat < 0,0,MinAbsLat)
    AbsLatRange <- MaxAbsLat - MinAbsLat
    AbsLatRange <- ifelse(AbsLatRange > 90,90,AbsLatRange)
    AbsLatRange <- ifelse(AbsLatRange == 0,1,AbsLatRange)
    NumContinents <- length(unique(df$Merged.Terrains))
    Endemic <- ifelse(NumContinents==1,1,0)
    NumRegions <- length(unique(df$Present.region))
    MaxDepth <- max(df$max_depth)
    MinDepth <- min(df$min_depth)
    MeanDepth <- mean(df$mean_depth)
    Mean.Loc.Depth <- mean(df$Loc.Depth)
    PropTrunc <- mean(df$Truncation)
    Pedicle.state <- max(df$Pedicle_state)
    Median.Comm.Rich <- median(df$Local.Genera)
    Prop.Cratonic <- mean(df$Cratonic)
    Prop.North <- mean(df$North.Hemisphere)
  }     
  else{
    eac=1
    gcd=10
    richness <- length(unique(df$Genus.Species))
    MaxLat <- max(df$Palaeolatitude) + range.buffer
    MinLat <- min(df$Palaeolatitude) - range.buffer
    MeanAbsLat <- round(mean(abs(df$Palaeolatitude))/bin,0)*bin
    LatRange <- 0
    MaxAbsLat <- max(abs(df$Palaeolatitude)) + range.buffer
    MinAbsLat <- min(abs(df$Palaeolatitude)) - range.buffer
    MaxAbsLat <- ifelse(MaxAbsLat > 90,90,MaxAbsLat)
    MinAbsLat <- ifelse(MinAbsLat < 0,0,MinAbsLat)
    AbsLatRange <- MaxAbsLat - MinAbsLat
    AbsLatRange <- ifelse(AbsLatRange > 90,90,AbsLatRange)
    AbsLatRange <- ifelse(AbsLatRange == 0,1,AbsLatRange)
    NumContinents <- length(unique(df$Merged.Terrains))
    Endemic <- ifelse(NumContinents==1,1,0)
    Single_occurrence <- 1
    NumRegions <- length(unique(df$Present.region))
    MaxDepth <- max(df$max_depth)
    MinDepth <- min(df$min_depth)
    MeanDepth <- mean(df$mean_depth)
    Mean.Loc.Depth <- mean(df$Loc.Depth)
    PropTrunc <- mean(df$Truncation)
    Pedicle.state <- max(df$Pedicle_state)
    Median.Comm.Rich <- median(df$Local.Genera)
    Prop.Cratonic <- mean(df$Cratonic)
    Prop.North <- mean(df$North.Hemisphere)
    
  }
  return(data.frame(eac,gcd,richness,MaxLat,MinLat,LatRange,MaxAbsLat,MinAbsLat,MeanAbsLat,AbsLatRange,NumContinents,Endemic,NumRegions,MaxDepth,MinDepth,MeanDepth,Single_occurrence,Pedicle.state,Median.Comm.Rich,Mean.Loc.Depth,Prop.Cratonic,Prop.North,PropTrunc))
}

byGenus <- ddply(Expanded.BrachData,.(Interval,Interval_Top,Interval_Bottom,Order,Family,Genus.Use),ranges,.progress="text")

## Extract first and last occurrences of genera from the processed dataset
FAD <- function(df) max(df$Interval_Top)
LAD <- function(df) min(df$Interval_Top)

byGenus.2 <- ddply(byGenus,.(Genus.Use),each(FAD,LAD),.progress="text")

byGenus.3 <- merge(byGenus,byGenus.2,by = "Genus.Use")

## indicate whether genus has first or last appearance in a given interval

Ex <- ifelse(byGenus.3$Interval_Top == byGenus.3$LAD,1,0)
Or <- ifelse(byGenus.3$Interval_Top == byGenus.3$FAD,1,0)
Age <- byGenus.3$FAD - byGenus.3$Interval_Top 
ExClass <- ifelse(Ex==1,"extinct","survive")
OrClass <- ifelse(Or==1,"origination","holdover")
byGenus.4 <- data.frame(byGenus.3,Ex,Or,Age,ExClass,OrClass)


## rename range parameters etc.
gcd <- byGenus.4$gcd
richness <- byGenus.4$richness
occupancy <- byGenus.4$eac
continents <- byGenus.4$NumContinents
endemic <- byGenus.4$Endemic
MaxLat <- byGenus.4$MaxAbsLat
MinLat <- byGenus.4$MinAbsLat
MinLatRel <- byGenus.4$MinLat
MaxLatRel <- byGenus.4$MaxLat
LatRangeRel <- MaxLatRel - MinLatRel
LatRange <- byGenus.4$LatRange
AbsLatRange <- byGenus.4$AbsLatRange
MaxLatRange.ever <- MaxLat-MinLat 
MeanAbsLat <-  byGenus.4$MeanAbsLat
MeanDepth <- byGenus.4$MeanDepth
Single_occurrence <- byGenus.4$Single_occurrence
Age <- byGenus.4$Age
Interval <- byGenus.4$Interval
Interval_Top <- byGenus.4$Interval_Top
Interval_Bottom <- byGenus.4$Interval_Bottom
MaxDepth <- byGenus.4$MaxDepth
MinDepth <- byGenus.4$MinDepth
DepthRange <- MaxDepth - MinDepth
Order <- as.factor(byGenus.4$Order)
Family <- as.factor(byGenus.4$Family)
Genus <- as.factor(byGenus.4$Genus.Use)
PropTrunc <- byGenus.4$PropTrunc
Pedicle.state <- byGenus.4$Pedicle.state
Median.Comm.Rich <- byGenus.4$Median.Comm.Rich
Mean.Loc.Depth <- byGenus.4$Mean.Loc.Depth
Prop.Cratonic <- byGenus.4$Prop.Cratonic
Ex <- byGenus.4$Ex
ExClass <- byGenus.4$ExClass
Or <- byGenus.4$Or
OrClass <- byGenus.4$OrClass
North.Hemisphere <- ifelse(byGenus.4$MaxLat > 0, 1,0)
Prop.North <- byGenus.4$Prop.North


final_data <- data.frame(Interval,Interval_Top,Interval_Bottom,Ex,Or,ExClass,OrClass,richness,occupancy,Single_occurrence,continents,endemic,MaxLat,MinLat,MinLatRel,MaxLatRel,LatRange,LatRangeRel,AbsLatRange,MeanAbsLat,MeanDepth,MaxDepth,MinDepth,DepthRange,Age,gcd,Order,Family,Genus,Pedicle.state,Median.Comm.Rich,Mean.Loc.Depth,MaxLatRange.ever,Prop.Cratonic,North.Hemisphere,Prop.North,PropTrunc)


## now determine possible lat ranges in Icehouse/Greenhouse world based on Katian 5 distributions

## bring in seasonality gradient based on modern world

lats <- c(-90,-80,-70,-60,-50,-40,-30,-20,-15,-10,0,10,15,20,30,40,50,60,70,80,90)
seasonality <-c(7,5.5,6,9.5,15,17,9.5,3,0,2,5,2,0,3,9.5,17,15,9.5,6,5.5,7)
seasonality <- seasonality +1
season <- approx(lats,seasonality, xout = seq(-90,90,by=1), method="linear", n=90)
season <- data.frame(season$x,season$y)
colnames(season) <- c("latitude","seasonality")

## multiply gradient by an arbitary factor tuned to match modern, use this to estimate summer max and winter min temperatures across the latitudinal gradient

seas.fac <- .25

## create mean annual temperature gradients for greenhouse and icehouse states based on Herrmann et al., 2004 models.  Greenhouse/Warm: 15X PAL CO2, high sea level, Ashgill paleogeography, Icehouse/Cool: 8X PAL CO2, low sea level, Ashgill paleogeography

latitude <- c(-90, -89, -88, -87, -86, -85, -84, -83, -82, -81, -80, -79, -78, -77, -76, -75, -74, -73, -72, -71, -70, -69, -68, -67, -66, -65, -64, -63, -62, -61, -60, -59, -58, -57, -56, -55, -54, -53, -52, -51, -50, -49, -48, -47, -46, -45, -44, -43, -42, -41, -40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90)
              
Ashgill.15X.HSL <- c(1, 1.1, 1.3, 1.4, 1.5, 1.7, 1.8, 1.9, 2.1, 2.2, 2.3, 2.5, 2.6, 2.7, 2.9, 3, 3.3, 3.5, 3.8, 4.1, 4.3, 4.6, 4.9, 5.1, 5.2, 5.3, 5.5, 5.6, 5.7, 5.9, 6, 6.3, 6.5, 6.8, 7.1, 7.3, 7.6, 7.9, 8.2, 8.6, 9, 9.4, 9.8, 10.2, 10.6, 11, 11.3, 11.5, 11.8, 12.1, 12.3, 12.6, 12.9, 13.7, 15, 16.3, 17.7, 19, 20.3, 21.7, 23, 23.5, 24.1, 24.6, 25.1, 25.7, 26.2, 26.7, 27.1, 27.4, 27.7, 27.9, 28.2, 28.5, 28.7, 29, 29.3, 29.5, 29.8, 30.1, 30.3, 30.6, 30.9, 30.7, 30.2, 29.7, 29.1, 28.6, 28.1, 27.5, 27, 27.5, 28.1, 28.6, 29.1, 29.7, 30.2, 30.7, 30.9, 30.6, 30.3, 30.1, 29.8, 29.5, 29.3, 29, 28.7, 28.5, 28.2, 27.9, 27.7, 27.4, 27.1, 26.7, 26.2, 25.7, 25.1, 24.6, 24.1, 23.5, 23, 21.7, 20.3, 19, 17.7, 16.3, 15, 13.7, 12.9, 12.6, 12.3, 12.1, 11.8, 11.5, 11.3, 11, 10.6, 10.2, 9.8, 9.4, 9, 8.6, 8.2, 7.9, 7.6, 7.3, 7.1, 6.8, 6.5, 6.3, 6, 5.9, 5.7, 5.6, 5.5, 5.3, 5.2, 5.1, 4.9, 4.6, 4.3, 4.1, 3.8, 3.5, 3.3, 3, 2.9, 2.7, 2.6, 2.5, 2.3, 2.2, 2.1, 1.9, 1.8, 1.7, 1.5, 1.4, 1.3, 1.1, 1)

Ashgill.8X.LSL <- c(-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1.5, -0.6, 0.3, 1.3, 2.2, 3.1, 4.1, 5, 5.3, 5.5, 5.8, 6.1, 6.3, 6.6, 6.9, 7.8, 9.4, 11, 12.6, 14.2, 15.8, 17.4, 19, 19.7, 20.3, 21, 21.7, 22.3, 23, 23.7, 24.2, 24.6, 25, 25.4, 25.8, 26.2, 26.6, 27, 26.9, 26.7, 26.6, 26.5, 26.3, 26.2, 26.1, 25.9, 25.8, 25.7, 25.5, 25.4, 25.3, 25.1, 25, 25.1, 25.3, 25.4, 25.5, 25.7, 25.8, 25.9, 26.1, 26.2, 26.3, 26.5, 26.6, 26.7, 26.9, 27, 26.6, 26.2, 25.8, 25.4, 25, 24.6, 24.2, 23.7, 23, 22.3, 21.7, 21, 20.3, 19.7, 19, 17.4, 15.8, 14.2, 12.6, 11, 9.4, 7.8, 6.9, 6.6, 6.3, 6.1, 5.8, 5.5, 5.3, 5, 4.1, 3.1, 2.2, 1.3, 0.3, -0.6, -1.5, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2)

temps <- data.frame(latitude,Ashgill.15X.HSL,Ashgill.8X.LSL)

temps <- merge(temps,season)

temps$Warm.SST <- temps$Ashgill.15X.HSL
temps$Cool.SST <- temps$Ashgill.8X.LSL
temps$Warm.SST.summer <- temps$Warm.SST + (temps$seasonality*seas.fac)
temps$Warm.SST.winter <- temps$Warm.SST - (temps$seasonality*seas.fac)
temps$Cool.SST.summer <- temps$Cool.SST + (temps$seasonality*seas.fac)
temps$Cool.SST.winter <- temps$Cool.SST - (temps$seasonality*seas.fac)
temps$Warm.SST.summer <- ifelse(temps$Warm.SST.summer > -2,temps$Warm.SST.summer ,-2)
temps$Warm.SST.winter <- ifelse(temps$Warm.SST.winter > -2,temps$Warm.SST.winter ,-2)
temps$Cool.SST.summer  <- ifelse(temps$Cool.SST.summer  > -2,temps$Cool.SST.summer ,-2)
temps$Cool.SST.winter  <- ifelse(temps$Cool.SST.winter  > -2,temps$Cool.SST.winter  ,-2)
temps$Cool.range <- temps$Cool.SST.summer - temps$Cool.SST.winter
temps$Warm.range <- temps$Warm.SST.summer - temps$Warm.SST.winter

diff <- temps$Warm.SST.summer - temps$Warm.SST.winter

temps2 <- temps[c("latitude","Warm.SST.summer","Warm.SST.winter","Cool.SST.summer","Cool.SST.winter")]
temps2 <- melt(temps2,id="latitude")
colnames(temps2) <- c("latitude","model","temp.")
temps2 <- drop.levels(temps2[(temps2$latitude >= 0),])
model <- levels(temps2$model)
states <- c("Icehouse summer SST","Icehouse winter SST","Greenhouse summer SST","Greenhouse winter SST")
states <- cbind(model,states)
temps2 <- merge(temps2,states)

## Make a temperature profile plot to compare icehouse and greenhouse models.

pdf("Temp.Profiles.pdf", width = 8, height = 4)
print(p <- ggplot(temps2,aes(latitude,temp.,colour=states,linetype=states)) + geom_line() + scale_colour_manual(values =c("red","red","blue","blue"),name=" ") + scale_linetype_manual(values =c(1,2,1,2),name=" ") + theme_bw() + ylab("mean sea surface temperature") + xlab("degrees absolute latitude") + scale_x_continuous(breaks=c(0,30,60,90)))
dev.off()

## get ranges of genera on each continent

K5.occs <- Expanded.BrachData[(Expanded.BrachData$Interval_Top == 444.7),]

Merged.Terrains<- function(df){
  # df <- final_data[23,]
  continents <- data.frame(unique(df$Merged.Terrains))
  return(continents)}
continent.ranges <- ddply(K5.occs,.(Genus.Use,Interval),each(Merged.Terrains),.progress="text")

Expanded.BrachData.A <- K5.occs
Expanded.BrachData.A$Palaeolatitude <- Expanded.BrachData.A$Palaeolatitude + range.buffer
Expanded.BrachData.B <- K5.occs
Expanded.BrachData.B$Palaeolatitude <- Expanded.BrachData.B$Palaeolatitude
Expanded.BrachData.C <- K5.occs
Expanded.BrachData.C$Palaeolatitude <- Expanded.BrachData.C$Palaeolatitude - range.buffer

Expanded.BrachData2 <- rbind(Expanded.BrachData.A,Expanded.BrachData.B,Expanded.BrachData.C)

Expanded.BrachData2 <- merge(Expanded.BrachData2,temps,by.x="Palaeolatitude",by.y="latitude")

## extract maximum and minimum temperatures on each continent/terrain. Do not allow temperature below -2 degrees C.

Expanded.BrachData2$Warm.SST.winter <- ifelse(Expanded.BrachData2$Warm.SST.winter < -2,-2,Expanded.BrachData2$Warm.SST.winter)
Expanded.BrachData2$Warm.SST.summer <- ifelse(Expanded.BrachData2$Warm.SST.summer < -2,-2,Expanded.BrachData2$Warm.SST.summer)
Expanded.BrachData2$Cool.SST.winter <- ifelse(Expanded.BrachData2$Cool.SST.winter < -2,-2,Expanded.BrachData2$Cool.SST.winter)
Expanded.BrachData2$Cool.SST.winter <- ifelse(Expanded.BrachData2$Cool.SST.winter < -2,-2,Expanded.BrachData2$Cool.SST.winter)

gen.ranges <- function(df) {
  #df <- Expanded.BrachData2[3,]
  GH.winter <- df$Warm.SST - df$Warm.SST*df$seasonality
  GH.gen.win.min <- min(df$Warm.SST - df$seasonality*seas.fac)
  GH.gen.win.max <- max(df$Warm.SST - df$seasonality*seas.fac)
  GH.gen.sum.min <- min(df$Warm.SST + df$seasonality*seas.fac)
  GH.gen.sum.max <- max(df$Warm.SST + df$seasonality*seas.fac)
  IH.gen.win.min <- min(df$Cool.SST - df$seasonality*seas.fac)
  IH.gen.win.max <- max(df$Cool.SST - df$seasonality*seas.fac)
  IH.gen.sum.min <- min(df$Cool.SST + df$seasonality*seas.fac)
  IH.gen.sum.max <- max(df$Cool.SST + df$seasonality*seas.fac)
  Max.Lat <- max(df$Palaeolatitude)
  Min.Lat <- min(df$Palaeolatitude)
  outs <- data.frame(GH.gen.win.min,GH.gen.win.max,GH.gen.sum.min,GH.gen.sum.max, IH.gen.win.min,IH.gen.win.max,IH.gen.sum.min,IH.gen.sum.max,Max.Lat,Min.Lat)
  return(outs)}
gen.data <- ddply(Expanded.BrachData2,.(Genus.Species,Genus.Use,Interval,Merged.Terrains),each(gen.ranges),.progress="text")

max.cont.lat <- function(df) max(df$Max.Palaeolatitude.Merged.450)
min.cont.lat <- function(df) min(df$Min.Palaeolatitude.Merged.450)
min.loc.lat <- function(df) {
  Katian5 <- df[(df$Interval == 12),]
  min(Katian5$Palaeolatitude)}
max.loc.lat <- function(df) {
  Katian5 <- df[(df$Interval == 12),]
  max(Katian5$Palaeolatitude)}

continent.lats <- ddply(Expanded.BrachData2,.(Merged.Terrains), each(min.cont.lat,max.cont.lat,min.loc.lat,max.loc.lat))
continent.lats <-(unique(do.call(data.frame,lapply(continent.lats, function(x) replace(x, is.infinite(x),0)))))

## choose widest interval to avoid blanks
continent.lats$min.cont.lat <- apply(continent.lats[,2:5],1,min)
continent.lats$max.cont.lat <- apply(continent.lats[,2:5],1,max)


gen.data2 <- merge(gen.data,continent.lats,by="Merged.Terrains",all.y=FALSE)

## set amount to add and subtract from thermal tolerance to test sensitivity ("extratemp").  Defaults to zero unless changed

extratemp <- 0
gen.data2$GH.winter.min <- ifelse(gen.data2$GH.gen.win.min-extratemp > -2,gen.data2$GH.gen.win.min-extratemp,-2)
gen.data2$GH.summer.max <- ifelse(gen.data2$GH.gen.sum.max+extratemp > -2,gen.data2$GH.gen.sum.max+extratemp,-2)
gen.data2$IH.winter.min <- ifelse(gen.data2$IH.gen.win.min-extratemp > -2,gen.data2$IH.gen.win.min-extratemp,-2)
gen.data2$IH.summer.max <- ifelse(gen.data2$IH.gen.sum.max+extratemp > -2,gen.data2$IH.gen.sum.max+extratemp,-2)

## now select the paleolat bins that fit within each genus's thermal tolerance range on each continent it occupies

range.change <- function(df) {
  
  temps$observed.range <-  ifelse(temps$latitude >= df$Min.Lat & temps$latitude <= df$Max.Lat,1,0) 
  temps$warming.sites.end <-  ifelse(temps$Warm.SST.summer  <= df$IH.summer.max & temps$Warm.SST.winter >= df$IH.winter.min,1,0)
  temps$warming.sites.start <- ifelse(temps$Cool.SST.summer  <= df$IH.summer.max & temps$Cool.SST.winter >= df$IH.winter.min,1,0)
  temps$cooling.sites.end <- ifelse(temps$Cool.SST.summer  <= df$GH.summer.max & temps$Cool.SST.winter >= df$GH.winter.min,1,0)
  temps$cooling.sites.start <- ifelse(temps$Warm.SST.summer  <= df$GH.summer.max & temps$Warm.SST.winter >= df$GH.winter.min,1,0)
  
  ## Make a linear vector for calculating paleolatitudinal dist.
  temps$linear.lat <- seq(1:181)
  
  ## Designate paleolatitudinal blocks for each variable
  out <- list()
  for (i in 1:182){
    breaks <- ifelse(temps$warming.sites.start[i] != temps$warming.sites.start[i-1],1,0)
    out[[i]] <- breaks
  }
  out <- do.call("rbind",out)
  outs2 <- cumsum(out)
  temps$block.warming.sites.start <- outs2
  temps$block.warming.sites.start <- c(0,0,temps$block.warming.sites.start[2:180])
  
  out <- list()
  for (i in 1:182){
    breaks <- ifelse(temps$warming.sites.end[i] != temps$warming.sites.end[i-1],1,0)
    out[[i]] <- breaks
  }
  out <- do.call("rbind",out)
  outs2 <- cumsum(out)
  temps$block.warming.sites.end <- outs2
  temps$block.warming.sites.end <- c(0,0,temps$block.warming.sites.end[2:180])
  
  out <- list()
  for (i in 1:182){
    breaks <- ifelse(temps$cooling.sites.start[i] != temps$cooling.sites.start[i-1],1,0)
    out[[i]] <- breaks
  }
  out <- do.call("rbind",out)
  outs2 <- cumsum(out)
  temps$block.cooling.sites.start <- outs2
  temps$block.cooling.sites.start <- c(0,0,temps$block.cooling.sites.start[2:180])
  
  out <- list()
  for (i in 1:182){
    breaks <- ifelse(temps$cooling.sites.end[i] != temps$cooling.sites.end[i-1],1,0)
    out[[i]] <- breaks
  }
  out <- do.call("rbind",out)
  outs2 <- cumsum(out)
  temps$block.cooling.sites.end <- outs2
  temps$block.cooling.sites.end <- c(0,0,temps$block.cooling.sites.end[2:180])
  
  
  ## Remove blocks that have zero latitdunal range (that is, there is no accessible habitat within the estimated thermal tolerance range of the population on the paleocontinent block)
  
  cooling.sites.end <- temps[(temps$cooling.sites.end == 1),]
  cooling.sites.start <- temps[(temps$cooling.sites.start == 1),]
  warming.sites.end <- temps[(temps$warming.sites.end == 1),]
  warming.sites.start <- temps[(temps$warming.sites.start == 1),]
  orig.sites <- temps[(temps$observed.range == 1),]
  
  ## Get mean paleolatitude of original range
  actual.range <-  mean(orig.sites$linear.lat)
  
  ## Make function to select the block closest to the starting lat, so that populations cannot disperse across areas that are ouside their themal range (for example, cannot shift from -20 degrees to 20 degrees if that involves crossing equatorial temperatures outside of the range they occupy)
  
  proximity <- function(df){
    meanlat <- mean(df$linear.lat)
    outs <- data.frame(meanlat)
    return(outs)}
  
  cooling.sites.end.dist <- ddply(cooling.sites.end,.(block.cooling.sites.end),each(proximity))
  cooling.sites.end.dist$diff <- abs(cooling.sites.end.dist$proximity - actual.range)
  mindiff <- min(cooling.sites.end.dist$diff)
  cooling.sites.end.dist <- cooling.sites.end.dist[(cooling.sites.end.dist$diff == mindiff),]
  block.num <- cooling.sites.end.dist$block.cooling.sites.end 
  cooling.sites.end <- cooling.sites.end[(cooling.sites.end$block.cooling.sites.end ==  block.num),]
  
  cooling.sites.start.dist <- ddply(cooling.sites.start,.(block.cooling.sites.start),each(proximity))
  cooling.sites.start.dist$diff <- abs(cooling.sites.start.dist$proximity - actual.range)
  mindiff <- min(cooling.sites.start.dist$diff)
  cooling.sites.start.dist <- cooling.sites.start.dist[(cooling.sites.start.dist$diff == mindiff),]
  block.num <- cooling.sites.start.dist$block.cooling.sites.start 
  cooling.sites.start <- cooling.sites.start[(cooling.sites.start$block.cooling.sites.start ==  block.num),]
  
  warming.sites.end.dist <- ddply(warming.sites.end,.(block.warming.sites.end),each(proximity))
  warming.sites.end.dist$diff <- abs(warming.sites.end.dist$proximity - actual.range)
  mindiff <- min(warming.sites.end.dist$diff)
  warming.sites.end.dist <- warming.sites.end.dist[(warming.sites.end.dist$diff == mindiff),]
  block.num <- warming.sites.end.dist$block.warming.sites.end 
  warming.sites.end <- warming.sites.end[(warming.sites.end$block.warming.sites.end ==  block.num),]
  
  warming.sites.start.dist <- ddply(warming.sites.start,.(block.warming.sites.start),each(proximity))
  warming.sites.start.dist$diff <- abs(warming.sites.start.dist$proximity - actual.range)
  mindiff <- min(warming.sites.start.dist$diff)
  warming.sites.start.dist <- warming.sites.start.dist[(warming.sites.start.dist$diff == mindiff),]
  block.num <- warming.sites.start.dist$block.warming.sites.start 
  warming.sites.start <- warming.sites.start[(warming.sites.start$block.warming.sites.start ==  block.num),]
  
  warming.nsites <- nrow(warming.sites.end)
  cooling.nsites <- nrow(cooling.sites.end)
  cooling.start.sites <- nrow(cooling.sites.start)
  warming.start.sites <- nrow(warming.sites.start)
  icehouse.maxlat <- ifelse(cooling.nsites > 0,max(abs(cooling.sites.end$latitude)),NA)
  icehouse.minlat <- ifelse(cooling.nsites > 0,min(abs(cooling.sites.end$latitude)),NA)
  icehouse.start.maxlat <- ifelse(cooling.start.sites > 0,max(abs(cooling.sites.start$latitude)),NA)
  icehouse.start.minlat <- ifelse(cooling.start.sites > 0,min(abs(cooling.sites.start$latitude)),NA)
  
  greenhouse.maxlat <- ifelse(warming.nsites > 0,max(abs(warming.sites.end$latitude)),NA)
  greenhouse.minlat <- ifelse(warming.nsites > 0,min(abs(warming.sites.end$latitude)),NA)
  greenhouse.start.maxlat <- ifelse(warming.start.sites > 0,max(abs(warming.sites.start$latitude)),NA)
  greenhouse.start.minlat <- ifelse(warming.start.sites > 0,min(abs(warming.sites.start$latitude)),NA)
  
  cool.inhabitable <- ifelse(cooling.nsites==0,0,1)
  warm.inhabitable <- ifelse(warming.nsites==0,0,1)
  projection.lat.diff <- warming.nsites - cooling.nsites
  
  IH.MaxTemp <- df$IH.summer.max
  IH.MinTemp <- df$IH.winter.min
  IH.TempRange <- IH.MaxTemp - IH.MinTemp 
  GH.MaxTemp <- df$GH.summer.max
  GH.MinTemp <- df$GH.winter.min
  GH.TempRange <- GH.MaxTemp - GH.MinTemp
  
  
  outs <- data.frame(warming.nsites,cooling.nsites,warming.start.sites,cooling.start.sites,cool.inhabitable,warm.inhabitable,icehouse.maxlat,icehouse.minlat,icehouse.start.maxlat,icehouse.start.minlat,greenhouse.maxlat,greenhouse.start.maxlat,greenhouse.minlat,greenhouse.start.minlat,projection.lat.diff,icehouse.maxlat,icehouse.start.maxlat,IH.MaxTemp,IH.MinTemp,IH.TempRange,GH.MaxTemp,GH.MinTemp,GH.TempRange)
  return(outs)
}

gen.data3 <- ddply(gen.data2,.(Genus.Species,Genus.Use,Interval,Merged.Terrains),each(range.change),.progress="text")

IH.Spec.Rich <- function(df){
  df$IH.Spec.LatRange <- df$icehouse.maxlat- df$icehouse.minlat
  df$IH.Spec.LatRange[is.na(df$IH.Spec.LatRange)] <- 0
  spec.ranges <- df[c("Genus.Species","IH.Spec.LatRange")]
  spec.ranges <- spec.ranges[(spec.ranges$IH.Spec.LatRange > 0),]
  IH.spec.rich <- length(unique(spec.ranges$Genus.Species))
  return(IH.spec.rich)
}
GH.Spec.Rich <- function(df){
  df$GH.Spec.LatRange <- df$greenhouse.maxlat- df$greenhouse.minlat
  df$GH.Spec.LatRange[is.na(df$GH.Spec.LatRange)] <- 0
  spec.ranges <- df[c("Genus.Species","GH.Spec.LatRange")]
  spec.ranges <- spec.ranges[(spec.ranges$GH.Spec.LatRange > 0),]
  GH.spec.rich <- length(unique(spec.ranges$Genus.Species))
  return(GH.spec.rich)
}

IH.Num.Inhabitable <- function(df) {
  df$IH.Spec.LatRange <- df$icehouse.maxlat- df$icehouse.minlat
  df$IH.Spec.LatRange[is.na(df$IH.Spec.LatRange)] <- 0
  IH.spec.conts <- df[c("Genus.Use","Merged.Terrains","IH.Spec.LatRange")]
  IH.spec.conts <- IH.spec.conts[(IH.spec.conts$IH.Spec.LatRange != 0),]
  IH.num.inhabitable <-  length(unique(IH.spec.conts$Merged.Terrains))
  return(IH.num.inhabitable)
}

GH.Num.Inhabitable <- function(df) {
  df$GH.Spec.LatRange <- df$greenhouse.maxlat- df$greenhouse.minlat
  df$GH.Spec.LatRange[is.na(df$GH.Spec.LatRange)] <- 0
  GH.spec.conts <- df[c("Genus.Use","Merged.Terrains","GH.Spec.LatRange")]
  GH.spec.conts <- GH.spec.conts[(GH.spec.conts$GH.Spec.LatRange != 0),]
  GH.num.inhabitable <-  length(unique(GH.spec.conts$Merged.Terrains))
  return(GH.num.inhabitable)
}

cooling.latrange <-function(df) {
  nsites <- length(na.omit(df$icehouse.maxlat)) 
  ifelse(nsites == 0,0, (max(na.omit(df$icehouse.maxlat)) - min(na.omit(df$icehouse.minlat))))}

cooling.start.latrange <-function(df) {
  nsites <- length(na.omit(df$icehouse.start.maxlat)) 
  ifelse(nsites == 0,0, (max(na.omit(df$icehouse.start.maxlat)) - min(na.omit(df$icehouse.start.minlat))))}

warming.latrange <-function(df) {
  nsites <- length(na.omit(df$greenhouse.maxlat)) 
  ifelse(nsites == 0,0, (max(na.omit(df$greenhouse.maxlat)) - min(na.omit(df$greenhouse.minlat))))}

IH.MaxTemp <- function(df) max(df$IH.MaxTemp)
IH.MinTemp <- function(df) min(df$IH.MaxTemp)
IH.TempRange <- function(df) max(df$IH.TempRange)
GH.MaxTemp <- function(df) max(df$GH.MaxTemp)
GH.MinTemp <- function(df) min(df$GH.MaxTemp)
GH.TempRange <- function(df) max(df$GH.TempRange)

warming.start.latrange <-function(df) {
  nsites <- length(na.omit(df$greenhouse.start.maxlat)) 
  ifelse(nsites == 0,0, (max(na.omit(df$greenhouse.start.maxlat)) - min(na.omit(df$greenhouse.start.minlat))))}

range.changes <- ddply(gen.data3,.(Interval,Genus.Use),each(IH.Spec.Rich,GH.Spec.Rich,IH.Num.Inhabitable,GH.Num.Inhabitable,cooling.latrange,cooling.start.latrange,warming.latrange,warming.start.latrange,IH.MaxTemp,IH.MinTemp,IH.TempRange,GH.MaxTemp,GH.MinTemp,GH.TempRange))

colnames(range.changes) <- c("Interval","Genus","IH Spec. Richness","GH Spec. Richness", "IH.Num.Inhabitable","GH.Num.Inhabitable","cooling.latrange","cooling.start.latrange","warming.latrange","warming.start.latrange","IH.MaxTemp","IH.MinTemp","IH.TempRange","GH.MaxTemp","GH.MinTemp","GH.TempRange")

final_data2 <- merge(final_data,range.changes,by=c("Genus","Interval"),all.x=TRUE)

## Write the processed dataframe as a csv file so that it can be opened in excel, if desired

write.csv(final_data2 ,"OS.brachs.processed.csv",row.names=FALSE)

## Rename the file for subsequent analyses

data <- final_data2 

data$stage_top <- data$Interval_Top

### make figure 2: proportional extinction of genera within depth:palaeolatitude bins.  Assumes genera are present throughout their depth and palaeolatitude ranges

lats <- seq(-90,90,by=5)
depths <- seq(1,6,by=1)
lat.depths <- merge(lats,depths,all=TRUE)
stages <- sort(rep(unique(data$stage_top),length(lat.depths[,1])))
lat.depths <- data.frame(lat.depths,stages)
colnames(lat.depths) <- c("Latitude","Depth","Interval_Top")

Subset.Genera <- function(df){
  genera <- data[(df$Latitude <= data$MaxLatRel  & df$Latitude >= data$MinLatRel  & df$Depth <= data$MaxDepth & df$Depth >= data$MinDepth  & data$Interval_Top == df$Interval_Top),] 
  genera
}
plot.gen <- ddply(lat.depths,.(Latitude,Depth,Interval_Top),each(Subset.Genera))

Prop.Ex <- function(df) mean(df$Ex)*100
Num.Gen <- function(df) length(df$Ex)

Gen.Ex <- ddply(plot.gen,.(Interval_Top,Latitude,Depth),each(Prop.Ex,Num.Gen))

## drop bins containing fewer than 3 genera
Gen.Ex$Prop.Ex <- ifelse(Gen.Ex$Num.Gen < 3,NA,Gen.Ex$Prop.Ex)

Stage <- c("Katian 1","Katian 2","Katian 3","Katian 4","Katian 5","Hirnantian","Rhuddanian","Aeronian 1","Aeronian 2")
Interval_Top <- c(452,448,446.5,445.5,444.7,443.4,439.0,438.0,436.0)

stages <- data.frame(Interval_Top,Stage)

Gen.Ex2 <- merge(Gen.Ex,stages,by="Interval_Top",all.x=FALSE)

Gen.Ex2$Stage <- factor(Gen.Ex2$Stage , levels=c("Katian 1","Katian 2","Katian 3", "Katian 4", "Katian 5", "Hirnantian", "Rhuddanian", "Aeronian 1","Aeronian 2"))

pdf("latitude.vs.depth.prop.extinction.pdf", width = 6, height = 7)
print(p <- ggplot(Gen.Ex2,aes(Latitude,Depth,fill=Prop.Ex)) + geom_tile(colour='black') + facet_wrap(~Stage,ncol=1) + scale_y_reverse(breaks =c(1,3,5)) + scale_fill_gradient2(low="dodgerblue",mid="chartreuse",high="orange",midpoint=30,"% Extinction",na.value="lightgray")+ theme_bw() + theme(strip.text.x = element_text(size = 9)) +theme(strip.background = element_blank()) + ylab("Depth (Brachiopod Assemblage Zone)") + xlab("Paleolatitude") + coord_cartesian(xlim=c(-90,40)) + scale_x_reverse(breaks=c(-60,-30,0,30)) + theme(legend.position = "top"))
dev.off()

##calculate lineage/mya extinction rates (using Foote, 2000 method) and make plot of extinction rate through time (figure 1 upper panel)

ExRate <- function(df) mean(df$Ex)
Foote.rate <- function(df) {
  df$Singleton <- ifelse(df$Ex==1 & df$Or ==1,1,0)
  df2 <-df[(df$Singleton == 0),]
  bt <- nrow(df2) - sum(df2$Ex)
  bL <- sum(df2$Ex)
  top <- mean(df$Interval_Top)
  bottom <- mean(df$Interval_Bottom)
  q <- -log(bt/(bt + bL))/(bottom-top)
  return(q)
}

exes <- ddply(data,.(Interval_Top,Interval_Bottom),each(ExRate,Foote.rate))

## define boundaries
stage_top <-c(456,452,448,447.5,447,446.5,446.2,445.8,445.5,445.3,445,444.7,444.05,443.4, 442.1,440.5,439.0,438.0,437.0,436.0,433.4,430.8)
stage_bottom <-c(460.9,456,452,448,447.5,447,446.5,446.2,445.8,445.5,445.3,445,444.7,444.05,443.4, 442.1,440.5,439.0,438.0,437.0,436.0,433.4)
bin <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

stage.tops <-c(456,452,448,446.5,445.5,444.7,443.4,439,436,430.2)
stage.bottoms <-c(460.5,456,452,448,446.5,445.5,444.7,443.4,439,436)
names <- c("San","K1","K2","K3","K4","K5","Hir.","Rhuddanian","Aeronian","Tel")
stages <- data.frame(stage.tops,stage.bottoms,names)
colnames(stages) <- c("stage.top","stage.bottom","stage.short")
stages$midpoint <- (stages$stage.top + stages$stage.bottom)/2

bases <- as.vector(stages$stage.bottom)

plot.data <- exes
bases <- as.vector(stages$stage.bottom)

plot.data <- plot.data[(plot.data$Interval_Top >= 436 & plot.data$Interval_Bottom <= 456),]
plot.data$midpoint <- (plot.data$Interval_Top + plot.data$Interval_Bottom)/2

Epochs <- data.frame(c(450,440),c("Late Ordovician","Early Silurian"))
colnames(Epochs) <- c("midpoint","epoch")

pdf("Extinction.Rates.pdf", width = 10, height = 3)
print(p <- ggplot(plot.data,aes(Interval_Top.new,Foote.rate)) + theme_bw() + scale_x_reverse() + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank()) + geom_hline(yintercept=2.6,size = .5, colour = "black")  + geom_line(data =  plot.data,aes(midpoint,Foote.rate),colour = "black") + geom_point(data=plot.data,aes(midpoint,Foote.rate),size=3,colour="black") + scale_x_reverse() + theme(panel.background = element_rect(colour="black", size =1, fill = NA)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank()) + ylab(expression("Per-capita extinction rate")) + xlab(expression("")) + theme(axis.title.x = element_text(size=16)) + theme(axis.title.y = element_text(size=16, angle = 90)) + theme(axis.text.x = element_text(size=12)) + theme(axis.text.y = element_text(size=12)) + geom_point(data=plot.data,aes(midpoint,Foote.rate),colour="black",size=2.75) + coord_cartesian(xlim = c(455.5,435.8),ylim = c(-.05,3.2)) + theme(legend.position="none") + geom_text(data = stages, aes(midpoint, 2.75, label=stage.short), size = 5, colour = "black")  + geom_hline(yintercept = -.25,colour = "black",size=.25) + geom_segment(data=stages,aes(x = stage.bottom,xend = stage.bottom,y=2.6,yend=2.9), size = .5, colour = "gray")  + geom_hline(yintercept=2.9,size = .5, colour = "black") + geom_text(data = Epochs, aes(midpoint, 3.05, label=epoch), size = 5, colour = "black") + geom_segment(aes(x = 443.4, xend=443.4,y=2.6,yend=3.2),colour="black"))
dev.off()

#make a column indicating whether each genus is a Lazarus (range-through) genus in the suceeding timeslice
data$stage <- factor(data$stage_top, levels=c(1:23))
data$stage2 <- factor(data$stage_top, levels=c(1:23))
data <- data[order(data$Genus, -data$stage_top),]
nexts <- data[c("Interval","Genus","stage_top")]
rows <- nrow(nexts)
nexts <- nexts[(2:nrow(nexts)),]
nexts[rows,] <- c(NA,NA,NA)
rownames(nexts) <- seq(nrow(nexts))
data$Gap.Next <- ifelse(data$Genus == nexts$Genus & data$Interval != (nexts$Interval -1),1,0)


## merge in clades as defined in the treatise
Order <- c("Athyridida","Atrypida","Billingsellida","Dictyonellida","Orthida","Orthotetida","Pentamerida","Productida","Protorthida","Rhynchonellida","Spiriferida","Strophomenida")
Clade <- c("Clade.3","Clade.3","Clade.1","Clade.4","Clade.2","Clade.1","Clade.3","Clade.1","Clade.2","Clade.3","Clade.3","Clade.1")
Clade.ID <- data.frame(Order,Clade)
data <- merge(data,Clade.ID,by="Order")

## calculate max occupancy for each interval

max.occ <- function(df) max(df$occupancy)
maxoccs <- ddply(data,.(Interval),each(max.occ))
data <- merge(data,maxoccs)

## round predictors to reduce false precision and so that number of categories is similar in all cases

round.fac <- 10

data$log.richness <- round(log(data$richness),0)
data$occupancy <- data$occupancy/data$max.occ
data$occupancy <- round((data$occupancy*round.fac),0)*round.fac
data$Prop.Cratonic <- round((data$Prop.Cratonic*round.fac),0)*round.fac
data$PropTrunc <- round((data$PropTrunc*round.fac),0)*round.fac
data$MaxLat <- ceiling(data$MaxLat/round.fac)*round.fac
data$MinLat <- floor(data$MinLat/round.fac)*round.fac
data$AbsLatRange <- round((data$AbsLatRange/round.fac),0)*round.fac
#data$AbsLatRange <- data$AbsLatRange + (runif(length(data$DepthRange))/10)
data$MeanAbsLat <- round(data$MeanAbsLat /round.fac,0) *round.fac
data$gcd <- round(data$gcd/1000,0)*1000
data$log.richness <- round(data$log.richness,1)
data$MeanDepth <- data$MeanDepth
data$DepthRange <- round(data$DepthRange,1)
data$MinDepth <- round(data$MinDepth,0)
data$MaxDepth <- round(data$MaxDepth,0)
data$Prop.North <- round((data$Prop.North*round.fac),0)*round.fac
data$cooling.latrange <- round(data$cooling.latrange/round.fac,0)*round.fac
data$cooling.latrange <- ifelse(data$cooling.latrange==91,90,data$cooling.latrange)
data$cooling.start.latrange <- round(data$cooling.start.latrange/round.fac,0)*round.fac
data$warming.latrange <- round(data$warming.latrange/round.fac,0)*round.fac
data$warming.latrange <- ifelse(data$warming.latrange==91,90,data$warming.latrange)
data$warming.start.latrange <- round(data$warming.start.latrange/round.fac,0)*round.fac
data$GH.IH.extinct <- ifelse(data$cooling.latrange==0,1,0)
data$IH.GH.extinct <- ifelse(data$warming.latrange==0,1,0)

## tabulate number of extinctions per timeslice
numex <- function(df) sum(df$Ex)
means <- ddply(data,.(stage_top),each(numex))
data <- merge(data,means)

## remove bins with too few extinctions
data <- data[(data$numex > 10),]

## make data frame with only variables to include in models
data2 <- data[c("stage_top","Ex","log.richness","gcd","Prop.Cratonic","PropTrunc","MaxDepth","MinDepth", "DepthRange","MaxLat","MinLat","AbsLatRange","occupancy","continents")]


## now use "gbm.step" function from dismo package to select "best" gbm model for each timeslice. This takes a few minutes to run

make.gbm <- function(df) {
  
  ##make weights proportional to class frequency
  obs <- length(df$Ex)
  exes <- sum(df$Ex)
  ExFreq <- (exes/obs)
  SurFreq <- (1-(exes/obs))
  MaxFreq <- max(ExFreq,SurFreq)
  ExWeight <- 1/(ExFreq/MaxFreq)
  SurWeight <- 1/(SurFreq/MaxFreq)
  train.weights <- ifelse(df$Ex==1,ExWeight,SurWeight)
  
  ## model for all analyses
  gbm.mod <- gbm.step(data=df , gbm.x = 3:14, gbm.y = 2,family = "bernoulli", tree.complexity = 5,learning.rate = 0.001, bag.fraction = .6,site.weights=train.weights,n.folds=20,n.trees=100,step.size=100,max.trees=10000,silent=FALSE,n.minobsinnode=5)
  
  return(gbm.mod)}

gbm.models <- dlply(data2,.(stage_top),failwith(NULL,make.gbm))

## make dataframe of partial dependences
partial.deps <- list()
index <- names(gbm.models)
for(i in index){
  model <- gbm.models[[i]]
  predictors <- model$var.names
  preds <- list()
  for(j in predictors){
    partials <- plot.gbm(model, i.var = j, n.trees = model$n.trees, continuous.resolution = 10, return.grid = TRUE, type="response") 
    partials$predictor <- rep(j,nrow(partials))
    colnames(partials) <- c("value","risk","var")
    preds[[j]] <- partials
  }
  preds <- do.call("rbind",preds)
  preds$interval <- i
  partial.deps[[i]] <- preds
}

partial.deps[["428"]] <- NULL
gbm.models[["428"]] <- NULL


partial.deps <- do.call("rbind",partial.deps)  

## make dataframe of model performance stats
perfstats <- list()
index <- names(gbm.models)
for (i in index){
  cv.correlation <- gbm.models[[i]]$cv.statistics$correlation.mean
  cv.correlation.se <- gbm.models[[i]]$cv.statistics$correlation.se
  cv.AUC <- gbm.models[[i]]$cv.statistics$discrimination.mean
  cv.AUC.se <- gbm.models[[i]]$cv.statistics$discrimination.se
  cv.deviance <- gbm.models[[i]]$cv.statistics$deviance.mean
  cv.deviance.se <- gbm.models[[i]]$cv.statistics$deviance.se
  cv.threshold <- gbm.models[[i]]$cv.statistics$cv.threshold
  cv.threshold.se <- gbm.models[[i]]$cv.statistics$cv.threshold.se
  interval <- i
  cv.correlation <- data.frame( cv.correlation, cv.correlation.se,cv.AUC,cv.AUC.se,cv.deviance, cv.deviance.se,cv.threshold,cv.threshold.se,interval)
  perfstats[[i]] <- cv.correlation
}
perfstats <- do.call("rbind",perfstats)

## make dataframe of relative influence
relinf <- list()
index <- names(gbm.models)
for (i in index){
  relative.influence <- data.frame(gbm.models[[i]]$contributions)
  relative.influence$interval <- rep(i,nrow(relative.influence))
  relinf[[i]] <- relative.influence
}
relative.influence <- do.call("rbind",relinf)


## make merged dataframe, multiple relinf by cv.correlation
plot.relinf <- merge(relative.influence,perfstats,by="interval")

## merge partial dependence and relative influence data.frames, make weighted partial dependence plot
all.plot.data <- merge(plot.relinf,partial.deps,by=c("interval","var"))
all.plot.data$rel.inf.class <- ifelse(all.plot.data$rel.inf > 5,2,1)
var.names <- levels(all.plot.data$var)

## assign interval and predictor names and reorder
interval2 <- c("Telychian 1","Aeronian 2","Aeronian 1", "Rhuddanian","Hirnantian","Katian 5","Katian 4","Katian 3","Katian 2","Katian 1")
interval <- c(433.4,436,438,439,443.4,444.7,445.5,446.5,448,452)
ints <- data.frame(interval,interval2)

all.plot.data <- merge(all.plot.data,ints,by="interval")

all.plot.data$interval2 <- factor(all.plot.data$interval2 , levels=c("Katian 1","Katian 2","Katian 3", "Katian 4","Katian 5", "Hirnantian", "Rhuddanian","Aeronian 1", "Aeronian 2", "Telychian 1"))

predictor <- all.plot.data$var
predictor <- revalue(predictor,c("AbsLatRange"="Paleolatitude Range","continents"="Num. Paleocontinents","DepthRange"="Depth Range","gcd" = "Great Circle Dist.","log.richness"="Log Richness","MaxDepth"="Maximum Depth","MaxLat"="Max. Paleolatitude","MinDepth"="Minimum Depth","MinLat"="Min. Paleolatitude","occupancy"="Grid Cell Occupancy","Prop.Cratonic"="% Cratonic","PropTrunc"="% Discontinuous"))
all.plot.data$predictor <- predictor

all.plot.data$predictor <- factor(all.plot.data$predictor, levels=c("Log Richness","Great Circle Dist.","Grid Cell Occupancy","Num. Paleocontinents","% Discontinuous", "% Cratonic","Min. Paleolatitude","Max. Paleolatitude","Paleolatitude Range","Minimum Depth","Maximum Depth","Depth Range"))

## rescale values to max

maxval <- function(df) max(df$value)
maxvals <- ddply(all.plot.data,.(predictor),each(maxval))
all.plot.data2 <- merge(all.plot.data,maxvals)
all.plot.data2$value2 <- all.plot.data2$value/all.plot.data2$maxval

all.plot.data2 <- all.plot.data2[(all.plot.data2$predictor != "Pedunculate" & all.plot.data2$predictor != "Mean Paleolatitude" & all.plot.data2$predictor != "Mean Depth"),]

## make partial dependence plot (Figure 1 lower panel)

pdf("Partial.Dependence.pdf", width = 10, height = 8)
print(p <- ggplot(all.plot.data2,aes(value2,risk,colour=rel.inf)) + geom_hline(yintercept=.5,colour="black") + geom_line(lineend = "round",size=1.5) + facet_grid(predictor~interval2,scales="free_x")  + theme_bw() + ylab("Marginal effect of predictor on extinction risk") + xlab("Predictor value (proportion of maximum)") + theme(strip.background = element_blank()) + theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=10,angle=0)) + scale_linetype_manual(values=c(3,1),guide=FALSE) + theme(axis.title.x = element_text(size=12,colour="black")) + theme(axis.title.y =  element_text(size=12, angle = 90,colour="black")) + scale_y_continuous(limits=c(.15,.8),breaks=c(.25,.5,.75))  + guides(colour=guide_legend(title="Relative influence of predictor on extinction risk:",override.aes = list(size=3))) + guides(size=guide_legend(title="% relative influence of predictor on extinction risk:")) + theme(legend.direction = "horizontal", legend.position = "bottom") + theme(legend.text=element_text(size=10)) +  theme(legend.title=element_text(size=12,face="plain")) + scale_colour_gradient2(low="black",high="orange",mid="red",midpoint=15) + theme(axis.text.x=element_text(size=8,angle=-90,hjust=0),axis.text.y=element_text(size=8))) 
dev.off()

## make partial dependence plot just for major predictors in Katian 5
all.plot.data3 <- all.plot.data[(all.plot.data$interval == "444.7"),]
all.plot.data3 <- all.plot.data3[(all.plot.data3$predictor == "% Cratonic" | all.plot.data3$predictor == "% Discontinuous" | all.plot.data3$predictor == "Paleolatitude Range" | all.plot.data3$predictor == "Minimum Depth"),]

pdf("Katian.5.Partial.Dependence.pdf", width = 9, height = 8)
print(p <- ggplot(all.plot.data3,aes(value,risk,colour=rel.inf)) + geom_hline(yintercept=.5,colour="black") + geom_line(lineend = "round",size=1.5) + facet_wrap(~predictor,scales="free_x")  + theme_linedraw() + ylab("Marginal effect of predictor on extinction risk") + xlab("Predictor value") + theme(strip.background = element_blank()) + theme(strip.text.x = element_text(size=10,colour="black")) + scale_linetype_manual(values=c(3,1),guide=FALSE) + theme(axis.title.x = element_text(size=12,colour="black")) + theme(axis.title.y =  element_text(size=12, angle = 90,colour="black")) + scale_y_continuous(limits=c(.15,.8),breaks=c(.25,.5,.75))  + guides(colour=guide_legend(title="Relative influence of predictor on extinction risk:",override.aes = list(size=3))) + guides(size=guide_legend(title="% relative influence of predictor on extinction risk:")) + theme(legend.direction = "horizontal", legend.position = "bottom") + theme(legend.text=element_text(size=10)) +  theme(legend.title=element_text(size=12,face="plain")) + scale_colour_gradient2(low="black",high="orange",mid="red",midpoint=15) + theme(axis.text.x=element_text(size=8,hjust=0),axis.text.y=element_text(size=8)))
dev.off()

##make Katian 5 relative influence plot
K5influence <- all.plot.data[, c("interval2","var", "predictor", "rel.inf")]
K5influence <- unique(K5influence[(K5influence$interval2 == "Katian 5"),])
K5influence <- K5influence[order(K5influence$rel.inf),]
names <- c(as.list(K5influence$predictor))
names <- c(do.call("cbind",names)) 
K5influence$predictor <- factor(K5influence$predictor,levels(K5influence$predictor)[names])

## make list of Katian5 factors to retain (only those with more influence than expected)
K5.med <- median(K5influence$rel.inf)
K5.retain  <- K5influence[(K5influence$rel.inf >= 100/12),]
K5.names <- c(as.list(as.character(K5.retain$var)))
K5.names <- c(do.call("cbind",K5.names))

K5.includes1 <- c("Ex",K5.names)
K5.includes2 <- c("Ex",K5.names,"GH.IH.extinct","IH.GH.extinct")

K5.data <- data[(data$Interval_Top == 444.7),]

K5.data1 <- K5.data[K5.includes1]
K5.data2 <- K5.data[K5.includes2]

## further cull the list of variables, include only those variables with conditional average p values of < .05 in the set of all possible models with deltaAICc values < 4

model.full <- glm(Ex~.,data = K5.data1, family="binomial",na.action="na.fail")
d.mod <- dredge(model.full,rank = "AICc",beta=FALSE)
avg.mod <- model.avg(d.mod, subset = delta < 4)
summary(avg.mod)

## now include these variables in three alternative models, two of which include predictions from climate change scenarios (greenhouse to icehouse and icehouse to greenhouse)

Cand.mod <- list()
Cand.mod[[1]] <- glm(Ex ~ AbsLatRange  + MinDepth + PropTrunc + Prop.Cratonic, data =K5.data2,family=binomial(link = "logit"))
Cand.mod[[2]] <- glm(Ex ~ AbsLatRange  + MinDepth + PropTrunc + Prop.Cratonic + GH.IH.extinct, data = K5.data2,family=binomial(link = "logit"))
Cand.mod[[3]] <- glm(Ex ~ AbsLatRange + MinDepth  + PropTrunc + Prop.Cratonic+ IH.GH.extinct, data = K5.data2,family=binomial(link = "logit"))

Modnames = c("Neutral","Greenhouse-Icehouse","Icehouse-Greenhouse")

mod.avg <- aictab(cand.set = Cand.mod, modnames = Modnames)

writeLines(capture.output(mod.avg),con="Table.1.txt")


## make classification tree (Figure 3) using the variables in the best model

K5.data3 <- K5.data2
K5.data3$Extinct <- ifelse(K5.data3$Ex == 1,"extinct","survive")
K5.data3$Prediction <- ifelse(K5.data3$GH.IH.extinct == 1,"extinct","survive")
K5tree <- rpart(Extinct ~ PropTrunc + Prop.Cratonic + Prediction + AbsLatRange + MinDepth, data=K5.data3, na.action = na.pass, method = "class",xval=20)

y <- ifelse(K5tree$frame$yval==1,K5tree$frame$dev/K5tree$frame$wt,1-K5tree$frame$dev/K5tree$frame$wt)

max <- 100
min <- 0
cols <- rainbow(99,end=.5,.65)[
  ifelse(y > .5, (y-min)  * (50-1)  / (y[1]-min) + 1,
         (y-min)  * (50-1)  / (y[1]-min) + 1)]

pdf("Classification.Tree.pdf", width = 7, height = 7)
prp(K5tree,extra=105,branch.type=5,box.col=cols,branch.col=cols)
dev.off()



