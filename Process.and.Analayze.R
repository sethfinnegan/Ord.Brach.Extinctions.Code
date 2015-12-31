## This script performs the analyses reported in the main text and reproduces (in rough form) figures 1 and 3. It begins by reading in local ranges of genera from the raw data file ("Local.Ranges.csv") and processessing these data to extract global first and last appearances of each genus and extimates of the geographic distribution of the genus in each timeslice (LINES ## through ###). It then uses the distribution of each genus on each paleocontinental block and modeled SST estimates from Herrmann et al. 2004 to estimate the thermal tolerance range of each genus on each paleoncontinental block and then determine whether it would continue to be able to access coastal habitats within its thermal tolerance ranges in a modeled climate state shift (greenhouse to icehouse and, for comparison, icehouse to greenhouse)(LINES ## through ###).  Paleocontinental blocks are defined based on paleolatitudinal reconstructions from Torsvik (TODO: get proper reference) as figured in Rasmussen and Harper 2011. Blocks that were joined or so close that the distance would not have presented a substantial dispersal barrier (for example Avalonia and Baltica) were considered as single blocks.  See supplemental materials for further information about blocks used. The data frame output by this script is used as input for analyses of extinction patterns within each timeslice (see Appendix S2: "Selectivity.Analyses.R")

## To run this script it is neccessary to modify line 17 to set to your local directory

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
setwd(dir="~/Dropbox/Rasmussen and Harper brach database")
BrachData <- read.csv("Local.Ranges.csv",header=TRUE,stringsAsFactors=FALSE)

## choose whether to use Genera or Subgenera
BrachData$Genus.Use <- BrachData$Operative_Genus.Subgenus

## remove genera first appearing in Llandovery to avoid edge effects
BrachData <- BrachData[(BrachData$Appearance != "Llandovery"),]

## choose binning factor for lat measures (determines how coarsely range metrics are binned)
bin <- 1

## remove unresolved occurrences, unnamed occurrences, and inarticulates
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
Interval_Tops <- c(457,453,449,448.5,448,447.5,447.2,446.8,446.5,446.3,446,445.6,444.6,443.7,442.1,440.5,439,438,437,436,433.4,430.8,428.2)
Combined_Intervals <- data.frame(Interval,Interval_Tops)
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

Loc.Gen <- ddply(Expanded.BrachData,.(Interval_Tops,Easting,Northing),each(Local.Genera,Loc.Depth,Deep.Only))
Expanded.BrachData <- merge(Expanded.BrachData,Loc.Gen,by=c("Interval_Tops","Easting","Northing"))
Expanded.BrachData$North.Hemisphere <- ifelse(Expanded.BrachData$Northing > 0,1,0)

#round lats
Expanded.BrachData$Northing <- round(Expanded.BrachData$Northing,0)
Expanded.BrachData$AbsNorthing <- abs(Expanded.BrachData$Northing)
Expanded.BrachData$Easting <- round(Expanded.BrachData$Easting,0)

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
  if(length(na.omit(df$Northing))>1){
    events <- data.frame(EID=1:length(df$Northing), X = as.numeric(as.vector(df$Easting)),Y=as.numeric(as.vector(df$Northing)))
    events <- as.EventData(na.omit(events), projection='UTM')
    fc <- findCells(events, polys)
    richness <- length(unique(df$Genus.Species))
    eac <-length(unique(subset(fc,select=c(PID,SID)))[,1])
    maxdist <- max(rdist.earth(events[,2:3]))
    gcd <-ifelse(maxdist > 10, maxdist, 10)
    MaxLat <- max(df$Northing) + range.buffer
    MinLat <- min(df$Northing) - range.buffer
    MeanAbsLat <- round(mean(abs(df$Northing))/bin,0)*bin
    LatRange <- MaxLat - MinLat
    Single_occurrence <- 0
    MaxAbsLat <- max(abs(df$Northing)) + range.buffer
    MinAbsLat <- min(abs(df$Northing)) - range.buffer
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
    MaxLat <- max(df$Northing) + range.buffer
    MinLat <- min(df$Northing) - range.buffer
    MeanAbsLat <- round(mean(abs(df$Northing))/bin,0)*bin
    LatRange <- 0
    MaxAbsLat <- max(abs(df$Northing)) + range.buffer
    MinAbsLat <- min(abs(df$Northing)) - range.buffer
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

byGenus <- ddply(Expanded.BrachData,.(Interval,Interval_Tops,Order,Family,Genus.Use),ranges,.progress="text")

## Extract first and last occurrences of genera from the processed dataset
FAD <- function(df) max(df$Interval_Tops)
LAD <- function(df) min(df$Interval_Tops)

byGenus.2 <- ddply(byGenus,.(Genus.Use),each(FAD,LAD),.progress="text")

byGenus.3 <- merge(byGenus,byGenus.2,by = "Genus.Use")

## indicate whether genus has first or last appearance in a given interval

Ex <- ifelse(byGenus.3$Interval_Tops == byGenus.3$LAD,1,0)
Or <- ifelse(byGenus.3$Interval_Tops == byGenus.3$FAD,1,0)
Age <- byGenus.3$FAD - byGenus.3$Interval_Tops 
ExClass <- ifelse(Ex==1,"extinct","survive")
OrClass <- ifelse(Or==1,"origination","holdover")
byGenus.6 <- data.frame(byGenus.3,Ex,Or,Age,ExClass,OrClass)

##TODO: check whether this is really required in final version
write.table(byGenus.6,"Processed_Brachiopod_Data.csv",sep=",")


## rename range parameters etc.
gcd <- byGenus.6$gcd
richness <- byGenus.6$richness
occupancy <- byGenus.6$eac
continents <- byGenus.6$NumContinents
endemic <- byGenus.6$Endemic
MaxLat <- byGenus.6$MaxAbsLat
MinLat <- byGenus.6$MinAbsLat
MinLatRel <- byGenus.6$MinLat
MaxLatRel <- byGenus.6$MaxLat
LatRangeRel <- MaxLatRel - MinLatRel
LatRange <- byGenus.6$LatRange
AbsLatRange <- byGenus.6$AbsLatRange
MaxLatRange.ever <- MaxLat-MinLat 
MeanAbsLat <-  byGenus.6$MeanAbsLat
MeanDepth <- byGenus.6$MeanDepth
Single_occurrence <- byGenus.6$Single_occurrence
Age <- byGenus.6$Age
Interval <- byGenus.6$Interval
Interval_Top <- byGenus.6$Interval_Tops
MaxDepth <- byGenus.6$MaxDepth
MinDepth <- byGenus.6$MinDepth
DepthRange <- MaxDepth - MinDepth
Order <- as.factor(byGenus.6$Order)
Family <- as.factor(byGenus.6$Family)
Genus <- as.factor(byGenus.6$Genus.Use)
PropTrunc <- byGenus.6$PropTrunc
Pedicle.state <- byGenus.6$Pedicle.state
Median.Comm.Rich <- byGenus.6$Median.Comm.Rich
Mean.Loc.Depth <- byGenus.6$Mean.Loc.Depth
Prop.Cratonic <- byGenus.6$Prop.Cratonic
Ex <- byGenus.6$Ex
ExClass <- byGenus.6$ExClass
Or <- byGenus.6$Or
OrClass <- byGenus.6$OrClass
North.Hemisphere <- ifelse(byGenus.6$MaxLat > 0, 1,0)
Prop.North <- byGenus.6$Prop.North


final_data <- data.frame(Interval,Interval_Top,Ex,Or,ExClass,OrClass,richness,occupancy,Single_occurrence,continents,endemic,MaxLat,MinLat,MinLatRel,MaxLatRel,LatRange,LatRangeRel,AbsLatRange,MeanAbsLat,MeanDepth,MaxDepth,MinDepth,DepthRange,Age,gcd,Order,Family,Genus,Pedicle.state,Median.Comm.Rich,Mean.Loc.Depth,MaxLatRange.ever,Prop.Cratonic,North.Hemisphere,Prop.North,PropTrunc)


##flag singletons:
Num_Intervals <- function(df) length(df$Genus)
Int.Num <- ddply(final_data,.(Genus),each(Num_Intervals))
final_data <- merge(final_data,Int.Num,by="Genus")

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

## read in mean annual temperature gradients for greenhouse and icehouse states based on Herrmann et al., 2004 models.  Greenhouse/Warm: 15X PAL CO2, high sea level, Ashgill paleogeography, Icehouse/Cool: 8X PAL CO2, low sea level, Ashgill paleogeography

temps <- read.csv("Herrmann.Temps.csv",header=TRUE)
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

temps2 <- data.frame(temps[,1],temps[,9:10],temps[,11:12])
temps2 <- melt(temps2,id="temps...1.")
colnames(temps2) <- c("latitude","model","temp.")
temps2 <- drop.levels(temps2[(temps2$latitude >= 0),])
model <- levels(temps2$model)
states <- c("Icehouse summer SST","Icehouse winter SST","Greenhouse summer SST","Greenhouse winter SST")
states <- cbind(model,states)
temps2 <- merge(temps2,states)

## Make a temperature profile plot to include in supplemental materials.

pdf("~/Dropbox/Rasmussen and Harper brach database/H&R paper drafts/Temp.Profiles.pdf", width = 8, height = 4)
p <- ggplot(temps2,aes(latitude,temp.,colour=states,linetype=states))
p + geom_line() + scale_colour_manual(values =c("red","red","blue","blue"),name=" ") + scale_linetype_manual(values =c(1,2,1,2),name=" ") + theme_bw() + ylab("mean sea surface temperature") + xlab("degrees absolute latitude") + scale_x_continuous(breaks=c(0,30,60,90))
dev.off()

## get ranges of genera on each continent
Merged.Terrains<- function(df){
  # df <- final_data[23,]
  continents <- data.frame(unique(df$Merged.Terrains))
  return(continents)}
continent.ranges <- ddply(Expanded.BrachData,.(Genus.Use,Interval),each(Merged.Terrains),.progress="text")

Expanded.BrachData.A <- Expanded.BrachData
Expanded.BrachData.A$Northing <- Expanded.BrachData.A$Northing + range.buffer
Expanded.BrachData.B <- Expanded.BrachData
Expanded.BrachData.B$Northing <- Expanded.BrachData.B$Northing
Expanded.BrachData.C <- Expanded.BrachData
Expanded.BrachData.C$Northing <- Expanded.BrachData.C$Northing - range.buffer

Expanded.BrachData2 <- rbind(Expanded.BrachData.A,Expanded.BrachData.B,Expanded.BrachData.C)

Expanded.BrachData2 <- merge(Expanded.BrachData2,temps,by.x="Northing",by.y="latitude")

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
  Max.Lat <- max(df$Northing)
  Min.Lat <- min(df$Northing)
  outs <- data.frame(GH.gen.win.min,GH.gen.win.max,GH.gen.sum.min,GH.gen.sum.max, IH.gen.win.min,IH.gen.win.max,IH.gen.sum.min,IH.gen.sum.max,Max.Lat,Min.Lat)
  return(outs)}
gen.data <- ddply(Expanded.BrachData2,.(Genus.Species,Genus.Use,Interval,Merged.Terrains),each(gen.ranges),.progress="text")

max.cont.lat <- function(df) max(df$Max.Northing.Merged.450)
min.cont.lat <- function(df) min(df$Min.Northing.Merged.450)
min.loc.lat <- function(df) {
  Katian5 <- df[(df$Interval == 12),]
  min(Katian5$Northing)}
max.loc.lat <- function(df) {
  Katian5 <- df[(df$Interval == 12),]
  max(Katian5$Northing)}

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

final_data2 <- merge(final_data,range.changes,by=c("Genus","Interval"))

## Write the processed dataframe (as a csv file) that will then be read in by "Selectivity.Analyses.Final.R". 

write.csv(final_data2,"OS.brachs.processed.csv",row.names=FALSE)

