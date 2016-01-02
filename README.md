# Ord.Brach.Extinctions.Code
Data and analyses for Finnegan et al., in prep., : determinants of brachiopod extinction risk during Late Ordovician Mass Extinction

You will need the following R packages:

gdata,
plyr,
miscTools,
reshape2,
ggplot2,
fields,
PBSmapping,
AICcmodavg,
dismo,
MuMIn,
gbm,
rpart,
rpart.plot

Neccessary data files are:

Local.Ranges.csv (local stratigraphic ranges of Late Ordovician and Early Silurian brachiopod species & genera along with ancillary information)

Herrmann.Temps.csv (meridional mean annual sea surface temperature gradients from Herrmann et al., 2004, Response of Late Ordovician paleoceanography to changes in sea level, continental drift, and atmospheric pCO2: potential causes for long-term cooling and glaciation)

To reproduce analyses download ZIP with all files and modify line 17 in "Process.and.Analyze.R" script to set working directory to your local folder.
