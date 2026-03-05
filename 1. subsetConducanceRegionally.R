## ------------------------------------------------------
## ------------------------------------------------------
##
## Biodiversity Data Science, Portugal
## Website: https://www.biodiversitydatascience.com
## Copyright: 2025
## Contact: jorgemfa@gmail.com
##
## ------------------------------------------------------

## set working directory

## ---------------------------
## Clean workspace

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ---------------------------
## Load libraries and functions

source("functions.R")

## ---------------------------
## Options

nCores <- 10
ea_crs <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
resultsFolder <- "../Results/passiveDispersalDirectedLogWeights"
subsetExtent <- c(-35,42.5,25,85)
subsetResultsFolder <- "../Results/ConnectivitySubset"

## ---------------------------

spatialGrid <- loadRData(paste0(resultsFolder,"/spatialGrid.RData"))$grid
spatialRaster <- loadRData(paste0(resultsFolder,"/spatialGrid.RData"))$shape
conductance <- loadRData(paste0(resultsFolder,"/conductance.RData"))

# Create a polygon from the WGS84 extent
bbox_poly <- rbind(c(subsetExtent[1], subsetExtent[3]), c(subsetExtent[2], subsetExtent[3]), 
                   c(subsetExtent[2], subsetExtent[4]), c(subsetExtent[1], subsetExtent[4]), 
                   c(subsetExtent[1], subsetExtent[3]))

# Project the polygon to Mollweide
subsetExtentProjected <- project(vect(bbox_poly, type="polygons", crs="+proj=longlat +datum=WGS84"), ea_crs)
subsetExtentProjected <- ext(subsetExtentProjected)
subsetExtentProjected <- as.vector(subsetExtentProjected)

# Subset cells within the projected extent
spatialGridSubset <- spatialGrid[spatialGrid$x >= subsetExtentProjected[1] & spatialGrid$x <= subsetExtentProjected[2] & spatialGrid$y >= subsetExtentProjected[3] & spatialGrid$y <= subsetExtentProjected[4], c("x","y","cell")]
cellsSubset <- spatialGridSubset[,"cell"]

spatialRaster[cellsSubset] <- 2
plot(spatialRaster)

spatialRasterSubset <- crop(spatialRaster, subsetExtentProjected)
plot(spatialRasterSubset)

conductance <- conductance[conductance$from %in% cellsSubset & conductance$to %in% cellsSubset,]

if( ! dir.exists(subsetResultsFolder) ) {
  dir.create(subsetResultsFolder)
  cat("Folder ",subsetResultsFolder," created","\n")
}

save(spatialGridSubset,file=paste0(subsetResultsFolder,"/spatialGrid.RData"))
save(spatialRaster,file=paste0(subsetResultsFolder,"/spatialRaster.RData"))
save(spatialRasterSubset,file=paste0(subsetResultsFolder,"/spatialRasterRegion.RData"))
save(conductance,file=paste0(subsetResultsFolder,"/conductance.RData"))
