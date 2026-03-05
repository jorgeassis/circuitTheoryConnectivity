## ------------------------------------------------------
## ------------------------------------------------------
##
## Biodiversity Data Science, Portugal
## Website: https://www.biodiversitydatascience.com
## Copyright: 2025
## Contact: jorgemfa@gmail.com
##
## ------------------------------------------------------
##
## Notes:
##   
## https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/07-1861.1
## https://www.rdocumentation.org/packages/gdistance/versions/1.6.4/topics/transition
## https://www.nature.com/articles/s43247-022-00569-5
## https://link.springer.com/article/10.1007/s10980-023-01690-2
##
## ---------------------------

## set working directory

setwd(".")

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

nCores <- 5

dataFolder <- "copernicusData"
dataAggFactor <- 1

ea_crs <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
wgs84_crs <- "EPSG:4326" 
extent <- c(-180,180,-90,90)
temporalWindow <- c("2010-01-01", "2020-12-31")

graphDirected <- TRUE
graphLogWeights <- TRUE
graphCollapseFun <- "mean"
disperalType <- "passive" # active passive

resultsFolder <- paste0("../Results/",disperalType,"Dispersal",ifelse(graphDirected,"Directed",paste0("UndirectedFun",graphCollapseFun)),ifelse(graphLogWeights,"LogWeights",""))

## ---------------------------
## ---------------------------

cm <- setupCopernicusEnvironment("jassis", "OwglJEYl_1")

## ------------------------------------------------------
## ------------------------------------------------------

if( ! dir.exists(resultsFolder) ) { dir.create(resultsFolder) }

## ----------------------------
## Define data structure

temporalSteps <- seq(as.Date(temporalWindow[1]),as.Date(temporalWindow[2]),by="month")
template <- rast(xmin = -180, xmax = 180,ymin = -90, ymax = 90,resolution = 1/12,crs = "EPSG:4326" )

copernicusData <- extractData2(
  cm = cm,
  variables = c("eastwardCurrentDirection", "northwardCurrentDirection"),
  lon = extent[1:2],
  lat = extent[3:4],
  startDate = as.character(temporalSteps)[1], endDate = as.character(temporalSteps)[1],
  depth = c(1), 
  frequency = "monthly",
  data_output = dataFolder )

shape <- rast(list.files(dataFolder, pattern = "_uo_", full.names = T)[1])
shape <- aggregate(shape, fact = dataAggFactor, fun = mean)
shape[!is.na(shape)] <- 1
names(shape) <- "ocean"

crs(shape) <- wgs84_crs
shape <- resample(shape, template, method = "near")
shape <- project(shape, ea_crs, method = "near")

oceanCells <- which(!is.na(values(shape)[,1]))
oceanCellsAdjacent <- adjacent(shape, cells = oceanCells, directions = 8, pairs = TRUE)
oceanCellsAdjacent <- as.data.frame(oceanCellsAdjacent)
oceanCellsAdjacent <- oceanCellsAdjacent[which(oceanCellsAdjacent[,2] %in% oceanCells),]

plot(shape)
unlink(paste0(dataFolder,"/*"), recursive = TRUE)

## ----------------------------
## Make spatial objects

spatialGrid <- list(shape=raster(shape), grid=data.frame(cell=oceanCells,x=xyFromCell(shape,oceanCells)[,1],y=xyFromCell(shape,oceanCells)[,2]))

# Make big matrix
if( ! dir.exists(paste0(resultsFolder,"/BM")) ) { dir.create(paste0(resultsFolder,"/BM")) }
file.remove( list.files( paste0(resultsFolder,"/BM") , pattern = "bigMemoryMatrix", full.names = TRUE) )
matrixBM <- big.matrix( nrow = nrow(oceanCellsAdjacent) , ncol = length(temporalSteps) , backingpath = paste0(resultsFolder,"/BM/"), backingfile = paste0("bigMemoryMatrix.bin"), descriptorfile = paste0("bigMemoryMatrix.desc"))

save(oceanCells,file=paste0(resultsFolder,"/oceanCells.RData"))
save(spatialGrid,file=paste0(resultsFolder,"/spatialGrid.RData"))
save(oceanCellsAdjacent,file=paste0(resultsFolder,"/oceanCellsAdjacent.RData"))

oceanCellsAdjacentVect <- paste(oceanCellsAdjacent$from, oceanCellsAdjacent$to)

## ----------------------------
## Determine conductance across temporal window

for( i in 1:length(temporalSteps) ) {
  
  cat("Processing month ",i," of ",length(temporalSteps),"\n")
  
  copernicusData <- extractData2(
    cm = cm,
    variables = c("eastwardCurrentDirection", "northwardCurrentDirection"),
    lon = extent[1:2],
    lat = extent[3:4],
    startDate = as.character(temporalSteps)[i], endDate = as.character(temporalSteps)[i],
    depth = c(1), 
    frequency = "monthly",
    data_output = dataFolder )
  
  uComponent <- rast(list.files(dataFolder, pattern = "_uo_", full.names = T))
  vComponent <- rast(list.files(dataFolder, pattern = "_vo_", full.names = T))
  
  uComponent <- aggregate(uComponent, fact = dataAggFactor, fun = mean)
  vComponent <- aggregate(vComponent, fact = dataAggFactor, fun = mean)

  speed <- sqrt(uComponent^2 + vComponent^2)
  direction <- (270 - (atan2(vComponent, uComponent) * 180 / pi)) %% 360

  crs(speed) <- wgs84_crs
  crs(direction) <- wgs84_crs
  
  speed <- resample(speed, template, method = "bilinear")
  speed <- project(speed, ea_crs, method = "bilinear")
  
  direction <- resample(direction, template, method = "bilinear")
  direction <- project(direction, ea_crs, method = "bilinear")
  
  currents <- c(direction, speed)
  names(currents) <- c("direction", "speed")
  
  Conductance <- flow.dispersion(x=stack(currents), type = disperalType, fun=cost.FMGS, output="transitionLayer")
  Conductance <- geoCorrection(Conductance, type="c", multpl=FALSE)
  
  transitionMatrix <- transitionMatrix(Conductance)
  edges <- summary(transitionMatrix)
  names(edges) <- c("from", "to", "weight")
  edges <- edges[which(edges$weight > 0),]
  
  idx <- match( paste(edges$from, edges$to) , oceanCellsAdjacentVect )
  matrixBM[idx[which(!is.na(idx))],i] <- edges$weight[which(!is.na(idx))]
  
  # ---------
  
  unlink(paste0(dataFolder,"/*"), recursive = TRUE)
  
}

## ----------------------------------------------------------
## ----------------------------------------------------------
## Aggregate conductance to average

oceanCells <- loadRData(paste0(resultsFolder,"/oceanCells.RData"))
oceanCellsAdjacent <- loadRData(paste0(resultsFolder,"/oceanCellsAdjacent.RData"))
oceanCellsAdjacentVect <- paste(oceanCellsAdjacent$from, oceanCellsAdjacent$to)

spatialGrid <- loadRData(paste0(resultsFolder,"/spatialGrid.RData"))
shape <- spatialGrid$shape

conductanceBM <- attach.big.matrix(paste0(resultsFolder,"/BM/bigMemoryMatrix.desc"))

rows <- nrow(conductanceBM)
cols <- ncol(conductanceBM)
chunk_size <- 10000
row_means <- numeric(rows)

for (i in 1:ceiling(rows / chunk_size)) {
  cat("Processing chunk ", i, " of ", ceiling(rows / chunk_size), "\n")
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, rows)
  chunk <- conductanceBM[start_row:end_row, ]
  chunk_means <- rowMeans(chunk)
  row_means[start_row:end_row] <- chunk_means
}

averageconductanceCells <- data.frame(oceanCellsAdjacent,vect=oceanCellsAdjacentVect,conductance=row_means)
save(averageconductanceCells,file=paste0(resultsFolder,"/conductanceCells.RData"))

## ------------------------------------
## ------------------------------------

spatialGrid <- loadRData(paste0(resultsFolder,"/spatialGrid.RData"))
shape <- spatialGrid$shape
oceanCells <- loadRData(paste0(resultsFolder,"/oceanCells.RData"))
averageconductanceCells <- loadRData(paste0(resultsFolder,"/conductanceCells.RData"))

# Get the last non-NA cell per row in a Mollweide raster

spatialGridDF <- spatialGrid$grid
oceanCellsAdjacentExtra <- data.frame()

for ( i in 1:length(unique(spatialGridDF$y)) ) {

  cat("Processing extra cell ",i," of ",length(unique(spatialGridDF$y)),"\n")
  
  y <- unique(spatialGridDF$y)[i]
  y_m1 <- unique(spatialGridDF$y)[i-1]
  y_p1 <- unique(spatialGridDF$y)[i+1]
  
  spatialGridDF.y <- spatialGridDF[spatialGridDF$y == y,]
  cell.i <- spatialGridDF.y[which.min(spatialGridDF.y$x),"cell"]
  cell.j <- spatialGridDF.y[which.max(spatialGridDF.y$x),"cell"]
  
  spatialGridDF.y.m_1 <- spatialGridDF[spatialGridDF$y == y_m1,]
  cell.j.y_m1 <- spatialGridDF.y.m_1[which.max(spatialGridDF.y.m_1$x),"cell"]
  if(length(cell.j.y_m1) == 0) { cell.j.y_m1 <- cell.j }
  
  spatialGridDF.y.p_1 <- spatialGridDF[spatialGridDF$y == y_p1,]
  cell.j.p_m1 <- spatialGridDF.y.p_1[which.max(spatialGridDF.y.p_1$x),"cell"]
  if(length(cell.j.p_m1) == 0) { cell.j.p_m1 <- cell.j }
  
  cond <- unlist(averageconductanceCells[ which(averageconductanceCells$to == cell.i), "conductance" ])
  
  oceanCellsAdjacentExtra <- rbind(oceanCellsAdjacentExtra, 
                                   data.frame(from=cell.i, 
                                              to=cell.j,
                                              conductance=max(cond)) ,
                                   data.frame(from=cell.i, 
                                              to=cell.j.y_m1,
                                              conductance=max(cond[cond != 0])) ,
                                   data.frame(from=cell.i, 
                                              to=cell.j.p_m1,
                                              conductance=max(cond[cond != 0]))
                                   )
}

oceanCellsAdjacentExtra <- oceanCellsAdjacentExtra[complete.cases(oceanCellsAdjacentExtra),]
oceanCellsAdjacentExtra <- oceanCellsAdjacentExtra[!duplicated(oceanCellsAdjacentExtra[,1:2]),]
oceanCellsAdjacentExtra <- oceanCellsAdjacentExtra[oceanCellsAdjacentExtra$from != oceanCellsAdjacentExtra$to,]
oceanCellsAdjacentExtra <- oceanCellsAdjacentExtra[oceanCellsAdjacentExtra$from %in% oceanCells,]
oceanCellsAdjacentExtra <- oceanCellsAdjacentExtra[oceanCellsAdjacentExtra$to %in% oceanCells,]

oceanCellsAdjacentExtra$vect <- paste(oceanCellsAdjacentExtra$from, oceanCellsAdjacentExtra$to)

averageconductanceCells <- rbind(averageconductanceCells, oceanCellsAdjacentExtra )
sum(duplicated(averageconductanceCells[,1:2]))

# -------

averageconductanceCells <- averageconductanceCells[,c("from","to","conductance")]
averageconductanceCells$from <- as.character(averageconductanceCells$from)
averageconductanceCells$to <- as.character(averageconductanceCells$to)
averageconductanceCells <- averageconductanceCells[averageconductanceCells$from != averageconductanceCells$to,]

averageconductanceCells$cost <- log( 1 / averageconductanceCells$conductance )
averageconductanceCells <- averageconductanceCells[averageconductanceCells$cost != Inf,]
save(averageconductanceCells,file=paste0(resultsFolder,"/conductance.RData"))

write.table(averageconductanceCells[,c("from","to","cost")], 
            file = paste0(resultsFolder,"/conductance.txt"), 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ------------------------------
# ------------------------------