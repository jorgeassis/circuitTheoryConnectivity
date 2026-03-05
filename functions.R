
library(rWind)
library(raster)
library(terra)
library(gdistance)
library(igraph)
library(foreach)
library(doParallel)
library(sf)
library(rnaturalearth)
library(dggridR)
library(ggplot2)
library(ggthemes)
library(exactextractr)
library(CopernicusDownloadR)
library(bigmemory)
library(clue)
library(data.table)
library(reticulate)

# -----------

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# -----------

extractData2 <- function (cm, variables, lon, lat, startDate, endDate, depth, frequency, csv_path = NULL, data_output = NULL) 
{
  if (is.null(csv_path)) {
    csv_path <- system.file("extdata", "CopernicusVariables&Datasets.csv", 
                            package = "CopernicusDownloadR")
  }
  if (!file.exists(csv_path)) {
    stop("The CSV file containing variable and dataset information was not found.")
  }
  data <- read.csv(csv_path)
  if (length(lon) == 1) {
    minLon <- lon
    maxLon <- lon
  }
  else if (length(lon) == 2) {
    minLon <- lon[1]
    maxLon <- lon[2]
  }
  else {
    stop("Longitude must be either a single value or a range of length 2.")
  }
  if (length(lat) == 1) {
    minLat <- lat
    maxLat <- lat
  }
  else if (length(lat) == 2) {
    minLat <- lat[1]
    maxLat <- lat[2]
  }
  else {
    stop("Latitude must be either a single value or a range of length 2.")
  }
  if (length(depth) == 1) {
    minDepth <- depth
    maxDepth <- depth
  }
  else if (length(depth) == 2) {
    minDepth <- depth[1]
    maxDepth <- depth[2]
  }
  else {
    stop("Depth must be either a single value or a range of length 2.")
  }
  metadata_list <- list()
  for (variable in variables) {
    variable_info <- data[data$variable == variable, ]
    if (nrow(variable_info) == 0) {
      stop(paste("The variable", variable, "does not exist in the dataset."))
    }
    if (frequency == "daily") {
      dataset_id <- variable_info$dataset_id_daily
    }
    else if (frequency == "monthly") {
      dataset_id <- variable_info$dataset_id_monthly
    }
    else {
      stop("Invalid frequency. Choose 'daily' or 'monthly'.")
    }
    variable_name_copernicus <- variable_info$variableNameCopernicus
    cat("Downloading data for variable:", variable, "\n")
    result <- cm$subset(dataset_id = dataset_id, start_datetime = startDate, 
                        end_datetime = endDate, variables = list(variable_name_copernicus), 
                        minimum_longitude = minLon, maximum_longitude = maxLon, 
                        minimum_latitude = minLat, maximum_latitude = maxLat, 
                        minimum_depth = minDepth, maximum_depth = maxDepth, 
                        output_directory = data_output)
    metadata_list[[variable]] <- list(variable = variable, 
                                      dataset_id = dataset_id, variableNameCopernicus = variable_name_copernicus, 
                                      startDate = startDate, endDate = endDate, minimum_longitude = minLon, 
                                      maximum_longitude = maxLon, minimum_latitude = minLat, 
                                      maximum_latitude = maxLat, minimum_depth = minDepth, 
                                      maximum_depth = maxDepth, resolution = ifelse(!is.na(variable_info$resolution), 
                                                                                    variable_info$resolution, "N/A"), doi = ifelse(!is.na(variable_info$doi), 
                                                                                                                                   variable_info$doi, "N/A"), download_date = Sys.Date())
  }
  metadata_df <- do.call(rbind, lapply(metadata_list, as.data.frame))
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  metadataFileName <- paste0("copernicus_metadata_", timestamp, 
                             ".csv")
  if (is.null(data_output)) {
    metadata_output <- file.path(getwd(), metadataFileName)
  }
  else {
    metadata_output <- file.path(data_output, metadataFileName)
  }
  write.csv(metadata_df, metadata_output, row.names = FALSE)
  cat("Data saved to:", metadata_output, "\n")
  return(result)
}

# ---------------------------------

membershipComplete <- function(membership,desiredN,fillValue=0) {
    if( sum(!desiredN %in% names(membership)) != 0 ) {
      missingNodes <- desiredN[which( ! desiredN %in% names(membership) )]
      vect <- rep(fillValue,length(missingNodes))
      names(vect) <- as.character(missingNodes)
      membership <- c(membership,vect)
      membership <- membership[sort(as.numeric(names(membership)), index.return=T)$ix]
    }
    return(membership)
}

# ---------------------------------

cluster_louvain_ensemble <- function(graph,repetitions=99,backendirectory=".",backendfile="temp_backend_file", nCores=1) {
    
  namesEdges <- V(graph)$name

  # ---
  
  file.remove( list.files( backendirectory , pattern = "temp_backend_file", full.names = TRUE) )
  matrixIteractionsBM <- big.matrix( nrow = length(namesEdges) , ncol = repetitions , backingpath = backendirectory, backingfile = paste0(backendfile,".bin"), descriptorfile = paste0(backendfile,".desc"))
  rm(matrixIteractionsBM)

  cl <- parallel::makeCluster(min(c(nCores,repetitions)))
  registerDoParallel(cl)

  parallelResults <- foreach(i = 1:repetitions, .combine = rbind, .export = "membershipComplete", .packages=c("igraph","bigmemory") ) %dopar% {
  
      matrixIteractionsBM <- attach.big.matrix(paste0(backendirectory,"/",backendfile,".desc"))

      communities <- cluster_louvain(graph, weights = E(graph)$conductance )
      membership <- communities$membership
      names(membership) <- communities$names
      membership <- membershipComplete(membership,namesEdges,fillValue=0)
      membership <- membership[sort(as.numeric(names(membership)), index.return=T)$ix]
      matrixIteractionsBM[,i] <- membership   

      rm(matrixIteractionsBM)
      return(NULL)
  
  }
  
  stopCluster(cl); rm(cl)
  closeAllConnections()
  gc(reset=TRUE)

  matrixIteractionsBM <- attach.big.matrix(paste0(backendirectory,"/",backendfile,".desc"))
  matrixIteractionsBM <- as.data.frame(as.matrix(matrixIteractionsBM))

  # Generate consensus

  N <- nrow(matrixIteractionsBM)
  K <- sapply(1:repetitions, function(x) { length(unique(matrixIteractionsBM[,x])) } )
  mainIter.i <- which(K == modal(K))[1]
  K <- modal(K)
  
  consensus <- matrix(0, nrow = N, ncol = repetitions)
  
  mainIter <- matrixIteractionsBM[, mainIter.i]
  mainIter <- as.numeric(factor(mainIter))
  for( i in 1:repetitions ) {
    iIter <- matrixIteractionsBM[, i]
    iIter <- as.numeric(factor(iIter))
    assigned.k <- 0
    for( j in unique(mainIter) ) {
      assigned.k <- assigned.k + 1
      alloc.j <- which(mainIter == j)
      modalAssignment <- as.numeric(names(sort(table(iIter[alloc.j]), decreasing = TRUE)[1]))
      consensus[which(iIter == modalAssignment), i] <- assigned.k
    }
  }
  
  stability <- sapply(1:N, function(x) { consensus.i <- consensus[x, ]; consensus.i <- consensus.i[consensus.i != 0]; sum(consensus.i == modal(consensus.i)[1]) / repetitions })
  membership <- sapply(1:N, function(x) { consensus.i <- consensus[x, ]; consensus.i <- consensus.i[consensus.i != 0]; modal(consensus.i, na.rm=T)[1] })
  membership <- as.numeric(factor(membership))
  names(membership) <- namesEdges

  if(  sum(membership == 0) != 0 ) {
      stop("Error: some nodes were not assigned to any cluster.")
  }

  meanStability <- mean(stability)
  sdStability <- sd(stability)
  minStability <- min(stability)
  maxStability <- max(stability)
    
  file.remove( list.files( backendirectory , pattern = "temp_backend_file", full.names = TRUE) )
  
  return(list(membership=membership,stability=stability,meanStability=meanStability,sdStability=sdStability,minStability=minStability,maxStability=maxStability))

}
    
