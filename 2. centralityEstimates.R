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

resultsFolder <- "../Results/ConnectivitySubset"
conductanceFile <- "conductance.txt"
conductanceFileRData <- "conductance.RData"

oceanCells <- loadRData(paste0(resultsFolder,"/spatialGrid.RData"))$cell
length(oceanCells)

# 70 000–100 000 Near-exact ranking (e approx 0.005); Publishable numbers (e approx 0.01)	15 000–25 000
nSamples <- 100000

## -----------

paste0(resultsFolder, "/conductance.txt")

## -----------------------------------------------------------------
## -----------------------------------------------------------------
## Betweenness centrality

conductance <- loadRData(paste0(resultsFolder, "/",conductanceFileRData))
write.table(conductance[,c("from","to","cost")], 
            file = paste0(resultsFolder,"/conductance.txt"), 
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

nk <- import("networkit")
graphio <- import("networkit.graphio")

# Set your thread count
nk$setNumberOfThreads(nCores)

# Read Graph and get Node Map
reader <- graphio$EdgeListReader(
  separator = " ",
  firstNode = 0L,
  continuous = FALSE,
  directed = TRUE
)
G <- reader$read(paste0(resultsFolder, "/",conductanceFile))
node_map <- reader$getNodeMap()

main <- import_main()
main$G <- G
main$node_map <- node_map

py_script <- paste0("
import networkit as nk
nk.setNumberOfThreads(",nCores,")
eb = nk.centrality.EstimateBetweenness(
    G,
    nSamples=",nSamples,",
    normalized=True,
    parallel=True
)
eb.run()

ranking = eb.ranking()
inverted_map = {v: k for k, v in node_map.items()}
output_filename = 'final_scores.csv'
with open(output_filename, 'w') as f:
    f.write('node,score\\n')
    for internal_id, score in ranking:
        original_id = inverted_map.get(internal_id, 'ID_NOT_FOUND')
        f.write(f'{original_id},{score}\\n')
")

# Execute the Python script block
py_run_string(py_script)

scores_df <- read.csv("final_scores.csv")
if(sum(is.na(match(oceanCells,scores_df$node))) > 0) {
  scores_df <- rbind(scores_df,data.frame(node=oceanCells[(is.na(match(oceanCells,scores_df$node)))], score=0))
}

centrality <- data.frame(cell=oceanCells,centrality=scores_df$score[match(oceanCells,scores_df$node)])
head(centrality)
save(centrality,file=paste0(resultsFolder,"/betweennessCentrality.RData"))

file.remove("final_scores.csv")
file.remove(paste0(resultsFolder,"/conductance.txt"))

## -----------------------------------------------------------------
## -----------------------------------------------------------------
## Eigenvector, Harmonic and Out-strenght Centrality

pairwiseConnectivity <- loadRData(paste0(resultsFolder, "/",conductanceFileRData))
g_cost <- graph_from_data_frame(pairwiseConnectivity[,c("from","to","cost")], directed = TRUE)
g_conductance <- graph_from_data_frame(pairwiseConnectivity[,c("from","to","conductance")], directed = TRUE)
epsilon <- .Machine$double.eps

# Harmonic Centrality
# Computes the sum of the reciprocal of distances from a node to all other nodes.
# Indicates how close a node is to all others, which can relate to connectivity importance.

nodesList <- split(V(g_cost), cut(seq_along(V(g_cost)), nCores*100, labels = FALSE))
Cluster <- makeCluster( nCores )
registerDoParallel( Cluster )
centrality <- foreach(nodesChunk = 1:length(nodesList), .combine = c, .packages = 'igraph') %dopar% {
  harmonic_centrality(g_cost, vids = nodesList[[nodesChunk]], weights = E(g_cost)$cost + epsilon, cutoff = 100, normalized=TRUE)
}
stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
closeAllConnections()

centrality <- unlist(centrality)
centrality <- data.frame(cell=oceanCells,centrality=centrality[match(oceanCells,as.numeric(names(centrality)))])
save(centrality,file=paste0(resultsFolder,"/harmonicCentrality.RData"))

# Closeness

Cluster <- makeCluster( nCores )
registerDoParallel( Cluster )
centrality <- foreach(nodesChunk = 1:length(nodesList), .combine = c, .packages = 'igraph') %dopar% {
  closeness(g_cost, vids = nodesList[[nodesChunk]], weights = E(g_cost)$cost, mode = "all", cutoff = 100, normalized=TRUE)
}
stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
closeAllConnections()

centrality <- unlist(centrality)
centrality <- data.frame(cell=oceanCells,centrality=centrality[match(oceanCells,as.numeric(names(centrality)))])
save(centrality,file=paste0(resultsFolder,"/closenessCentrality.RData"))

# Out-Strength Centrality
# Sum of weights of edges connected to a node.

# Cluster <- makeCluster( nCores )
# registerDoParallel( Cluster )
# centrality <- foreach(nodesChunk = nodesList, .combine = c, .packages = 'igraph') %dopar% {
#   strength(g_conductance, vids = unlist(nodesChunk), weights = E(g_conductance)$conductance, mode = "out")
# }
# stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
# closeAllConnections()

centrality <- strength(g_conductance, weights = E(g_conductance)$conductance, mode = "out")
centrality <- data.frame(cell=oceanCells,centrality=centrality[match(oceanCells,as.numeric(names(centrality)))])
save(centrality,file=paste0(resultsFolder,"/strengthOutCentrality.RData"))

# Eigenvector
# Measures the influence of a node based on the importance of its neighbors.
# Nodes connected to high-scoring nodes receive higher scores themselves.
# Interpretation: Focuses on nodes connected to other influential nodes, which may differ from betweenness centrality.
# Identify nodes that are important because they are connected to other important nodes.

centrality <- eigen_centrality(g_conductance, directed = graphDirected , weights = E(g_conductance)$conductance )$vector
centrality <- data.frame(cell=oceanCells,centrality=centrality[match(oceanCells,as.numeric(names(centrality)))])
save(centrality,file=paste0(resultsFolder,"/eigenvectorCentrality.RData"))
