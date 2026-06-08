# Load igraph to perform graph-based analysis.
library(igraph)

# Convert Symbols to Ensembl IDs.
library("org.Hs.eg.db")

# Read in the filtered networks.
indir <- NULL
outdir <- NULL
dir.create(outdir)
adipose <- read.table(paste0(indir, "/subcutaneous_adipose.sig"))[,c("V1", "V6")]
aorta <- read.table(paste0(indir, "/aorta.sig"))[,c("V1", "V6")]
lung <- read.table(paste0(indir, "/lung.sig"))[,c("V1", "V6")]
muscle <- read.table(paste0(indir, "/skeletal_muscle.sig"))[,c("V1", "V6")]
skin <- read.table(paste0(indir, "/skin.sig"))[,c("V1", "V6")]

# Select genes from papers to focus on. PMC311426 for muscle.
# MB = myoglobin, MYH2 = myosin heavy chain 2a, MYL2 = myosin light chain 2,
# DES = desmin, TNNC2 = fast skeletal troponin C, TNNC1 = EST (slow twitch sk troponin 1),
# ENO3 = β enolase, MYL3 = fast myosin alkali light chain 3, TTN = titin, 
# TPM1 = α-tropomyosin, TCAP = telethonin, MYL1 = myosin alkali light chain 1f
# TPM3,4 = tropomyosinm, COX5B = cytochrome c oxidase 5b, COX5A = cytochrome c oxidase 5a
# MYL2 =  EST (myosin regulatory light chain 2), CKMT2 = creatine kinase, sarcomeric mitochondrial,
# α tubulin = TUBA1A, TUBA1B, TUBA4A, TUBA1C, TUBA3C, TUBA8, TUBA3D, TUBA3E, TUBA4B
muscle_dev_genes <- c("MB", "MYH2", "MYL2", "DES", "TNNC1", "TNNC2", "ENO3", "MYL3",
                      "TTN", "TPM1", "TCAP", "MYL1", "TPM3", "TPM4", "COX5B",
                      "COX5A", "CKMT2", "TUBA1A", "TUBA1B", "TUBA4A", "TUBA1C", 
                      "TUBA3C", "TUBA8", "TUBA3D", "TUBA3E", "TUBA4B")
muscle_dev_ensembl <- mapIds(org.Hs.eg.db, keys = muscle_dev_genes, keytype = "SYMBOL", column="ENSEMBL")
adipogenic_genes <- c("PPARG", "FASN", "SREBF1", "SCD", "CEBPA", "ADIPOQ", "FABP4")
adipogenic_ensembl <- mapIds(org.Hs.eg.db, keys = adipogenic_genes, keytype = "SYMBOL", column="ENSEMBL")
seeds_ensembl <- unique(c(muscle_dev_ensembl, adipogenic_ensembl))
seeds <- unique(c(muscle_dev_genes, adipogenic_genes))

# Run SP in the skeletal muscle.
getAllShortestPaths <- function(network, seeds){
  graph <- igraph::graph_from_edgelist(as.matrix(network), directed = FALSE)
  sps <- do.call(rbind, lapply(seeds, function(seed1){
    edgepath <- igraph::shortest_paths(graph, from = seed1, to=setdiff(seeds, seed1), output="epath")
    retval <- do.call(rbind, lapply(edgepath$epath, function(p) {
      retval2 <- NULL
      if(length(p) > 0){
        ends <- igraph::ends(graph, p)
        retval2 <- data.frame(source = ends[, 1], target = ends[, 2], stringsAsFactors = FALSE)
      }
      return(retval2)
    }))
    cat(".")
    return(retval)
  }))
  str(sps)
  return(sps)
}
subgraphSkin <- getAllShortestPaths(skin, seeds_ensembl)
subgraphAorta <- getAllShortestPaths(aorta, seeds_ensembl)
subgraphLung <- getAllShortestPaths(lung, seeds_ensembl)
subgraphMuscle <- getAllShortestPaths(muscle, seeds_ensembl)
subgraphAdipose <- getAllShortestPaths(adipose, seeds_ensembl)
write.csv(subgraphSkin, paste0(outdir, "skin.csv"))
write.csv(subgraphAorta, paste0(outdir, "aorta.csv"))
write.csv(subgraphLung, paste0(outdir, "lung.csv"))
write.csv(subgraphMuscle, paste0(outdir, "muscle.csv"))
write.csv(subgraphAdipose, paste0(outdir, "adipose.csv"))

# Adipogenic markers from 36647068
geneColorMapping <- data.frame(gene = c(muscle_dev_ensembl, adipogenic_ensembl), 
                               color = c(rep("hotpink", length(muscle_dev_genes)), 
                                           rep("goldenrod1", length(adipogenic_genes)))) 

# Get the exclusive networks.
getExclusive <- function(network, otherNetworks, file){
  netNames <- paste(network[,1], network[,2], sep = "__")
  otherNetNames <- unlist(lapply(otherNetworks, function(net){
    return(paste(net[,1], net[,2], sep = "__"))
  }))
  str(netNames)
  str(otherNetNames)
  exclusiveNames <- setdiff(netNames, otherNetNames)
  str(exclusiveNames)
  whichExclusiveNames <- unlist(lapply(exclusiveNames, function(name){
    return(which(netNames == name)[1])
  }))
  exclusiveNet <- network[whichExclusiveNames,]
  colnames(exclusiveNet) <- c("tf", "gene")
  str(exclusiveNet)
  write.csv(exclusiveNet, file)
  return(exclusiveNet)
}
adiposeOnly <- getExclusive(subgraphAdipose, list(subgraphSkin, subgraphAorta, subgraphLung, subgraphMuscle),
                            paste0(outdir, "adiposeOnly.csv"))
muscleOnly <- getExclusive(subgraphMuscle, list(subgraphSkin, subgraphAorta, subgraphLung, subgraphAdipose),
                           paste0(outdir, "muscleOnly.csv"))
skinOnly <- getExclusive(subgraphSkin, list(subgraphMuscle, subgraphAorta, subgraphLung, subgraphAdipose),
                         paste0(outdir, "skinOnly.csv"))
lungOnly <- getExclusive(subgraphLung, list(subgraphMuscle, subgraphAorta, subgraphSkin, subgraphAdipose),
                         paste0(outdir, "lungOnly.csv"))
aortaOnly <- getExclusive(subgraphAorta, list(subgraphMuscle, subgraphLung, subgraphSkin, subgraphAdipose),
                          paste0(outdir, "aortaOnly.csv"))

# Plot
pdf(paste0(outdir, "/skeletal_muscle_only.pdf"))
PlotNetwork(network = muscleOnly,genesOfInterest = seeds_ensembl,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

pdf(paste0(outdir, "/adipose_only.pdf"))
PlotNetwork(network = adiposeOnly,genesOfInterest = seeds_ensembl,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

pdf(paste0(outdir, "/skin_only.pdf"))
PlotNetwork(network = skinOnly,genesOfInterest = seeds_ensembl,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

pdf(paste0(outdir, "/lung_only.pdf"))
PlotNetwork(network = lungOnly,genesOfInterest = seeds_ensembl,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

pdf(paste0(outdir, "/aorta_only.pdf"))
PlotNetwork(network = aortaOnly,genesOfInterest = seeds_ensembl,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

# Check the average number of TFs connecting a pair of genes.
TFCountAllPairs <- function(genes, network){
  return(unlist(lapply(genes, function(i){
    return(unlist(lapply(setdiff(genes, i), function(j){
      tfsConnectedToi <- unique(network[which(network$gene == i), "tf"])
      tfsConnectedToj <- unique(network[which(network$gene == j), "tf"])
      return(length(intersect(tfsConnectedToi, tfsConnectedToj)))
    })))
  })))
}
tfCountAllPairsMuscleMuscle <- TFCountAllPairs(muscle_dev_ensembl, muscleOnly)
print(mean(tfCountAllPairsMuscleMuscle))
tfCountAllPairsMuscleAdipose <- TFCountAllPairs(adipogenic_ensembl, muscleOnly)
print(mean(tfCountAllPairsMuscleAdipose))
tfCountAllPairsAdiposeMuscle <- TFCountAllPairs(muscle_dev_ensembl, adiposeOnly)
print(mean(tfCountAllPairsAdiposeMuscle))
tfCountAllPairsAdiposeAdipose <- TFCountAllPairs(adipogenic_ensembl, adiposeOnly)
print(mean(tfCountAllPairsAdiposeAdipose))
tfCountAllPairsSkinMuscle <- TFCountAllPairs(muscle_dev_ensembl, skinOnly)
print(mean(tfCountAllPairsSkinMuscle))
tfCountAllPairsSkinAdipose <- TFCountAllPairs(adipogenic_ensembl, skinOnly)
print(mean(tfCountAllPairsSkinAdipose))
tfCountAllPairsLungMuscle <- TFCountAllPairs(muscle_dev_ensembl, lungOnly)
print(mean(tfCountAllPairsLungMuscle))
tfCountAllPairsLungAdipose <- TFCountAllPairs(adipogenic_ensembl, lungOnly)
print(mean(tfCountAllPairsLungAdipose))
tfCountAllPairsAortaMuscle <- TFCountAllPairs(muscle_dev_ensembl, aortaOnly)
print(mean(tfCountAllPairsAortaMuscle))
tfCountAllPairsAortaAdipose <- TFCountAllPairs(adipogenic_ensembl, aortaOnly)
print(mean(tfCountAllPairsAortaAdipose))