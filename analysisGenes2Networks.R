library(igraph)

# Read in the results from Genes2Networks.
indir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/BLOBFISH_results/Genes2Networks_results/"
outdir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/BLOBFISH_results/Genes2Networks_results/"
skinSubnet <- read.table(paste0(outdir, "/skinOut.sig"))[,c("V1", "V6")]
skeletalMuscleSubnet <- read.table(paste0(outdir, "/muscleOut.sig"))[,c("V1", "V6")]
lungSubnet <- read.table(paste0(outdir, "/lungOut.sig"))[,c("V1", "V6")] 
aortaSubnet <- read.table(paste0(outdir, "/aortaOut.sig"))[,c("V1", "V6")]
adiposeSubnet <- read.table(paste0(outdir, "/adiposeOut.sig"))[,c("V1", "V6")]
colnames(skinSubnet) <- colnames(skeletalMuscleSubnet) <- colnames(aortaSubnet) <- colnames(adiposeSubnet) <- colnames(lungSubnet) <- c("tf", "gene")

# Genes
muscle_dev_genes <- c("MB", "MYH2", "MYL2", "DES", "TNNC1", "TNNC2", "ENO3", "MYL3",
                      "TTN", "TPM1", "TCAP", "MYL1", "TPM3", "TPM4", "COX5B",
                      "COX5A", "CKMT2", "TUBA1A", "TUBA1B", "TUBA4A", "TUBA1C", 
                      "TUBA3C", "TUBA8", "TUBA3D", "TUBA3E", "TUBA4B")
muscle_dev_ensembl <- mapIds(org.Hs.eg.db, keys = muscle_dev_genes, keytype = "SYMBOL", column="ENSEMBL")
adipogenic_genes <- c("PPARG", "FASN", "SREBF1", "SCD", "CEBPA", "ADIPOQ", "FABP4")
adipogenic_ensembl <- mapIds(org.Hs.eg.db, keys = adipogenic_genes, keytype = "SYMBOL", column="ENSEMBL")

# Adipogenic markers from 36647068
geneColorMapping <- data.frame(gene = c(muscle_dev_genes, adipogenic_genes), 
                               color = c(rep("hotpink", length(muscle_dev_genes)), 
                                           rep("goldenrod1", length(adipogenic_genes)))) 

# Convert genes to symbol.
convertToSymbol <- function(subnet){
  subnet$gene <- unlist(lapply(subnet$gene, function(gene){
    retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
    if(length(retval) == 0){
      retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    }
    return(retval[1])
  }))
  return(subnet)
}
adiposeSymbol <- convertToSymbol(adiposeSubnet)
muscleSymbol <- convertToSymbol(skeletalMuscleSubnet)
lungSymbol <- convertToSymbol(lungSubnet)
aortaSymbol <- convertToSymbol(aortaSubnet)
skinSymbol <- convertToSymbol(skinSubnet)

# Remove duplicate edges. Keep only the first (order doesn't matter here because
# we don't have weights).
removeDups <- function(network){
  netNames <- paste(network[,1], network[,2], sep = "__")
  str(netNames)
  netNameCounts <- table(netNames)
  nameWhichDup <- names(netNameCounts)[which(netNameCounts > 1)]
  print(nameWhichDup)
  toRemove <- unlist(lapply(nameWhichDup, function(name){
    whichName <- which(netNames == name)
    whichNotFirst <- whichName[2:lenght(whichName)]
    return(whichNotFirst)
  }))
  newNet <- network[setdiff(1:nrow(network), toRemove),]
}
adiposeSubnetNew <- removeDups(adiposeSymbol)
muscleSubnetNew <- removeDups(muscleSymbol)
lungSubnetNew <- removeDups(lungSymbol)
aortaSubnetNew <- removeDups(aortaSymbol)
skinSubnetNew <- removeDups(skinSymbol)

# Get exclusive networks.
getExclusive <- function(network, otherNetworks){
  netNames <- paste(network[,1], network[,2], sep = "__")
  otherNetNames <- unlist(lapply(otherNetworks, function(net){
    return(paste(net[,1], net[,2], sep = "__"))
  }))
  str(netNames)
  str(otherNetNames)
  exclusiveNames <- setdiff(netNames, otherNetNames)
  str(exclusiveNames)
  exclusiveNet <- network[which(netNames %in% exclusiveNames),]
  return(exclusiveNet)
}
adiposeOnly <- getExclusive(adiposeSubnetNew, list(muscleSubnetNew, aortaSubnetNew, skinSubnetNew, lungSubnetNew))
muscleOnly <- getExclusive(muscleSubnetNew, list(adiposeSubnetNew, aortaSubnetNew, skinSubnetNew, lungSubnetNew))
skinOnly <- getExclusive(skinSubnetNew, list(adiposeSubnetNew, aortaSubnetNew, muscleSubnetNew, lungSubnetNew))
lungOnly <- getExclusive(lungSubnetNew, list(adiposeSubnetNew, skinSubnetNew, aortaSubnetNew, muscleSubnetNew))
aortaOnly <- getExclusive(aortaSubnetNew, list(adiposeSubnetNew, skinSubnetNew, muscleSubnetNew, lungSubnetNew))

# Plot exclusive networks.
pdf(paste0(outdir, "/adipose_Genes2Networks.pdf"))
PlotNetwork(network = adiposeOnly,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
dev.off()

pdf(paste0(outdir, "/muscle_Genes2Networks.pdf"))
PlotNetwork(network = muscleOnly,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
dev.off()

pdf(paste0(outdir, "/skin_Genes2Networks.pdf"))
PlotNetwork(network = skinOnly,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
dev.off()

pdf(paste0(outdir, "/aorta_Genes2Networks.pdf"))
PlotNetwork(network = aortaOnly,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
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
tfCountAllPairsMuscleMuscle <- TFCountAllPairs(muscle_dev_genes, muscleOnly)
print(mean(tfCountAllPairsMuscleMuscle))
tfCountAllPairsMuscleAdipose <- TFCountAllPairs(adipogenic_genes, muscleOnly)
print(mean(tfCountAllPairsMuscleAdipose))
tfCountAllPairsAdiposeMuscle <- TFCountAllPairs(muscle_dev_genes, adiposeOnly)
print(mean(tfCountAllPairsAdiposeMuscle))
tfCountAllPairsAdiposeAdipose <- TFCountAllPairs(adipogenic_genes, adiposeOnly)
print(mean(tfCountAllPairsAdiposeAdipose))
tfCountAllPairsSkinMuscle <- TFCountAllPairs(muscle_dev_genes, skinOnly)
print(mean(tfCountAllPairsSkinMuscle))
tfCountAllPairsSkinAdipose <- TFCountAllPairs(adipogenic_genes, skinOnly)
print(mean(tfCountAllPairsSkinAdipose))
tfCountAllPairsLungMuscle <- TFCountAllPairs(muscle_dev_genes, lungOnly)
print(mean(tfCountAllPairsLungMuscle))
tfCountAllPairsLungAdipose <- TFCountAllPairs(adipogenic_genes, lungOnly)
print(mean(tfCountAllPairsLungAdipose))
tfCountAllPairsAortaMuscle <- TFCountAllPairs(muscle_dev_genes, aortaOnly)
print(mean(tfCountAllPairsAortaMuscle))
tfCountAllPairsAortaAdipose <- TFCountAllPairs(adipogenic_genes, aortaOnly)
print(mean(tfCountAllPairsAortaAdipose))

# Compare to BLOBFISH results.
blobfishDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/GRAND_BLOBFISH_results/"
lungBlobfish <- read.csv(paste0(blobfishDir, "lungSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)
adiposeBlobfish <- read.csv(paste0(blobfishDir, "adiposeSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)
aortaBlobfish <- read.csv(paste0(blobfishDir, "aortaSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)
muscleBlobfish <- read.csv(paste0(blobfishDir, "skeletalMuscleSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)
skinBlobfish <- read.csv(paste0(blobfishDir, "skinSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)

compareNets <- function(net1, net2){
  netNames1 <- paste(net1[,1], net1[,2], sep = "__")
  netNames2 <- paste(net2[,1], net2[,2], sep = "__")
  str(setdiff(netNames1, netNames2))
  str(setdiff(netNames2, netNames1))
}
compareNets(lungBlobfish, lungSubnet)
compareNets(adiposeBlobfish, adiposeSubnet)
compareNets(aortaBlobfish, aortaSubnet)
compareNets(muscleBlobfish, skeletalMuscleSubnet)
compareNets(skinBlobfish, skinSubnet)