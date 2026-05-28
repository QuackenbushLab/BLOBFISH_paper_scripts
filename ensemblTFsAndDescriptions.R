# Convert Symbols to Ensembl IDs.
library("org.Hs.eg.db")
library("AnnotationDbi")

# Set paths.
blobfishDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/GRAND_BLOBFISH_results/"
motifFile <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Roaa/Null_PANDA_prior/tissues_motif.txt"

# Load BLOBFISH networks
lungBlobfish <- read.csv(paste0(blobfishDir, "lungSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)
adiposeBlobfish <- read.csv(paste0(blobfishDir, "adiposeSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)
aortaBlobfish <- read.csv(paste0(blobfishDir, "aortaSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)
muscleBlobfish <- read.csv(paste0(blobfishDir, "skeletalMuscleSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)
skinBlobfish <- read.csv(paste0(blobfishDir, "skinSubnet_05_updatedNull_moreGenes.csv"), row.names = 1)

# Load BLOBFISH exclusive networks.
lungBlobfishX <- read.csv(paste0(blobfishDir, "lung_only_updated.csv"), row.names = 1)
adiposeBlobfishX <- read.csv(paste0(blobfishDir, "subcutaneous_adipose_only_updated.csv"), row.names = 1)
aortaBlobfishX <- read.csv(paste0(blobfishDir, "aorta_only_updated.csv"), row.names = 1)
muscleBlobfishX <- read.csv(paste0(blobfishDir, "skeletal_muscle_only_updated.csv"), row.names = 1)
skinBlobfishX <- read.csv(paste0(blobfishDir, "skin_only_updated.csv"), row.names = 1)

# Load motif prior.
motif <- read.table(motifFile)
motif <- motif[which(motif$V3 == 1),]

# Get Ensembl IDs for all TFs, updating old names.
tfsToInput <- uniqueTFs <- unique(motif$V1)
uniqueTFsModMapping <- data.frame(symbolsUsed = c("ARNTL", "PRKRIR", "T"), commonSymbols = c("BMAL1", "THAP12", "TBXT"))
tfsToInput[which(tfsToInput %in% uniqueTFsModMapping$symbolsUsed)] <- 
  unlist(lapply(tfsToInput[which(tfsToInput %in% uniqueTFsModMapping$symbolsUsed)], function(oldSym){
  return(uniqueTFsModMapping[which(uniqueTFsModMapping$symbolsUsed == oldSym), "commonSymbols"])
}))
tfs_ensembl <- mapIds(org.Hs.eg.db, keys = tfsToInput, keytype = "SYMBOL", column="ENSEMBL")
tfs_descname <- mapIds(org.Hs.eg.db, keys = tfsToInput, keytype = "SYMBOL", column="GENENAME")
names(tfs_ensembl)[which(names(tfs_ensembl) %in% uniqueTFsModMapping$commonSymbols)] <- 
  unlist(lapply(names(tfs_ensembl)[which(names(tfs_ensembl) %in% uniqueTFsModMapping$commonSymbols)], function(newSym){
    return(uniqueTFsModMapping[which(uniqueTFsModMapping$commonSymbols == newSym), "symbolsUsed"])
  }))
names(tfs_descname)[which(names(tfs_descname) %in% uniqueTFsModMapping$commonSymbols)] <- 
  unlist(lapply(names(tfs_descname)[which(names(tfs_descname) %in% uniqueTFsModMapping$commonSymbols)], function(newSym){
    return(uniqueTFsModMapping[which(uniqueTFsModMapping$commonSymbols == newSym), "symbolsUsed"])
  }))

# For all BLOBFISH networks, match data.
addTfInfo <- function(network){
  network$tfEnsembl <- tfs_ensembl[network$tf]
  network$tfDesc <- tfs_descname[network$tf]
  str(network)
  return(network)
}
aortaTf <- addTfInfo(aortaBlobfishX)
adiposeTf <- addTfInfo(adiposeBlobfishX)
lungTf <- addTfInfo(lungBlobfishX)
muscleTf <- addTfInfo(muscleBlobfishX)
skinTf <- addTfInfo(skinBlobfishX)
aortaTfX <- addTfInfo(aortaBlobfishX)
adiposeTfX <- addTfInfo(adiposeBlobfishX)
lungTfX <- addTfInfo(lungBlobfishX)
muscleTfX <- addTfInfo(muscleBlobfishX)
skinTfX <- addTfInfo(skinBlobfishX)

# Get descriptions for all genes.
getDesc <- function(network, isEnsembl){
  uniqueGenes <- unique(network$gene)
  gene_descname <- NULL
  if(isEnsembl == FALSE){
    gene_descname <- mapIds(org.Hs.eg.db, keys = uniqueGenes, keytype = "SYMBOL", column="GENENAME")
  }else{
    gene_descname <- mapIds(org.Hs.eg.db, keys = uniqueGenes, keytype = "ENSEMBL", column="GENENAME")
  }
  network$geneDesc <- gene_descname[network$gene]
  str(network)
  return(network)
}
adiposeAll <- getDesc(adiposeTf, TRUE)
aortaAll <- getDesc(aortaTf, TRUE)
lungAll <- getDesc(lungTf, TRUE)
muscleAll <- getDesc(muscleTf, TRUE)
skinAll <- getDesc(skinTf, TRUE)
adiposeAllX <- getDesc(adiposeTfX, FALSE)
aortaAllX <- getDesc(aortaTfX, FALSE)
lungAllX <- getDesc(lungTfX, FALSE)
muscleAllX <- getDesc(muscleTfX, FALSE)
skinAllX <- getDesc(skinTfX, FALSE)

# Save results.
write.csv(adiposeAll, paste0(blobfishDir, "adiposeComplete.csv"))
write.csv(aortaAll, paste0(blobfishDir, "aortaComplete.csv"))
write.csv(lungAll, paste0(blobfishDir, "lungComplete.csv"))
write.csv(muscleAll, paste0(blobfishDir, "muscleComplete.csv"))
write.csv(skinAll, paste0(blobfishDir, "skinComplete.csv"))
write.csv(adiposeAllX, paste0(blobfishDir, "adiposeCompleteX.csv"))
write.csv(aortaAllX, paste0(blobfishDir, "aortaCompleteX.csv"))
write.csv(lungAllX, paste0(blobfishDir, "lungCompleteX.csv"))
write.csv(muscleAllX, paste0(blobfishDir, "muscleCompleteX.csv"))
write.csv(skinAllX, paste0(blobfishDir, "skinCompleteX.csv"))