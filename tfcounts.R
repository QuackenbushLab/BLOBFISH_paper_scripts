# Read in the files.
dir_randoutputs <- NULL
randset_file <- NULL
randsets <- readRDS(randset_file)

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

# Get the TF counts.
tissues <- c("adipose", "aorta", "lung", "muscle", "skin")
tissueSpecificMeans <- list()
tissueSpecificMeans[["adipose"]] <- list()
tissueSpecificMeans[["aorta"]] <- list()
tissueSpecificMeans[["lung"]] <- list()
tissueSpecificMeans[["muscle"]] <- list()
tissueSpecificMeans[["skin"]] <- list()

for(randsetIndex in 1:100){
   blobfishes <- list()
   blobfishes[["adipose"]] <- readRDS(paste0(dir_randoutputs, "/adipose_", randsetIndex, ".RDS"))
   blobfishes[["aorta"]] <- readRDS(paste0(dir_randoutputs, "/aorta_", randsetIndex, ".RDS"))
   blobfishes[["lung"]] <- readRDS(paste0(dir_randoutputs, "/lung_", randsetIndex, ".RDS"))
   blobfishes[["muscle"]] <- readRDS(paste0(dir_randoutputs, "/muscle_", randsetIndex, ".RDS"))
   blobfishes[["skin"]] <- readRDS(paste0(dir_randoutputs, "/skin_", randsetIndex, ".RDS"))
   
   for(tissue in tissues){
      allOtherBlobfishes <- unlist(lapply(setdiff(tissues, tissue), function(tis){return(rownames(blobfishes[[tis]]))}))
      blobfish <- blobfishes[[tissue]][setdiff(rownames(blobfishes[[tissue]]), allOtherBlobfishes),]
      tfCounts <- TFCountAllPairs(genes = randsets[[randsetIndex]], network = blobfish)
      tissueSpecificMeans[[tissue]][[randsetIndex]] <- mean(tfCounts)
   }
}

for(tissue in tissues){
   print(tissue)
   print(mean(unlist(tissueSpecificMeans[[tissue]])))
   str(unlist(tissueSpecificMeans[[tissue]]))
   print(sd(unlist(tissueSpecificMeans[[tissue]])))
}
