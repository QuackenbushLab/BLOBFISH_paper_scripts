# Load the network zoo.
library("netZooR")

# Read in the files.
dir_input <- NULL
output_file <- NULL
skin <- readRDS(paste0(outdir, "skin.RDS"))
subcutaneous_adipose <- readRDS(paste0(dir_input, "adipose.RDS"))
skeletal_muscle <- readRDS(paste0(dir_input, "muscle.RDS"))
lung <- readRDS(paste0(dir_input, "lung.RDS"))
print("read all files")

# Generate random sets.
num_rand_sets <- 100
num_genes <- 33
sharedGenes <- unique(skeletal_muscle[[1]]$gene)
for(tissue in list(skeletal_muscle, subcutaneous_adipose, skin, lung, aorta)){
  for(i in 1:length(tissue)){
    sharedGenes <- intersect(sharedGenes, unique(tissue[[i]]$gene))
  }
}
str(sharedGenes)
randsets <- lapply(1:num_rand_sets, function(i){
  return(sample(sharedGenes, num_genes))
})
str(randsets)
saveRDS(randsets, output_file)
