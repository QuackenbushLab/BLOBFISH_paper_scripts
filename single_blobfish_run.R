# Load the network zoo.
library("netZooR")

# Read in the files.
args <- commandArgs(trailingOnly = TRUE)
cellType <- args[1]
randsetIndex <- as.numeric(args[2])
pvalfile <- args[3]
outdir <- NULL
null_file <- NULL
tissue_dir <- NULL
randset_file <- NULL
tissue <- readRDS(paste0(tissue_dir, cellType, ".RDS"))
null <- readRDS(null_file)
randsets <- readRDS(randset_file)

# Run BLOBFISH on each randomization, for each tissue.
str(cellType)
str(pvalfile)
str(randsets[[randsetIndex]])
blobfish <- netZooR::RunBLOBFISH(geneSet = randsets[[randsetIndex]],
                networks = tissue, alpha = 0.05, hopConstraint = 2,
                nullDistribution = null, pValueFile = pvalfile,
                verbose = TRUE, loadPValues = TRUE,
                doFDRAdjustment = TRUE)
saveRDS(blobfish, paste0(outdir, "/", cellType, "_", randsetIndex, ".RDS"))
