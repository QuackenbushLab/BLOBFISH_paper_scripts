# Load the network zoo.
library("netZooR")

# Add PPI file and motif file paths and the directory where you wish to save the
# null model.
ppiFilePath <- NULL
motifFilePath <- NULL
nullFilePath <- NULL

# Generate and save null model.
null <- GenerateNullPANDADistribution(ppiFile = ppiFilePath, motifFile = motifFilePath)
saveRDS(null, nullFilePath)
