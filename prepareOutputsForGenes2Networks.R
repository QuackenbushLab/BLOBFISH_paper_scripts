# Load the package to flatten the matrix.
library(reshape2)

# Load the network zoo (BLOBFISH package).
library(netZooR)

# Load igraph to perform graph-based analysis.
library(igraph)

# Convert Symbols to Ensembl IDs.
library("org.Hs.eg.db")

# Load the package to flatten the matrix.
library(reshape2)

# Read in the p-value files. We will use these to filter the networks for input to Genes2Networks.
indir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/GRAND_BLOBFISH_results/"
outdir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/BLOBFISH_results/Genes2Networks_results/"
dir.create(outdir)
skin <- readRDS(paste0(indir, "/skinPval.RDS"))
subcutaneous_adipose <- readRDS(paste0(indir, "/adiposePval.RDS"))
skeletal_muscle <- readRDS(paste0(indir, "/skeletalMusclePval.RDS"))
lung <- readRDS(paste0(indir, "/lungPval.RDS"))
aorta <- readRDS(paste0(indir, "/aortaPval.RDS"))

# Save the sig files.
writeSigFile <- function(tissue, sigFile){
 
  # Build edges.
  sigEdges <- tissue[which(tissue < 0.05)] 
  source <- unlist(lapply(names(sigEdges), function(edge){return(strsplit(edge, "__")[[1]][1])}))
  target <- unlist(lapply(names(sigEdges), function(edge){return(strsplit(edge, "__")[[1]][2])}))
  edges <- data.frame(source = source, target = target, weight = rep(1, length(source)))
  edges <- data.frame(
    source = source,
    c2 = "NA",
    c3 = "NA",
    source_type = "protein",
    source_location = "NA",
    target = target,
    c7 = "NA",
    c8 = "NA",
    target_type = "protein",
    target_location = "NA",
    interaction = "interaction",
    database = "skin",
    pmid = paste0("SKIN_", seq_len(length(source)))
  )
  str(edges)
  
  # Write it.
  write.table(edges, file = sigFile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
writeSigFile(skin, paste0(outdir, "/skin.sig"))
writeSigFile(subcutaneous_adipose, paste0(outdir, "/subcutaneous_adipose.sig"))
writeSigFile(skeletal_muscle, paste0(outdir, "/skeletal_muscle.sig"))
writeSigFile(lung, paste0(outdir, "/lung.sig"))
writeSigFile(aorta, paste0(outdir, "/aorta.sig"))

# Select genes from papers to focus on. PMC311426 for muscle.
# MB = myoglobin, MYH2 = myosin heavy chain 2a, MYL2 = myosin light chain 2,
# DES = desmin, TNNC2 = fast skeletal troponin C, TNNC1 = EST (slow twitch sk troponin 1),
# ENO3 = β enolase, MYL3 = fast myosin alkali light chain 3, TTN = titin, 
# TPM1 = α-tropomyosin, TCAP = telethonin, MYL1 = myosin alkali light chain 1f
# TPM3,4 = tropomyosinm, COX5B = cytochrome c oxidase 5b, COX5A = cytochrome c oxidase 5a
# MYL2 =  EST (myosin regulatory light chain 2), CKMT2 = creatine kinase, sarcomeric mitochondrial,
# α tubulin = TUBA1A, TUBA1B, TUBA4A, TUBA1C, TUBA3C, TUBA8, TUBA3D, TUBA3E, TUBA4B
muscle_dev_genes <- c("MB", "MYH2", "MYL2", "DES", "TNNC1", "TNNC2", "ENO3", "MYL3",
                      "TTN", "TPM1", "TCAP", "MYL1", "TPM3", "TPM4", "COX5B", "MYL2",
                      "COX5A", "CKMT2", "TUBA1A", "TUBA1B", "TUBA4A", "TUBA1C", 
                      "TUBA3C", "TUBA8", "TUBA3D", "TUBA3E", "TUBA4B")
muscle_dev_ensembl <- mapIds(org.Hs.eg.db, keys = muscle_dev_genes, keytype = "SYMBOL", column="ENSEMBL")
adipogenic_genes <- c("PPARG", "FASN", "SREBF1", "SCD", "CEBPA", "ADIPOQ", "FABP4")
adipogenic_ensembl <- mapIds(org.Hs.eg.db, keys = adipogenic_genes, keytype = "SYMBOL", column="ENSEMBL")
genes <- unique(c(unname(muscle_dev_ensembl), unname(adipogenic_ensembl)))

# Write the genes.
writeLines(toupper(genes), paste0(outdir, "/inputGenesGenes2Network.txt"))