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

# Load the package to run pathway analysis.
library(fgsea)

# Modify the directories as needed. These include the directories containing
# networks from each tissue type (downloaded from GRAND), expression for
# each tissue type (downloaded from GRAND), and the output directory.
dir_subcutaneous_adipose <- NULL
dir_expression <- NULL
dir_skeletal_muscle <- NULL
dir_skin <- NULL
dir_lung <- NULL
dir_aorta <- NULL
outdir <- NULL
gmt_pathway_file <- NULL
null_file <- NULL
null_output_file <- NULL

dir.create(outdir)
skin <- lapply(list.files(dir_skin), function(file){
 adj <- read.csv(paste0(dir_skin, "/", file))
 flatadj <- reshape2::melt(adj)
 colnames(flatadj) <- c("tf", "gene", "score")
 flatadj$gene <- as.character(flatadj$gene)
 str(flatadj)
 return(flatadj)
})
saveRDS(skin, paste0(outdir, "/skin.RDS"))
subcutaneous_adipose <- lapply(list.files(dir_subcutaneous_adipose), function(file){
  adj <- read.csv(paste0(dir_subcutaneous_adipose, "/", file))
  flatadj <- reshape2::melt(adj)
  colnames(flatadj) <- c("tf", "gene", "score")
  flatadj$gene <- as.character(flatadj$gene)
  str(flatadj)
  return(flatadj)
})
saveRDS(subcutaneous_adipose, paste0(outdir, "/adipose.RDS"))
skeletal_muscle <- lapply(list.files(dir_skeletal_muscle), function(file){
  adj <- read.csv(paste0(dir_skeletal_muscle, "/", file))
  flatadj <- reshape2::melt(adj)
  colnames(flatadj) <- c("tf", "gene", "score")
  flatadj$gene <- as.character(flatadj$gene)
  str(flatadj)
  return(flatadj)
})
saveRDS(skeletal_muscle, paste0(outdir, "/muscle.RDS"))
lung <- lapply(list.files(dir_lung), function(file){
  adj <- read.csv(paste0(dir_lung, "/", file))
  flatadj <- reshape2::melt(adj)
  colnames(flatadj) <- c("tf", "gene", "score")
  flatadj$gene <- as.character(flatadj$gene)
  str(flatadj)
  return(flatadj)
})
saveRDS(lung, paste0(outdir, "/lung.RDS"))
aorta <- lapply(list.files(dir_aorta), function(file){
  adj <- read.csv(paste0(dir_aorta, "/", file))
  flatadj <- reshape2::melt(adj)
  colnames(flatadj) <- c("tf", "gene", "score")
  flatadj$gene <- as.character(flatadj$gene)
  str(flatadj)
  return(flatadj)
})
saveRDS(aorta, paste0(outdir, "/aorta.RDS"))

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

# Load null distribution.
pandas<-readRDS(null_file)
skin <- readRDS(paste0(outdir, "/skin.RDS"))
skeletal_muscle <- readRDS(paste0(outdir, "/muscle.RDS"))
subcutaneous_adipose <- readRDS(paste0(outdir, "/adipose.RDS"))
lung <- readRDS(paste0(outdir, "/lung.RDS"))
aorta <- readRDS(paste0(outdir, "/aorta.RDS"))
str("read files")

# Run BLOBFISH in the skeletal muscle.
null <- sample(pandas, 1000)
saveRDS(null, null_output_file)
skeletalMuscleSubnet <- RunBLOBFISH(geneSet = unname(adipogenicmuscle_dev_ensembl),
                                    networks = skeletal_muscle, alpha = 0.05, hopConstraint = 2,
                                    nullDistribution = null, pValueFile = paste0(outdir, "/skeletalMusclePval.RDS"),
                                    verbose = TRUE,
                                    doFDRAdjustment = TRUE)
write.csv(skeletalMuscleSubnet, paste0(outdir, "/skeletalMuscleSubnet_05_updatedNull.csv"))
skeletalMuscleSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                                    networks = skeletal_muscle, alpha = 0.05, hopConstraint = 2,
                                    nullDistribution = null, pValueFile = paste0(outdir, "/skeletalMusclePval.RDS"),
                                    verbose = TRUE, loadPValues = TRUE,
                                    doFDRAdjustment = TRUE)
write.csv(skeletalMuscleSubnet, paste0(outdir, "/skeletalMuscleSubnet_05_updatedNull_moreGenes.csv"))

# Run BLOBFISH in the subcutaneous adipose.
adiposeSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                                    networks = subcutaneous_adipose, alpha = 0.05, hopConstraint = 2,
                                    nullDistribution = null, pValueFile = paste0(outdir, "/adiposePval.RDS"),
                                    verbose = TRUE, loadPValues = TRUE,
                                    doFDRAdjustment = TRUE)
write.csv(adiposeSubnet, paste0(outdir, "/adiposeSubnet_05_updatedNull_moreGenes.csv"))
skinSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                                    networks = skin, alpha = 0.05, hopConstraint = 2,
                                    nullDistribution = null, pValueFile = paste0(outdir, "/skinPval.RDS"),
                                    verbose = TRUE, loadPValues = TRUE,
                                    doFDRAdjustment = TRUE)
write.csv(skinSubnet, paste0(outdir, "/skinSubnet_05_updatedNull_moreGenes.csv"))
lungSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                          networks = lung, alpha = 0.05, hopConstraint = 2,
                          nullDistribution = null, pValueFile = paste0(outdir, "/lungPval.RDS"),
                          verbose = TRUE, #loadPValues = TRUE,
                          doFDRAdjustment = TRUE)
write.csv(lungSubnet, paste0(outdir, "/lungSubnet_05_updatedNull_moreGenes.csv"))
str("lung done")
aortaSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                           networks = aorta, alpha = 0.05, hopConstraint = 2,
                           nullDistribution = null, pValueFile = paste0(outdir, "/aortaPval.RDS"),
                           verbose = TRUE, #loadPValues = TRUE,
                           doFDRAdjustment = TRUE)
write.csv(aortaSubnet, paste0(outdir, "/aortaSubnet_05_updatedNull_moreGenes.csv"))
str("aorta done")

# Adipogenic markers from 36647068
geneColorMapping <- data.frame(gene = c(muscle_dev_genes, adipogenic_genes), 
                               color = c(rep("hotpink", length(muscle_dev_genes)), 
                                           rep("goldenrod1", length(adipogenic_genes)))) 

# Plot subcutaneous adipose.
adiposeSubnet$gene <- unlist(lapply(adiposeSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
pdf(paste0(outdir, "/subcutaneous_adipose_subnet_updatedNull.pdf"))
PlotNetwork(network = adiposeSubnet,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
dev.off()

# Plot skeletal muscle.
skeletalMuscleSubnet$gene <- unlist(lapply(skeletalMuscleSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
pdf(paste0(outdir, "/skeletal_muscle_subnet_updateNull.pdf"),)
PlotNetwork(network = skeletalMuscleSubnet,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
dev.off()

# Plot skin.
skinSubnet$gene <- unlist(lapply(skinSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
pdf(paste0(outdir, "/skin_subnet_updateNull.pdf"),)
PlotNetwork(network = skinSubnet, genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
dev.off()

# Plot lung.
lungSubnet$gene <- unlist(lapply(lungSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
pdf(paste0(outdir, "/lung_subnet_updateNull.pdf"),)
PlotNetwork(network = lungSubnet, genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
dev.off()

# Plot aorta.
aortaSubnet$gene <- unlist(lapply(aortaSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
pdf(paste0(outdir, "/aorta_subnet_updateNull.pdf"),)
PlotNetwork(network = aortaSubnet, genesOfInterest = c(muscle_dev_genes, adipogenic_genes),vertexLabels = NULL,
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.50, geneColorMapping = geneColorMapping)
dev.off()

# Get the exclusive networks.
skeletal_muscle_only <- skeletalMuscleSubnet[setdiff(rownames(skeletalMuscleSubnet), 
                                                c(rownames(adiposeSubnet),
                                                  rownames(skinSubnet), rownames(lungSubnet),
                                                  rownames(aortaSubnet))),]
write.csv(skeletal_muscle_only, paste0(outdir, "/skeletal_muscle_only_updated.csv"))
pdf(paste0(outdir, "/skeletal_muscle_only_updated_withlabels.pdf"))
PlotNetwork(network = skeletal_muscle_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, skeletal_muscle_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

subcutaneous_adipose_only <- adiposeSubnet[setdiff(rownames(adiposeSubnet), 
                                                c(rownames(skeletalMuscleSubnet),
                                                  rownames(skinSubnet), rownames(lungSubnet),
                                                  rownames(aortaSubnet))),]
write.csv(subcutaneous_adipose_only, paste0(outdir, "/subcutaneous_adipose_only_updated.csv"))
pdf(paste0(outdir, "/adipose_only_updated_withlabels.pdf"))
PlotNetwork(network = subcutaneous_adipose_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, subcutaneous_adipose_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

skin_only <- skinSubnet[setdiff(rownames(skinSubnet), 
                                                          c(rownames(skeletalMuscleSubnet),
                                                            rownames(adiposeSubnet),
                                                            rownames(lungSubnet), rownames(aortaSubnet))),]
write.csv(skin_only, paste0(outdir, "/skin_only_updated.csv"))
pdf(paste0(outdir, "/skin_only_updated_withlabels.pdf"))
PlotNetwork(network = skin_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, skin_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

lung_only <- lungSubnet[setdiff(rownames(lungSubnet), 
                          c(rownames(skeletalMuscleSubnet),
                            rownames(adiposeSubnet),
                            rownames(skinSubnet), rownames(aortaSubnet))),]
write.csv(lung_only, paste0(outdir, "/lung_only_updated.csv"))
pdf(paste0(outdir, "/lung_only_updated_withlabels.pdf"))
PlotNetwork(network = lung_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, lung_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

aorta_only <- aortaSubnet[setdiff(rownames(aortaSubnet), 
                          c(rownames(skeletalMuscleSubnet),
                            rownames(adiposeSubnet),
                            rownames(skinSubnet), rownames(lungSubnet))),]
write.csv(aorta_only, paste0(outdir, "/aorta_only_updated.csv"))
pdf(paste0(outdir, "/aorta_only_updated_withlabels.pdf"))
PlotNetwork(network = aorta_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, aorta_only$tf),
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
tfCountAllPairsMuscleMuscle <- TFCountAllPairs(muscle_dev_genes, skeletal_muscle_only)
print(mean(tfCountAllPairsMuscleMuscle))
tfCountAllPairsMuscleAdipose <- TFCountAllPairs(adipogenic_genes, skeletal_muscle_only)
print(mean(tfCountAllPairsMuscleAdipose))
tfCountAllPairsAdiposeMuscle <- TFCountAllPairs(muscle_dev_genes, subcutaneous_adipose_only)
print(mean(tfCountAllPairsAdiposeMuscle))
tfCountAllPairsAdiposeAdipose <- TFCountAllPairs(adipogenic_genes, subcutaneous_adipose_only)
print(mean(tfCountAllPairsAdiposeAdipose))
tfCountAllPairsSkinMuscle <- TFCountAllPairs(muscle_dev_genes, skin_only)
print(mean(tfCountAllPairsSkinMuscle))
tfCountAllPairsSkinAdipose <- TFCountAllPairs(adipogenic_genes, skin_only)
print(mean(tfCountAllPairsSkinAdipose))
tfCountAllPairsLungMuscle <- TFCountAllPairs(muscle_dev_genes, lung_only)
print(mean(tfCountAllPairsLungMuscle))
tfCountAllPairsLungAdipose <- TFCountAllPairs(adipogenic_genes, lung_only)
print(mean(tfCountAllPairsLungAdipose))
tfCountAllPairsAortaMuscle <- TFCountAllPairs(muscle_dev_genes, aorta_only)
print(mean(tfCountAllPairsAortaMuscle))
tfCountAllPairsAortaAdipose <- TFCountAllPairs(adipogenic_genes, aorta_only)
print(mean(tfCountAllPairsAortaAdipose))

# Do pathway analysis.
library(fgsea)
pathways = fgsea::gmtPathways(gmt_pathway_file)
inputMuscleOnly <- c(table(skeletal_muscle_only$tf), table(skeletal_muscle_only$gene))
pathwayResultMuscle <- fgsea::fgsea(pathways = pathways, stats = inputMuscleOnly,
                              scoreType = "pos")
inputAdiposeOnly <- c(table(subcutaneous_adipose_only$tf), table(subcutaneous_adipose_only$gene))
pathwayResultAdipose <- fgsea::fgsea(pathways = pathways, stats = inputAdiposeOnly,
                              scoreType = "pos")
skinOnly <- c(table(skin_only$tf), table(skin_only$gene))
pathwayResultSkin <- fgsea::fgsea(pathways = pathways, stats = skinOnly,
                                    scoreType = "pos")
lungOnly <- c(table(lung_only$tf), table(lung_only$gene))
pathwayResultLung <- fgsea::fgsea(pathways = pathways, stats = lungOnly,
                                  scoreType = "pos")
aortaOnly <- c(table(aorta_only$tf), table(aorta_only$gene))
pathwayResultAorta <- fgsea::fgsea(pathways = pathways, stats = aortaOnly,
                                  scoreType = "pos")
write.csv(pathwayResultMuscle[, 1:7], paste0(outdir, "/muscle_only_pathways.csv"))
write.csv(pathwayResultAdipose[, 1:7], paste0(outdir, "/adipose_only_pathways.csv"))
write.csv(pathwayResultSkin[, 1:7], paste0(outdir, "/skin_only_pathways.csv"))
write.csv(pathwayResultLung[, 1:7], paste0(outdir, "/lung_only_pathways.csv"))
write.csv(pathwayResultAorta[, 1:7], paste0(outdir, "/aorta_only_pathways.csv"))

# Do pathway analysis for differentially expressed genes only. Modify paths as needed.
expressionAdipose <- read.csv(paste0(dir_expression, "/Adipose_Subcutaneous.csv"), row.names = 1)
expressionAdiposeSubset <- expressionAdipose[c(muscle_dev_ensembl, adipogenic_ensembl),
                                             make.names(c("GTEX-11P82-1726-SM-5Q5AT", "GTEX-11TT1-2426-SM-5EQMK",
                                                          "GTEX-12126-0426-SM-5Q5AP", "GTEX-12C56-1626-SM-5FQUO",
                                                          "GTEX-144GM-1926-SM-5LUAN", "GTEX-O5YT-0226-SM-32PK5",
                                                          "GTEX-QEG5-0326-SM-2S1PB", "GTEX-QESD-1526-SM-2S1QT",
                                                          "GTEX-S4Q7-1626-SM-3K2AE", "GTEX-SNMC-1326-SM-2XCFK",
                                                          "GTEX-T5JC-0526-SM-32PM7", "GTEX-WFG8-2326-SM-3GILF",
                                                          "GTEX-WFG8-2326-SM-5S2UE", "GTEX-XQ3S-1626-SM-4WAYN",
                                                          "GTEX-YJ8A-0526-SM-5IFHT", "GTEX-ZTX8-1526-SM-5N9GI",
                                                          "GTEX-ZVTK-0526-SM-5GZWT", "GTEX-ZZ64-1626-SM-5E43W"))]
expressionMuscle <- read.csv(paste0(dir_expression, "/Skeletal_Muscle.csv"), row.names = 1)
expressionMuscleSubset <- expressionMuscle[c(muscle_dev_ensembl, adipogenic_ensembl),
                                           make.names(c("GTEX-11NSD-2026-SM-5HL5U", "GTEX-11P82-1826-SM-5PNYJ",
                                                        "GTEX-11TT1-2326-SM-5GU6N", "GTEX-12126-0526-SM-5PNW4",
                                                        "GTEX-12C56-1926-SM-5FQUG", "GTEX-13RTK-1626-SM-5QGQS",
                                                        "GTEX-144GM-2026-SM-5LU3D", "GTEX-O5YT-1626-SM-32PK6",
                                                        "GTEX-PW2O-1726-SM-2S1OO", "GTEX-QEG5-0426-SM-2I5GJ",
                                                        "GTEX-QESD-1626-SM-2S1RB", "GTEX-R55E-0526-SM-2TC6B",
                                                        "GTEX-S4Q7-1526-SM-3K2AG", "GTEX-SNMC-1426-SM-2XCFM",
                                                        "GTEX-T5JC-0626-SM-3NMA6", "GTEX-WFG7-2226-SM-3GIKP",
                                                        "GTEX-WFG8-2426-SM-3GILL", "GTEX-WHSE-0526-SM-4SOJ6",
                                                        "GTEX-WHSE-0526-SM-5GZZD", "GTEX-WHSE-0526-SM-5HL79",
                                                        "GTEX-WOFM-1326-SM-3MJFR", "GTEX-XPT6-2026-SM-4B64V",
                                                        "GTEX-XQ3S-0426-SM-4BOOA", "GTEX-ZA64-2026-SM-5PNXT",
                                                        "GTEX-ZTX8-1626-SM-51MRY", "GTEX-ZVTK-0626-SM-5GU6Y",
                                                        "GTEX-ZZ64-1526-SM-5E43K"))]
expressionLung <- read.csv(paste0(dir_expression, "/Lung.csv"), row.names = 1)
expressionLungSubset <- expressionLung[c(muscle_dev_ensembl, adipogenic_ensembl),
                                       make.names(c("GTEX-12WSE-0826-SM-5S2VL", "GTEX-144GM-0126-SM-5Q5AX",
                                                    "GTEX-O5YT-0526-SM-32PK8", "GTEX-PW2O-0526-SM-2I3DX",
                                                    "GTEX-QEG5-1126-SM-2I5GH", "GTEX-QEG5-1126-SM-33HC2",
                                                    "GTEX-QEG5-1126-SM-5SI88", "GTEX-QESD-0626-SM-2I5G4",
                                                    "GTEX-S4Q7-0426-SM-3K2BJ", "GTEX-T5JC-0826-SM-32PMC",
                                                    "GTEX-WFG7-0526-SM-3GIKI", "GTEX-WFG8-0926-SM-3GIKJ",
                                                    "GTEX-WOFM-0126-SM-3MJFE", "GTEX-XQ3S-0926-SM-4BOPI",
                                                    "GTEX-ZA64-0326-SM-5HL8T"))]
expressionSkin <- read.csv(paste0(dir_expression, "/Skin.csv"), row.names = 1)
expressionSkinSubset <- expressionSkin[c(muscle_dev_ensembl, adipogenic_ensembl),
                                       make.names(c("GTEX-11NSD-0626-SM-5A5LU", "GTEX-11NSD-2226-SM-5986V",
                                                    "GTEX-11P82-1426-SM-5BC5L", "GTEX-11P82-1626-SM-59879",
                                                    "GTEX-11TT1-1826-SM-5EQLV", "GTEX-11TT1-2526-SM-5EGIB",
                                                    "GTEX-12126-0326-SM-5PNW3", "GTEX-12126-0826-SM-5FQTZ",
                                                    "GTEX-12C56-1526-SM-5FQUQ", "GTEX-12C56-1726-SM-5EQ61",
                                                    "GTEX-O5YT-0126-SM-48TBW", "GTEX-QEG5-0226-SM-2I5GI",
                                                    "GTEX-QESD-1426-SM-2S1R9", "GTEX-R55E-0126-SM-2TC5Y",
                                                    "GTEX-SNMC-1226-SM-2XCFP", "GTEX-T5JC-0426-SM-32PLO",
                                                    "GTEX-WFG7-2026-SM-3GIL7", "GTEX-WFG7-2026-SM-5SI7P",
                                                    "GTEX-WFG8-2126-SM-3GIKQ", "GTEX-WFG8-2226-SM-3GIL9",
                                                    "GTEX-WHSE-0126-SM-3NMBT", "GTEX-XQ3S-1426-SM-4BOPR",
                                                    "GTEX-XQ3S-1526-SM-4BOOC", "GTEX-YJ8A-0426-SM-5IFID",
                                                    "GTEX-ZA64-1026-SM-5HL8R", "GTEX-ZTX8-1326-SM-4YCEL",
                                                    "GTEX-ZTX8-1426-SM-5DUVO", "GTEX-ZVTK-0226-SM-51MRG",
                                                    "GTEX-ZVTK-0426-SM-51MRS", "GTEX-ZZ64-1026-SM-5GZXG",
                                                    "GTEX-ZZ64-1726-SM-5GZYB"))]
expressionAorta <- read.csv(paste0(dir_expression, "/Artery_Aorta.csv"), row.names = 1)
expressionAortaSubset <- expressionAorta[c(muscle_dev_ensembl, adipogenic_ensembl),
                                         make.names(c("GTEX-11P82-0126-SM-5HL72", "GTEX-11TT1-1226-SM-5Q5AV",
                                                      "GTEX-12C56-0226-SM-5N9FB", "GTEX-13RTK-0126-SM-5RQHQ",
                                                      "GTEX-O5YT-0426-SM-3MJHD", "GTEX-PW2O-0426-SM-48TCC",
                                                      "GTEX-QEG5-1226-SM-447AR", "GTEX-QESD-0426-SM-4R1JZ",
                                                      "GTEX-R55E-1126-SM-48FDZ", "GTEX-T5JC-0326-SM-4DM5C",
                                                      "GTEX-WFG7-0926-SM-4LMK7", "GTEX-WFG8-0826-SM-4LVN5",
                                                      "GTEX-XPT6-0326-SM-4B66V", "GTEX-XQ3S-0126-SM-4BOO9",
                                                      "GTEX-ZA64-0126-SM-5HL8Z"))]
# Formula for log fold change.
calculateLogFoldChange <- function(expressionCase, expressionControl){
  return(rowMeans(log2(expressionCase)) - rowMeans(log2(expressionControl)))
}

# Calculate log fold changes.
logFCAdipose <- calculateLogFoldChange(expressionAdiposeSubset, 
                                       do.call(cbind, list(expressionMuscleSubset,
                                                           expressionAortaSubset,
                                                           expressionSkinSubset,
                                                           expressionLungSubset)))
write.csv(logFCAdipose, paste0(dir_expression, "/adiposeLogFC.csv"))
logFCMuscle <- calculateLogFoldChange(expressionMuscleSubset, 
                                      do.call(cbind, list(expressionAdiposeSubset,
                                                          expressionAortaSubset,
                                                          expressionSkinSubset,
                                                          expressionLungSubset)))
write.csv(logFCMuscle, paste0(dir_expression, "/muscleLogFC.csv"))
logFCAorta <- calculateLogFoldChange(expressionAortaSubset, 
                                     do.call(cbind, list(expressionAdiposeSubset,
                                                         expressionMuscleSubset,
                                                         expressionSkinSubset,
                                                         expressionLungSubset)))
write.csv(logFCAorta, paste0(dir_expression, "/aortaLogFC.csv"))
logFCLung <- calculateLogFoldChange(expressionLungSubset, 
                                    do.call(cbind, list(expressionAdiposeSubset,
                                                        expressionMuscleSubset,
                                                        expressionSkinSubset,
                                                        expressionAortaSubset)))
write.csv(logFCLung, paste0(dir_expression, "/lungLogFC.csv"))
logFCSkin <- calculateLogFoldChange(expressionSkinSubset, 
                                    do.call(cbind, list(expressionAdiposeSubset,
                                                        expressionMuscleSubset,
                                                        expressionLungSubset,
                                                        expressionAortaSubset)))
write.csv(logFCSkin, paste0(dir_expression, "/skinLogFC.csv"))


# Do pathway analysis.
doPathwayAnalysis <- function(logFC, file){
  names(logFC) <- c(muscle_dev_genes, adipogenic_genes)
  str(logFC)
  res <- fgsea::fgsea(pathways = pathways, stats = logFC,
                      scoreType = "pos", maxSize = 500)
  str(res)
  write.csv(res[,c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")], 
            file)
}
doPathwayAnalysis(logFCAdipose,  paste0(dir_expression, "/pathwayResultAdipose.csv"))
doPathwayAnalysis(logFCMuscle, paste0(dir_expression, "/pathwayResultMuscle.csv"))
doPathwayAnalysis(logFCAorta, paste0(dir_expression, "/pathwayResultAorta.csv"))
doPathwayAnalysis(logFCLung, paste0(dir_expression, "/pathwayResultLung.csv"))
doPathwayAnalysis(logFCSkin, paste0(dir_expression, "/pathwayResultSkin.csv"))