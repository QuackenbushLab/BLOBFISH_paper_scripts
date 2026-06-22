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

# Read in the files from the eQTL directory.
eqtlDir <- NULL
dummyPDir <- NULL
refFile <- NULL
ref <- as.data.frame(fread(refFile))
dir.create(dummyPDir)

# Function to build BLOBFISH formatted networks
formatBlobfish <- function(rds, pvals){
  # Extract edges.
  data <- readRDS(rds)$edges
  colnames(data) <- c("tf", "gene", "score")
  
  # Remove version from ENSG IDs.
  data$gene <- unlist(lapply(data$gene, function(gene){
    return(strsplit(gene, split = ".", fixed = TRUE)[[1]][1])
  }))
  
  # Add row names.
  rownames(data) <- paste(data$tf, data$gene, sep = "__")
  
  # Create dummy p-values.
  dummyP <- rep(0.001, nrow(data))
  saveRDS(dummyP, pvals)
  
  # Return the data.
  return(data)
}
adipose <- formatBlobfish(paste0(eqtlDir, "/network_qscores_Adipose_Subcutaneous.Rds"),
                       paste0(dummyPDir, "/dummyPvals_Adipose_Subcutaneous.Rds"))
aorta <- formatBlobfish(paste0(eqtlDir, "/network_qscores_Artery_Aorta.Rds"),
                     paste0(dummyPDir, "/dummyPvals_Artery_Aorta.Rds"))
muscle <- formatBlobfish(paste0(eqtlDir, "/network_qscores_Muscle_Skeletal.Rds"),
                        paste0(dummyPDir, "/dummyPvals_Muscle.Rds"))
lung <- formatBlobfish(paste0(eqtlDir, "/network_qscores_Lung.Rds"),
                         paste0(dummyPDir, "/dummyPvals_Lung.Rds"))
skin <- formatBlobfish(paste0(eqtlDir, "/network_qscores_Skin_Sun_Exposed_Lower_leg.Rds"),
                       paste0(dummyPDir, "/dummyPvals_Skin.Rds"))

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

# Run BLOBFISH.
dummyNull <- c(1,2,3)
adiposeSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                                    networks = adipose, alpha = 0.05, hopConstraint = 2, loadPValues = TRUE,
                                    nullDistribution = dummyNull, pValueFile = paste0(dummyPDir, "/dummyPvals_Adipose_Subcutaneous.Rds"),
                                    verbose = TRUE)
write.csv(adiposeSubnet, paste0(eqtlDir, "/adiposeSubnet_eQTL.csv"))
muscleSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                             networks = muscle, alpha = 0.05, hopConstraint = 2, loadPValues = TRUE,
                             nullDistribution = dummyNull, pValueFile = paste0(dummyPDir, "/dummyPvals_Muscle.Rds"),
                             verbose = TRUE)
write.csv(muscleSubnet, paste0(eqtlDir, "/muscleSubnet_eQTL.csv"))
lungSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                            networks = lung, alpha = 0.05, hopConstraint = 2, loadPValues = TRUE,
                            nullDistribution = dummyNull, pValueFile = paste0(dummyPDir, "/dummyPvals_Lung.Rds"),
                            verbose = TRUE)
write.csv(lungSubnet, paste0(eqtlDir, "/lungSubnet_eQTL.csv"))
skinSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                          networks = skin, alpha = 0.05, hopConstraint = 2, loadPValues = TRUE,
                          nullDistribution = dummyNull, pValueFile = paste0(dummyPDir, "/dummyPvals_Skin.Rds"),
                          verbose = TRUE)
write.csv(skinSubnet, paste0(eqtlDir, "/skinSubnet_eQTL.csv"))
aortaSubnet <- RunBLOBFISH(geneSet = c(unname(adipogenic_ensembl),unname(muscle_dev_ensembl)),
                          networks = aorta, alpha = 0.05, hopConstraint = 2, loadPValues = TRUE,
                          nullDistribution = dummyNull, pValueFile = paste0(dummyPDir, "/dummyPvals_Artery_Aorta.Rds"),
                          verbose = TRUE)
write.csv(aortaSubnet, paste0(eqtlDir, "/aortaSubnet_eQTL.csv"))


# Adipogenic markers from 36647068
geneColorMapping <- data.frame(gene = c(muscle_dev_genes, adipogenic_genes), 
                               color = c(rep("hotpink", length(muscle_dev_genes)), 
                                           rep("goldenrod1", length(adipogenic_genes)))) 

# Get the exclusive networks.
muscleSubnet$gene <- unlist(lapply(muscleSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
skeletal_muscle_only <- muscleSubnet[setdiff(rownames(muscleSubnet), 
                                                c(rownames(adiposeSubnet),
                                                  rownames(skinSubnet), rownames(lungSubnet),
                                                  rownames(aortaSubnet))),]
write.csv(skeletal_muscle_only, paste0(eqtlDir, "/skeletal_muscle_only_eQTL.csv"))
pdf(paste0(eqtlDir, "/skeletal_muscle_only_updated_withlabels.pdf"))
PlotNetwork(network = skeletal_muscle_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, skeletal_muscle_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

adiposeSubnet$gene <- unlist(lapply(adiposeSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
subcutaneous_adipose_only <- adiposeSubnet[setdiff(rownames(adiposeSubnet), 
                                                c(rownames(muscleSubnet),
                                                  rownames(skinSubnet), rownames(lungSubnet),
                                                  rownames(aortaSubnet))),]
write.csv(subcutaneous_adipose_only, paste0(eqtlDir, "/subcutaneous_adipose_only_eQTL.csv"))
pdf(paste0(eqtlDir, "/adipose_only_updated_withlabels.pdf"))
PlotNetwork(network = subcutaneous_adipose_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, subcutaneous_adipose_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

skinSubnet$gene <- unlist(lapply(skinSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
skin_only <- skinSubnet[setdiff(rownames(skinSubnet), 
                                                          c(rownames(muscleSubnet),
                                                            rownames(adiposeSubnet),
                                                            rownames(lungSubnet), rownames(aortaSubnet))),]
write.csv(skin_only, paste0(eqtlDir, "/skin_only_eQTL.csv"))
pdf(paste0(eqtlDir, "/skin_only_updated_withlabels.pdf"))
PlotNetwork(network = skin_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, skin_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

lungSubnet$gene <- unlist(lapply(lungSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
lung_only <- lungSubnet[setdiff(rownames(lungSubnet), 
                          c(rownames(muscleSubnet),
                            rownames(adiposeSubnet),
                            rownames(skinSubnet), rownames(aortaSubnet))),]
write.csv(lung_only, paste0(eqtlDir, "/lung_only_eQTL.csv"))
pdf(paste0(eqtlDir, "/lung_only_updated_withlabels.pdf"))
PlotNetwork(network = lung_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, lung_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

aortaSubnet$gene <- unlist(lapply(aortaSubnet$gene, function(gene){
  retval <- names(muscle_dev_ensembl)[which(muscle_dev_ensembl == gene)]
  if(length(retval) == 0){
    retval <- names(adipogenic_ensembl)[which(adipogenic_ensembl == gene)]
    str(retval)
  }
  return(retval[1])
}))
aorta_only <- aortaSubnet[setdiff(rownames(aortaSubnet), 
                          c(rownames(muscleSubnet),
                            rownames(adiposeSubnet),
                            rownames(skinSubnet), rownames(lungSubnet))),]
write.csv(aorta_only, paste0(eqtlDir, "/aorta_only_eQTL.csv"))
pdf(paste0(eqtlDir, "/aorta_only_updated_withlabels.pdf"))
PlotNetwork(network = aorta_only,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, aorta_only$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 0.2, geneColorMapping = geneColorMapping)
dev.off()

# Map the SNPs to genes.
mapSNPs <- function(eqtls, ref){
  # Separate out the components of the ID.
  chrom <- unlist(lapply(eqtls$tf, function(snp){
    return(strsplit(snp, "_")[[1]][1])
  }))
  pos <- as.numeric(unlist(lapply(eqtls$tf, function(snp){
    return(strsplit(snp, "_")[[1]][2])
  })))
  snpInfo <- data.frame(CHROM = chrom, POS = pos)
  snpGene <- do.call(rbind, lapply(1:nrow(snpInfo), function(i){
    
    # Subset chromosome.
    chrName <- as.character(snpInfo[i,"CHROM"])
    refChrom <- ref[which(ref$V1 == chrName),]
    
    # Find starting points where SNP position is greater.
    refGreater <- refChrom[which(refChrom$V2 <= snpInfo[i,"POS"]),]
    
    # Find ending points where SNP position is less.
    refLess <- refGreater[which(refGreater$V3 >= snpInfo[i,"POS"]),]
    
    # Order of importance.
    refImportant <- refLess
    ooi <- c("CDS", "start_codon", "stop_codon", "Selenocysteine", "exon", "gene", "transcript",
             "UTR")
    if(nrow(refImportant) > 0){
      found <- FALSE
      j = 1
      while(found == FALSE && j < length(ooi)){
        if(ooi[j] %in% refLess$V8){
          refImportant <- refLess[which(refLess$V8 == ooi[j]),]
          found <- TRUE
        }else{
          j <- j + 1
        }
      }
    }else{
      # If the variant doesn't map to a gene, then map it to its closest gene.
      closestStart <- refChrom[which.min(abs(snpInfo[i,"POS"] - refChrom$V2)),]
      closestEnd <- refChrom[which.min(abs(refChrom$V3 - snpInfo[i,"POS"])),]
      refImportant <- closestStart
      if(abs(snpInfo[i,"POS"] - closestStart$V2) < abs(closestEnd$V3 - snpInfo[i,"POS"])){
        refImportant <- closestEnd
      }
    }
    snpDetails <- refImportant[,ncol(refImportant)]
    if(nrow(refImportant) > 1){
      refImportant <- refImportant[1,]
    }
    
    refImportant <- refImportant[,c(2:4, 6, 8)]
    colnames(refImportant) <- c("geneStart", "geneEnd", "name", "strand", "region")

    # Paste together the SNP and gene information.
    snpGeneInfo <- do.call(cbind, list(eqtls[i,], snpInfo[i,], refImportant))
    splitDetails <- strsplit(snpDetails, "; ")[[1]]
    geneType <- splitDetails[startsWith(splitDetails, "gene_type")]
    geneTypeSplit <- strsplit(geneType, split = "\"", fixed = TRUE)[[1]][2]
    geneName <- splitDetails[startsWith(splitDetails, "gene_name")]
    geneNameSplit <- strsplit(geneName, split = "\"", fixed = TRUE)[[1]][2]
    snpGeneInfo$geneName <- geneNameSplit
    snpGeneInfo$geneType <- geneTypeSplit
    return(snpGeneInfo)
  }))
  str(snpGene)
  return(snpGene)
}
aortaSNPs <- mapSNPs(aorta_only, ref)
adiposeSNPs <- mapSNPs(subcutaneous_adipose_only, ref)
muscleSNPs <- mapSNPs(skeletal_muscle_only, ref)
skinSNPs <- mapSNPs(skin_only, ref)
lungSNPs <- mapSNPs(lung_only, ref)

# Save the results.
write.csv(aortaSNPs, paste0(eqtlDir, "/aorta_only_snpLabels.csv"))
write.csv(adiposeSNPs, paste0(eqtlDir, "/adipose_only_snpLabels.csv"))
write.csv(muscleSNPs, paste0(eqtlDir, "/muscle_only_snpLabels.csv"))
write.csv(skinSNPs, paste0(eqtlDir, "/skin_only_snpLabels.csv"))
write.csv(lungSNPs, paste0(eqtlDir, "/lung_only_snpLabels.csv"))

# Consolidate SNP networks by gene.
adiposeConsol <- adiposeSNPs[,c("geneName", "gene")]
colnames(adiposeConsol)[1] <- "tf"
pdf(paste0(eqtlDir, "/adipose_only_consolidated_withlabels.pdf"))
PlotNetwork(network = adiposeConsol,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, adiposeConsol$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 1, geneColorMapping = geneColorMapping,
            nodeSize = 3)
dev.off()

muscleConsol <- muscleSNPs[,c("geneName", "gene")]
colnames(muscleConsol)[1] <- "tf"
pdf(paste0(eqtlDir, "/muscle_only_consolidated_withlabels.pdf"))
PlotNetwork(network = muscleConsol,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, muscleConsol$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 1, geneColorMapping = geneColorMapping,
            nodeSize = 3)
dev.off()

aortaConsol <- aortaSNPs[,c("geneName", "gene")]
colnames(aortaConsol)[1] <- "tf"
pdf(paste0(eqtlDir, "/aorta_only_consolidated_withlabels.pdf"))
PlotNetwork(network = aortaConsol,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, aortaConsol$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 1, geneColorMapping = geneColorMapping,
            nodeSize = 3)
dev.off()

lungConsol <- lungSNPs[,c("geneName", "gene")]
colnames(lungConsol)[1] <- "tf"
pdf(paste0(eqtlDir, "/lung_only_consolidated_withlabels.pdf"))
PlotNetwork(network = lungConsol,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, lungConsol$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 1, geneColorMapping = geneColorMapping,
            nodeSize = 3)
dev.off()

skinConsol <- skinSNPs[,c("geneName", "gene")]
colnames(skinConsol)[1] <- "tf"
pdf(paste0(eqtlDir, "/skin_only_consolidated_withlabels.pdf"))
PlotNetwork(network = skinConsol,genesOfInterest = c(muscle_dev_genes, adipogenic_genes),
            vertexLabels = c(muscle_dev_genes, adipogenic_genes, skinConsol$tf),
            layoutBipartite = FALSE,edgeWidth = 1, vertexLabelSize = 1, geneColorMapping = geneColorMapping,
            nodeSize = 3)
dev.off()

# Analyze why the adipogenic genes did not overlap by looking at the original eQTLs.
adiposeAdipose <- adipose[which(adipose$gene %in% adipogenic_ensembl),]
aortaAdipose <- aortaAdipose[which(aorta$gene %in% adipogenic_ensembl),]
muscleAdipose <- muscle[which(muscle$gene %in% adipogenic_ensembl),]
lungAdipose <- lung[which(lung$gene %in% adipogenic_ensembl),]
skinAdipose <- skin[which(skin$gene %in% adipogenic_ensembl),]
adiposeAdiposeOnly <- adiposeAdipose[setdiff(rownames(adiposeAdipose), c(rownames(lungAdipose), rownames(muscleAdipose), rownames(skinAdipose), rownames(aortaAdipose)))]
adiposeAdiposeOnlyGenes <- adiposeAdiposeOnly
adiposeAdiposeOnlyGenes$gene <- unlist(lapply(adiposeAdiposeOnly$gene, function(gene){
  return(adipogenic_genes[which(adipogenic_ensembl == gene)])
}))

PlotNetwork(network = adiposeAdiposeOnlyGenes, genesOfInterest = adipogenic_genes, 
            vertexLabels = adipogenic_genes, layoutBipartite = FALSE, edgeWidth = 1, vertexLabelSize = 0.5,
            geneColorMapping = geneColorMapping, nodeSize = 3)
