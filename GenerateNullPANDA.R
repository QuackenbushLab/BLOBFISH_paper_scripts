setwd("/home/ubuntu/netZooR")
roxygen2::roxygenize()

null <- GenerateNullPANDADistribution(ppiFile = "/home/ubuntu/ppi.txt", motifFile = "/home/ubuntu/motif.txt")
saveRDS(null, "/home/ubuntu/nullPANDA.RDS")
