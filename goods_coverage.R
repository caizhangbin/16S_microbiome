

library(QsRutils)

# load your OTU table as a community data matrix
# make sure that samples are rows and OTUs are columns
com <- read.csv("/rarefied.removed.low.reads_OTU.csv", row.names = 1)

# calculate Good's coverage
coverage <- goods(com)

write.csv(coverage, "coverage.csv", row.names = TRUE)


