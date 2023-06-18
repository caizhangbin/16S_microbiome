library(tidyverse)


original_data <- read.csv("rarefied.removed.low.reads_OTU.csv", row.names = 1)
Group <- original_data$Group 

otu_table <- original_data%>%
  select(-Group)

# Calculate relative abundance
rel_abund <- t(t(otu_table) / rowSums(otu_table))

# Convert relative abundance to percentage
rel_abund_pct <- rel_abund * 100

rel_abund_group <- cbind(Group, rel_abund_pct)

write.csv(rel_abund_group, "rel_abund.csv", row.names = TRUE)

