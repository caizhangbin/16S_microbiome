library(vegan)
library(dplyr)

# Load OTU data
original_data <- read_csv("/removed.low.reads_OTU.csv") 

groups <- original_data$Group 

otu_data <- original_data %>%
  select(-Group)


# Remove rows with missing values
otu_data <- na.omit(otu_data)

# Convert otu_data to a numeric matrix
otu_matrix <- as.matrix(otu_data)

# Set the minimum number of sequences
min_seqs <- min(rowSums(otu_matrix))

# Set the number of permutations
nperm <- 1000

# Initialize a matrix to store the sum of the rarefied OTU tables
sum_rarefied_otu_tables <- matrix(0, nrow = nrow(otu_data), ncol = ncol(otu_data))

# Perform the rarefaction and calculate the sum of the rarefied OTU tables
for (i in 1:nperm) {
  sum_rarefied_otu_tables <- sum_rarefied_otu_tables + rrarefy(otu_data, min_seqs)
}

# Calculate the average rarefied OTU table
avg_rarefied_otu_table <- sum_rarefied_otu_tables / nperm

rarefied_table_with_info <- cbind(groups, avg_rarefied_otu_table)

write.csv(rarefied_table_with_info, "rarefied_OTU.csv", row.names = FALSE)

