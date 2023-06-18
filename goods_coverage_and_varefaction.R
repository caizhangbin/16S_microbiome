library(vegan)
library(tidyverse)
library(dplyr)
library(QsRutils)

shared<-read_csv("/Users/zhangbincai/Desktop/study/projects/sherbrook innoculation/second round 16s/analysis/removed.low.reads_OTU.csv") %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  #summarize(totall= sum(value))%>%
  #arrange(totall) %>%
  group_by(name)%>%
  mutate(total = sum(value)) %>% 
  filter(total!=0)%>%
  ungroup()%>%
  select(-total)

rand <- shared %>%
  uncount (value) %>%
  mutate(name = sample (name)) %>%
  count (Group, name, name ="value")

sampling_coverage <- shared %>%
  group_by(Group) %>%
  summarise(n_seqs = sum(value))

sampling_coverage %>%
  ggplot(aes( x =n_seqs)) +
  geom_histogram(binwidth = 500) +
  coord_cartesian(xlim = c(0, 5000))


sampling_coverage %>%
  ggplot(aes(x=1, y=n_seqs)) +
  geom_jitter() +
  scale_y_log10() +
  geom_hline(yintercept = 3210)

sampling_coverage %>%
  arrange(n_seqs)%>%
  ggplot(aes(x=1:nrow(.), y=n_seqs)) +
  geom_line() +
  coord_cartesian(xlim = c(0, 50), ylim = c(0,5000))

sampling_coverage %>%
  arrange(n_seqs)%>%
  print(n=60)

# calculate Good's coverage
com <- read.csv("/Users/zhangbincai/rarefied_OTU.csv") %>%
  select(starts_with("Otu"))
coverage <- goods(com)

write.csv(coverage, "coverage.csv", row.names = TRUE)

shared %>%
  group_by(Group)%>%
  summarise(n_seqs = sum(value),
            n_sings = sum(value ==1),
            goods = 100*(1- n_sings/n_seqs)) %>%
  #filter(n_seqs >3210)%>%
  ggplot(aes(x=n_seqs, y=goods))+
  geom_point()

shared %>%
  group_by(Group)%>%
  summarise(n_seqs = sum(value))%>%
  filter(n_seqs >3210)%>%
  ungroup()%>%
  select(-numOtus, -Group)

original_data <-read.table("/Users/zhangbincai/Desktop/study/projects/sherbrook innoculation/second round 16s/analysis/stability.trim.contigs.good.unique.good.filter.unique.precluster.abund.pick.pick.opti_mcc.shared",sep = "\t", header = TRUE)
labels <- original_data$label
numOtus <- original_data$numOtus
groups <- original_data$Group 

# calculate row sums of OTU counts
otu_sums <- original_data %>%
  select(starts_with("Otu")) %>%
  rowSums()
# add row sums to data frame
original_data$otu_sums <- otu_sums

# filter out samples with less than 3210 reads
filtered_data <- original_data %>%
  filter(otu_sums >= 3210)

otu_data <- filtered_data %>%
  select(-label, -numOtus, -Group)

# Convert otu_data to a numeric matrix
otu_matrix <- as.matrix(otu_data)

otu_df <- as.data.frame(otu_matrix)

min_seqs <- 3210


# initialize a matrix to store the sum of the rarefied OTU tables
sum_rarefied_otu_tables <- matrix(0, nrow = nrow(otu_df), ncol = ncol(otu_df))

nperm <- 1000


# Perform the rarefaction and calculate the sum of the rarefied OTU tables
for (i in 1:nperm) {
  sum_rarefied_otu_tables <- sum_rarefied_otu_tables + rrarefy(otu_df, min_seqs)
}
# Calculate the average rarefied OTU table
avg_rarefied_otu_table <- sum_rarefied_otu_tables / nperm

rarefied_table_with_info <- cbind(labels, numOtus, groups, avg_rarefied_otu_table)

write.csv(rarefied_table_with_info, "rarefied_OTU.csv", row.names = FALSE)

refaried_shareshared<-read_csv("/Users/zhangbincai/rarefied_OTU.csv") %>%
  select(-label, -numOtus) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  #summarize(totall= sum(value))%>%
  #arrange(totall) %>%
  group_by(name)%>%
  mutate(total = sum(value)) %>% 
  filter(total!=0)%>%
  ungroup()%>%
  select(-total) 

rand <- shared %>%
  uncount (value) %>%
  mutate(name = sample (name)) %>%
  count (Group, name, name ="value")
  
  
