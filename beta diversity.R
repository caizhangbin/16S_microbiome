library(vegan)
library(tidyverse)
library(ggplot2)

metadata <- read.csv("/Users/zhangbincai/Desktop/study/projects/sherbrook_innoculation/second_round_16s/analysis/D7.csv")
shared <- read.csv("/Users/zhangbincai/Desktop/study/projects/sherbrook_innoculation/second_round_16s/analysis/rarefied.D7.OTU.csv") %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  #summarize(totall= sum(value))%>%
  #arrange(totall) %>%
  group_by(name)%>%
  mutate(total = sum(value)) %>% 
  filter(total!=0)%>%
  ungroup()%>%
  select(-total) %>%
  pivot_wider(Group)


rownames(shared) <-shared$Group

group_names <- shared$Group
shared <- shared %>% select(-Group)
shared <- as.matrix(shared)
rownames(shared) <- group_names

shared <- shared[,-1]
shared <-as.matrix(shared)

set.seed(19910712)
dist <- vegdist(shared, method = "bray")

set.seed(1)
nmds <- metaMDS(dist) %>%
  scores() %>%
  as_tibble(rownames="Group")



nmds <- metaMDS(dist) %>%
  scores() %>%
  `rownames<-`(., value = rownames(shared)) %>%
  as_tibble(rownames = "Group")

metadata_ndmds <- inner_join(metadata, nmds)

adonis2 (dist ~ metadata $ Treatment)

# Run the adonis2 function
adonis2_result <- adonis2(dist ~ metadata$Treatment)

str(adonis2_result)

# Extract the p-value
p_value <- adonis2_result$`Pr(>F)`[1]

# Create the plot title
plot_title <- paste("Adonis2 p-value:", round(p_value, 3))

# Create the plot
metadata_ndmds %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(fill = Treatment), shape = 21, color = "black")+
  stat_ellipse(geom = "polygon", aes(x0 = mean(NMDS1), y0 = mean(NMDS2), a = sd(NMDS1)*2, b = sd(NMDS2)*2, angle = 45, fill = Treatment), alpha = 0.5, scale = 2, color = NA) +
  scale_fill_manual(values = c("#3399ff", "#cc0066")) +
  ggtitle(plot_title)




