library(vegan)
library(ggplot2)

# Read in your OTU table and metadata
otu <- read.csv("rarefied.D3.OTU.csv", row.names = 1)
meta <- read.csv("D3.csv")

# Calculate alpha diversity indices
alpha_div <- data.frame(
  #Chao1 = estimateR(otu),
  Shannon = diversity(otu, index = "shannon"),
  Simpson = diversity(otu, index = "simpson")
)

# Combine alpha diversity with metadata
alpha_div_meta <- cbind(meta, alpha_div)

# Calculate p-values using Mann-Whitney test
#p_chao1 <- wilcox.test(Chao1 ~ Treatment, data = alpha_div_meta)
p_shannon <- wilcox.test(Shannon ~ Treatment, data = alpha_div_meta)
p_simpson <- wilcox.test(Simpson ~ Treatment, data = alpha_div_meta)



# Create box plots for each index
#p1 <- ggplot(alpha_div_meta, aes(x = Treatment, y = Chao1)) +
#geom_boxplot(aes(fill = Treatment)) +
#geom_jitter(aes(color = Treatment), width = 0.2) +
#scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
#scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
#labs(title = paste("Chao1 Index (p-value =", round(p_chao1$p.value, digits = 3), ")"))

p2 <- ggplot(alpha_div_meta, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot(alpha=0.4) +
  geom_jitter(show.legend = FALSE, color = 'black', shape = 21, width = 0.2) +
  scale_fill_manual(values = c("#3399ff", "#cc0066")) +
  #scale_color_manual(values = c("#6699cc", "#ff99ff", "#619CFF")) +
  labs(title = paste("Shannon Index (p-value =", round(p_shannon$p.value, digits = 3), ")"))


p3 <- ggplot(alpha_div_meta, aes(x = Treatment, y = Simpson, fill = Treatment)) +
  geom_boxplot(alpha=0.4) +
  geom_jitter(show.legend = FALSE, color = 'black', shape = 21, width = 0.2) +
  scale_fill_manual(values = c("#3399ff", "#cc0066")) +
  labs(title = paste("Simpson Index (p-value =", round(p_simpson$p.value, digits = 3), ")"))

# Print p-values
#cat("Chao1 p-value:", p_chao1$p.value, "\n")
cat("Shannon p-value:", p_shannon$p.value, "\n")
cat("Simpson p-value:", p_simpson$p.value)

