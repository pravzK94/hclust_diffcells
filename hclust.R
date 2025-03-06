set.seed(1234)
getwd()

library(readxl)
library(dplyr)
library(dendextend)

#Loading the data
data_3 <- as.data.frame(read.csv(file = "data_3.csv", sep = ",", header= T, row.names=1))
data_2 <- read_excel("data_2.xlsx")

#Filtering the cells with maturation stage > 10 and differentiated genes
ms10 <- as.numeric(data_3[3, ]) > 10
genes <- data_2$`Fig. 1 (Differentiation and layer identity genes)`

filter_data <- data_3[, ms10]
datasub <- as.data.frame(t(filter_data))
datasub <- datasub %>% select(all_of(genes))
#Saving the file in CSV
write.csv(datasub, file = "subset.csv", row.names = T)

#'hclust' function on the cells (in rows)
cells_dmatrix <- dist(datasub, method = "euclidean")
cells_clust <- hclust(cells_dmatrix, method = "ward.D2")

# Plotting parameters with Cell-Infomap cluster colored labels
dend_cells <- as.dendrogram(cells_clust)
labels(dend_cells) <- rownames(datasub)
infomap_labels <- as.numeric(as.factor(as.character(filter_data[4, ])))
dend_cells <- color_labels(dend_cells, k = length(unique(infomap_labels)))

#Plotting and saving in pdf
pdf("hclust_cells.pdf", width = 30, height = 15)
plot(dend_cells, type = "rectangle", cex.main = 4,
     main = "Hierarchical clustering of cells with Maturation stage above 10")
dev.off()