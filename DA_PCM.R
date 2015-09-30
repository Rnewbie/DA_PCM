library(readxl)
library(caret)
library(paran)
library(cowplot)

setwd("/Volumes/SAW SIMEON/DA_PCM")
diene <- read_excel("Diene.xlsx")
dienophile <- read_excel("Dienophile.xlsx")
### pca analysis for diene
diene_df <- diene[, 4:ncol(diene)]
#diene_df <- diene_df[ , -nearZeroVar(diene_df)]
descriptors <- c("MeanAbs", "RBN", "nCIC", "nHDon", "nHAcc", "TPSA(Tot)", "Energy", "Dipole",
                 "MW", "HOMO", "LUMO", "GAP", "ALOGP")
diene_df <- diene_df[, descriptors]
plot <- list(13)
for (i in descriptors) {
  data <- diene_df[, i]
  plot <- ggplot(data = diene_df, aes(data)) + geom_histogram(col = "black", fill = "blue") + xlab(i) +
    ylab("Frequency") + geom_vline(data = diene_df, aes(xintercept = mean(data), na.rm = TRUE),
                                   colour = "black", linetype = "dashed", size = 1) +
    theme(
      axis.text.x = element_text(colour = "black", face = "bold"),
      axis.text.y = element_text(colour = "black", face = "bold"),
      panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)
    )
  print(plot)
}





paran(diene_df)
pca <- prcomp(diene_df, center = TRUE, scale. = TRUE, retx = TRUE)
summary(pca)
scores <- pca$x[,1:5]
loadings <- pca$rotation[,1:5]
km <- kmeans(scores, center=3, nstart=5)
ggdata <- data.frame(scores, Cluster=km$cluster)
#ggdata <- cbind(compoundname, ggdata)
### paper numbering
library(grid)
set.seed(23)
x <- ggplot(ggdata, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(fill=factor(Cluster)), size=5, shape=20, pch = 21, alpha = 0.8) +
  ggtitle("") +
  stat_ellipse(aes(fill=factor(Cluster)), colour = "black", 
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  #geom_text(aes(label=compoundnumber), size=7, hjust=0.5, vjust= 1.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    #plot.title = element_text(size=20, face="bold", colour="black", vjust = -1, hjust=-0.21),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) +
  coord_cartesian(ylim = c(-5, 5), xlim = c(-10, 10))
x

km <- kmeans(loadings, center=2, nstart=5)
ggdata <- data.frame(loadings, Cluster=km$cluster)
### paper numbering
library(grid)
set.seed(23)
a <- ggplot(ggdata, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(fill=factor(Cluster)), size=5, shape=20, pch = 21, alpha = 0.8) +
  ggtitle("") +
  stat_ellipse(aes(fill=factor(Cluster)), colour = "black", 
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  #geom_text(aes(label=labelcompound), size=7, hjust=0.5, vjust= 1.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    #plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) +
  coord_cartesian(ylim = c(-1, 1), xlim = c(-1, 1))


### pca analysis for dienophile
dienophile_df <- dienophile[, 4:ncol(dienophile)]
descriptors <- c("MeanAbs", "RBN", "nCIC", "nHDon", "nHAcc", "TPSA(Tot)", "Energy", "Dipole",
                 "MW", "HOMO", "LUMO", "GAP", "ALOGP")
dienophile_df <- dienophile_df[, descriptors]
pca_dienophile <- prcomp(dienophile_df, center = TRUE, scale. = TRUE, retx = TRUE)
paran(dienophile_df)
summary(pca_dienophile)
scores_dienophile <- pca_dienophile$x[,1:5]
loadings_dienophile <- pca_dienophile$rotation[,1:5]

km <- kmeans(scores_dienophile, center=3, nstart=5)
ggdata <- data.frame(scores_dienophile, Cluster=km$cluster)
### paper numbering
library(grid)
set.seed(2300)
y <- ggplot(ggdata, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(fill=factor(Cluster)), size=5, shape=20, pch = 21, alpha = 0.8) +
  ggtitle(" ") +
  stat_ellipse(aes(fill=factor(Cluster)), colour = "black", 
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  #geom_text(aes(label=proteinname), size=7, hjust=0.5, vjust= 1.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    # plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) +
  coord_cartesian(ylim = c(-5, 5), xlim = c(-4, 4))
set.seed(200)
km <- kmeans(loadings_dienophile, center=2, nstart=5)
ggdata <- data.frame(loadings_dienophile, Cluster=km$cluster)
labelprotein <- rownames(loadings)
ggdata <- cbind(labelprotein, ggdata)
### paper numbering
library(grid)
set.seed(23)
b <- ggplot(ggdata, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(fill=factor(Cluster)), size=5, shape=20, pch = 21, alpha = 0.8) +
  ggtitle(" ") +
  stat_ellipse(aes(fill=factor(Cluster)), colour = "black", 
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  #geom_text(aes(label=labelprotein), size=7, hjust=0.5, vjust= 1.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    #plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) +
  coord_cartesian(ylim = c(-0.8, 1), xlim = c(-0.8, 0.8))

plot_grid(x, y,a, b,  labels = c("(a)", "(b)", "(c)", "(d)"), ncol = 2, label_size = 20, scale = 1)


