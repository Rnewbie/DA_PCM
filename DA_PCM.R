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

### preparing for the PCM dataset with three descriptors
diene <- read_excel("Diene.xlsx")
dienophile <- read_excel("Dienophile.xlsx")
## subset for diene
diene_df <- diene[, 4:ncol(diene)]
#diene_df <- diene_df[ , -nearZeroVar(diene_df)]
descriptors <- c("HOMO", "LUMO", "GAP")
diene_df <- diene_df[, descriptors]
dienophile_df <- dienophile[, 4:ncol(dienophile)]
dienophile_df <- dienophile_df[, descriptors]
Product <- diene$Product
Product <- as.factor(Product)
dieneXdienophile <- getCPI(diene_df, dienophile_df, type = "tensorprod")
dieneXdienophile <- as.data.frame(dieneXdienophile)
dfDiene <- names(data.frame(diene_df[, 1:3]))
dfDienophile <- names(data.frame(dienophile_df[, 1:3]))
dieneNamecross <- rep(dfDiene, each = 3)
dienophileNamecross <- rep(dfDienophile, times = 3)
label <- paste(dieneNamecross, dienophileNamecross, sep = "_")
colnames(dieneXdienophile) <- label
dieneXdienophile <- as.data.frame(dieneXdienophile)
### diene self scorss terms. 
dieneXdiene <- getCPI(diene_df, diene_df, type = "tensorprod")
dieneName2 <- rep(dfDiene, times = 3)
dieneName1 <- rep(dfDiene, each = 3)
label_diene <- paste(dieneName1, dieneName2, sep = "_")
colnames(dieneXdiene) <- label_diene
dieneXdiene <- as.data.frame(dieneXdiene)
index <- seq(1, 9, by =  4)
dieneselfcross <- dieneXdiene[, -index]
transposedIndexed_diene <- t(dieneselfcross)
index1 <- which(duplicated(transposedIndexed_diene))
removed_duplicated_diene <- transposedIndexed_diene[-index1, ]
dieneXdiene <- t(removed_duplicated_diene)
dieneXdiene <- as.data.frame(dieneXdiene)

dienophileXdienophile <- getCPI(dienophile_df, dienophile_df, type = "tensorprod")
dienophileName2 <- rep(dfDienophile, times = 3)
dienophileName1 <- rep(dfDienophile, each = 3)
label <- paste(dienophileName1, dienophileName2, sep = "_")
colnames(dienophileXdienophile) <- label
dienophileXdienophile <- as.data.frame(dienophileXdienophile)
index <- seq(1, 9, by =  4)
dienophileselfcross <- dienophileXdienophile[, -index]
transposedIndexed_dienophile <- t(dienophileselfcross)
index1 <- which(duplicated(transposedIndexed_dienophile))
removed_duplicated_dienophile <- transposedIndexed_dienophile[-index1, ]
dienophileXdienophile <- t(removed_duplicated_dienophile)
dienophileXdienophile <- as.data.frame(dienophileXdienophile)

diene <- diene_df
dienophile <- dienophile_df
diene_dienophile <- cbind(diene, dienophile)
diene_dienophile_dieneXdienophile_block_scale <- cbind(diene, dienophile,
                                           dieneXdienophile) * (1/sqrt(length(diene)+
                                                                         length(dienophile)+
                                                                         length(dieneXdienophile)))
diene_dienophile_dieneXdiene_block_scale <- cbind(diene, dienophile,
                                                  dieneXdiene) * (1/sqrt(length(diene)+
                                                                           length(dienophile)+
                                                                           length(dieneXdiene)))
diene_dienophile_dienophileXdienophile_block_scale <- cbind(diene, dienophile,
                                                            dienophileXdienophile) * (1/sqrt(length(diene)+
                                                                                               length(dienophile)+
                                                                                               length(dienophileXdienophile)))
diene_dienophile_dieneXdienophile_dieneXdiene_block_scale <- cbind(diene,
                                                                   dienophile,
                                                                   dieneXdienophile,
                                                                   dieneXdiene) * (1/sqrt(length(diene)+
                                                                                            length(dienophile)+
                                                                                            length(dieneXdienophile)+
                                                                                            length(dieneXdiene)))
diene_dienophile_dieneXdienophile_dienophileXdienophile_block_scale <- cbind(diene,
                                                                              dienophile,
                                                                              dieneXdienophile,
                                                                             dienophileXdienophile) * (1/sqrt(length(diene)+
                                                                                                                length(dienophile)+
                                                                                                                length(dieneXdienophile)+
                                                                                                                length(dienophileXdienophile)))
diene_dienophile_dieneXdiene_dienophileXdienophile_block_scale <- cbind(diene,
                                                                        dienophile,
                                                                        dieneXdiene,
                                                                        dienophileXdienophile) * (1/sqrt(length(diene)+
                                                                                                           length(dienophile)+
                                                                                                           length(dieneXdiene)+
                                                                                                           length(dienophileXdienophile)))
diene_dienophile_dieneXdienophile_dieneXdiene_dienophileXdienophile_block_scale <- cbind(diene,
                                                                                         dienophile,
                                                                                         dieneXdienophile,
                                                                                         dieneXdiene,
                                                                                         dienophileXdienophile) * (1/sqrt(length(diene)+
                                                                                                                            length(dienophile)+
                                                                                                                            length(dieneXdienophile)+
                                                                                                                            length(dieneXdiene)+
                                                                                                                            length(dienophileXdienophile)))
diene <- cbind(Product, diene)
dienophile <- cbind(Product,dienophile)
dieneXdienophile <- cbind(Product, dieneXdienophile)
dieneXdiene <- cbind(Product, dieneXdiene)
dienophileXdienophile <- cbind(Product, dienophileXdienophile)
diene_dienophile <- cbind(Product, diene_dienophile)
diene_dienophile_dieneXdienophile <- cbind(Product, diene_dienophile_dieneXdienophile_block_scale)
diene_dienophile_dieneXdiene <- cbind(Product, diene_dienophile_dieneXdiene_block_scale)
diene_dienophile_dienophileXdienophile <- cbind(Product, diene_dienophile_dienophileXdienophile_block_scale)
diene_dienophile_dieneXdienophile_dieneXdiene <- cbind(Product, 
      diene_dienophile_dieneXdienophile_dieneXdiene_block_scale)
diene_dienophile_dieneXdienophile_dienophileXdienophile <- cbind(Product, 
      diene_dienophile_dieneXdienophile_dienophileXdienophile_block_scale)
diene_dienophile_dieneXdiene_dienophileXdienophile <- cbind(Product, 
                                                            diene_dienophile_dieneXdiene_dienophileXdienophile_block_scale)
diene_dienophile_dieneXdienophile_dieneXdiene_dienophileXdienophile <- cbind(Product, 
      diene_dienophile_dieneXdienophile_dieneXdiene_dienophileXdienophile_block_scale)



#### training results using J48
J48_training <- function(x, Product){
  if (Product == "Meta") {
    library(parallel)
    library(doSNOW)
    cl <- makeCluster(2)
    registerDoSNOW(cl)
    
  ok <- list(100)
  ok <- foreach(i = 1:100) %dopar% { 
    in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
    train <- x[in_train, ]
    test <- x[-in_train, ]
    model_train <- RWeka::J48(Product~., data = train)
    actual <- train$Product
    prediction <- predict(model_train, train)
    results <- caret::confusionMatrix(prediction, actual)
    results <- results$table
    results <- table(prediction, actual)
    results <- as.numeric(results)
    ok[[i]] <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
    #Ortho <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
    #Para <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
  }
}  else if (Product == "Ortho") {
    cl <- makeCluster(2)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      model_train <- RWeka::J48(Product~., data = train)
      actual <- train$Product
      prediction <- predict(model_train, train)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      #Meta <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
      ok[[i]] <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
      #Para <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
    } 
}  else if (Product == "Para") {
  cl <- makeCluster(2)
  registerDoSNOW(cl)
  
  ok <- list(100)
  ok <- foreach(i = 1:100) %dopar% { 
    in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
        train <- x[in_train, ]
        test <- x[-in_train, ]
        model_train <- RWeka::J48(Product~., data = train)
        actual <- train$Product
        prediction <- predict(model_train, train)
        results <- caret::confusionMatrix(prediction, actual)
        results <- results$table
        results <- table(prediction, actual)
        results <- as.numeric(results)
        #Meta <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
        #Ortho <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
        ok[[i]] <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
      }
  return(ok)
  stopCluster(cl)
} }

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}


results_training_Meta <- function(x) {
  yes <- J48_training(x, Product = "Meta")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}
  
results_training_Ortho <- function(x) {
  yes <- J48_training(x, Product = "Ortho")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}


results_training_Para <- function(x) {
  yes <- J48_training(x, Product = "Para")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                         "MCC_Mean", "MCC_SD")
  return(results_all)
}

J48_training_all <- function(x) {
  meta <- results_training_Meta(x)
  ortho <- results_training_Ortho(x)
  para <- results_training_Para(x)
  result_all <- cbind(meta, ortho, para)
  total <- apply(result_all, 1, mean)
  result_all_mean <- cbind(result_all, total)
  colnames(result_all_mean) <- c("Meta", "Ortho", "Para", "Overall")
  return(result_all_mean)
}

input <- list(diene = diene,
              dienophile = dienophile,
              dieneXdienophile = dieneXdienophile,
              dieneXdiene = dieneXdiene,
              dienophileXdienophile = dienophileXdienophile,
              diene_dienophile = diene_dienophile,
              diene_dienophile_dieneXdienophile = diene_dienophile_dieneXdienophile,
              diene_dienophile_dieneXdiene = diene_dienophile_dieneXdiene,
              diene_dienophile_dienophileXdienophile = diene_dienophile_dienophileXdienophile,
              diene_dienophile_dieneXdienophile_dieneXdiene = diene_dienophile_dieneXdienophile_dieneXdiene,
              diene_dienophile_dieneXdienophile_dienophileXdienophile = diene_dienophile_dieneXdienophile_dienophileXdienophile,
              diene_dienophile_dieneXdiene_dienophileXdienophile = diene_dienophile_dieneXdiene_dienophileXdienophile,
              diene_dienophile_dieneXdienophile_dieneXdiene_dienophileXdienophile = diene_dienophile_dieneXdienophile_dieneXdiene_dienophileXdienophile)


#### 10fold  fold cross validation

#### training results using J48
J48_10_CV <- function(x, Product){
  if (Product == "Meta") {
    library(parallel)
    library(doSNOW)
    cl <- makeCluster(2)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
      myData <- x[in_train, ]
      test <- x[-in_train, ]
      k = 10
      index <- sample(1:k, nrow(myData), replace = TRUE)
      folds <- 1:k
      myRes <- data.frame()
      for (j in 1:k)
        training <- subset(myData, index %in% folds[-j])
      testing <- subset(myData, index %in% c(j))
      model_train <- RWeka::J48(Product~., data = training)
      actual <- testing$Product
      prediction <- predict(model_train, testing)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      results <- rbind(myRes, results)
      ok[[i]] <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
      #Ortho <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
      #Para <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
    }
  }  else if (Product == "Ortho") {
    cl <- makeCluster(2)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
      myData <- x[in_train, ]
      test <- x[-in_train, ]
      k = 10
      index <- sample(1:k, nrow(myData), replace = TRUE)
      folds <- 1:k
      myRes <- data.frame()
      for (j in 1:k)
        training <- subset(myData, index %in% folds[-j])
      testing <- subset(myData, index %in% c(j))
      model_train <- RWeka::J48(Product~., data = training)
      actual <- testing$Product
      prediction <- predict(model_train, testing)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      results <- rbind(myRes, results)
      #Meta <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
      ok[[i]] <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
      #Para <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
    } 
  }  else if (Product == "Para") {
    cl <- makeCluster(2)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
      myData <- x[in_train, ]
      test <- x[-in_train, ]
      k = 10
      index <- sample(1:k, nrow(myData), replace = TRUE)
      folds <- 1:k
      myRes <- data.frame()
      for (j in 1:k)
        training <- subset(myData, index %in% folds[-j])
      testing <- subset(myData, index %in% c(j))
      model_train <- RWeka::J48(Product~., data = training)
      actual <- testing$Product
      prediction <- predict(model_train, testing)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      results <- rbind(myRes, results)
      #Meta <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
      #Ortho <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
      ok[[i]] <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
    }
    return(ok)
    stopCluster(cl)
  } }

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}


results_CV_Meta <- function(x) {
  yes <- J48_10_CV(x, Product = "Meta")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

results_CV_Ortho <- function(x) {
  yes <- J48_10_CV(x, Product = "Ortho")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}


results_CV_Para <- function(x) {
  yes <- J48_10_CV(x, Product = "Para")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

J48_CV_all <- function(x) {
  meta <- results_CV_Meta(x)
  ortho <- results_CV_Ortho(x)
  para <- results_CV_Para(x)
  result_all <- cbind(meta, ortho, para)
  total <- apply(result_all, 1, mean)
  result_all_mean <- cbind(result_all, total)
  colnames(result_all_mean) <- c("Meta", "Ortho", "Para", "Overall")
  return(result_all_mean)
}




#### training results using J48
J48_testing <- function(x, Product){
  if (Product == "Meta") {
    library(parallel)
    library(doSNOW)
    cl <- makeCluster(2)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      model_train <- RWeka::J48(Product~., data = train)
      actual <- test$Product
      prediction <- predict(model_train, test)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      ok[[i]] <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
      #Ortho <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
      #Para <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
    }
  }  else if (Product == "Ortho") {
    cl <- makeCluster(2)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      model_train <- RWeka::J48(Product~., data = train)
      actual <- test$Product
      prediction <- predict(model_train, test)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      #Meta <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
      ok[[i]] <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
      #Para <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
    } 
  }  else if (Product == "Para") {
    cl <- makeCluster(2)
    registerDoSNOW(cl)
    
    ok <- list(100)
    ok <- foreach(i = 1:100) %dopar% { 
      in_train <- caret::createDataPartition(x$Product, p = 0.80, list = FALSE)
      train <- x[in_train, ]
      test <- x[-in_train, ]
      model_train <- RWeka::J48(Product~., data = train)
      actual <- test$Product
      prediction <- predict(model_train, test)
      results <- caret::confusionMatrix(prediction, actual)
      results <- results$table
      results <- table(prediction, actual)
      results <- as.numeric(results)
      #Meta <- cbind(results[[1]], (results[[2]] + results[[3]]), (results[[4]] + results[[7]]), (results[[5]] + results[[9]]))
      #Ortho <- cbind(results[5], (results[2] + results[8]), (results[4] + results[6]), (results[1] + results[9]))
      ok[[i]] <- cbind(results[9], (results[3] + results[6]), (results[4] + results[6]), (results[1] + results[5]))
    }
    return(ok)
    stopCluster(cl)
  } }

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}


results_testing_Meta <- function(x) {
  yes <- J48_testing(x, Product = "Meta")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

results_testing_Ortho <- function(x) {
  yes <- J48_testing(x, Product = "Ortho")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}


results_testing_Para <- function(x) {
  yes <- J48_testing(x, Product = "Para")
  great <- data.frame(yes)
  TP <- seq(from = 1, to = 400, by = 4)
  FN <- seq(from = 2, to = 400, by = 4)
  FP <- seq(from = 3, to = 400, by = 4)
  TN <- seq(from = 4, to = 400, by = 4)
  results <- mapply(c, great[TP], great[FN], great[FP], great[TN])
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}

J48_testing_all <- function(x) {
  meta <- results_testing_Meta(x)
  ortho <- results_testing_Ortho(x)
  para <- results_testing_Para(x)
  result_all <- cbind(meta, ortho, para)
  total <- apply(result_all, 1, mean)
  result_all_mean <- cbind(result_all, total)
  colnames(result_all_mean) <- c("Meta", "Ortho", "Para", "Overall")
  return(result_all_mean)
}

### training
results_J48_training <- lapply(input, function(x) {
  models <- J48_training_all(x)
  return(models)
})
### 10 CV
results_J48_CV <- lapply(input, function(x) {
  models <- J48_CV_all(x)
  return(models)
})
### testing
results_J48_testing <- lapply(input, function(x) {
  models <- xJ48_testing_all(x)
  return(models)
})



