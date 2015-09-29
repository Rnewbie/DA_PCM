library(readxl)
Diene <- read_excel("Diene.xlsx")
Dienophile <- read_excel("Dienophile.xlsx")
Product <- Diene$Product

Diene <- Diene[, c("MeanAbs", "RBN","ALOGP", "nCIC", "nHDon", "nHAcc", "TPSA(Tot)", "Energy", "Dipole", "MW",  "HOMO", "LUMO", "GAP")]
Dienophile <- Dienophile[, c("MeanAbs","ALOGP", "RBN", "nCIC", "nHDon", "nHAcc", "TPSA(Tot)", "Energy", "Dipole", "MW",  "HOMO", "LUMO", "GAP")]

results_Diene <- data.frame(t(sapply(Diene, function(cl) list (Mean_of_Diene = round(mean(cl, na.rm = TRUE), digits = 3),
                                                                   SD_Diene = round(sd(cl, na.rm  = TRUE), digits = 3)))))
results_Dienophile <- data.frame(t(sapply(Dienophile, function(cl) list(Mean_of_Dienophile= round(mean(cl, na.rm = TRUE), digits = 3),
                                                                    SD_Dienophile = round(sd(cl, na.rm = TRUE), digits = 3)))))
results <- cbind(results_Diene, results_Dienophile)
Descriptor_Name <- row.names(results)
df <- cbind(Descriptor_Name, results)
my.df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
write.csv(my.df, file = "EDA_Substructure.csv", row.names = FALSE)
print(df, row.names = FALSE)

####latext table
library(xtable)
table <- xtable(results, row.names = FALSE)
print(table)

normality_test_Diene <- t(sapply(Diene, function(x)  shapiro.test(x)[c("statistic", "p.value",
                                                                        "method")]))

df <- normality_test_Diene[, 1:2]
df <- round(df, digits = 3)
my.df$significant_level <- ifelse(my.df$p.value < 0.05, "Significant", "Unsignificant")
label <- row.names(df)
normality <- cbind(label, my.df)
table <- xtable(normality, row.names = FALSE)
print(table)

### plot 
library(readxl)
Diene <- read_excel("Diene.xlsx")
Dienophile <- read_excel("Dienophile.xlsx")
Product <- Diene$Product
Product <- as.factor(Product)

data <- rbind(Diene, Dienophile)
data <- cbind(Product, data)

Mann_Whitney <- t(sapply(data[-1], function(x) unlist(pairwise.wilcox.test(x, data$Product, 
                                                                           alt = "two.sided",
                                                                           conf.int = TRUE,
                                                                           conf.level = 0.95,
                                                                           exact = FALSE)[c("statistic.W", "p.value", 
                                                                                            "method", "conf.int",
                                                                                            "conf.level")])))

colnames(data)[7] <-  "TPSA"
top3 <- data[, c("Product", "HOMO", "LUMO", "GAP")]
fit <- randomForest(Product~., data = data)
J48 <- J48(Product~., data = data)




a <- ggplot(data, aes(x = HOMO, fill = Product)) +
  geom_histogram(binwidth= 1, alpha = 0.5, colour = "black", position = "dodge") +
  theme(
    legend.position=("none"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  #labs(y = expression( Predicted~(pIC[50])~ activity)) +
  #labs(x = expression( Exprimental~(pIC[50])~ activity))
  labs(y =  "Frequency") +
coord_cartesian(ylim = c(0, 120), xlim = c(-2, 2))

b <- ggplot(data, aes(x = LUMO, fill = Product)) +
  geom_histogram(binwidth= 1, alpha = 0.5, colour = "black", position = "dodge") +
  theme(
    legend.position=("none"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) + 
  #labs(y = expression( Predicted~(pIC[50])~ activity)) +
  #labs(x = expression( Exprimental~(pIC[50])~ activity))
  labs(y =  "Frequency") +
coord_cartesian(ylim = c(0, 120), xlim = c(-2, 2))

c <- ggplot(data, aes(x = GAP, fill = Product)) +
  geom_histogram(binwidth= 1, alpha = 0.5, colour = "black", position = "dodge") +
  theme(
    legend.position=("none"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  #labs(y = expression( Predicted~(pIC[50])~ activity)) +
  #labs(x = expression( Exprimental~(pIC[50])~ activity))
  labs(y =  "Frequency") +
  coord_cartesian(ylim = c(0, 120), xlim = c(-2, 2))
library(cowplot)
plot_grid(a, b, c, labels = c(" ", " ", " "), ncol =3, label_size = 20)
### Diene VS Dienophile
Diene$Label <- "Diene"
Dienophile$Label <- "Dienophile"
data <- rbind(Diene, Dienophile)
d <- ggplot(data, aes(x = HOMO, fill = Label)) +
  geom_histogram(binwidth= 1, alpha = 0.5, colour = "black", position = "dodge") +
  theme(
    legend.position=("none"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  #labs(y = expression( Predicted~(pIC[50])~ activity)) +
  #labs(x = expression( Exprimental~(pIC[50])~ activity))
  labs(y =  "Frequency") +
  coord_cartesian(ylim = c(0, 120), xlim = c(-2, 2))


e <- ggplot(data, aes(x = LUMO, fill = Label)) +
  geom_histogram(binwidth= 1, alpha = 0.5, colour = "black", position = "dodge") +
  theme(
    legend.position=("none"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  #labs(y = expression( Predicted~(pIC[50])~ activity)) +
  #labs(x = expression( Exprimental~(pIC[50])~ activity))
  labs(y =  "Frequency") +
  coord_cartesian(ylim = c(0, 120), xlim = c(-2, 2))

f <- ggplot(data, aes(x = GAP, fill = Label)) +
  geom_histogram(binwidth= 1, alpha = 0.5, colour = "black", position = "dodge") +
  theme(
    legend.position=("none"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  #labs(y = expression( Predicted~(pIC[50])~ activity)) +
  #labs(x = expression( Exprimental~(pIC[50])~ activity))
  labs(y =  "Frequency") +
  coord_cartesian(ylim = c(0, 120), xlim = c(-2, 2))

### preparing the results
Diene <- Diene[, -nearZeroVar(Diene)]
Dienophile <- Dienophile[, -nearZeroVar(Dienophile)]
DieneXDienophile <- getCPI(Diene, Dienophile, type = "tensorprod")
DieneXDienophile <- as.data.frame(DieneXDienophile)

Label 