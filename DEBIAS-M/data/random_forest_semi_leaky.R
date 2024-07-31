#!/usr/bin/env/Rscript 

#### Goal: preform a random forest with 17 datasets as testing data and the remaining one as the training data ####
# Do this 15 times and generate 15 ROC curves
# Rscript randomForest.R 3 > run3_out.txt

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(BSDA))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
# library(stats)
# library(tibble)
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
file = args[1]
permutations = as.numeric(args[2])
ID = as.numeric(args[3])

set.seed(100)
DEBIAS_data <- read.csv(file, colClasses = c("Phenotype" = "factor"))


phenos <- read.csv("./csv_files/phenotypes.csv") # for naming the graphs


auc.df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(auc.df) <-   c("Study_ID", "AUC", "Permutation", "DEBIAS")


png(paste("./output/post_DEBIAS-M_RF_lognorm_ROC_", ID, ".png", sep=""))#, height = 24, width = 24)


training <- filter(DEBIAS_data, Study_ID != ID)
training <- training[, c(2, 4:length(training))]
testing <- filter(DEBIAS_data, Study_ID == ID)
testing <- testing[, c(2, 4:length(testing))]

RF_fit <- randomForest(Phenotype~., method = "class", data = training)

RF_pred <- predict(RF_fit, testing, type = "prob")

rf_roc <- roc(testing[,1], RF_pred[,1])
auc.df[nrow(auc.df) + 1,] = c(ID, auc(rf_roc), FALSE, TRUE)

p <- plot(rf_roc, add = FALSE, col = "red", print.auc = TRUE)
phen <- filter(phenos, ID == ID)
phen <- phen[1,2]
title(paste("Training without: ", ID, " (", phen, ")",sep=""), line = + 2.5)

for(j in 1:permutations){
  training$Phenotype <- sample(training$Phenotype)
  testing$Phenotype <- sample(testing$Phenotype)
  RF_fit <- randomForest(Phenotype~., method = "class", data = training)
  RF_pred <- predict(RF_fit, testing, type = "prob")
  rf_roc <- roc(testing[,1], RF_pred[,1])
  p <- plot(rf_roc, print.auc=FALSE, add = TRUE)
  auc.df[nrow(auc.df) + 1,] = c(ID, auc(rf_roc), TRUE, TRUE)
  print(paste(j, " done post debias", sep=""))
}
p

dev.off()  

png(paste("./output/post_DEBIAS-M_RF_lognorm_histogram_", ID, ".png", sep=""))

a <- auc.df[auc.df$Permutation == 0,]$AUC
samp <- auc.df[auc.df$Permutation == 1,]$AUC
z = (a-mean(samp))/(sd(samp)/sqrt(1))
for.pval = pnorm(z, lower.tail = FALSE)

g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
  geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
  labs(title = paste("Training without: ", ID, " (", phen, ")", sep=""), y = "count") +
  annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) 


print("done")  

print(g)

dev.off()  


print(paste(ID, " done."))


if(file.exists("./csv_files/AUCs.csv")){
  all <- read.csv("./csv_files/AUCs.csv")
  all <- rbind(all, auc.df)
  write.csv(all, "./csv_files/AUCs.csv", row.names = FALSE)
}else{(write.csv(auc.df, "./csv_files/semi_leaky_AUCs.csv", row.names = FALSE))}

write.csv(auc.df, paste("./csv_files/post_DEBIAS_",ID,"_AUCs.csv",sep=""), row.names = FALSE)

print("done writing to csv")


