#!/usr/bin/env/Rscript 

#### Goal: preform a random forest with 17 datasets as testing data and the remaining one as the training data ####
# Do this 18 times and generate 18 ROC curves
# Rscript randomForest.R 3 > run3_out.txt

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(BSDA))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
# library(stats)
# library(tibble)
suppressPackageStartupMessages(library(ggplot2))

set.seed(100)
DEBIAS_data <- read.csv("./csv_files/debiased_lognorm.csv", colClasses = c("Phenotype" = "factor"))
lognorm_out <- read.csv("./csv_files/log_DEBIAS-M.csv", colClasses = c("Phenotype" = "factor")) # after lognorm before DEBIAS-M
print("lognorm_out before subset")
print(colnames(lognorm_out)[1:7])
print(dim(lognorm_out))
lognorm_out <- lognorm_out[,c(1:3,6:length(lognorm_out))]
print("lognorm_out after subset")
print(colnames(lognorm_out)[1:7])
print(dim(lognorm_out))
phenos <- read.csv("./csv_files/phenotypes.csv") # for naming the graphs
# phenos <- phenos[grep("autism", phenos$phenotype),] 


# IDs <- distinct(data, Study_ID)$Study_ID
# IDs <- IDs[IDs %in% phenos$ID]

#### log normalization ####

# # LOGNORM
# lognorm <- function(table)
# {
#   avg <- sum(rowSums(table))/nrow(table)
#   table <- sweep(table,1,rowSums(table),"/")
#   table <- log10(table*avg + 1)
#   return(table)
# }
# 
# lognorm_out <- lognorm(data[4:length(data)])
# lognorm_out <- add_column(lognorm_out, Study_ID=data$Study_ID, .before = colnames(lognorm_out)[1])
# lognorm_out <- add_column(lognorm_out, Phenotype=data$Phenotype, .before = colnames(lognorm_out)[1])
# lognorm_out <- add_column(lognorm_out, sample_name=data$sample_name, .before = colnames(lognorm_out)[1])
# 
# lognorm_out <- lognorm_out %>% filter(Study_ID != 13241) %>% filter(Study_ID != 14812) %>% filter(Study_ID != 11712)
# IDs <- distinct(lognorm_out, Study_ID)$Study_ID
# 
# #### Get lognorm_out ready for DEBIAS-M ####
# 
# df <- lognorm_out
# df$ID = 0
# for(i in 1:length(IDs)){
#   
#   # print(distinct(df, ID))
#   # print(i)
#   df$ID[df$Study_ID == IDs[i]] <- i-1
# }
# 
# df$case <- case_when(
#                    df$Phenotype == 1 ~ TRUE,
#                    df$Phenotype == 0 ~ FALSE,
#                    )
# df <- relocate(df, ID, .after = Study_ID)
# df <- relocate(df, case, .after = Study_ID)
# 
# write.csv(df, "./csv_files/log_DEBIAS-M.csv", row.names =FALSE)

#### permutations testing as a loop ####

# c <- DEBIAS_data[sample(nrow(DEBIAS_data), 1000), 1:200]
# c$Phenotype <- sample(c$Phenotype)
# c %>% group_by(Phenotype) %>% count(Study_ID) %>% print(n = 30)
# IDs <- distinct(c, Study_ID)$Study_ID

# lognorm
IDs <- distinct(lognorm_out, Study_ID)$Study_ID

auc.df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(auc.df) <-   c("Study_ID", "AUC", "Permutation", "DEBIAS")


args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])
permutations = as.numeric(args[2])


png(paste("./output/pre_DEBIAS-M_RF_lognorm_ROC_", IDs[i], ".png", sep=""))#, height = 24, width = 24)
# par(mar=c(3,3,1,0))#, mfrow=c(2,2))


training <- filter(lognorm_out, Study_ID != IDs[i])
training <- training[, c(2, 4:length(training))]
testing <- filter(lognorm_out, Study_ID == IDs[i])
testing <- testing[, c(2, 4:length(testing))]

# training <- filter(c, Study_ID != IDs[i])
# training <- training[, c(2, 4:length(training))]
# testing <- filter(c, Study_ID == IDs[i])
# testing %>% group_by(Phenotype) %>% count(Study_ID) %>% print(n = 30)
# testing <- testing[, c(2, 4:length(testing))]

RF_fit <- randomForest(Phenotype~., method = "class", data = training)
    
RF_pred <- predict(RF_fit, testing, type = "prob")
      
rf_roc <- roc(testing[,1], RF_pred[,1])
auc.df[nrow(auc.df) + 1,] = c(IDs[i], auc(rf_roc), FALSE, FALSE)

p <- plot(rf_roc, add = FALSE, col = "red", print.auc = TRUE)
phen <- filter(phenos, ID == IDs[i])
phen <- phen[1,2]
title(paste("Training without: ", IDs[i], " (", phen, ")",sep=""), line = + 2.5)
      
 for(j in 1:permutations){
     training$Phenotype <- sample(training$Phenotype)
     testing$Phenotype <- sample(testing$Phenotype)
     RF_fit <- randomForest(Phenotype~., method = "class", data = training)
     RF_pred <- predict(RF_fit, testing, type = "prob")
     rf_roc <- roc(testing[,1], RF_pred[,1])
     p <- plot(rf_roc, print.auc=FALSE, add = TRUE)
     auc.df[nrow(auc.df) + 1,] = c(IDs[i], auc(rf_roc), TRUE, FALSE)
     print(paste(j, " done", sep=""))
     }
p

dev.off()  

png(paste("./output/pre_DEBIAS-M_RF_lognorm_histogram_", IDs[i], ".png", sep=""))

a <- auc.df[auc.df$Permutation == 0,]$AUC
samp <- auc.df[auc.df$Permutation == 1,]$AUC
z = (a-mean(samp))/(sd(samp)/sqrt(1))
for.pval = pnorm(z, lower.tail = FALSE)

for.pval <- z.test(samp, mu = a, sigma.x = sd(samp), alternative = "less")$p.value
g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
  geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
  labs(title = paste("Training without: ", IDs[i], " (", phen, ")", sep=""), y = "count") +
  annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) 


print("done")  

print(g)

dev.off()  


print(paste(i, "of", length(IDs), " done.", IDs[i]))


if(file.exists("./csv_files/AUCs.csv")){
  all <- read.csv("./csv_files/AUCs.csv")
  all <- rbind(all, auc.df)
  write.csv(all, "./csv_files/AUCs.csv", row.names = FALSE)
}else{(write.csv(auc.df, "./csv_files/AUCs.csv", row.names = FALSE))}

write.csv(auc.df, paste("./csv_files/pre_DEBIAS_",IDs[i],"_AUCs.csv",sep=""), row.names = FALSE)

print("done writing to csv")








#### post-debias ####
IDs <- distinct(DEBIAS_data, Study_ID)$Study_ID

auc.df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(auc.df) <-   c("Study_ID", "AUC", "Permutation", "DEBIAS")


png(paste("./output/post_DEBIAS-M_RF_lognorm_ROC_", IDs[i], ".png", sep=""))#, height = 24, width = 24)


training <- filter(DEBIAS_data, Study_ID != IDs[i])
training <- training[, c(2, 4:length(training))]
testing <- filter(DEBIAS_data, Study_ID == IDs[i])
testing <- testing[, c(2, 4:length(testing))]

RF_fit <- randomForest(Phenotype~., method = "class", data = training)

RF_pred <- predict(RF_fit, testing, type = "prob")

rf_roc <- roc(testing[,1], RF_pred[,1])
auc.df[nrow(auc.df) + 1,] = c(IDs[i], auc(rf_roc), FALSE, TRUE)

p <- plot(rf_roc, add = FALSE, col = "red", print.auc = TRUE)
phen <- filter(phenos, ID == IDs[i])
phen <- phen[1,2]
title(paste("Training without: ", IDs[i], " (", phen, ")",sep=""), line = + 2.5)

for(j in 1:permutations){
  training$Phenotype <- sample(training$Phenotype)
  testing$Phenotype <- sample(testing$Phenotype)
  RF_fit <- randomForest(Phenotype~., method = "class", data = training)
  RF_pred <- predict(RF_fit, testing, type = "prob")
  rf_roc <- roc(testing[,1], RF_pred[,1])
  p <- plot(rf_roc, print.auc=FALSE, add = TRUE)
  auc.df[nrow(auc.df) + 1,] = c(IDs[i], auc(rf_roc), TRUE, TRUE)
  print(paste(j, " done post debias", sep=""))
}
p

dev.off()  

png(paste("./output/post_DEBIAS-M_RF_lognorm_histogram_", IDs[i], ".png", sep=""))

a <- auc.df[auc.df$Permutation == 0,]$AUC
samp <- auc.df[auc.df$Permutation == 1,]$AUC
z = (a-mean(samp))/(sd(samp)/sqrt(1))
for.pval = pnorm(z, lower.tail = FALSE)

g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
  geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
  labs(title = paste("Training without: ", IDs[i], " (", phen, ")", sep=""), y = "count") +
  annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) 


print("done")  

print(g)

dev.off()  


print(paste(i, "of", length(IDs), " done.", IDs[i]))


if(file.exists("./csv_files/AUCs.csv")){
  all <- read.csv("./csv_files/AUCs.csv")
  all <- rbind(all, auc.df)
  write.csv(all, "./csv_files/AUCs.csv", row.names = FALSE)
}else{(write.csv(auc.df, "./csv_files/AUCs.csv", row.names = FALSE))}

write.csv(auc.df, paste("./csv_files/post_DEBIAS_",IDs[i],"_AUCs.csv",sep=""), row.names = FALSE)

print("done writing to csv")


