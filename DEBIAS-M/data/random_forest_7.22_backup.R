#### Goal: preform a random forest with 17 datasets as testing data and the remaining one as the training data ####
# Do this 18 times and generate 18 ROC curves
library(dplyr)
library(randomForest)
library(pROC)
library(stats)
library(tibble)
library(ggplot2)
library(patchwork)
library(gridGraphics)
library(cowplot)
set.seed(100)
data <- read.csv("./csv_files/DEBIAS-M.csv", colClasses = c("Phenotype" = "factor")) # before lognorm and DEBIAS-M
DEBIAS_data <- read.csv("./csv_files/debiased_lognorm.csv", colClasses = c("Phenotype" = "factor"))
phenos <- read.csv("./csv_files/phenotypes.csv") # for naming the graphs
# phenos <- phenos[grep("autism", phenos$phenotype),] 


IDs <- distinct(data, Study_ID)$Study_ID
# IDs <- IDs[IDs %in% phenos$ID]

#### log normalization ####

# LOGNORM
lognorm <- function(table)
{
  avg <- sum(rowSums(table))/nrow(table)
  table <- sweep(table,1,rowSums(table),"/")
  table <- log10(table*avg + 1)
  return(table)
}

lognorm_out <- lognorm(data[4:length(data)])
lognorm_out <- add_column(lognorm_out, Study_ID=data$Study_ID, .before = colnames(lognorm_out)[1])
lognorm_out <- add_column(lognorm_out, Phenotype=data$Phenotype, .before = colnames(lognorm_out)[1])
lognorm_out <- add_column(lognorm_out, sample_name=data$sample_name, .before = colnames(lognorm_out)[1])

lognorm_out <- lognorm_out %>% filter(Study_ID != 13241) %>% filter(Study_ID != 14812) %>% filter(Study_ID != 11712)
IDs <- distinct(lognorm_out, Study_ID)$Study_ID

#### Get lognorm_out ready for DEBIAS-M ####

df <- lognorm_out
df$ID = 0
for(i in 1:length(IDs)){
  
  # print(distinct(df, ID))
  # print(i)
  df$ID[df$Study_ID == IDs[i]] <- i-1
}

df$case <- case_when(
                   df$Phenotype == 1 ~ TRUE,
                   df$Phenotype == 0 ~ FALSE,
                   )
df <- relocate(df, ID, .after = Study_ID)
df <- relocate(df, case, .after = Study_ID)

write.csv(df, "./csv_files/log_DEBIAS-M.csv", row.names =FALSE)

#### permutations testing as a loop ####



# IDs <- distinct(data, Study_ID)$Study_ID
IDs <- distinct(DEBIAS_data, Study_ID)$Study_ID

# c <- DEBIAS_data[sample(nrow(DEBIAS_data), 300), ]
# c$Phenotype <- sample(c$Phenotype)
# c %>% group_by(Phenotype) %>% count(Study_ID) %>% print(n = 29)

args <- commandArgs(trailingOnly = TRUE)
pdf("./output/test_permutated/post_DEBIAS-M_RF_lognorm.pdf", height = 4, width = 8)
par(mar=c(3,3,1,0), mfrow=c(2,2))

for(i in 1:length(IDs)){
  training <- filter(DEBIAS_data, Study_ID != IDs[i])
  training <- training[, c(2, 4:length(training))]
  testing <- filter(DEBIAS_data, Study_ID == IDs[i])
  testing <- testing[, c(2, 4:length(testing))]
  # 
  # training <- filter(c, Study_ID != IDs[i])
  # training <- training[, c(2, 4:length(training))]
  # testing <- filter(c, Study_ID == IDs[i])
  # testing <- testing[, c(2, 4:length(testing))]

  
  RF_fit <- randomForest(Phenotype~., method = "class", data = training)
      
  RF_pred <- predict(RF_fit, testing, type = "prob")
        
  rf_roc <- roc(testing[,1], RF_pred[,1])
  p <- plot(rf_roc, add = FALSE, col = "red", print.auc = TRUE)
  phen <- filter(phenos, ID == IDs[i])
  phen <- phen[1,2]
  title(paste("Training without:", IDs[i], "(", phen, ")"), line = - 0.5)
   # permutation       
   for(j in 1:4){
       training$Phenotype <- sample(training$Phenotype)
       testing$Phenotype <- sample(testing$Phenotype)
       RF_fit <- randomForest(Phenotype~., method = "class", data = training)
       RF_pred <- predict(RF_fit, testing, type = "prob")
       rf_roc <- roc(testing[,1], RF_pred[,1])
       p <- plot(rf_roc, print.auc=FALSE, add = TRUE)
  
       }
  p
          
  print(paste(i, "of", length(IDs), " done."))
}

dev.off()



#### Dummy Data for troubleshooting ####

auc.df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(auc.df) <-   c("Study_ID", "AUC", "Permutation")

roc.df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(roc.df) <-   c("sensitivities", "specificities", "Permutation")

c <- DEBIAS_data[sample(nrow(DEBIAS_data), 300), ]
c$Phenotype <- sample(c$Phenotype)
c %>% group_by(Phenotype) %>% count(Study_ID) %>% print(n = 29)
IDs <- distinct(c, Study_ID)$Study_ID


pdf("./output/test2.pdf", height = 4, width = 8)
# par(mar=c(3,3,1,0), mfrow=c(2,2))


for(i in 1:length(IDs[1:2])){

  training <- filter(c, Study_ID != IDs[i])
  training <- training[, c(2, 4:length(training))]
  testing <- filter(c, Study_ID == IDs[i])
  testing <- testing[, c(2, 4:length(testing))]

  
  RF_fit <- randomForest(Phenotype~., method = "class", data = training)
  
  RF_pred <- predict(RF_fit, testing, type = "prob")
  
  rf_roc <- roc(testing[,1], RF_pred[,1])
  auc.df[nrow(auc.df) + 1,] = c(IDs[i], auc(rf_roc), FALSE)

  df <- data.frame(sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities, Permutation = 0)
  roc.df <- rbind(roc.df, df)
  phen <- filter(phenos, ID == IDs[i])
  phen <- phen[1,2]
  
  a <- round(auc(rf_roc), 2)
  # p <- ggplot() + ggroc(rf_roc, color="red") + 
  #   annotate("text", x=0.25, y=0.15, label = round(auc(rf_roc), 2)) + 
  #   labs(title = paste("Training without:", IDs[i], "(", phen, ")"))
  
  
  for(j in 1:2){
    training$Phenotype <- sample(training$Phenotype)
    testing$Phenotype <- sample(testing$Phenotype)
    RF_fit <- randomForest(Phenotype~., method = "class", data = training)
    RF_pred <- predict(RF_fit, testing, type = "prob")
    rf_roc <- roc(testing[,1], RF_pred[,1])
    
    df <- data.frame(sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities, Permutation = j)
    roc.df <- rbind(roc.df, df)
    
    auc.df[nrow(auc.df) + 1,] = c(IDs[i], auc(rf_roc), TRUE)
    
    # p <- plot(rf_roc, print.auc=FALSE, add = TRUE)
  }
  
  p <- ggplot() + 
    # geom_line(roc.df, aes(x = specificity, y = sensitivity, group=Permutations)) +
    geom_line(data = filter(roc.df, Permutation == 0), aes(x = specificities, y = sensitivities, group = Permutation), color = "red") +
    geom_line(data = filter(roc.df, Permutation != 0), aes(x = specificities, y = sensitivities))
  
  print(p)
  # current <- auc.df %>% filter(Study_ID == IDs[i])

  # g <- ggplot() + geom_histogram(data = filter(current, Permutation == TRUE), aes(x = AUC), bins = 10) +
  #   geom_vline(filter(current, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") 
  # 
  # plot_grid(p, g)
  # p
  # rm(p)  
  print(paste(i, "of", length(IDs[1:4]), " done."))
}

dev.off()






df_1 = data.frame(x = c(1:100), y = rnorm(100))
df_2 = data.frame(x = c(1:100), y = rnorm(100))
df_3 = data.frame(x = c(1:100), y = rnorm(100))
df_4 = data.frame(x = c(1:100), y = rnorm(100))

all <- df_1
all$y2 <- df_2[['y']]
all$y3 <- df_3[['y']]
all$y4 <- df_4[['y']]
all$x <- rep(c(0,1)) 
all$x <- as.factor(all$x)


pdf("test.pdf", height = 4, width = 12)
par(mar=c(3,3,1,0), mfrow=c(2,2))

for(i in 2:4){
  dat <- all[,i:i+1]
  RF_test <- randomForest(x~dat, method = "class", data = all, importance = TRUE)
  roc_test <- roc(all$x, RF_test$votes[,2]) # same as rf_roc
  plot(roc_test, print.auc=TRUE)
  title(paste("test plot"), line = - 0.5)
  
}

dev.off()

df1 <- data.frame(e1 = sort(runif(5, 0.05, 0.25)),
                  e2 = sort(runif(5, 0.05, 0.25)),
                  e3 = sort(runif(5, 0.05, 0.25)),
                  t1 = sort(runif(5, 1, 100)),
                  t2 = sort(runif(5, 1, 100)),
                  t3 = sort(runif(5, 1, 100))
)
### reshape this to give a column indicating group
df2 <- with(df1,
            as.data.frame(cbind( c(t1, t2, t3),
                                 c(e1, e2, e3),
                                 rep(seq(3), each=5) )
            ))
colnames(df2) <- c("temp","activity","team")

